//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Bryan Clark, bclark@Princeton.edu, Princeton University
//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Paul R. C. Kent, kentpr@ornl.gov, Oak Ridge National Laboratory
//                    Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Peter W. Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////

#include "QMCHamiltonian.h"
#include "Particle/DistanceTable.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCWaveFunctions/Fermion/MultiSlaterDetTableMethod.h"
#include "Utilities/TimerManager.h"
#include "BareKineticEnergy.h"
#include "Containers/MinimalContainers/RecordArray.hpp"
#include "type_traits/ConvertToReal.h"
#include "CPU/math.hpp"

namespace qmcplusplus
{
struct QMCHamiltonian::QMCHamiltonianMultiWalkerResource : public Resource
{
  QMCHamiltonianMultiWalkerResource() : Resource("QMCHamiltonian") {}
  // the listeners represet the connection of a particular crowds estimators to the crowds lead QMCHamiltonian.
  // So you can not clone them.
  std::unique_ptr<Resource> makeClone() const override
  {
    return std::make_unique<QMCHamiltonianMultiWalkerResource>(*this);
  }
  std::vector<ListenerVector<RealType>> kinetic_listeners_;
  std::vector<ListenerVector<RealType>> potential_listeners_;
  std::vector<ListenerVector<RealType>> ion_kinetic_listeners_;
  std::vector<ListenerVector<RealType>> ion_potential_listeners_;
};

/** constructor
*/
QMCHamiltonian::QMCHamiltonian(const std::string& aname)
    : myIndex(0),
      numCollectables(0),
      myName(aname),
      hasPhysicalNLPP_(false),
      l2_ptr(nullptr),
      ham_timer_(createGlobalTimer("Hamiltonian:" + aname + "::evaluate", timer_level_medium)),
      eval_vals_derivs_timer_(createGlobalTimer("Hamiltonian:" + aname + "::ValueParamDerivs", timer_level_medium)),
      eval_ion_derivs_fast_timer_(
          createGlobalTimer("Hamiltonian:" + aname + ":::evaluateIonDerivsFast", timer_level_medium))
#if !defined(REMOVE_TRACEMANAGER)
      ,
      streaming_position(false),
      id_sample(nullptr),
      pid_sample(nullptr),
      step_sample(nullptr),
      gen_sample(nullptr),
      age_sample(nullptr),
      mult_sample(nullptr),
      weight_sample(nullptr),
      position_sample(nullptr)
#endif
{}

///// copy constructor is distable by declaring it as private
//QMCHamiltonian::QMCHamiltonian(const QMCHamiltonian& qh) {}

/** destructor
 */
QMCHamiltonian::~QMCHamiltonian()
{
  //@todo clean up H and auxH
}

bool QMCHamiltonian::get(std::ostream& os) const
{
  for (int i = 0; i < H.size(); i++)
  {
    os << "  " << std::setw(16) << std::left << H[i]->getName();
    H[i]->get(os);
    os << "\n";
  }
  return true;
}

/** add a new Hamiltonian the the list of Hamiltonians.
 * @param h an operator
 * @param aname name of h
 * @param physical if true, a physical operator
 */
void QMCHamiltonian::addOperator(std::unique_ptr<OperatorBase>&& h, const std::string& aname, bool physical)
{
  //change UpdateMode[PHYSICAL] of h so that cloning can be done correctly
  h->getUpdateMode()[OperatorBase::PHYSICAL] = physical;
  if (physical)
  {
    for (int i = 0; i < H.size(); ++i)
    {
      if (H[i]->getName() == aname)
      {
        app_warning() << "QMCHamiltonian::addOperator cannot " << aname << ". The name is already used" << std::endl;
        return;
      }
    }
    app_log() << "  QMCHamiltonian::addOperator " << aname << " to H, physical Hamiltonian " << std::endl;
    h->setName(aname);
    H.push_back(std::move(h));
    std::string tname = "Hamiltonian:" + aname;
    my_timers_.push_back(createGlobalTimer(tname, timer_level_fine));
  }
  else
  {
    //ignore timers for now
    for (int i = 0; i < auxH.size(); ++i)
    {
      if (auxH[i]->getName() == aname)
      {
        app_warning() << "QMCHamiltonian::addOperator cannot " << aname << ". The name is already used" << std::endl;
        return;
      }
    }
    app_log() << "  QMCHamiltonian::addOperator " << aname << " to auxH " << std::endl;
    h->setName(aname);
    auxH.push_back(std::move(h));
  }

  //assign save NLPP if found
  //  name is fixed in ECPotentialBuilder::put()
  if (physical && (aname == "NonLocalECP" || aname == "SOECP"))
    hasPhysicalNLPP_ = true;

  //save L2 potential if found
  //  name is fixed in ECPotentialBuilder::put()
  if (aname == "L2")
  {
    if (l2_ptr == nullptr)
    {
      l2_ptr = physical ? dynamic_cast<L2Potential*>(H.back().get()) : dynamic_cast<L2Potential*>(auxH.back().get());
    }
    else
    {
      APP_ABORT("QMCHamiltonian::addOperator l2_ptr is supposed to be null. Something went wrong!");
    }
  }
}


void QMCHamiltonian::addOperatorType(const std::string& name, const std::string& type)
{
  app_log() << "QMCHamiltonian::addOperatorType added type " << type << " named " << name << std::endl;
  operator_types[name] = type;
}


const std::string& QMCHamiltonian::getOperatorType(const std::string& name)
{
  std::map<std::string, std::string>::iterator type = operator_types.find(name);
  if (type == operator_types.end())
    APP_ABORT("QMCHamiltonian::getOperatorType\n  operator type not found for name " + name);
  return type->second;
}

///** remove a named Hamiltonian from the list
// *@param aname the name of the Hamiltonian
// *@return true, if the request hamiltonian exists and is removed.
// */
//bool
//QMCHamiltonian::remove(const std::string& aname)
//{
//  return false;
//}

void QMCHamiltonian::updateSource(ParticleSet& s)
{
  for (int i = 0; i < H.size(); ++i)
    H[i]->updateSource(s);
  for (int i = 0; i < auxH.size(); ++i)
    auxH[i]->updateSource(s);
}

/** add a number of properties to the ParticleSet
 * @param P ParticleSet to which multiple columns to be added
 *
 * QMCHamiltonian can add any number of properties to a ParticleSet.
 * Hindex contains the index map to the ParticleSet::PropertyList.
 * This enables assigning the properties evaluated by each OperatorBase
 * object to the correct property column.
 */
//void
//QMCHamiltonian::addObservables(PropertySetType& plist)
//{
//  //first add properties to Observables
//  Observables.clear();
//  for(int i=0; i<H.size(); ++i) H[i]->addObservables(Observables);
//  for(int i=0; i<auxH.size(); ++i) auxH[i]->addObservables(Observables);
//
//  myIndex=plist.add(Observables.Names[0]);
//  for(int i=1; i<Observables.size(); ++i) plist.add(Observables.Names[i]);
//
//  app_log() << "  QMCHamiltonian::add2WalkerProperty added " << Observables.size() << " data to PropertyList" << std::endl;
//  app_log() << "    starting Index = " << myIndex << std::endl;
//}
//
int QMCHamiltonian::addObservables(ParticleSet& P)
{
  //first add properties to Observables
  Observables.clear();
  //ParticleSet::mcObservables (large data, e.g. density) are accumulated while evaluations
  P.Collectables.clear();
  P.Collectables.rewind();
  for (int i = 0; i < H.size(); ++i)
    H[i]->addObservables(Observables, P.Collectables);
  for (int i = 0; i < auxH.size(); ++i)
    auxH[i]->addObservables(Observables, P.Collectables);
  myIndex = P.PropertyList.add(Observables.Names[0]);
  for (int i = 1; i < Observables.size(); ++i)
    P.PropertyList.add(Observables.Names[i]);
  numCollectables = P.Collectables.size();
  app_log() << "\n  QMCHamiltonian::add2WalkerProperty added"
            << "\n    " << Observables.size() << " to P::PropertyList "
            << "\n    " << P.Collectables.size() << " to P::Collectables "
            << "\n    starting Index of the observables in P::PropertyList = " << myIndex << std::endl;
  return Observables.size();
}

void QMCHamiltonian::resetObservables(int start, int ncollects)
{
  Observables.clear();
  BufferType collectables;
  collectables.rewind();
  for (int i = 0; i < H.size(); ++i)
    H[i]->addObservables(Observables, collectables);
  for (int i = 0; i < auxH.size(); ++i)
    auxH[i]->addObservables(Observables, collectables);
  if (collectables.size() != ncollects)
  {
    APP_ABORT("  QMCHamiltonian::resetObservables numCollectables != ncollects");
  }
  myIndex         = start;
  numCollectables = ncollects;
}

void QMCHamiltonian::registerObservables(std::vector<ObservableHelper>& h5desc, hdf_archive& file) const
{
  for (int i = 0; i < H.size(); ++i)
    H[i]->registerObservables(h5desc, file);
  for (int i = 0; i < auxH.size(); ++i)
    auxH[i]->registerObservables(h5desc, file);
}

void QMCHamiltonian::registerCollectables(std::vector<ObservableHelper>& h5desc, hdf_archive& file) const
{
  //The physical operators cannot add to collectables
  for (int i = 0; i < auxH.size(); ++i)
    auxH[i]->registerCollectables(h5desc, file);
}

void QMCHamiltonian::mw_registerKineticListener(QMCHamiltonian& ham_leader, ListenerVector<RealType> listener)
{
  // This creates a state replication burder of unknown scope when operators are cloned.
  ham_leader.mw_res_handle_.getResource().kinetic_listeners_.push_back(listener);
}

void QMCHamiltonian::mw_registerLocalEnergyListener(QMCHamiltonian& ham_leader, ListenerVector<RealType> listener)
{
  // This creates a state replication burder of unknown scope when operators are cloned.
  // A local energy listener listens to both the kinetic operator and all involved in the potential.
  ham_leader.mw_res_handle_.getResource().kinetic_listeners_.push_back(listener);
  ham_leader.mw_res_handle_.getResource().potential_listeners_.push_back(listener);
}

void QMCHamiltonian::mw_registerLocalPotentialListener(QMCHamiltonian& ham_leader, ListenerVector<RealType> listener)
{
  // This creates a state replication burder of unknown scope when operators are cloned.
  ham_leader.mw_res_handle_.getResource().potential_listeners_.push_back(listener);
}

void QMCHamiltonian::mw_registerLocalIonPotentialListener(QMCHamiltonian& ham_leader, ListenerVector<RealType> listener)
{
  // This creates a state replication burder of unknown scope when operators are cloned.
  ham_leader.mw_res_handle_.getResource().ion_potential_listeners_.push_back(listener);
}

void QMCHamiltonian::informOperatorsOfListener()
{
  for (int i = 0; i < H.size(); ++i)
    H[i]->informOfPerParticleListener();
}

#if !defined(REMOVE_TRACEMANAGER)
void QMCHamiltonian::initialize_traces(TraceManager& tm, ParticleSet& P)
{
  static bool first_init = true;
  bool trace_log         = first_init && tm.verbose && omp_get_thread_num() == 0;
  if (trace_log)
    app_log() << "\n  Hamiltonian is initializing traces" << std::endl;

  //fill std::string vectors for combined trace quantities
  std::vector<std::string> Eloc;
  std::vector<std::string> Vloc;
  std::vector<std::string> Vq, Vc, Vqq, Vqc, Vcc;
  for (int i = 0; i < H.size(); ++i)
    Eloc.push_back(H[i]->getName());
  for (int i = 1; i < H.size(); ++i)
    Vloc.push_back(H[i]->getName());

  // These contributions are based on the potential energy components.
  // Loop starts at one to skip the kinetic energy component.
  for (int i = 1; i < H.size(); ++i)
  {
    OperatorBase& h = *H[i];
    if (h.isQuantum())
      Vq.push_back(h.getName());
    else if (h.isClassical())
      Vc.push_back(h.getName());
    else if (h.isQuantumQuantum())
      Vqq.push_back(h.getName());
    else if (h.isQuantumClassical())
      Vqc.push_back(h.getName());
    else if (h.isClassicalClassical())
      Vcc.push_back(h.getName());
    else if (omp_get_thread_num() == 0)
      app_log() << "  warning: potential named " << h.getName()
                << " has not been classified according to its quantum domain (q,c,qq,qc,cc)\n    estimators depending "
                   "on this classification may not function properly"
                << std::endl;
  }


  //make trace quantities available
  request.contribute_scalar("id", true);           //default trace quantity
  request.contribute_scalar("parent_id", true);    //default trace quantity
  request.contribute_scalar("step", true);         //default trace quantity
  request.contribute_scalar("generation", true);   //default trace quantity
  request.contribute_scalar("age", true);          //default trace quantity
  request.contribute_scalar("multiplicity", true); //default trace quantity
  request.contribute_scalar("weight", true);       //default trace quantity
  request.contribute_array("position");
  for (int i = 0; i < H.size(); ++i)
    H[i]->contributeTraceQuantities();
  for (int i = 0; i < auxH.size(); ++i)
    auxH[i]->contributeTraceQuantities();


  //note availability of combined quantities
  request.contribute_combined("LocalEnergy", Eloc, true);
  request.contribute_combined("LocalPotential", Vloc, true, true);
  if (Vq.size() > 0)
    request.contribute_combined("Vq", Vq, true, true);
  if (Vc.size() > 0)
    request.contribute_combined("Vc", Vc, true, true);
  if (Vqq.size() > 0)
    request.contribute_combined("Vqq", Vqq, true, true);
  if (Vqc.size() > 0)
    request.contribute_combined("Vqc", Vqc, true, true);
  if (Vcc.size() > 0)
    request.contribute_combined("Vcc", Vcc, true, true);


  //collect trace requests
  std::vector<TraceRequest*> requests;
  //  Hamiltonian request (id, step, weight, positions)
  requests.push_back(&request);
  //  requests from Hamiltonian components
  for (int i = 0; i < H.size(); ++i)
    requests.push_back(&H[i]->getRequest());
  //  requests from other observables
  for (int i = 0; i < auxH.size(); ++i)
    requests.push_back(&auxH[i]->getRequest());

  //collect trace quantity availability/requests from contributors/requestors
  for (int i = 0; i < requests.size(); ++i)
    tm.request.incorporate(*requests[i]);

  //balance requests with availability, mark quantities as streaming/writing
  tm.request.determine_stream_write();

  //relay updated streaming information to all contributors/requestors
  for (int i = 0; i < requests.size(); ++i)
    tm.request.relay_stream_info(*requests[i]);

  //set streaming/writing traces in general
  tm.update_status();

  // setup traces, if any quantities should be streaming

  // tracing
  bool tracing = request.streaming();
  if (tracing != tm.streaming_traces)
    APP_ABORT("QMCHamiltonian::initialize_traces  trace request failed to initialize properly");
  if (!tracing)
  {
    // Empty. Do not log if nothing will be done

    if (trace_log)
      app_log() << "    no traces streaming" << std::endl;
  }
  else
  {
    if (trace_log)
      app_log() << "    traces streaming" << std::endl;
    //checkout trace quantities
    //(requested sources checkout arrays to place samples in for streaming)
    //  checkout walker trace quantities
    streaming_position = request.streaming_array("position");
    if (request.streaming_default_scalars)
    {
      id_sample     = tm.checkout_int<1>("id");
      pid_sample    = tm.checkout_int<1>("parent_id");
      step_sample   = tm.checkout_int<1>("step");
      gen_sample    = tm.checkout_int<1>("generation");
      age_sample    = tm.checkout_int<1>("age");
      mult_sample   = tm.checkout_int<1>("multiplicity");
      weight_sample = tm.checkout_real<1>("weight");
    }
    if (streaming_position)
      position_sample = tm.checkout_real<2>("position", P, DIM);
    //  checkout observable trace quantities
    for (int i = 0; i < H.size(); ++i)
    {
      if (trace_log)
        app_log() << "    OperatorBase::checkoutTraceQuantities  " << H[i]->getName() << std::endl;
      H[i]->checkoutTraceQuantities(tm);
    }
    for (int i = 0; i < auxH.size(); ++i)
    {
      if (trace_log)
        app_log() << "    OperatorBase::checkoutTraceQuantities  " << auxH[i]->getName() << std::endl;
      auxH[i]->checkoutTraceQuantities(tm);
    }
    //setup combined traces that depend on H information
    //  LocalEnergy, LocalPotential, Vq, Vc, Vqq, Vqc, Vcc
    if (Vloc.size() > 0 && request.streaming("LocalPotential"))
      tm.make_combined_trace("LocalPotential", Vloc);
    if (Eloc.size() > 0 && request.streaming("LocalEnergy"))
      tm.make_combined_trace("LocalEnergy", Eloc);
    if (Vq.size() > 0 && request.streaming("Vq"))
      tm.make_combined_trace("Vq", Eloc);
    if (Vc.size() > 0 && request.streaming("Vc"))
      tm.make_combined_trace("Vc", Eloc);
    if (Vqq.size() > 0 && request.streaming("Vqq"))
      tm.make_combined_trace("Vqq", Eloc);
    if (Vqc.size() > 0 && request.streaming("Vqc"))
      tm.make_combined_trace("Vqc", Eloc);
    if (Vcc.size() > 0 && request.streaming("Vcc"))
      tm.make_combined_trace("Vcc", Eloc);

    //all trace samples have been created ( streaming instances)
    //  mark the ones that will be writing also
    tm.screen_writes();

    //observables that depend on traces check them out
    if (trace_log)
      app_log() << "\n  Hamiltonian is fulfilling trace requests from observables" << std::endl;
    for (int i = 0; i < auxH.size(); ++i)
    {
      if (trace_log)
        app_log() << "    OperatorBase::getRequiredTraces  " << auxH[i]->getName() << std::endl;
      auxH[i]->getRequiredTraces(tm);
    }
    //report

    //write traces status to the log
    if (trace_log)
      tm.user_report();

    first_init = false;
  }
}


void QMCHamiltonian::collect_walker_traces(Walker_t& walker, int step)
{
  if (request.streaming_default_scalars)
  {
    (*id_sample)(0)     = walker.getWalkerID();
    (*pid_sample)(0)    = walker.getParentID();
    (*step_sample)(0)   = step;
    (*gen_sample)(0)    = walker.Generation;
    (*age_sample)(0)    = walker.Age;
    (*mult_sample)(0)   = walker.Multiplicity;
    (*weight_sample)(0) = walker.Weight;
  }
  if (streaming_position)
    for (int i = 0; i < walker.R.size(); ++i)
      for (int d = 0; d < DIM; ++d)
        (*position_sample)(i, d) = walker.R[i][d];
}


void QMCHamiltonian::finalize_traces()
{
  if (request.streaming_default_scalars)
  {
    delete id_sample;
    delete pid_sample;
    delete step_sample;
    delete gen_sample;
    delete age_sample;
    delete mult_sample;
    delete weight_sample;
  }
  if (streaming_position)
    delete position_sample;
  if (request.streaming())
  {
    for (int i = 0; i < H.size(); ++i)
      H[i]->deleteTraceQuantities();
    for (int i = 0; i < auxH.size(); ++i)
      auxH[i]->deleteTraceQuantities();
  }
  streaming_position = false;
  request.reset();
}
#endif

/** Evaluate all the Hamiltonians for the N-particle  configuration
 *@param P input configuration containing N particles
 *@return the local energy
 */
QMCHamiltonian::FullPrecRealType QMCHamiltonian::evaluate(ParticleSet& P)
{
  ScopedTimer local_timer(ham_timer_);
  LocalEnergy = 0.0;
  for (int i = 0; i < H.size(); ++i)
  {
    ScopedTimer h_timer(my_timers_[i]);
    H[i]->evaluate(P);
    updateComponent(*H[i], *this, P);
#if !defined(REMOVE_TRACEMANAGER)
    H[i]->collectScalarTraces();
#endif
  }
  updateKinetic(*this, P);
  return LocalEnergy;
}

QMCHamiltonian::FullPrecRealType QMCHamiltonian::evaluateDeterministic(ParticleSet& P)
{
  ScopedTimer local_timer(ham_timer_);
  LocalEnergy = 0.0;
  for (int i = 0; i < H.size(); ++i)
  {
    ScopedTimer h_timer(my_timers_[i]);
    H[i]->evaluateDeterministic(P);
    updateComponent(*H[i], *this, P);
#if !defined(REMOVE_TRACEMANAGER)
    H[i]->collectScalarTraces();
#endif
  }
  updateKinetic(*this, P);
  return LocalEnergy;
}
void QMCHamiltonian::updateComponent(OperatorBase& op, QMCHamiltonian& ham, ParticleSet& pset)
{
  // It's much better to be able to see where this is coming from.  It is caught just fine.
  if (qmcplusplus::isnan(op.getValue()))
  {
    std::ostringstream msg;
    msg << "QMCHamiltonian::updateComponent component " << op.getName() << " returns NaN." << std::endl;
    pset.print(msg);
    throw std::runtime_error(msg.str());
  }
  // The following is a ridiculous breach of encapsulation.
  ham.LocalEnergy += op.getValue();
  op.setObservables(ham.Observables);
  op.setParticlePropertyList(pset.PropertyList, ham.myIndex);
}

void QMCHamiltonian::updateKinetic(QMCHamiltonian& ham, ParticleSet& pset)
{
  ham.KineticEnergy                     = ham.H[0]->getValue();
  pset.PropertyList[WP::LOCALENERGY]    = ham.LocalEnergy;
  pset.PropertyList[WP::LOCALPOTENTIAL] = ham.LocalEnergy - ham.KineticEnergy;
}

std::vector<QMCHamiltonian::FullPrecRealType> QMCHamiltonian::mw_evaluate(
    const RefVectorWithLeader<QMCHamiltonian>& ham_list,
    const RefVectorWithLeader<TrialWaveFunction>& wf_list,
    const RefVectorWithLeader<ParticleSet>& p_list)
{
  auto& ham_leader = ham_list.getLeader();
  ScopedTimer local_timer(ham_leader.ham_timer_);
  for (QMCHamiltonian& ham : ham_list)
    ham.LocalEnergy = 0.0;

  const int num_ham_operators = ham_leader.H.size();

  // It is an invariant of the class that H[0] be the kinetic energy operator.
  // This is enforced by HamiltonianFactory's constuctor and not QMCHamiltonians
  int kinetic_index = 0;
  {
    ScopedTimer h_timer(ham_leader.my_timers_[kinetic_index]);
    const auto HC_list(extract_HC_list(ham_list, kinetic_index));
    if (ham_leader.mw_res_handle_.getResource().kinetic_listeners_.size() > 0)
      ham_leader.H[kinetic_index]
          ->mw_evaluatePerParticle(HC_list, wf_list, p_list, ham_leader.mw_res_handle_.getResource().kinetic_listeners_,
                                   ham_leader.mw_res_handle_.getResource().ion_kinetic_listeners_);
    else
      ham_leader.H[kinetic_index]->mw_evaluate(HC_list, wf_list, p_list);
    for (int iw = 0; iw < ham_list.size(); iw++)
      updateComponent(HC_list[iw], ham_list[iw], p_list[iw]);
  }

  for (int i_ham_op = 1; i_ham_op < num_ham_operators; ++i_ham_op)
  {
    ScopedTimer h_timer(ham_leader.my_timers_[i_ham_op]);
    const auto HC_list(extract_HC_list(ham_list, i_ham_op));

    if (ham_leader.mw_res_handle_.getResource().potential_listeners_.size() > 0)
      ham_leader.H[i_ham_op]->mw_evaluatePerParticle(HC_list, wf_list, p_list,
                                                     ham_leader.mw_res_handle_.getResource().potential_listeners_,
                                                     ham_leader.mw_res_handle_.getResource().ion_potential_listeners_);
    else
      ham_leader.H[i_ham_op]->mw_evaluate(HC_list, wf_list, p_list);
    for (int iw = 0; iw < ham_list.size(); iw++)
      updateComponent(HC_list[iw], ham_list[iw], p_list[iw]);
  }

  for (int iw = 0; iw < ham_list.size(); iw++)
    updateKinetic(ham_list[iw], p_list[iw]);

  std::vector<FullPrecRealType> local_energies(ham_list.size(), 0.0);
  for (int iw = 0; iw < ham_list.size(); ++iw)
    local_energies[iw] = ham_list[iw].getLocalEnergy();

  return local_energies;
}

QMCHamiltonian::FullPrecRealType QMCHamiltonian::evaluateValueAndDerivatives(ParticleSet& P,
                                                                             const opt_variables_type& optvars,
                                                                             Vector<ValueType>& dlogpsi,
                                                                             Vector<ValueType>& dhpsioverpsi)
{
  // The first componennt must be BareKineticEnergy for both handling KineticEnergy and dlogpsi computation
  // by calling TWF::evaluateDerivatives inside BareKineticEnergy::evaluateValueAndDerivatives
  assert(dynamic_cast<BareKineticEnergy*>(H[0].get()) &&
         "BUG: The first componennt in Hamiltonian must be BareKineticEnergy.");
  ScopedTimer local_timer(eval_vals_derivs_timer_);

  {
    ScopedTimer h_timer(my_timers_[0]);
    LocalEnergy = KineticEnergy = H[0]->evaluateValueAndDerivatives(P, optvars, dlogpsi, dhpsioverpsi);
  }

  for (int i = 1; i < H.size(); ++i)
  {
    ScopedTimer h_timer(my_timers_[i]);
    LocalEnergy += H[i]->evaluateValueAndDerivatives(P, optvars, dlogpsi, dhpsioverpsi);
  }
  return LocalEnergy;
}

std::vector<QMCHamiltonian::FullPrecRealType> QMCHamiltonian::mw_evaluateValueAndDerivatives(
    const RefVectorWithLeader<QMCHamiltonian>& ham_list,
    const RefVectorWithLeader<TrialWaveFunction>& wf_list,
    const RefVectorWithLeader<ParticleSet>& p_list,
    const opt_variables_type& optvars,
    RecordArray<ValueType>& dlogpsi,
    RecordArray<ValueType>& dhpsioverpsi)
{
  std::vector<FullPrecRealType> local_energies(ham_list.size(), 0.0);
  for (int iw = 0; iw < ham_list.size(); iw++)
    ham_list[iw].LocalEnergy = 0.0;

  if (ham_list.size() > 0)
  {
    auto& ham_leader            = ham_list.getLeader();
    const int num_ham_operators = ham_leader.H.size();
    for (int i_ham_op = 0; i_ham_op < num_ham_operators; ++i_ham_op)
    {
      ScopedTimer local_timer(ham_leader.my_timers_[i_ham_op]);
      const auto HC_list(extract_HC_list(ham_list, i_ham_op));

      ham_leader.H[i_ham_op]->mw_evaluateWithParameterDerivatives(HC_list, p_list, optvars, dlogpsi, dhpsioverpsi);

      for (int iw = 0; iw < ham_list.size(); iw++)
        updateComponent(HC_list[iw], ham_list[iw], p_list[iw]);
    }

    for (int iw = 0; iw < ham_list.size(); iw++)
      updateKinetic(ham_list[iw], p_list[iw]);

    for (int iw = 0; iw < ham_list.size(); ++iw)
      local_energies[iw] = ham_list[iw].getLocalEnergy();
  }

  return local_energies;
}

QMCHamiltonian::FullPrecRealType QMCHamiltonian::evaluateVariableEnergy(ParticleSet& P, bool free_nlpp)
{
  RealType nlpp = 0.0;
  RealType ke   = H[0]->evaluate(P);
  if (free_nlpp)
    for (int i = 1; i < H.size(); ++i)
    {
      if (H[i]->isNonLocal())
        nlpp += H[i]->evaluate(P);
    }
  return ke + nlpp;
}

void QMCHamiltonian::auxHevaluate(ParticleSet& P)
{
  for (int i = 0; i < auxH.size(); ++i)
  {
    RealType sink = auxH[i]->evaluate(P);
    auxH[i]->setObservables(Observables);
#if !defined(REMOVE_TRACEMANAGER)
    auxH[i]->collectScalarTraces();
#endif
    auxH[i]->setParticlePropertyList(P.PropertyList, myIndex);
    //H[i]->setParticlePropertyList(P.PropertyList,myIndex);
  }
}

///This is more efficient. Only calculate auxH elements if moves are accepted.
void QMCHamiltonian::auxHevaluate(ParticleSet& P, Walker_t& ThisWalker)
{
#if !defined(REMOVE_TRACEMANAGER)
  collect_walker_traces(ThisWalker, P.current_step);
#endif
  for (int i = 0; i < auxH.size(); ++i)
  {
    auxH[i]->setHistories(ThisWalker);
    RealType sink = auxH[i]->evaluate(P);
    auxH[i]->setObservables(Observables);
#if !defined(REMOVE_TRACEMANAGER)
    auxH[i]->collectScalarTraces();
#endif
    auxH[i]->setParticlePropertyList(P.PropertyList, myIndex);
  }
}
///Evaluate properties only.
void QMCHamiltonian::auxHevaluate(ParticleSet& P, Walker_t& ThisWalker, bool do_properties, bool do_collectables)
{
#if !defined(REMOVE_TRACEMANAGER)
  collect_walker_traces(ThisWalker, P.current_step);
#endif
  for (int i = 0; i < auxH.size(); ++i)
  {
    bool is_property    = !(auxH[i]->getMode(OperatorBase::COLLECTABLE));
    bool is_collectable = (auxH[i]->getMode(OperatorBase::COLLECTABLE));
    if ((is_property && do_properties) || (is_collectable && do_collectables))
    {
      auxH[i]->setHistories(ThisWalker);
      RealType sink = auxH[i]->evaluate(P);
      auxH[i]->setObservables(Observables);
#if !defined(REMOVE_TRACEMANAGER)
      auxH[i]->collectScalarTraces();
#endif
      auxH[i]->setParticlePropertyList(P.PropertyList, myIndex);
    }
  }
}

/** Looks like a hack see DMCBatched.cpp and DMC.cpp weight is used like temporary flag
 *  from DMC.
 */
void QMCHamiltonian::rejectedMove(ParticleSet& P, Walker_t& ThisWalker)
{
  // weight should be 0 from DMC
  //   since other traced properties will be invalid
  //   (they will be from the walker moved before this one)
#if !defined(REMOVE_TRACEMANAGER)
  collect_walker_traces(ThisWalker, P.current_step);
#endif
  //   ThisWalker.rejectedMove();
  for (int i = 0; i < auxH.size(); ++i)
  {
    auxH[i]->setHistories(ThisWalker);
    RealType sink = auxH[i]->rejectedMove(P);
  }
}

QMCHamiltonian::FullPrecRealType QMCHamiltonian::evaluateWithToperator(ParticleSet& P)
{
  ScopedTimer local_timer(ham_timer_);
  LocalEnergy = 0.0;
  for (int i = 0; i < H.size(); ++i)
  {
    ScopedTimer h_timer(my_timers_[i]);
    H[i]->evaluateWithToperator(P);
    updateComponent(*H[i], *this, P);
#if !defined(REMOVE_TRACEMANAGER)
    H[i]->collectScalarTraces();
#endif
  }
  updateKinetic(*this, P);
  return LocalEnergy;
}

std::vector<QMCHamiltonian::FullPrecRealType> QMCHamiltonian::mw_evaluateWithToperator(
    const RefVectorWithLeader<QMCHamiltonian>& ham_list,
    const RefVectorWithLeader<TrialWaveFunction>& wf_list,
    const RefVectorWithLeader<ParticleSet>& p_list)
{
  for (QMCHamiltonian& ham : ham_list)
    ham.LocalEnergy = 0.0;

  auto& ham_leader            = ham_list.getLeader();
  const int num_ham_operators = ham_leader.H.size();

  const int kinetic_index = 0;
  {
    ScopedTimer h_timer(ham_leader.my_timers_[kinetic_index]);
    const auto HC_list(extract_HC_list(ham_list, kinetic_index));
    if (ham_leader.mw_res_handle_.getResource().kinetic_listeners_.size() > 0)
      ham_leader.H[kinetic_index]
          ->mw_evaluatePerParticleWithToperator(HC_list, wf_list, p_list,
                                                ham_leader.mw_res_handle_.getResource().kinetic_listeners_,
                                                ham_leader.mw_res_handle_.getResource().ion_kinetic_listeners_);
    else
      ham_leader.H[kinetic_index]->mw_evaluateWithToperator(HC_list, wf_list, p_list);
    for (int iw = 0; iw < ham_list.size(); iw++)
      updateComponent(HC_list[iw], ham_list[iw], p_list[iw]);
  }

  for (int i_ham_op = 1; i_ham_op < num_ham_operators; ++i_ham_op)
  {
    ScopedTimer local_timer(ham_leader.my_timers_[i_ham_op]);
    const auto HC_list(extract_HC_list(ham_list, i_ham_op));
    if (ham_leader.mw_res_handle_.getResource().potential_listeners_.size() > 0)
      ham_leader.H[i_ham_op]
          ->mw_evaluatePerParticleWithToperator(HC_list, wf_list, p_list,
                                                ham_leader.mw_res_handle_.getResource().potential_listeners_,
                                                ham_leader.mw_res_handle_.getResource().ion_potential_listeners_);
    else
      ham_leader.H[i_ham_op]->mw_evaluateWithToperator(HC_list, wf_list, p_list);
    for (int iw = 0; iw < ham_list.size(); ++iw)
      updateComponent(HC_list[iw], ham_list[iw], p_list[iw]);
  }

  for (int iw = 0; iw < ham_list.size(); iw++)
    updateKinetic(ham_list[iw], p_list[iw]);

  std::vector<FullPrecRealType> local_energies(ham_list.size());
  for (int iw = 0; iw < ham_list.size(); ++iw)
    local_energies[iw] = ham_list[iw].getLocalEnergy();

  return local_energies;
}
void QMCHamiltonian::evaluateElecGrad(ParticleSet& P,
                                      TrialWaveFunction& psi,
                                      ParticleSet::ParticlePos& Egrad,
                                      RealType delta)
{
  int nelec = P.getTotalNum();
  RealType ep(0.0);
  RealType em(0.0);
  RealType e0(0.0);
  for (int iel = 0; iel < nelec; iel++)
  {
    for (int dim = 0; dim < OHMMS_DIM; dim++)
    {
      RealType r0 = P.R[iel][dim];
      ep          = 0;
      em          = 0;
      //Plus
      RealType rp   = r0 + delta;
      P.R[iel][dim] = rp;
      P.update();
      psi.evaluateLog(P);
      ep = evaluateDeterministic(P);

      //minus
      RealType rm   = r0 - delta;
      P.R[iel][dim] = rm;
      P.update();
      psi.evaluateLog(P);
      em = evaluateDeterministic(P);

      Egrad[iel][dim] = (ep - em) / (2.0 * delta);
      P.R[iel][dim]   = r0;
      P.update();
      psi.evaluateLog(P);
    }
  }
}

void QMCHamiltonian::evaluateIonDerivs(ParticleSet& P,
                                       ParticleSet& ions,
                                       TrialWaveFunction& psi,
                                       ParticleSet::ParticlePos& hf_term,
                                       ParticleSet::ParticlePos& pulay_terms,
                                       ParticleSet::ParticlePos& wf_grad)
{
  ParticleSet::ParticleGradient wfgradraw_(ions.getTotalNum());
  wfgradraw_ = 0.0;

  for (int i = 0; i < H.size(); ++i)
    H[i]->evaluateIonDerivs(P, ions, psi, hf_term, pulay_terms);

  for (int iat = 0; iat < ions.getTotalNum(); iat++)
  {
    wfgradraw_[iat] = psi.evalGradSource(P, ions, iat);
    convertToReal(wfgradraw_[iat], wf_grad[iat]);
  }
}

QMCHamiltonian::FullPrecRealType QMCHamiltonian::getEnsembleAverage()
{
  FullPrecRealType sum = 0.0;
  for (int i = 0; i < H.size(); i++)
    sum += H[i]->getEnsembleAverage();
  return sum;
}

/** return pointer to the QMCHamtiltonian with the name
 *@param aname the name of Hamiltonian
 *@return the pointer to the named term.
 *
 * If not found, return 0
 */
OperatorBase* QMCHamiltonian::getHamiltonian(const std::string& aname)
{
  for (int i = 0; i < H.size(); ++i)
    if (H[i]->getName() == aname)
      return H[i].get();
  for (int i = 0; i < auxH.size(); ++i)
    if (auxH[i]->getName() == aname)
      return auxH[i].get();
  return nullptr;
}

RefVector<OperatorBase> QMCHamiltonian::getTWFDependentComponents()
{
  RefVector<OperatorBase> components;
  for (int i = 0; i < H.size(); i++)
    if (H[i]->dependsOnWaveFunction())
      components.push_back(*H[i]);
  return components;
}

void QMCHamiltonian::resetTargetParticleSet(ParticleSet& P)
{
  for (int i = 0; i < H.size(); i++)
    H[i]->resetTargetParticleSet(P);
  for (int i = 0; i < auxH.size(); i++)
    auxH[i]->resetTargetParticleSet(P);
}

void QMCHamiltonian::setRandomGenerator(RandomBase<FullPrecRealType>* rng)
{
  for (int i = 0; i < H.size(); i++)
    H[i]->setRandomGenerator(rng);
  for (int i = 0; i < auxH.size(); i++)
    auxH[i]->setRandomGenerator(rng);
}

int QMCHamiltonian::makeNonLocalMoves(ParticleSet& P, NonLocalTOperator& move_op)
{
  int num_moves = 0;
  for (int i = 0; i < H.size(); ++i)
    num_moves += H[i]->makeNonLocalMovesPbyP(P, move_op);
  return num_moves;
}


std::vector<int> QMCHamiltonian::mw_makeNonLocalMoves(const RefVectorWithLeader<QMCHamiltonian>& ham_list,
                                                      const RefVectorWithLeader<TrialWaveFunction>& wf_list,
                                                      const RefVectorWithLeader<ParticleSet>& p_list,
                                                      NonLocalTOperator& move_op)
{
  std::vector<int> num_accepts(ham_list.size(), 0);
  for (int iw = 0; iw < ham_list.size(); ++iw)
    num_accepts[iw] = ham_list[iw].makeNonLocalMoves(p_list[iw], move_op);
  return num_accepts;
}

void QMCHamiltonian::createResource(ResourceCollection& collection) const
{
  auto resource_index = collection.addResource(std::make_unique<QMCHamiltonianMultiWalkerResource>());
  for (int i = 0; i < H.size(); ++i)
    H[i]->createResource(collection);
}

void QMCHamiltonian::acquireResource(ResourceCollection& collection,
                                     const RefVectorWithLeader<QMCHamiltonian>& ham_list)
{
  auto& ham_leader          = ham_list.getLeader();
  ham_leader.mw_res_handle_ = collection.lendResource<QMCHamiltonianMultiWalkerResource>();
  for (int i_ham_op = 0; i_ham_op < ham_leader.H.size(); ++i_ham_op)
  {
    const auto HC_list(extract_HC_list(ham_list, i_ham_op));
    ham_leader.H[i_ham_op]->acquireResource(collection, HC_list);
  }
}

void QMCHamiltonian::releaseResource(ResourceCollection& collection,
                                     const RefVectorWithLeader<QMCHamiltonian>& ham_list)
{
  auto& ham_leader = ham_list.getLeader();
  collection.takebackResource(ham_leader.mw_res_handle_);
  for (int i_ham_op = 0; i_ham_op < ham_leader.H.size(); ++i_ham_op)
  {
    const auto HC_list(extract_HC_list(ham_list, i_ham_op));
    ham_leader.H[i_ham_op]->releaseResource(collection, HC_list);
  }
}

std::unique_ptr<QMCHamiltonian> QMCHamiltonian::makeClone(ParticleSet& qp, TrialWaveFunction& psi) const
{
  auto myclone = std::make_unique<QMCHamiltonian>(myName);
  for (int i = 0; i < H.size(); ++i)
    H[i]->add2Hamiltonian(qp, psi, *myclone);
  for (int i = 0; i < auxH.size(); ++i)
    auxH[i]->add2Hamiltonian(qp, psi, *myclone);
  //sync indices
  myclone->resetObservables(myIndex, numCollectables);
  //Hamiltonian needs to make sure qp.Collectables are the same as defined by the original Hamiltonian
  if (numCollectables)
  {
    qp.Collectables.clear();
    qp.Collectables.resize(numCollectables);
  }
  //Assume tau is correct for the Kinetic energy operator and assign to the rest of the clones
  //Return_t tau = H[0]->Tau;
  //myclone->setTau(tau);
  return myclone;
}

RefVectorWithLeader<OperatorBase> QMCHamiltonian::extract_HC_list(const RefVectorWithLeader<QMCHamiltonian>& ham_list,
                                                                  int id)
{
  RefVectorWithLeader<OperatorBase> HC_list(*ham_list.getLeader().H[id]);
  HC_list.reserve(ham_list.size());
  for (QMCHamiltonian& H : ham_list)
    HC_list.push_back(*(H.H[id]));
  return HC_list;
}

void QMCHamiltonian::evaluateIonDerivsFast(ParticleSet& P,
                                           ParticleSet& ions,
                                           TrialWaveFunction& psi_in,
                                           TWFFastDerivWrapper& psi_wrapper_in,
                                           ParticleSet::ParticlePos& dEdR,
                                           ParticleSet::ParticlePos& wf_grad)
{
  ScopedTimer local_timer(eval_ion_derivs_fast_timer_);
  P.update();
  //resize everything;
  const int ngroups = psi_wrapper_in.numGroups();

  std::vector<ValueMatrix> X_;       // [Nptcl, Nptcl] auxiliary X matrix (Minv.B.Minv)
  std::vector<ValueMatrix> Minv_;    // [Nptcl, Nptcl] inverse slater matrix over GS orbs
  std::vector<ValueMatrix> B_;       // [Nptcl, Norb] B matrix (Op(M)) over all orbs
  std::vector<ValueMatrix> B_gs_;    // [Nptcl, Nocc] subset of B over GS orbs
  std::vector<ValueMatrix> M_;       // [Nptcl, Norb] slater matrix over all orbs
  std::vector<ValueMatrix> M_gs_;    // [Nptcl, Nocc] subset of M over GS orbs
  std::vector<ValueMatrix> Minv_B_;  // [Nptcl, Norb] Minv.B
  std::vector<ValueMatrix> Minv_Mv_; // [Nptcl, Nvirt] Minv.M over virtual orbs

  std::vector<std::vector<ValueMatrix>> dM_;      // Derivative of slater matrix
  std::vector<std::vector<ValueMatrix>> dB_;      // Derivative of B matrix
  std::vector<std::vector<ValueMatrix>> dM_gs_;   // subset of dM over GS orbs
  std::vector<std::vector<ValueMatrix>> dB_gs_;   // subset of dB over GS orbs
  std::vector<std::vector<ValueMatrix>> Minv_dM_; // Minv.dM
  std::vector<std::vector<ValueMatrix>> Minv_dB_; // Minv.dB

  {
    M_.resize(ngroups);
    M_gs_.resize(ngroups);
    X_.resize(ngroups);
    B_.resize(ngroups);
    B_gs_.resize(ngroups);
    Minv_.resize(ngroups);
    Minv_B_.resize(ngroups);
    Minv_Mv_.resize(ngroups);

    for (int gid = 0; gid < ngroups; gid++)
    {
      const int sid    = psi_wrapper_in.getTWFGroupIndex(gid);
      const int norbs  = psi_wrapper_in.numOrbitals(sid);
      const int first  = P.first(gid);
      const int last   = P.last(gid);
      const int nptcls = last - first;
      const int nvirt  = norbs - nptcls;

      M_[sid].resize(nptcls, norbs);
      B_[sid].resize(nptcls, norbs);
      Minv_B_[sid].resize(nptcls, norbs);
      Minv_Mv_[sid].resize(nptcls, nvirt);

      M_gs_[sid].resize(nptcls, nptcls);
      Minv_[sid].resize(nptcls, nptcls);
      B_gs_[sid].resize(nptcls, nptcls);
      X_[sid].resize(nptcls, nptcls);
    }

    dM_.resize(OHMMS_DIM);
    dB_.resize(OHMMS_DIM);
    dM_gs_.resize(OHMMS_DIM);
    dB_gs_.resize(OHMMS_DIM);
    Minv_dM_.resize(OHMMS_DIM);
    Minv_dB_.resize(OHMMS_DIM);

    for (int idim = 0; idim < OHMMS_DIM; idim++)
    {
      dM_[idim].resize(ngroups);
      dB_[idim].resize(ngroups);
      dM_gs_[idim].resize(ngroups);
      dB_gs_[idim].resize(ngroups);
      Minv_dM_[idim].resize(ngroups);
      Minv_dB_[idim].resize(ngroups);

      for (int gid = 0; gid < ngroups; gid++)
      {
        const int sid    = psi_wrapper_in.getTWFGroupIndex(gid);
        const int norbs  = psi_wrapper_in.numOrbitals(sid);
        const int first  = P.first(gid);
        const int last   = P.last(gid);
        const int nptcls = last - first;

        dM_[idim][sid].resize(nptcls, norbs);
        dB_[idim][sid].resize(nptcls, norbs);
        dM_gs_[idim][sid].resize(nptcls, nptcls);
        dB_gs_[idim][sid].resize(nptcls, nptcls);
        Minv_dM_[idim][sid].resize(nptcls, norbs);
        Minv_dB_[idim][sid].resize(nptcls, norbs);
      }
    }
    psi_wrapper_in.wipeMatrices(M_);
    psi_wrapper_in.wipeMatrices(M_gs_);
    psi_wrapper_in.wipeMatrices(X_);
    psi_wrapper_in.wipeMatrices(B_);
    psi_wrapper_in.wipeMatrices(Minv_);
    psi_wrapper_in.wipeMatrices(B_gs_);
    psi_wrapper_in.wipeMatrices(Minv_B_);
    psi_wrapper_in.wipeMatrices(Minv_Mv_);

    for (int idim = 0; idim < OHMMS_DIM; idim++)
    {
      psi_wrapper_in.wipeMatrices(dM_[idim]);
      psi_wrapper_in.wipeMatrices(dB_[idim]);
      psi_wrapper_in.wipeMatrices(dM_gs_[idim]);
      psi_wrapper_in.wipeMatrices(dB_gs_[idim]);
      psi_wrapper_in.wipeMatrices(Minv_dM_[idim]);
      psi_wrapper_in.wipeMatrices(Minv_dB_[idim]);
    }
  }
  ParticleSet::ParticleGradient wfgradraw_(ions.getTotalNum());
  ParticleSet::ParticleGradient pulay_(ions.getTotalNum());
  ParticleSet::ParticleGradient hf_(ions.getTotalNum());
  ParticleSet::ParticleGradient dedr_complex(ions.getTotalNum());
  ParticleSet::ParticlePos pulayterms_(ions.getTotalNum());
  ParticleSet::ParticlePos hfdiag_(ions.getTotalNum());
  wfgradraw_ = 0.0;

  {
    psi_wrapper_in.getM(P, M_);
    psi_wrapper_in.getGSMatrices(M_, M_gs_);
    psi_wrapper_in.invertMatrices(M_gs_, Minv_);
  }

  //Build B-matrices.  Only for non-diagonal observables right now.
  for (int i = 0; i < H.size(); ++i)
    if (H[i]->dependsOnWaveFunction())
      H[i]->evaluateOneBodyOpMatrix(P, psi_wrapper_in, B_);
    else
      H[i]->evaluateIonDerivs(P, ions, psi_in, hfdiag_, pulayterms_);


  {
    psi_wrapper_in.getGSMatrices(B_, B_gs_);
    // X_ is now built by buildIntermediates; could go back to buildX depending on future refactoring
    // psi_wrapper_in.buildX(Minv_, B_gs_, X_);
    psi_wrapper_in.buildIntermediates(Minv_, B_, M_, X_, Minv_B_, Minv_Mv_);
  }


  // ===== Initialize some multidet quantities =====

  // values for MultiDiracDet i, excited det j (also include GS det at j==0)
  std::vector<Vector<ValueType>> fvals_O;     // (O D[i][j]/D[i][j])
  std::vector<Vector<ValueType>> fvals_dmu;   // d/dmu(log(D[i][j])
  std::vector<Vector<ValueType>> fvals_dmu_O; // d/dmu(O D[i][j]/D[i][j])

  // same order as Dets in msd; index of associated SPOset in psi_wrapper_in.sposets_
  std::vector<int> mdd_spo_ids;
  std::vector<const WaveFunctionComponent*> mdd_list;

  if (psi_wrapper_in.hasMultiSlaterDet())
  {
    const auto& msd = static_cast<const MultiSlaterDetTableMethod&>(psi_wrapper_in.getMultiSlaterDet());

    auto n_mdd = msd.getDetSize();
    fvals_O.resize(n_mdd);
    fvals_dmu.resize(n_mdd);
    fvals_dmu_O.resize(n_mdd);

    for (size_t i_mdd = 0; i_mdd < n_mdd; i_mdd++)
    {
      const MultiDiracDeterminant& multidiracdet_i = msd.getDet(i_mdd);
      mdd_list.push_back(static_cast<const WaveFunctionComponent*>(&multidiracdet_i));
      // particle group id for this multidiracdet
      const int gid = P.getGroupID(multidiracdet_i.getFirstIndex());
      // SPOSet location in psi_wrapper_in.sposets_ for this particle group
      const int sid = psi_wrapper_in.getTWFGroupIndex(gid);
      mdd_spo_ids.push_back(sid);
      fvals_O[i_mdd].resize(multidiracdet_i.getNumDets());
      fvals_dmu[i_mdd].resize(multidiracdet_i.getNumDets());
      fvals_dmu_O[i_mdd].resize(multidiracdet_i.getNumDets());
    }
    psi_wrapper_in.wipeVectors(fvals_O);

    // compute (OD/D) for all excited DiracDets D
    psi_wrapper_in.computeMDDerivatives_Obs(Minv_Mv_, Minv_B_, mdd_spo_ids, mdd_list, fvals_O);
  }


  //And now we compute the 3N force derivatives.  3 at a time for each atom.
  for (int iat = 0; iat < ions.getTotalNum(); iat++)
  {
    //The total wavefunction derivative has two contributions.  One from determinantal piece,
    //One from the Jastrow.  Jastrow is easy, so we evaluate it here, then add on the
    //determinantal piece at the end of this block.

    wfgradraw_[iat] = psi_wrapper_in.evaluateJastrowGradSource(P, ions, iat);
    for (int idim = 0; idim < OHMMS_DIM; idim++)
    {
      psi_wrapper_in.wipeMatrices(dM_[idim]);
      psi_wrapper_in.wipeMatrices(dB_[idim]);
      psi_wrapper_in.wipeMatrices(dM_gs_[idim]);
      psi_wrapper_in.wipeMatrices(dB_gs_[idim]);
      psi_wrapper_in.wipeMatrices(Minv_dM_[idim]);
      psi_wrapper_in.wipeMatrices(Minv_dB_[idim]);
    }

    {
      //ion derivative of slater matrix.
      psi_wrapper_in.getIonGradM(P, ions, iat, dM_);
    }

    for (int i = 0; i < H.size(); ++i)
      if (H[i]->dependsOnWaveFunction())
        H[i]->evaluateOneBodyOpMatrixForceDeriv(P, ions, psi_wrapper_in, iat, dB_);


    psi_wrapper_in.buildIntermediates_dmu(Minv_, dB_, dM_, Minv_dB_, Minv_dM_);

    for (int idim = 0; idim < OHMMS_DIM; idim++)
    {
      psi_wrapper_in.getGSMatrices(dB_[idim], dB_gs_[idim]);
      psi_wrapper_in.getGSMatrices(dM_[idim], dM_gs_[idim]);

      psi_wrapper_in.wipeVectors(fvals_dmu);
      psi_wrapper_in.wipeVectors(fvals_dmu_O);

      // d_mu(OPsi/Psi)
      ValueType fval_dmu_O = 0.0;
      // d_mu(log Psi)
      ValueType fval_dmu = 0.0;

      // this is OD/D terms, which does not depend on d_mu, so will be same for every ionid/dim
      // we don't use it here, but computeMDDerivatives_total returns it (and needs it internally anyway, so it adds no additional cost)
      ValueType fval_O = 0.0;


      /// NOTE: depending on which matrix products we are precomputing, we may be able to save a factor of Nptcl here
      ///       (Minv.dB - X.dM) == (Minv_dB - Minv_B.Minv_dM)
      ///       not worth doing the extra gemm if we only need the trace here,
      ///       but if we're already doing the gemm for the multidet terms, then we can use that here
      fval_dmu_O = psi_wrapper_in.computeGSDerivative(Minv_, X_, dM_gs_[idim], dB_gs_[idim]);

      /// TODO: can save a factor of nptcl here because we already have Minv_dM
      fval_dmu = psi_wrapper_in.trAB(Minv_, dM_gs_[idim]);

      dedr_complex[iat][idim] = fval_dmu_O;
      wfgradraw_[iat][idim] += fval_dmu;

      if (psi_wrapper_in.hasMultiSlaterDet())
      {
        const auto& msd = static_cast<const MultiSlaterDetTableMethod&>(psi_wrapper_in.getMultiSlaterDet());
        // compute d_mu(OD/D) and d_mu(logD) for all excited DiracDets D
        psi_wrapper_in.computeMDDerivatives_dmu(Minv_Mv_, Minv_B_, Minv_dM_[idim], Minv_dB_[idim], mdd_spo_ids,
                                                mdd_list, fvals_dmu_O, fvals_dmu);

        // compute {d_mu(O Psi/Psi), d_mu(log Psi), (O Psi/Psi)} (don't need (O Psi/Psi) here, but we get it for free)
        std::tie(fval_dmu_O, fval_dmu, fval_O) =
            psi_wrapper_in.computeMDDerivatives_total(mdd_list, fvals_dmu_O, fvals_O, fvals_dmu);

        dedr_complex[iat][idim] += fval_dmu_O; // multidet part of d_mu(OPsi/Psi)
        wfgradraw_[iat][idim] += fval_dmu;     // multidet part of d_mu(log(Psi))
      }
    }
    convertToReal(dedr_complex[iat], dEdR[iat]);
    convertToReal(wfgradraw_[iat], wf_grad[iat]);
  }
  dEdR += hfdiag_;
}
} // namespace qmcplusplus
