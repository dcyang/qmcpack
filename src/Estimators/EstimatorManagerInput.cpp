//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2025 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File refactored from: EstimatorManagerNew.cpp
//////////////////////////////////////////////////////////////////////////////////////

#include "EstimatorManagerInput.h"
#include "EstimatorInputDelegates.h"
#include <algorithm>
#include <variant>
#include "MagnetizationDensityInput.h"
#include "ModernStringUtils.hpp"

namespace qmcplusplus
{
EstimatorManagerInput::EstimatorManagerInput(xmlNodePtr cur) { readXML(cur); }

EstimatorManagerInput::EstimatorManagerInput(std::initializer_list<EstimatorManagerInput> emil)
{
  // \todo the following code fusing two vectors can be written in more consice way with std::back_inserter. See history.
  // Right now needs to use a clumsy way to make intel classic compiler 19.1.x happy when using gcc 9 on some OS distros.
  size_t est_offset        = 0;
  size_t scalar_est_offset = 0;
  for (const EstimatorManagerInput& emi : emil)
  {
    est_offset += emi.estimator_inputs_.size();
    scalar_est_offset += emi.scalar_estimator_inputs_.size();
  }
  estimator_inputs_.resize(est_offset);
  scalar_estimator_inputs_.resize(scalar_est_offset);

  est_offset        = 0;
  scalar_est_offset = 0;
  for (const EstimatorManagerInput& emi : emil)
  {
    std::copy(emi.estimator_inputs_.begin(), emi.estimator_inputs_.end(), estimator_inputs_.begin() + est_offset);
    est_offset += emi.estimator_inputs_.size();
    std::copy(emi.scalar_estimator_inputs_.begin(), emi.scalar_estimator_inputs_.end(),
              scalar_estimator_inputs_.begin() + scalar_est_offset);
    scalar_est_offset += emi.scalar_estimator_inputs_.size();
  }
}

void EstimatorManagerInput::readXML(xmlNodePtr cur)
{
  const std::string error_tag{"EstimatorManager input:"};
  std::string cur_name{lowerCase(castXMLCharToChar(cur->name))};
  xmlNodePtr child;
  if (cur_name == "estimators")
    child = cur->xmlChildrenNode;
  else
    child = cur; // the case when 'estimator's are not encapsulated by a 'estimators' node
  while (child != nullptr)
  {
    std::string cname{lowerCase(castXMLCharToChar(child->name))};
    if (cname == "estimator")
    {
      std::string atype(lowerCase(getXMLAttributeValue(child, "type")));
      std::string aname(lowerCase(getXMLAttributeValue(child, "name")));
      if (atype.empty() && !aname.empty())
        atype = aname;
      if (aname.empty() && !atype.empty())
        aname = atype;
      if (atype == "localenergy" || atype == "elocal")
        appendScalarEstimatorInput<LocalEnergyInput>(child);
      else if (atype == "cslocalenergy")
      {
        appendScalarEstimatorInput<CSLocalEnergyInput>(child);
        app_warning() << "CSLocalEnergyEstimator support is at best experimental with batch drivers" << std::endl;
      }
      else if (atype == "rmc")
      {
        appendScalarEstimatorInput<RMCLocalEnergyInput>(child);
        app_warning() << "RMCLocalEnergyEstimator support is at best experimental with batch drivers" << std::endl;
      }
      else if (atype == lowerCase(OneBodyDensityMatricesInput::type_tag))
        appendEstimatorInput<OneBodyDensityMatricesInput>(child);
      else if (atype == lowerCase(SpinDensityInput::type_tag))
        appendEstimatorInput<SpinDensityInput>(child);
      else if (atype == lowerCase(MomentumDistributionInput::type_tag))
        appendEstimatorInput<MomentumDistributionInput>(child);
      else if (atype == lowerCase(SelfHealingOverlapInput::type_tag))
        appendEstimatorInput<SelfHealingOverlapInput>(child);
      else if (atype == lowerCase(PerParticleHamiltonianLoggerInput::type_tag))
        appendEstimatorInput<PerParticleHamiltonianLoggerInput>(child);
      else if (atype == lowerCase(MagnetizationDensityInput::type_tag))
        appendEstimatorInput<MagnetizationDensityInput>(child);
      else if (atype == lowerCase(EnergyDensityInput::type_tag))
        appendEstimatorInput<EnergyDensityInput>(child);
      else if (atype == "structurefactor")
        appendEstimatorInput<StructureFactorInput>(child);
      else
        throw UniformCommunicateError(error_tag + "unparsable <estimator> node, name: " + aname + " type: " + atype +
                                      " in Estimators input.");
    }
    else if (cname != "text")
    {
      std::string atype(lowerCase(getXMLAttributeValue(child, "type")));
      std::string aname(lowerCase(getXMLAttributeValue(child, "name")));
      throw UniformCommunicateError(error_tag + "<estimators> can only contain <estimator> nodes");
    }

    if (cur_name == "estimators")
      child = child->next;
    else
    {
      app_summary() << "<estimator> nodes not contained in <estimators>...</estimators> is a deprecated input xml idiom"
                    << std::endl;
      break;
    }
  }
}

template<typename T>
std::vector<int> EstimatorManagerInput::getEstimatorTypeIndexes() const
{
  std::vector<int> type_indexes;
  for (int i = 0; i < estimator_inputs_.size(); ++i)
  {
    if (std::holds_alternative<T>(estimator_inputs_[i]))
      type_indexes.push_back(i);
  }
  return type_indexes;
}


void EstimatorManagerInput::append(const EstimatorInput& ei) { estimator_inputs_.emplace_back(ei); }
void EstimatorManagerInput::append(const ScalarEstimatorInput& sei) { scalar_estimator_inputs_.emplace_back(sei); }

template std::vector<int> EstimatorManagerInput::getEstimatorTypeIndexes<EnergyDensityInput>() const;
template std::vector<int> EstimatorManagerInput::getEstimatorTypeIndexes<SpinDensityInput>() const;
template std::vector<int> EstimatorManagerInput::getEstimatorTypeIndexes<MomentumDistributionInput>() const;
template std::vector<int> EstimatorManagerInput::getEstimatorTypeIndexes<OneBodyDensityMatricesInput>() const;
template std::vector<int> EstimatorManagerInput::getEstimatorTypeIndexes<SelfHealingOverlapInput>() const;
template std::vector<int> EstimatorManagerInput::getEstimatorTypeIndexes<MagnetizationDensityInput>() const;
template std::vector<int> EstimatorManagerInput::getEstimatorTypeIndexes<PerParticleHamiltonianLoggerInput>() const;
template std::vector<int> EstimatorManagerInput::getEstimatorTypeIndexes<StructureFactorInput>() const;


} // namespace qmcplusplus
