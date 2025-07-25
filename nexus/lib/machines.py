##################################################################
##  (c) Copyright 2015-  by Jaron T. Krogel                     ##
##################################################################


#====================================================================#
#  machines.py                                                       #
#    Representations of local machine environments including         #
#    workstations and supercomputers and the jobs that will be       #
#    executed on them.                                               #
#                                                                    #
#  Content summary:                                                  #
#    Job                                                             #
#      Class to represent a generic simulation job.                  #
#                                                                    #
#    Machine                                                         #
#      Represents a generic machine.                                 #
#      Base class for specific machine types.                        #
#                                                                    #
#    Workstation                                                     #
#      Represents a workstation with a fixed number of cores.        #
#                                                                    #
#    InteractiveCluster                                              #
#      Represents a supercomputer in interactive mode.               #
#      Similar to a workstation with many cores.                     #
#                                                                    #
#    Supercomputer                                                   #
#      Represents a generic supercomputer with a batch queue.        #
#      Base class for specific supercomputers.                       #
#      See Jaguar, Kraken, Golub, OIC5, Hopper, Edison, BlueWatersXE,#
#        BlueWatersXK, Titan, EOS, Vesta, Cetus, Mira, Lonestar,     #
#        Matisse, Komodo, and Amos                                   #
#                                                                    #
#    cpu_count                                                       #
#      Function to return the number of cores on the local machine.  #
#                                                                    #
#    Options                                                         #
#      Class representing command line options for a simulation job, #
#      including arguments to the simulation executable,             #
#      the run launcher (aprun/mpirun, etc), and the job submitter   #
#      (e.g. qsub).                                                  #
#                                                                    #
#====================================================================#


import os
import time
#from multiprocessing import cpu_count
from socket import gethostname
from subprocess import Popen,PIPE
from numpy import array,mod,floor,ceil,round,log,empty
from generic import obj
from developer import DevBase,to_str
from nexus_base import NexusCore,nexus_core
from execute import execute
from debug import *

import re,subprocess

import importlib.util
import importlib.machinery

def our_load_source(modname, filename):
    """" Replacement for the deprecated imp.load_source function"""
    loader = importlib.machinery.SourceFileLoader(modname, filename)
    spec = importlib.util.spec_from_file_location(modname, filename, loader=loader)
    module = importlib.util.module_from_spec(spec)
    # The module is always executed and not cached in sys.modules.
    # Uncomment the following line to cache the module.
    # sys.modules[module.__name__] = module
    loader.exec_module(module)
    return module

def cpu_count():
    """ Number of virtual or physical CPUs on this system, i.e.
    user/real as output by time(1) when called with an optimally scaling
    userspace-only program"""

    # Python 2.6+
    try:
        import multiprocessing
        return multiprocessing.cpu_count()
    except (ImportError,NotImplementedError):
        None
    #end try

    # POSIX
    try:
        res = int(os.sysconf('SC_NPROCESSORS_ONLN'))

        if res > 0:
            return res
        #end if
    except (AttributeError,ValueError):
        None
    #end try
#end def cpu_count



class Options(DevBase):
    def __init__(self,**kwargs):
        self.add(**kwargs)
    #end def __init__


    def add(self,**kwargs):
        self.transfer_from(kwargs)
    #end def add


    def read(self,options):
        if isinstance(options,(dict,obj)):
            self.add(**options)
        elif isinstance(options,list):
            for o in options:
                if not isinstance(o,str):
                    self.error('each option must be a string')
                #end if
                self[str(len(self))] = o
            #end for
        elif isinstance(options,str):
            self[str(len(self))] = options
        else:
            self.error('invalid type provided to Options')
        #end if
    #end def read


    def write(self):
        s = ''
        for k in sorted(self.keys()):
            s += ' '+str(self[k])
        #end for
        return s
    #end def write
#end class Options



job_defaults_assign = obj(
    name               = 'jobname',
    type               = None,
    directory          = None,
    subdir             = None,
    app_name           = None, # name of/path to application
    app_command        = None, # command used to launch application
    app_props          = None,
    full_command       = None, # custom command including e.g. mpirun
    outfile            = None,
    errfile            = None,
    env                = None,
    user_env           = True, # import user environment
    presub             = '',   # shell text executed just prior to submission
    postsub            = '',   # shell text executed just after submission
    queue              = None,
    bundled_jobs       = None,
    relative           = False,
    cores              = None, # number of cores for the job
    nodes              = None, # number of nodes for the job
    threads            = 1,    # number of openmp threads for the job
    hyperthreads       = None,
    ppn                = None,
    gpus               = None, # number of gpus per node
    serial             = False, # run job serially, no mpi
    local              = False, # run job locally, no queue submission
    days               = 0,
    hours              = 0,
    minutes            = 0,
    seconds            = 0,
    subfile            = None,
    grains             = None,
    procs              = None,
    processes          = None,
    processes_per_proc = None,
    processes_per_node = None,
    account            = None,
    email              = None,
    constraint         = None, # slurm specific, Cori
    core_spec          = None, # slurm specific, Cori
    switches           = None, # slurm specific, SuperMUC-NG
    alloc_flags        = None, # lsf specific, Summit
    filesystems        = None, # pbs specific
    qos                = None,
    group_list         = None,
    default_cpus_per_task = False, # optionally bypass typical nexus processing for supermucng
    ntasks_per_core    = None,
    cpus_per_task      = None,
    template           = None,
    )

    # these are not assigned directly
job_defaults_nonassign = obj(
    fake               = False,
    app                = None, # name of/path to application
    machine            = None,
    options            = None,
    app_flags          = None,
    app_options        = None,
    run_options        = None,
    sub_options        = None,
    skip_machine       = False,
    )

job_defaults = obj(job_defaults_assign,job_defaults_nonassign)

class Job(NexusCore):

    machine = None #default machine if none is specified in settings

    states = obj(
        none      = 0,
        waiting   = 1,
        submitted = 2,
        running   = 3,
        finished  = 4
        )
    state_names = states.inverse()

    job_count = 0


    @staticmethod
    def restore_default_settings():
        Job.machine = None
    #end def restore_default_settings


    @staticmethod
    def generate_jobid():
        Job.job_count += 1
        return Job.job_count
    #end def generate_jobid


    @classmethod
    def zero_time(cls):
        time = obj(days=0,hours=0,minutes=0,seconds=0)
        return time
    #end def zero_time


    def __init__(self,**kwargs):
        # rewrap keyword arguments
        kw = obj(**kwargs)

        # save information used to initialize job object
        self.init_info = kw.copy()

        # set defaults
        kw.set_optional(**job_defaults)

        # extract keywords not assigned
        app          = kw.delete('app')
        machine      = kw.delete('machine')
        options      = kw.delete('options')
        app_flags    = kw.delete('app_flags')
        app_options  = kw.delete('app_options')
        run_options  = kw.delete('run_options')
        sub_options  = kw.delete('sub_options')
        env          = kw.delete('env')
        fake         = kw.delete('fake')
        skip_machine = kw.delete('skip_machine')

        # assign keywords
        self.set(**kw)

        # assign fake job
        self.fake_job           = fake

        # check template
        if self.template is not None and not isinstance(self.template,str):
            self.error('template must be a string\nReceived type: {}'.format(self.template.__class__.__name__))
        #end if

        # initialize other internal variables
        self.app_options        = Options()
        self.run_options        = Options()
        self.sub_options        = Options()
        self.env                = None
        self.internal_id        = None
        self.system_id          = None
        self.tot_cores          = None
        self.identifier         = None
        self.submitted          = False
        self.status             = self.states.none
        self.crashed            = False
        self.overtime           = False
        self.successful         = False
        self.finished           = False

        self.user_app_command = self.app_command is not None

        if app is not None:
            self.app_name = app
        #end if
        if app_options is not None:
            self.app_options.read(app_options)
        #end if
        if run_options is not None:
            self.run_options.read(run_options)
        #end if
        if sub_options is not None:
            self.sub_options.read(sub_options)
        #end if
        if app_flags is not None:
            self.app_options.read(app_flags)
        #end if
        if options is not None:
            self.run_options.read(options)
        #end if

        if self.app_props is None:
            self.app_props = []
        #end if

        if self.serial:
            self.cores = 1
            self.nodes = None
        #end if

        if skip_machine:
            self.machine = machine
        else:
            if machine is not None:
                self.machine = machine
            #end if
            if env is not None:
                self.set_environment(**env)
            #end if
            #check that the machine exists and have it complete the job info
            self.process()

            machine = self.get_machine()
            self.batch_mode = machine.in_batch_mode()

            if self.bundled_jobs is not None and not machine.batch_capable:
                self.error('running batched/bundled jobs on {0} is either not possible or not yet implemented, sorry.'.format(machine.name))
            #end if
        #end if

        self.normalize_time()

    #end def __init__


    def get_machine(self):
        return Machine.get(self.machine)
    #end def get_machine


    def process(self,machine=None):
        if self.template is not None:
            return
        #end if
        if machine is None:
            machine = self.get_machine()
        #end if
        machine.process_job(self)
    #end def process


    # test needed
    def process_options(self,machine=None):
        if self.template is not None:
            return
        #end if
        if machine is None:
            machine = self.get_machine()
        #end if
        machine.process_job_options(self)
    #end def process_options


    # test needed
    def initialize(self,sim):
        self.set_id()
        self.identifier = sim.identifier
        machine = self.get_machine()
        if machine.prefixed_output:
            sim.outfile = sim.identifier + machine.outfile_extension
            sim.errfile = sim.identifier + machine.errfile_extension
        #end if
        if self.directory is None:
            self.directory = sim.locdir
            self.abs_dir   = os.path.abspath(sim.locdir)
        elif self.abs_dir is None:
            self.abs_dir = os.path.abspath(self.directory)
        #end if
        if self.subdir is None:
            if machine.local_directory!=None:
                self.subdir = os.path.join(machine.local_directory,nexus_core.runs,sim.path)
                self.abs_subdir = self.subdir
            else:
                self.subdir = self.directory
                self.abs_subdir = self.abs_dir
            #end if
        #end if
        if self.app_name is None:
            app_name = sim.app_name
        else:
            app_name = self.app_name
        #end if
        if app_name!=None and not '/' in app_name:
            ads = machine.app_directories
            ad  = machine.app_directory
            new_app_name = None
            if ads!=None and self.app_name in ads:
                new_app_name = os.path.join(ads[self.app_name],app_name)
            elif ad!=None:
                new_app_name = os.path.join(ad,app_name)
            #end if
            if new_app_name!=None:
                self.app_name = new_app_name
            #end if
        #end if
        sim.set_app_name(app_name)
        self.set(
            name    = sim.identifier,
            simid   = sim.simid,
            outfile = sim.outfile,
            errfile = sim.errfile
            )
        if self.app_command is None:
            self.app_command = sim.app_command()
        #end if
        if self.app_props==None:
            self.app_props   = list(sim.app_props)
        #end if
        # ensure job is processed properly by this initialization stage
        self.process()
    #end def initialize


    # test needed
    def renew_app_command(self,sim):
        if not self.user_app_command:
            self.app_command = sim.app_command()
        #end if
    #end def renew_app_command


    def set_id(self):
        self.internal_id = Job.generate_jobid()
    #end def set_id


    # remove?
    def set_processes(self):
        if self.processes is None:
            self.error('processes should have been set before now\ncontact the developers and have them fix this','Developer')
            self.processes = int(ceil(float(self.cores)/self.threads))
        #end if
    #end def set_processes


    def set_environment(self,limited_env=False,clear_env=False,**env):
        machine = self.get_machine()
        if isinstance(machine,Supercomputer):
            limited_env = True
        #end if
        if self.env is None:
            self.env = os.environ.copy()
            if limited_env:
                self.env.clear()
            #end if
        #end if
        if clear_env:
            self.env.clear()
        #end if
        for n,v in env.items():
            self.env[n]=str(v)
        #end for
    #end def set_environment


    def divert_out_err(self):
        self.identifier += '_divert'
    #end def divert_out_err


    def get_time(self):
        time = obj(
            days = self.days,
            hours = self.hours,
            minutes = self.minutes,
            seconds = self.seconds
            )
        return time
    #end def get_time


    def max_time(self,time):
        t  = time.seconds + 60*(time.minutes+60*(time.hours+24*time.days))
        ts = self.seconds + 60*(self.minutes+60*(self.hours+24*self.days))
        if ts>t:
            time.days    = self.days
            time.hours   = self.hours
            time.minutes = self.minutes
            time.seconds = self.seconds
        #end if
        return time
    #end def max_time


    def serial_only(self):
        return 'serial' in self.app_props and len(self.app_props)==1
    #end if


    # remove?
    def determine_end_status(self,status):
        if not nexus_core.generate_only:
            self.successful = False # not really implemented yet
        #end if
    #end def determine_end_status


    # test needed
    def write(self,file=False):
        machine = self.get_machine()
        return machine.write_job(self,file=file)
    #end def write


    # test needed
    def submit(self):
        machine = self.get_machine()
        machine.add_job(self)
        self.submitted = True
    #end def submit


    # test needed
    def reenter_queue(self):
        machine = self.get_machine()
        machine.requeue_job(self)
    #end def reenter_queue


    def run_command(self,launcher=None,redirect=False,serial=False):
        if self.template is not None:
            return ''
        #end if
        machine = self.get_machine()
        if launcher is None:
            launcher = machine.app_launcher
        #end if
        c = ''
        if self.bundled_jobs is None:
            if self.full_command is not None:
                c = self.full_command
            else:
                if self.app_command is None:
                    self.error('app_command has not been provided')
                #end if
                if launcher=='runjob':
                    separator = ' : '
                else:
                    separator = ' '
                #end if
                if self.serial and self.processes==1:
                    c = ''
                else:
                    c = launcher + self.run_options.write() + separator
                #end if
                c+=self.app_command+self.app_options.write()
                if redirect:
                    c+=' >'+self.outfile+' 2>'+self.errfile
                    if not serial:
                        c+='&'
                    #end if
                elif machine.redirect_output and self.outfile is not None:
                    c+=' >'+self.outfile+' 2>'+self.errfile
                #end if
            #end if
        elif machine.special_bundling:
            c+='\n'
            c+=machine.specialized_bundle_commands(self,launcher,serial)
        elif self.relative:
            cdir = self.abs_subdir
            c+='\n'
            for job in self.bundled_jobs:
                c+='\ncd '+os.path.relpath(job.abs_subdir,cdir)+'\n'
                c+=job.run_command(launcher,redirect=True,serial=serial)+'\n'
                cdir = job.abs_subdir
            #end for
            c+='\nwait\n'
        else:
            c+='\n'
            for job in self.bundled_jobs:
                c+='\ncd '+job.abs_subdir+'\n'
                c+=job.run_command(launcher,redirect=True,serial=serial)+'\n'
            #end for
            c+='\nwait\n'
        #end if
        return c
    #end def run_command


    def pbs_walltime(self):
        walltime=\
            str(int(self.hours   )).zfill(2)+':'\
            +str(int(self.minutes)).zfill(2)+':'\
            +str(int(self.seconds)).zfill(2)
        if self.days!=0:
            walltime = str(self.days)+':'+walltime
        #end if
        return walltime
    #end def pbs_walltime


    def sbatch_walltime(self):
        walltime=\
            str(int(24*self.days+self.hours)).zfill(2)+':'\
            +str(int(self.minutes)).zfill(2)+':'\
            +str(int(self.seconds)).zfill(2)
        return walltime
    #end def sbatch_walltime


    def ll_walltime(self):
        walltime=\
            str(int(24*self.days+self.hours)).zfill(2)+':'\
            +str(int(self.minutes)).zfill(2)+':'\
            +str(int(self.seconds)).zfill(2)
        return walltime
    #end def ll_walltime


    def lsf_walltime(self):
        walltime=\
            str(int(24*self.days+self.hours)).zfill(2)+':'\
            +str(int(self.minutes)).zfill(2)
        return walltime
    #end def lsf_walltime


    def normalize_time(self):
        t = self.total_seconds()
        d = int(t/(24*3600))
        t -= d*24*3600
        h = int(t/3600)
        t -= h*3600
        m = int(t/60)
        t -= m*60
        s = int(t)
        self.days    = d
        self.hours   = h
        self.minutes = m
        self.seconds = s
    #end def normalize_time


    def total_seconds(self):
        return self.seconds+60*(self.minutes+60*(self.hours+24*self.days))
    #end def total_seconds


    def total_minutes(self):
        return int(self.total_seconds()/60)
    #end def total_minutes


    def total_hours(self):
        return int(self.total_seconds()/3600)
    #end def total_hours


    def total_days(self):
        return int(self.total_seconds()/(24*3600))
    #end def total_days


    def clone(self):
        job = self.copy()
        job.set_id()
        return job
    #end def clone


    def serial_clone(self):
        kw = self.init_info.copy()
        kw.serial=True
        return Job(**kw)
    #end def serial_clone


    def split_nodes(self,n):
        run_options = self.run_options
        if not isinstance(n,int):
            self.error('cannot split job by nodes\nrequested split value must be an integer\nreceived type: {0}\nwith value: {1}'.format(n.__class__.__name__,n))
        elif n<1 or n>=self.nodes:
            self.error('cannot split job by nodes\nrequested split must be in the range [1,{0})\nrequested split: {1}'.format(self.nodes,n))
        #end if
        m = self.get_machine()
        if m.app_launcher=='srun':
            self.error('splitting jobs by nodes is not currently supported on machine "{0}" (SLURM)'.format(m.name))
        #end if
        job1 = self.clone()
        job2 = self.clone()
        job1.nodes = n
        job2.nodes = self.nodes - n
        m.process_job(job1)
        m.process_job(job2)
        return job1,job2
    #end def split_nodes
#end class Job




class Machine(NexusCore):

    machines = obj()

    modes = obj(
        none        = 0,
        interactive = 1,
        batch       = 2
        )
    mode = modes.none

    batch_capable       = False
    requires_account    = False
    executable_subfile  = False
    redirect_output     = False
    query_with_username = False
    special_bundling    = False

    prefixed_output    = False
    outfile_extension  = None
    errfile_extension  = None

    allow_warnings = True
    queue_configs = None

    @staticmethod
    def get_hostname():
        hostname = gethostname()
        if '.' in hostname:
            machine_name = hostname.split('.')[0]
        else:
            machine_name = hostname
        #end if
        return machine_name.lower()
    #end def get_hostname


    @staticmethod
    def exists(machine_name):
        return machine_name in Machine.machines
    #end def exists


    @staticmethod
    def is_unique(machine):
        return id(machine)==id(Machine.machines[machine.name])
    #end def is_unique


    @staticmethod
    def add(machine):
        if not isinstance(machine,Machine):
            Machine.class_error('attempted to add non-machine instance')
        #end if
        if not 'name' in machine:
            Machine.class_error('attempted to add a machine without a name')
        #end if
        name = machine.name
        if not name in Machine.machines:
            Machine.machines[name] = machine
        else:
            Machine.class_error('attempted to create machine {0}, but it already exists'.format(name))
        #end if
    #end def add


    @staticmethod
    def get(machine_name):
        if isinstance(machine_name,str):
            machine_name = machine_name.lower()
        else:
            Machine.class_error('machine name must be a string, you provided a '+machine_name.__class__.__name__)
        #end if
        if Machine.exists(machine_name):
            machine = Machine.machines[machine_name]
        else:
            machs = sorted(Machine.machines.keys())
            Machine.class_error('attempted to get machine '+machine_name+', but it is unknown\nknown options are '+str(machs))
        #end if
        return machine
    #end def get


    def warn(self,*args,**kwargs):
        if Machine.allow_warnings:
            NexusCore.warn(self,*args,**kwargs)
        #end if
    #end def warn


    def validate(self):
        if Machine.exists(self.name):
            if not Machine.is_unique(self):
                self.error('duplicate instance of machine '+self.name+' encountered\n  this is either a developer error, or you have created a duplicate machine')
            #end if
        else:
            self.error('machine {0} id {1} was created without calling Machine.__init__() and is therefore invalid'.format(self.name,id(self)))
        #end if
    #end def validate


    def in_batch_mode(self):
        return self.mode==self.modes.batch
    #end def in_batch_mode


    def query_queue(self):
        self.not_implemented()
    #end def query_queue

    def submit_jobs(self):
        self.not_implemented()
    #end def submit_jobs

    # update all job information, must be idempotent
    def process_job(self,job):
        self.not_implemented()
    #end def process_job

    def process_job_options(self,job):
        self.not_implemented()
    #end def process_job_options

    def write_job(self,job,file=False):
        self.not_implemented()
    #end def write_job

    def submit_job(self,job):
        self.not_implemented()
    #end def submit_job

    def specialized_bundle_commands(self,job,launcher,serial):
        self.not_implemented()
    #end def specialized_bundle_commands

    def __init__(self,name,queue_size=0):
        self.name = name
        self.queue_size = queue_size
        self.processes = obj()
        self.jobs = obj()
        self.waiting = set()
        self.running = set()
        self.finished= set()

        #user defined variables
        self.account         = None
        self.user            = None
        self.local_directory = None
        self.app_directory   = None
        self.app_directories = None

        if not isinstance(name,str):
            self.error('machine name must be a string\nyou provided '+str(name))
        #end if

        Machine.add(self)
    #end def __init__


    def restore_default_settings(self):
        self.account         = None
        self.user            = None
        self.local_directory = None
        self.app_directory   = None
        self.app_directories = None
    #end def restore_default_settings


    def add_job(self,job):
        if isinstance(job,Job):
            if job.template is None:
                self.process_job(job)
            #end if
            self.write_job(job)
            jid = job.internal_id
            self.jobs[jid] = job
            job.status = job.states.waiting
            self.waiting.add(jid)
            #self.write_job_states('add_job')
        else:
            self.error('add_job received non-Job instance '+job.__class__.__name__)
        #end if
    #end def add_job


    def requeue_job(self,job):
        None
    #end def requeue_job


    allowed_user_info = set(['account','local_directory','app_directory','app_directories'])
    def incorporate_user_info(self,infoin):
        info = obj(**infoin)
        vars = set(info.keys())
        invalid = vars-self.allowed_user_info
        if len(invalid)>0:
            self.error('invalid inputs encountered in incorporate_user_info\nallowed inputs: {0}\n  invalid inputs: {1}'.format(list(self.allowed_user_info),list(invalid)))
        #end if
        if 'app_directories' in info:
            ad = info.app_directories
            if not isinstance(ad,dict) and not isinstance(ad,obj):
                self.error('app_directories must be of type dict or obj\nyou provided '+ad.__class__.__name__)
            #end if
        #end if
        self.transfer_from(info)
    #end def incorporate_user_info
#end class Machine




class Workstation(Machine):

    mode = Machine.modes.interactive

    batch_capable = False

    def __init__(self,
                 name                = 'workstation',
                 cores               = None,
                 app_launcher        = 'mpirun',
                 process_granularity = 1
                 ):
        Machine.__init__(self,name)
        self.app_launcher = app_launcher
        if cores==None:
            self.cores = cpu_count()
        else:
            self.cores = cores
        #end if
        self.queue_size = cores
        self.process_granularity = process_granularity
    #end def __init__


    def process_job(self,job):
        if job.serial_only():
            job.cores=1
        elif job.cores==None:
            if job.processes!=None:
                job.cores = job.processes*job.threads
            else:
                job.cores = self.cores
            #end if
        #end if
        job.processes = max(1,int(floor(float(job.cores)/job.threads)))
        grains = int(ceil(float(job.cores)/self.process_granularity))
        if abs(grains-1-float(job.cores)/self.process_granularity)<1e-6:
            grains-=1
        #end if
        job.grains = grains
        job.cores = grains*self.process_granularity

        self.process_job_options(job)
    #end def process_job


    def process_job_options(self,job):
        job.run_options.add(np='-np '+str(job.processes))
    #end def process_job_options


    def write_job_states(self,title=''):
        self.log(title,n=2)
        n=3
        self.log('{0} {1} {2} job states'.format(self.__class__.__name__,self.name,id(self)),n=n )
        self.log('processes',n=n+1)
        for process in self.processes:
            job = process.job
            self.log('{0:>4} {1:>10} {2:>4} {3}'.format(job.internal_id,job.name,job.simid,job.directory),n=n+2)
        #end for
        self.log('jobs',n=n+1)
        jobids = list(self.jobs.keys())
        jobids.sort()
        for jobid in jobids:
            job = self.jobs[jobid]
            self.log('{0:>4} {1:>10} {2:>4} {3}'.format(job.internal_id,job.name,job.simid,job.directory),n=n+2)
        #end for
        self.log('waiting',n=n+1)
        jobids = list(self.waiting)
        jobids.sort()
        for jobid in jobids:
            job = self.jobs[jobid]
            self.log('{0:>4} {1:>10} {2:>4} {3}'.format(job.internal_id,job.name,job.simid,job.directory),n=n+2)
        #end for
        self.log('running',n=n+1)
        jobids = list(self.running)
        jobids.sort()
        for jobid in jobids:
            job = self.jobs[jobid]
            self.log('{0:>4} {1:>10} {2:>4} {3}'.format(job.internal_id,job.name,job.simid,job.directory),n=n+2)
        #end for
        self.log('finished',n=n+1)
        jobids = list(self.finished)
        jobids.sort()
        for jobid in jobids:
            job = self.jobs[jobid]
            self.log('{0:>4} {1:>10} {2:>4} {3}'.format(job.internal_id,job.name,job.simid,job.directory),n=n+2)
        #end for
        self.log('end job states',n=1)
    #end def write_job_states


    def query_queue(self):
        #self.write_job_states('query queue')
        self.validate()
        done = []
        for pid,process in self.processes.items():
            if nexus_core.generate_only or not nexus_core.monitor:
                qpid,status = pid,0
            else:
                qpid,status = os.waitpid(pid,os.WNOHANG)
            #end if
            if pid==qpid:
                job = process.job
                job.status = job.states.finished
                job.finished = True
                job.determine_end_status(status)
                iid = job.internal_id
                self.running.remove(iid)
                self.finished.add(iid)
                done.append(pid)
                if not nexus_core.generate_only:
                    job.out.close()
                    job.err.close()
                #end if
            #end if
        #end for
        for pid in done:
            del self.processes[pid]
        #end for
    #end def query_queue


    def submit_jobs(self):
        cores_used = 0
        for process in self.processes:
            cores_used += process.job.cores
        #end for
        cores_available = self.cores-cores_used

        core_req = []
        job_req  = []
        for iid in self.waiting:
            job = self.jobs[iid]
            job_req.append(job)
            core_req.append(job.cores)
        #end for
        core_req = array(core_req,dtype=int)

        # The following line does not work correctly under Numpy 1.10 or greater.
        # It should create an ndarray of Job objects from a list of Job objects.
        # Instead it creates nested ndarray's with inner type of bool
        #job_req  = array(job_req ,dtype=object)

        job_req_tmp = job_req
        job_req  = empty(len(job_req_tmp) ,dtype=object)
        for idx,job in enumerate(job_req_tmp):
            job_req[idx] = job
        #end for

        order    = core_req.argsort()
        job_req  = job_req[order]

        for job in job_req:
            if job.cores>self.cores and not nexus_core.generate_only:
                self.error('job '+str(job.internal_id)+' is too large to run on this machine\ncores requested: '+str(job.cores)+'\nmachine cores: '+str(self.cores))
            #end if
            if job.cores<=cores_available:
                iid = job.internal_id
                self.waiting.remove(iid)
                self.running.add(iid)
                self.submit_job(job)
                cores_available-=job.cores
            elif job.cores>self.cores:
                self.error('job requested more cores than are present on '+self.name+'\ncores requested: {0}\ncores present: {1}'.format(job.cores,self.cores))
            else:
                break
            #end if
        #end for
    #end def submit_jobs


    def job_command(self,job,pad=None):
        command = 'export OMP_NUM_THREADS='+str(job.threads)+'\n'
        if len(job.presub)>0:
            command += job.presub+'\n'
        #end if
        if job.serial is not None:
            command += job.run_command(self.app_launcher,serial=job.serial)
        else:
            command += job.run_command(self.app_launcher)
        #end if
        if len(job.postsub)>0:
            command += job.postsub+'\n'
        #end if
        if pad!=None:
            command = ('\n'+command).replace('\n','\n  '+pad)
        #end if
        return command
    #end def job_command


    def write_job(self,job,file=False):
        c = self.job_command(job)
        return c
    #end def write_job


    def submit_job(self,job):
        pad = self.enter(job.directory,msg=job.simid)
        command = self.job_command(job,pad=pad)
        job.status = job.states.running
        process = obj()
        process.job = job
        if nexus_core.generate_only:
            self.log(pad+'Would have executed:  '+command)
            job.system_id = job.internal_id
        else:
            if nexus_core.monitor:
                self.log(pad+'Executing:  '+command)
                job.out = open(job.outfile,'w')
                job.err = open(job.errfile,'w')
                p = Popen(command,env=job.env,stdout=job.out,stderr=job.err,shell=True)
                process.popen = p
                job.system_id = p.pid
            else:
                command+=' >'+job.outfile+' 2>'+job.errfile+'&'
                self.log(pad+'Executing:  '+command)
                os.system(command)
                job.system_id = job.internal_id
            #end if
        #end if
        self.processes[job.system_id] = process
        self.leave()
    #end def submit_job
#end class Workstation




# test needed
class InteractiveCluster(Workstation):

    def __init__(self,*args,**kwargs):
        if len(args)==0 or not isinstance(args[0],Supercomputer):
            self.init_from_args(*args,**kwargs)
        else:
            super = args[0]
            cores = args[1]
            self.init_from_supercomputer(super,cores)
        #end if
        Machine.__init__(self,self.name,self.queue_size)
    #end def __init__


    def init_from_args(self,
                       name                = 'icluster',
                       nodes               = None,
                       procs_per_node      = None,
                       cores_per_proc      = None,
                       process_granularity = None,
                       ram_per_node        = None,
                       app_launcher        = None
                       ):
        self.name           = name
        self.nodes          = nodes
        self.procs_per_node = procs_per_node
        self.cores_per_proc = cores_per_proc
        self.process_granularity = process_granularity
        self.ram_per_node   = ram_per_node
        self.app_launcher   = app_launcher

        self.cores_per_node = self.cores_per_proc*self.procs_per_node
        if process_granularity is None:
            self.process_granularity = self.cores_per_node
        #end if

        self.procs = self.procs_per_node*self.nodes
        self.cores = self.cores_per_proc*self.procs
        self.ram   = self.ram_per_node*self.nodes

        self.queue_size = self.cores
    #end def init_from_args


    def init_from_supercomputer(self,super,cores):
        nodes = cores//super.cores_per_node
        if cores-nodes*super.cores_per_node!=0:
            self.error('interactive cores corresponds to a fractional number of nodes\n  cores '+str(cores)+'\n  cores per node '+str(super.cores_per_node))
        #end if
        self.init_from_args(super.name+'_interactive',nodes,super.procs_per_node,
                            super.cores_per_proc,super.cores_per_node,
                            super.ram_per_node,super.app_launcher)
    #end def init_from_supercomputer


    def process_job(self,job):
        job.cores = min(job.cores,self.cores)

        Workstation.process_job(self,job)

        job.nodes = job.grains
        job.procs = job.nodes*self.procs_per_node

        if mod(job.processes,job.nodes)!=0:
            job.processes_per_node = None
        else:
            job.processes_per_node = job.processes//job.nodes
        #end if

        if mod(job.processes,job.procs)!=0:
            job.processes_per_proc = None
        else:
            job.processes_per_proc = job.processes//job.procs
        #end if
    #end def process_job
#end class InteractiveCluster




class Supercomputer(Machine):
    mode = Machine.modes.batch
    name = 'supercomputer'

    batch_capable = False #only set to true for specific machines

    aprun_options = set(['n','d'])

    required_inputs = [
        'nodes',
        'procs_per_node',
        'cores_per_proc',
        'ram_per_node',
        'queue_size',
        'app_launcher',
        'sub_launcher',
        'queue_querier',
        'job_remover'
        ]

    def __init__(self,
                 nodes          = None,
                 procs_per_node = None,
                 cores_per_proc = None,
                 ram_per_node   = None,
                 queue_size     = 0,
                 app_launcher   = None,
                 sub_launcher   = None,
                 queue_querier  = None,
                 job_remover    = None,
                 name           = None,
                 ):
        if name is None:
            if self.name is not None:
                name = self.name
            else:
                name = self.__class__.__name__.lower()
            #end if
        #end if
        Machine.__init__(self,name)
        self.nodes          = nodes          #  # of nodes
        self.procs_per_node = procs_per_node #  # of processors/sockets on a node
        self.cores_per_proc = cores_per_proc #  # of cores on a processor/socket
        self.ram_per_node   = ram_per_node
        self.queue_size     = queue_size
        self.app_launcher   = app_launcher
        self.sub_launcher   = sub_launcher
        self.queue_querier  = queue_querier
        self.job_remover    = job_remover

        for var in Supercomputer.required_inputs:
            if self[var] is None:
                self.error('input variable '+var+' is required to initialize Supercomputer object.')
            #end if
        #end for

        self.cores_per_node = self.cores_per_proc*self.procs_per_node

        self.procs = self.procs_per_node*self.nodes
        self.cores = self.cores_per_proc*self.procs
        self.ram   = self.ram_per_node*self.nodes

        # 'complete' is the only actively used status so far
        #   At least one queue state should correspond to 'complete',
        #   though even this is not strictly necessary.
        #   In general if a pid is missing from the queue,
        #   that job is assumed to be complete.

        self.system_queue = obj()
        if self.queue_querier=='qstat':
            self.job_states=dict(E = 'exiting',
                                 H = 'held',
                                 Q = 'queued',
                                 R = 'running',
                                 S = 'suspended',
                                 T = 'transferring',
                                 W = 'waiting',
                                 C = 'complete', 
                                 F = 'complete',
                                 B = 'has_subjob',
                                 M = 'moved_to_another_server',
                                 U = 'user_suspended',
                                 X = 'subjob_completed',
                                 )
        elif self.queue_querier=='qstata':
            #already gives status as queued, running, etc.
            None
        elif  self.queue_querier=='squeue':
            self.job_states=dict(CG = 'exiting',
                                 TO = 'timeout',
                                 CA = 'failed',
                                 F = 'failed',
                                 NF = 'node_fail',
                                 PD = 'waiting',
                                 R = 'running',
                                 S = 'suspended',
                                 CD = 'complete',
                                 RD = 'held',
                                 BF = 'failed',
                                 CF = 'configuring',
                                 DL = 'deadline',
                                 OOM= 'out_of_memory',
                                 PR = 'preempted',
                                 RF = 'requeue_fed',
                                 RH = 'requeue_hold',
                                 RQ = 'requeued',
                                 RS = 'resizing',
                                 RV = 'revoked',
                                 SI = 'signaling',
                                 SE = 'special_exit',
                                 SO = 'stage_out',
                                 ST = 'stopped',
                                 )
        elif self.queue_querier=='sacct':
            self.job_states=dict(CANCELLED = 'failed',  #long form
                                 COMPLETED = 'complete',
                                 COMPLETING = 'exiting',
                                 CONFIGURING = 'waiting',
                                 FAILED = 'failed',
                                 PREMEEMPTED = 'failed',
                                 PENDING = 'waiting',
                                 NODE_FAIL = 'failed',
                                 RESIZING = 'resizing',
                                 RUNNING = 'running',
                                 SUSPENDED = 'suspended',
                                 TIMEOUT = 'failed',
                                 CA = 'failed',        #short form
                                 CD = 'complete',
                                 CG = 'exiting',
                                 CF = 'waiting',
                                 F = 'failed',
                                 PR = 'failed',
                                 PD = 'waiting',
                                 NF = 'failed',
                                 RS = 'resizing',
                                 R = 'running',
                                 S = 'suspended',
                                 TO = 'failed'
                                 )
        elif self.queue_querier=='llq':
            self.job_states=dict(I  = 'idle',
                                 NQ = 'not_queued',
                                 H  = 'user_hold',
                                 S  = 'system_hold',
                                 HS = 'user_system_hold',
                                 D  = 'deferred',
                                 R  = 'running',
                                 P  = 'pending',
                                 ST = 'starting',
                                 C  = 'complete',
                                 CA = 'canceled',
                                 E  = 'preempted',
                                 EP = 'preempt_pending',
                                 MP = 'resume_pending',
                                 )
        elif self.queue_querier=='bjobs':
            self.job_states=dict(PEND  = 'pending',
                                 RUN   = 'running',
                                 DONE  = 'complete',
                                 EXIT  = 'failed',
                                 PSUSP = 'suspended',
                                 USUSP = 'suspended',
                                 SSUSP = 'suspended',
                                 )
        elif self.queue_querier=='test_query':
            None
        else:
            self.error('ability to query queue with '+self.queue_querier+' has not yet been implemented')
        #end if

    #end def __init__


    # test needed
    def interactive_representation(self,cores):
        return InteractiveCluster(self,cores)
    #end def interactive_representation


    # test needed
    def requeue_job(self,job):
        if isinstance(job,Job):
            jid = job.internal_id
            pid = job.system_id
            if pid is None:
                self.error('job {0} does not have a process id issued by the scheduler'.format(jid))
            #end if
            self.process_job(job)
            self.jobs[jid] = job
            job.status = job.states.running
            self.running.add(jid)
            process = obj(job=job)
            self.processes[pid] = process
        else:
            self.error('requeue_job received non-Job instance '+job.__class__.__name__)
        #end if
    #end def requeue_job


    def process_job(self,job):
        if job.fake_job:
            return
        #end if

        self.pre_process_job(job)

        job.subfile = job.name+'.'+self.sub_launcher+'.in'
        no_cores = job.cores is None
        no_nodes = job.nodes is None
        if no_cores and no_nodes:
            self.error('job did not specify cores or nodes\nAt least one must be provided')
        elif no_cores:
            job.cores = self.cores_per_node*job.nodes
        elif no_nodes:
            job.nodes = int(ceil(float(job.cores)/self.cores_per_node))
            if abs(job.nodes-1-float(job.cores)/self.cores_per_node)<1e-6:
                job.nodes-=1
            #end if
        else:
            job.cores = min(job.cores,job.nodes*self.cores_per_node)
        #end if
        if job.processes_per_node is not None:
            job.processes = job.nodes*job.processes_per_node
        else:
            job.processes = max(1,int(float(job.cores)/job.threads))
        #end if
        job.tot_cores = job.nodes*self.cores_per_node
        job.procs = job.nodes*self.procs_per_node

        if mod(job.processes,job.nodes)!=0:
            job.processes_per_node = None
        else:
            job.processes_per_node = job.processes//job.nodes
        #end if

        if mod(job.processes,job.procs)!=0:
            job.processes_per_proc = None
        else:
            job.processes_per_proc = job.processes//job.procs
        #end if

        if job.ppn is None:
            job.ppn = self.cores_per_node
        #end if

        if job.account is None:
            if self.account is not None:
                job.account = self.account
            elif self.requires_account:
                self.error('account not specified for job on '+self.name)
            #end if
        #end if
        self.post_process_job(job)

        job.set_environment(OMP_NUM_THREADS=job.threads)

        self.process_job_options(job)
    #end def process_job


    def process_job_options(self,job):
        launcher = self.app_launcher
        if launcher=='mpirun':
            job.run_options.add(np='-np '+str(job.processes))
        elif launcher=='mpiexec':
            job.run_options.add(n='-n '+str(job.processes))
        elif launcher=='aprun':
            if 'n' in self.aprun_options:
                job.run_options.add(n='-n '+str(job.processes))
            #end if
            if 'd' in self.aprun_options and job.threads>1:
                job.run_options.add(d='-d '+str(job.threads))
            #end if
            if 'N' in self.aprun_options and job.processes_per_node is not None:
                job.run_options.add(N='-N '+str(job.processes_per_node))
            #end if
            if 'S' in self.aprun_options and job.processes_per_proc is not None:
                job.run_options.add(S='-S '+str(job.processes_per_proc))
            #end if
        elif launcher=='runjob':
            #bypass setup_environment
            if job.env is not None:
                envs='--envs'
                for name,value in job.env.items():
                    envs+=' {0}={1}'.format(name,value)
                #end for
                job.env = None
            elif 'envs' in job.run_options:
                envs = job.run_options.envs
            else:
                self.error('failed to set env options for runjob')
            #end if
            job.run_options.add(
                np       = '--np '+str(job.processes),
                p        = '-p '+str(job.processes_per_node),
                xlocargs = '$LOCARGS',
                verbose  = '--verbose=INFO',
                envs     = envs
                )
        elif launcher=='srun':  # Amos contribution from Ryan McAvoy
            None
        elif launcher=='ibrun': # Lonestar contribution from Paul Young
            job.run_options.add(
	        np	= '-n '+str(job.processes),
	        p	= '-o '+str(0),
	        )
        elif launcher=='jsrun': # Summit
            None # Summit class takes care of this in post_process_job
        elif launcher=='lrun': # Lassen
            None # Lassen class takes care of this in post_process_job
        else:
            self.error(launcher+' is not yet implemented as an application launcher')
        #end if
    #end def process_job_options


    def pre_process_job(self,job):
        None
    #end def pre_process_job


    def post_process_job(self,job):
        None
    #end def post_process_job


    def query_queue(self,out=None):
        self.system_queue.clear()
        if self.query_with_username and self.user is None:
            self.error('querying queue on machine "{}" requires user name\nplease provide username via the "user" keyword in settings'.format(self.name))
        #end if
        if self.queue_querier=='qstat':
            if out is None:
                out,err,rc = execute('qstat -a')
            #end if
            lines = out.splitlines()
            for line in lines:
                tokens=line.split()
                if len(tokens)>0:
                    if '.' in tokens[0]:
                        spid = tokens[0].split('.')[0]
                    else:
                        spid = tokens[0]
                    #endif
                    if spid.isdigit() and len(tokens)==11:
                        pid = int(spid)
                        jid,uname,queue,jname,sessid,nodes,tasks,mem,rtime,status,etime = tokens
                        if status in self.job_states:
                            self.system_queue[pid] = self.job_states[status]
                        else:
                            self.error('job state '+status+' is unrecognized')
                        #end if
                    #end if
                #end if
            #end for
        elif self.queue_querier=='qstata':
            if out is None:
                out,err,rc = execute('qstat')
            #end if
            lines = out.splitlines()
            for line in lines:
                tokens=line.split()
                if len(tokens)>0:
                    if '.' in tokens[0]:
                        spid = tokens[0].split('.')[0]
                    else:
                        spid = tokens[0]
                    #endif
                    if spid.isdigit() and len(tokens)==6:
                        pid = int(spid)
                        jid,uname,wtime,nodes,status,loc = tokens
                        self.system_queue[pid] = status
                    #end if
                #end if
            #end for
        elif self.queue_querier=='squeue': # contributed by Ryan McAvoy
            if out is None:
                if isinstance(self.user,bool) and self.user==False:
                    extra = ''
                elif self.user is not None:
                    extra = ' -u {}'.format(self.user)
                else:
                    extra = ' --user=$USER'
                #end if
                out,err,rc = execute('squeue'+extra)
            #end if
            lines = out.splitlines()
            for line in lines:
                tokens=line.split()
                if len(tokens)>0:
                    if '.' in tokens[0]:
                        spid = tokens[0].split('.')[0]
                    else:
                        spid = tokens[0]
                    #endif
                    if spid.isdigit():
                        pid = int(spid)
                        status = None
                        jid,loc,name,uname,status,wtime,nodes,reason = tokens[:8]
                        if status is not None:
                            if status in self.job_states:
                                self.system_queue[pid] = self.job_states[status]
                            else:
                                self.error('job state '+status+' is unrecognized')
                            #end if
                        #end if
                    #end if
                #end if
            #end for
        elif self.queue_querier=='sacct': # contributed by Ryan McAvoy
            if out is None:
                out,err,rc = execute('sacct')
            #end if
            lines = out.splitlines()
            for line in lines:
                tokens=line.split()
                if len(tokens)>0:
                    if '.' in tokens[0]:
                        spid = tokens[0].split('.')[0]
                    else:
                        spid = tokens[0]
                    #endif
                    if spid.isdigit() and len(tokens)==6:  #if account is empty, only 6 tokens.

                        pid = int(spid)
                        jid,name,loc,cores,status,exit_code = tokens
                        status = status.split('+')[0]  ## get rid of '+' in the end
                        if status in self.job_states:
                            self.system_queue[pid] = self.job_states[status]
                        else:
                            self.error('job state '+status+' is unrecognized')
                        #end if
                    elif spid.isdigit() and len(tokens)==7:

                        pid = int(spid)
                        jid,name,loc,uname,cores,status,exit_code = tokens
                        status = status.split('+')[0]  ## get rid of '+' in the end
                        if status in self.job_states:
                            self.system_queue[pid] = self.job_states[status]
                        else:
                            self.error('job state '+status+' is unrecognized')
                        #end if
                    #end if
                #end if
            #end for
        elif self.queue_querier=='llq':
            if out is None:
                out,err,rc = execute('sacct')
            #end if
            lines = out.splitlines()
            for line in lines:
                tokens=line.split()
                if len(tokens)>0:
                    if '.' in tokens[0]:
                        spid = tokens[0].split('.')[1]
                    else:
                        spid = tokens[0]
                    #endif
                    if spid.isdigit() and (len(tokens)==7 or len(tokens)==8):
                        pid = int(spid)
                        if len(tokens)==7:
                            jid,owner,subdate,subtime,status,pri,class_ = tokens
                        elif len(tokens)==8:
                            jid,owner,subdate,subtime,status,pri,class_,running_on = tokens
                        #end if
                        if status in self.job_states:
                            self.system_queue[pid] = self.job_states[status]
                        else:
                            self.error('job state '+status+' is unrecognized')
                        #end if
                    #end if
                #end if
            #end for
        elif self.queue_querier=='bjobs':
            if out is None:
                out,err,rc = execute('bjobs')
            #end if
            lines = out.splitlines()
            for line in lines:
                tokens=line.split()
                if len(tokens)>0:
                    spid = tokens[0]
                    if spid.isdigit() and len(tokens)==8:
                        pid = int(spid)
                        jid,uname,status,slots,queue,start,finish,jname = tokens
                        if status in self.job_states:
                            self.system_queue[pid] = self.job_states[status]
                        else:
                            self.error('job state '+status+' is unrecognized')
                        #end if
                    #end if
                #end if
            #end for
        elif self.queue_querier=='test_query': # for testing
            # pretend that all jobs have finished
            for pid in self.processes.keys():
                self.system_queue[pid] = 'complete'
            #end for
        else:
            self.error('ability to query queue with '+self.queue_querier+' has not yet been implemented')
        #end if
        done = []
        for pid,process in self.processes.items():
            if not pid in self.system_queue or self.system_queue[pid]=='complete' or nexus_core.generate_only:
                job = process.job
                job.status = job.states.finished
                job.finished = True
                iid = job.internal_id
                self.running.remove(iid)
                self.finished.add(iid)
                done.append(pid)
            #end if
        #end for
        for pid in done:
            del self.processes[pid]
        #end for
        return self.system_queue
    #end def query_queue


    def submit_jobs(self):
        nprocesses_running = len(self.processes)
        queue_slots_available = self.queue_size-nprocesses_running
        remove = []
        for iid in self.waiting:
            if queue_slots_available>0:
                remove.append(iid)
                self.running.add(iid)
                job = self.jobs[iid]
                self.submit_job(job)
                queue_slots_available -= 1
            else:
                break
            #end if
        #end for
        for iid in remove:
            self.waiting.remove(iid)
        #end for
    #end def submit_jobs


    def submit_job(self,job):
        pad = self.enter(job.directory,msg=job.internal_id)
        if job.subfile==None:
            self.error('submission file not specified for job')
        elif not os.path.exists(job.subfile):
            self.error('job submission file was not written prior to submission\n  submission file: '+os.path.join(job.directory,job.subfile))
        #end if
        command = self.sub_command(job)
        if nexus_core.generate_only:
            self.log(pad+'Would have executed:  '+command)
            job.status = job.states.running
            process = obj()
            process.job = job
            self.processes[job.internal_id] = process
        else:
            self.log(pad+'Executing:  '+command)
            job.status = job.states.running
            process = obj()
            process.job = job
            out,err,rc = execute(command)
            output=out+'\n'+err
            pid = self.read_process_id(output)
            if pid is None:
                self.error('process id could not be determined from submission output\n  output:\n'+output)
            else:
                self.log(pad+'  pid: {0}'.format(pid))
            #end if
            #pid = 'fakepid_'+str(job.internal_id)
            job.system_id = pid
            self.processes[pid] = process
        #end if
        self.leave()
    #end def submit_job


    def sub_command(self,job):
        return self.sub_launcher+job.sub_options.write()+' '+job.subfile
    #end def sub_command


    def remove_job(self,job):
        if self.job_remover=='qdel':
            command = 'qdel '+str(job.system_id)
        elif self.job_remover=='scancel':
            command = 'scancel '+str(job.system_id)
        else:
            self.error('ability to remove job using '+self.job_remover+' has not yet been implemented')
        #endif
        os.system(command)
    #end def remove_job


    def setup_environment(self,job):
        env = ''
        if job.env is not None:
            for name,val in job.env.items():
                env +='export {0}={1}\n'.format(name,val)
            #end for
        #end if
        return env
    #end def setup_environment


    def write_job(self,job,file=False):
        job.subfile = job.name+'.'+self.sub_launcher+'.in'
        if job.template is None:
            env = self.setup_environment(job)
            command = job.run_command(self.app_launcher,serial=job.serial)

            c = self.write_job_header(job)+'\n'
            if len(job.presub)>0:
                c+=job.presub+'\n'
            #end if
            c+=env
            c+=command+'\n'
            if len(job.postsub)>0:
                c+=job.postsub+'\n'
            #end if
        else:
            c = job.template
        #end if
        if file:
            filepath = os.path.join(job.directory,job.subfile)
            fobj = open(filepath,'w')
            fobj.write(c)
            fobj.close()
            if self.executable_subfile:
                os.system('chmod +x '+filepath)
            #end if
        #end if
        return c
    #end def write_job


    def write_job_header(self,job):
        self.not_implemented()
    #end def write_job_header

    @staticmethod
    def walltime_to_seconds(walltime_str):
        """
        Convert walltime string to total seconds
        Handles formats: 'dd:hh:mm:ss', 'hh:mm:ss', 'mm:ss', 'seconds'
        """
        try:
            parts = walltime_str.split(':')
            seconds = 0
            if len(parts) == 4:    # dd:hh:mm:ss
                seconds += int(parts[0]) * 24 * 3600  # days
                seconds += int(parts[1]) * 3600       # hours
                seconds += int(parts[2]) * 60         # minutes
                seconds += int(parts[3])              # seconds
            elif len(parts) == 3:  # hh:mm:ss
                seconds += int(parts[0]) * 3600       # hours
                seconds += int(parts[1]) * 60         # minutes
                seconds += int(parts[2])              # seconds
            elif len(parts) == 2:  # mm:ss
                seconds += int(parts[0]) * 60         # minutes
                seconds += int(parts[1])              # seconds
            elif len(parts) == 1:  # seconds only
                seconds += int(parts[0])
            else:
                raise ValueError(f"Invalid walltime format: {walltime_str}, accepted formats: 'dd:hh:mm:ss', 'hh:mm:ss', 'mm:ss', 'seconds'")
            return seconds
        except (ValueError, IndexError) as e:
            raise ValueError(f"Failed to parse walltime '{walltime_str}': {str(e)}")
        #end try
    #end def walltime_to_seconds
    @ staticmethod
    def seconds_to_walltime(seconds):
        """
        Convert total seconds to walltime string
        Handles formats: 'dd:hh:mm:ss', 'hh:mm:ss', 'mm:ss', 'seconds'
        """
        days = seconds // 86400
        hours = (seconds % 86400) // 3600
        minutes = (seconds % 3600) // 60
        seconds = seconds % 60
        return f"{days}:{hours:02d}:{minutes:02d}:{seconds:02d}"
    #end def seconds_to_walltime

    def validate_queue_config(self, job):
        """
        Validate job against queue configuration constraints
        Returns True if valid, raises an error with detailed message if invalid

        Queue config format:
        queue_configs = {
            'default': 'queue_name',
            'queue_name': {
                'min_nodes'     : minimum number of nodes,
                'max_nodes'     : maximum number of nodes,
                'max_walltime'      : maximum walltime in hours,
                'cores_per_node': cores per node,
                'ram_per_node'  : RAM per node in GB,
                'constraints'   : {
                    'knl': {
                        'max_nodes': knl specific max nodes,
                        'max_time': knl specific max time,
                        ...
                    },
                    'haswell': {...},
                    'gpu': {...},
                    'cpu': {...}
                }
            }
        }
        """
        if self.queue_configs is None:
            return True
        #end if
        
        # Use only warnings here, as smaller computers may not have a queue
        if job.queue is None:
            if 'default' in self.queue_configs:
                job.queue = self.queue_configs['default']
                self.warn('No default queue is specified. Using default queue {}'.format(job.queue))
            else:
                # No queue or default queue is specified
                self.warn('No queue or default queue is specified.')
                return True
            #end if
        #end if

        # Check if queue exists
        if job.queue not in self.queue_configs:
            # Queue is defined but config is not available
            self.warn('Queue "{}" is not available. Available queues: {}'.format(
                job.queue, list(self.queue_configs.keys())))
            return False
        else:
            # Queue is defined and config is available
            config = self.queue_configs[job.queue]
        #end if

        
        errors = []

        # Get constraint-specific config if applicable
        constraint_config = None
        if hasattr(job, 'constraint') and job.constraint is not None:
            if 'constraints' in config and job.constraint in config['constraints']:
                constraint_config = config['constraints'][job.constraint]
            else:
                valid_constraints = list(config.get('constraints', {}).keys())
                if valid_constraints:
                    errors.append('Invalid constraint "{}". Valid constraints for queue {}: {}'.format(
                        job.constraint, job.queue, valid_constraints))
                #end if
            #end if
        #end if

        # Use constraint specific config if available, otherwise use base config
        active_config = constraint_config if constraint_config is not None else config

        # Validate nodes
        if 'min_nodes' in active_config and job.nodes < active_config['min_nodes']:
            errors.append('Number of nodes ({}) is below minimum ({}) for queue {}'.format(
                job.nodes, active_config['min_nodes'], job.queue))
        #end if
        if 'max_nodes' in active_config and job.nodes > active_config['max_nodes']:
            errors.append('Number of nodes ({}) exceeds maximum ({}) for queue {}'.format(
                job.nodes, active_config['max_nodes'], job.queue))
        #end if

        # Validate cores per node
        if 'cores_per_node' in active_config and hasattr(job, 'processes_per_node'):
            if job.processes_per_node > active_config['cores_per_node']:
                errors.append('Processes per node ({}) exceeds cores per node ({}) for queue {}'.format(
                    job.processes_per_node, active_config['cores_per_node'], job.queue))
            #end if
        #end if

        # Validate memory
        if 'ram_per_node' in active_config and hasattr(job, 'ram_per_node'):
            if job.ram_per_node > active_config['ram_per_node']:
                errors.append('RAM per node ({} GB) exceeds maximum ({} GB) for queue {}'.format(
                    job.ram_per_node, active_config['ram_per_node'], job.queue))
            #end if
        #end if

        # Validate walltime
        if 'max_time' in active_config:
            job.total_hours = job.days*24 + job.hours + job.minutes/60.0 + job.seconds/3600.0
            if job.total_hours > active_config['max_time']:
                errors.append('Walltime ({} hours) exceeds maximum ({} hours) for queue {}'.format(
                    job.total_hours, active_config['max_time'], job.queue))
            #end if
        #end if

        # Validate GPU requirements if specified
        if 'gpus_per_node' in active_config and hasattr(job, 'gpus_per_node'):
            if job.gpus_per_node > active_config['gpus_per_node']:
                errors.append('GPUs per node ({}) exceeds maximum ({}) for queue {}'.format(
                    job.gpus_per_node, active_config['gpus_per_node'], job.queue))
            #end if
        #end if

        # Report all validation errors
        if errors:
            self.error('Queue validation failed for queue {}:\n  {}'.format(
                job.queue, '\n  '.join(errors)))
            return False
        #end if

        return True
    #end def validate_queue_config

    def read_process_id(self,output):
        pid = None
        lines = output.splitlines()
        if self.sub_launcher=='llsubmit': # specialization for load leveler (SuperMUC)
            for line in lines:
                if 'llsubmit: The job' in line and '"' in line:
                    spid = line.split('"')[1].split('.')[1].strip()
                    if spid.isdigit():
                        pid = int(spid)
                        break
                    #end if
                #end if
            #end for
        else:  # most other machines follow the pattern below
            for line in lines:
                ls = line.strip()
                if ls.isdigit():
                    pid = int(ls)
                    break
                elif '.' in line:
                    spid = line.split('.')[0]
                    if spid.isdigit():
                        pid = int(spid)
                        break
                    #end if
                elif ' ' in line: # specialized for Amos?
                    spid = line.split(' ')[-1]
                    if spid.isdigit():
                        pid = int(spid)
                        break
                    #end if
                #end if
            #end for
        #end if
        return pid
    #end def read_process_id

#end class Supercomputer



# Load local class for local cluster's setting from ~/.nexus/local.py
# The following is an example of this file (see machines.py for more examples):
'''
from machines import Supercomputer
class Clustername(Supercomputer):
    name = 'clustername'
    requires_account = False
    batch_capable    = True
    def write_job_header(self,job):
        if job.queue is None:
            job.queue='batch'
        #end if
        c= '#!/bin/bash\n'
        c+='#PBS -l nodes={0}:ppn={1}\n'.format(job.nodes,job.ppn)
        c+='#PBS -l walltime='+job.pbs_walltime()+'\n'
        c+='#PBS -N '+job.name +'\n'
        c+='#PBS -o '+job.outfile+'\n'
        c+='#PBS -e '+job.errfile+'\n'
        c+='#PBS -V\n'
        c+='#PBS -q '+job.queue+'\n'
        c+=' \n '
        return c
    #end def write_job_header
#end class Clustername
#            nodes sockets cores ram qslots  qlaunch  qsubmit     qstatus   qdelete
Clustername(      4,   1,    16,   24,    4, 'mpirun',     'qsub',   'qstat',    'qdel')
'''
try:
    our_load_source('*',os.path.expanduser('~/.nexus/local_machines.py'))
except IOError:
    pass
except:
    raise
#end try

#Decommissioned
class Kraken(Supercomputer):

    name = 'kraken'

    requires_account = True

    def write_job_header(self,job):
        c='#!/bin/bash\n'
        c+='#PBS -A '+str(job.account)+'\n'
        c+='#PBS -N '+str(job.name)+'\n'
        c+='#PBS -l walltime='+job.pbs_walltime()+'\n'
        c+='#PBS -l size='+str(job.tot_cores)+'\n'
        c+='#PBS -o '+job.outfile+'\n'
        c+='#PBS -e '+job.errfile+'\n'
        if job.user_env:
            c+='#PBS -V\n'
        #end if
        c+='''
cd $PBS_O_WORKDIR
export MPICH_PTL_SEND_CREDITS=-1
export MPICH_MAX_SHORT_MSG_SIZE=1024
export MPICH_PTL_UNEX_EVENTS=800000
export MPICH_UNEX_BUFFER_SIZE=16M
export MPI_MSGS_PER_PROC=32768
'''
        return c
    #end def write_job_header
#end class Kraken


#Decommissioned
class Jaguar(Supercomputer):
    name = 'jaguar'

    requires_account = True

    def write_job_header(self,job):
        if job.queue is None:
            job.queue='batch'
        #end if
        c='#!/bin/bash\n'
        c+='#PBS -A '+str(job.account)+'\n'
        c+='#PBS -q '+job.queue+'\n'
        c+='#PBS -N '+str(job.name)+'\n'
        c+='#PBS -l walltime='+job.pbs_walltime()+'\n'
        c+='#PBS -l size='+str(job.tot_cores)+'\n'
        c+='#PBS -l gres=widow2%widow3\n'
        c+='#PBS -o '+job.outfile+'\n'
        c+='#PBS -e '+job.errfile+'\n'
        if job.user_env:
            c+='#PBS -V\n'
        #end if
        c+='''
echo $PBS_O_WORKDIR
cd $PBS_O_WORKDIR
export MPICH_PTL_SEND_CREDITS=-1
export MPICH_MAX_SHORT_MSG_SIZE=1024
export MPICH_PTL_UNEX_EVENTS=800000
export MPICH_UNEX_BUFFER_SIZE=16M
export MPI_MSGS_PER_PROC=32768
'''
        return c
    #end def write_job_header
#end class Jaguar



#Unknown
class Golub(Supercomputer):
    name = 'golub'
    def write_job_header(self,job):
        if job.queue is None:
            job.queue='secondary'
        #end if
        c=''
        c+='#PBS -q '+job.queue+'\n'
        c+='#PBS -N '+job.name+'\n'
        c+='#PBS -l nodes={0}:ppn={1}\n'.format(job.nodes,job.ppn)
        c+='#PBS -l walltime='+job.pbs_walltime()+'\n'
        c+='#PBS -e '+job.errfile+'\n'
        c+='#PBS -o '+job.outfile+'\n'
        if job.user_env:
            c+='#PBS -V\n'
        #end if
        c+='''
cd ${PBS_O_WORKDIR}

'''
        return c
    #end def write_job_header
#end class Golub



# Decommissioned
class OIC5(Supercomputer):

    name = 'oic5'
    batch_capable = True

    def write_job_header(self,job):
        if job.queue is None:
            job.queue = 'mstqmc13q'
        #end if

        ppn = job.processes_per_node
        #ppn = 32/job.threads
        #if ppn*job.threads!=32:
        #    self.error('ppn is not being set properly for OIC5\n  perhaps the number of threads requested does not evenly divide the 32 cores\n  you requested {0} threads'.format(job.threads))
        ##end if

        c='#!/bin/bash\n'
        c+='#PBS -q '+job.queue+'\n'
        c+='#PBS -N '+str(job.name)+'\n'
        c+='#PBS -l walltime='+job.pbs_walltime()+'\n'
        c+='#PBS -l nodes={0}:ppn={1}\n'.format(job.nodes,ppn)
        c+='#PBS -W x=\"NACCESSPOLICY:SINGLEJOB\"\n'
        c+='#PBS -o '+job.outfile+'\n'
        c+='#PBS -e '+job.errfile+'\n'
        if job.user_env:
            c+='#PBS -V\n'
        #end if
        c+='''
echo $PBS_O_WORKDIR
cd $PBS_O_WORKDIR
'''
        return c
    #end def write_job_header


    def read_process_id(self,output):
        pid = None
        lines = output.splitlines()
        for line in lines:
            if 'oic.ornl.gov' in line:
                spid = line.split('.')[0]
                if spid.isdigit():
                    pid = int(spid)
                #end if
            #end if
        #end for
        return pid
    #end def read_process_id
#end class OIC5




class NerscMachine(Supercomputer):
    batch_capable = True

    def write_job_header(self,job):
        if job.queue is None:
            job.queue = 'regular'
        #end if
        c='#!/bin/bash\n'
        c+='#SBATCH -p '+job.queue+'\n'
        c+='#SBATCH -J '+str(job.name)+'\n'
        c+='#SBATCH -t '+job.sbatch_walltime()+'\n'
        c+='#SBATCH -N '+str(job.nodes)+'\n'
        c+='#SBATCH --ntasks-per-node={0}\n'.format(job.processes_per_node)
        c+='#SBATCH --cpus-per-task={0}\n'.format(job.threads)
        c+='#SBATCH -o '+job.outfile+'\n'
        c+='#SBATCH -e '+job.errfile+'\n'
        if job.user_env:
            c+='#SBATCH --export=ALL\n'   # equiv to PBS -V
        else:
            c+='#SBATCH --export=NONE\n'
        #end if
        c+='''
echo $SLURM_SUBMIT_DIR
cd $SLURM_SUBMIT_DIR
'''
        return c
    #end def write_job_header
#end class NerscMachine

#Decommissioned
class Cori(NerscMachine):
    name = 'cori'

    def pre_process_job(self,job):
        if job.queue is None:
            job.queue = 'regular'
        #end if
        if job.constraint is None:
            job.constraint = 'knl'
        #end if
        # account for dual nature of Cori
        if 'knl' in job.constraint:
            self.nodes          = 9688
            self.procs_per_node = 1
            self.cores_per_node = 68
            self.ram_per_node   = 96
        elif 'haswell' in job.constraint:
            self.nodes = 2388
            self.procs_per_node = 2
            self.cores_per_node = 32
            self.ram_per_node   = 128
        elif 'amd' in job.constraint:
            self.nodes = 20
            self.procs_per_node = 2
            self.cores_per_node = 32
            self.ram_per_node   = 2048
        else:
            self.error('SLURM input "constraint" must contain either "knl", "haswell", or "amd" on Cori\nyou provided: {0}'.format(job.constraint))
        #end if
        if job.core_spec is not None:
            self.cores_per_node -= job.core_spec
        #end if
    #end def pre_process_job

    def write_job_header(self,job):
        self.pre_process_job(job) # sync machine view with job
        if 'knl' in job.constraint:
            hyperthreads   = 4
        elif 'haswell' in job.constraint:
            hyperthreads   = 2
        elif 'amd' in job.constraint:
            hyperthreads   = 2
        else:
            self.error('SLURM input "constraint" must contain either "knl", "haswell" or "amd" on Cori\nyou provided: {0}'.format(job.constraint))
        #end if
        cpus_per_task = int(floor(float(self.cores_per_node)/job.processes_per_node))*hyperthreads
        c='#!/bin/bash\n'
        if job.account is not None:
            c+= '#SBATCH -A '+job.account+'\n'
        #end if
        c+='#SBATCH -p '+job.queue+'\n'
        c+='#SBATCH -C '+str(job.constraint)+'\n'
        c+='#SBATCH -J '+str(job.name)+'\n'
        c+='#SBATCH -t '+job.sbatch_walltime()+'\n'
        c+='#SBATCH -N '+str(job.nodes)+'\n'
        c+='#SBATCH --tasks-per-node={0}\n'.format(job.processes_per_node)
        c+='#SBATCH --cpus-per-task={0}\n'.format(cpus_per_task)
        c+='#SBATCH -o '+job.outfile+'\n'
        c+='#SBATCH -e '+job.errfile+'\n'
        if job.user_env:
            c+='#SBATCH --export=ALL\n'   # equiv to PBS -V
        else:
            c+='#SBATCH --export=NONE\n'
        #end if
        c+='''
echo $SLURM_SUBMIT_DIR
cd $SLURM_SUBMIT_DIR
'''
        if job.threads>1:
            c+='''
export OMP_PROC_BIND=true
export OMP_PLACES=threads
'''
        #end if
        return c
    #end def write_job_header
#end class Cori



# Active
class Perlmutter(NerscMachine):
    name = 'perlmutter'

    def pre_process_job(self,job):
        # Set default queue and node type
        if job.queue is None:
            job.queue = 'regular'
        #end if
        if job.constraint is None:
            job.constraint = 'cpu'
        #end if
        # Account for dual nature of Perlmutter
        if 'cpu' in job.constraint:
            self.nodes          = 3072
            self.procs_per_node = 2
            self.cores_per_node = 128
            self.ram_per_node   = 512
        elif 'gpu' in job.constraint:
            self.nodes          = 1536
            self.procs_per_node = 1
            self.cores_per_node = 64
            self.ram_per_node   = 256
            self.gpus_per_node  = 4
        else:
            self.error('SLURM input "constraint" must contain either "cpu" or "gpu" on Perlmutter\nyou provided: {0}'.format(job.constraint))
        #end if
    #end def pre_process_job

    def write_job_header(self,job):
        self.pre_process_job(job) # sync machine view with job

        # Check if the user gave reasonable processes_per_node
        if 'cpu' in job.constraint:
            if job.processes_per_node > self.cores_per_node:
                self.error('processes_per_node can not be greater than logical CPUs per node (256)\nyou provided: {0}'.format(job.processes_per_node))
            #end if
        elif 'gpu' in job.constraint:
            if job.processes_per_node > self.gpus_per_node:
                self.error('processes_per_node can not be greater than GPUs per node (4)\nyou provided: {0}'.format(job.processes_per_node))
            #end if
            # Also check if the user forgot to include '_g' in the account name for GPU jobs
            if ('_g' in job.account) == False:
                job.account = job.account + '_g'
            #end if
        #end if

        # Check if the user gave reasonable queue inputs
        if job.queue == 'debug':
            base_partition = 1
            max_partition = 8
            max_time = 0.5
        elif job.queue == 'regular':
            base_partition = 1
            max_partition = self.nodes
            max_time = 12
        elif job.queue == 'preempt':
            base_partition = 1
            max_partition = 128
            max_time = 24
        elif job.queue == 'overrun':
            base_partition = 1
            max_partition = self.nodes
            max_time = 12
        else:
            self.error('The requested queue is not implemented.')
        #end if
        job.total_hours = job.days*24 + job.hours + job.minutes/60.0 + job.seconds/3600.0
        if job.total_hours > max_time:
            self.error('The maximum runtime on {0} queue should not be more than {1} hours\n  you requested: {2} hours'.format(job.queue,max_time,job.total_hours))
        #end if
        if job.nodes<base_partition:
            self.error('The number of nodes on {0} queue should not be less than {1}\n  you requested: {2}'.format(job.queue,base_partition,job.nodes))
        elif job.nodes>max_partition:
            self.error('The number of nodes on {0} queue should not be more than {1}\n  you requested: {2}'.format(job.queue,max_partition,job.nodes))
        #end if

        # Use the user cpus_per_task if specified. If not specified, then use available cpus for each process
        if job.cpus_per_task is not None:
            cpus_per_task = job.cpus_per_task
        else:
            hyperthreads = 2 # Both CPU and GPU nodes use the same AMD EPYC 7763 (Milan) CPUs
            cpus_per_task = int(floor(float(self.cores_per_node)/job.processes_per_node))*hyperthreads
        #end if

        c='#!/bin/bash\n'
        if job.account is not None:
            c+= '#SBATCH -A '+job.account+'\n'
        #end if
        c+='#SBATCH -C '+str(job.constraint)+'\n'
        c+='#SBATCH -q '+job.queue+'\n'
        c+='#SBATCH -t '+job.sbatch_walltime()+'\n'
        c+='#SBATCH -N '+str(job.nodes)+'\n'
        c+='#SBATCH --ntasks-per-node={0}\n'.format(job.processes_per_node)
        c+='#SBATCH -c '+str(cpus_per_task)+'\n'
        c+='#SBATCH -J '+str(job.name)+'\n'
        c+='#SBATCH -o '+job.outfile+'\n'
        c+='#SBATCH -e '+job.errfile+'\n'

        if job.email is not None:
            c+='#SBATCH --mail-type=ALL\n'
            c+='#SBATCH --mail-user={0}\n'.format(job.email)

        if 'gpu' in job.constraint:
            gpus_per_task = int(floor(float(self.gpus_per_node)/job.processes_per_node))
            c+='#SBATCH --gpus-per-task={0}\n'.format(gpus_per_task)
        #end if

        if job.user_env:
            c+='#SBATCH --export=ALL\n'   # equiv to PBS -V
        else:
            c+='#SBATCH --export=NONE\n'
        #end if
        c+='''
echo $SLURM_SUBMIT_DIR
cd $SLURM_SUBMIT_DIR
'''
        if (job.threads>1) and ('cpu' in job.constraint):
            c+='''
export OMP_PROC_BIND=true
export OMP_PLACES=threads
'''
        #end if
        if 'gpu' in job.constraint:
            c+='''
export SLURM_CPU_BIND="cores"
'''
        #end if
        return c
    #end def write_job_header
#end class Perlmutter



# Active
class BlueWatersXK(Supercomputer):

    name = 'bluewaters_xk'
    requires_account = False
    batch_capable    = True

    def write_job_header(self,job):
        c='#!/bin/bash\n'
        c+='#PBS -N '+str(job.name)+'\n'
        c+='#PBS -l walltime='+job.pbs_walltime()+'\n'
        c+='#PBS -l nodes={0}:ppn={1}:xk\n'.format(job.nodes,job.ppn)
        c+='#PBS -o '+job.outfile+'\n'
        c+='#PBS -e '+job.errfile+'\n'
        if job.user_env:
            c+='#PBS -V\n'
        #end if
        c+='''
echo $PBS_O_WORKDIR
cd $PBS_O_WORKDIR
'''
        return c
    #end def write_job_header
#end class BlueWatersXK



# Active
class BlueWatersXE(Supercomputer):

    name = 'bluewaters_xe'
    requires_account = False
    batch_capable    = True

    def write_job_header(self,job):
        c='#!/bin/bash\n'
        c+='#PBS -N '+str(job.name)+'\n'
        c+='#PBS -l walltime='+job.pbs_walltime()+'\n'
        c+='#PBS -l nodes={0}:ppn={1}:xe\n'.format(job.nodes,job.ppn)
        c+='#PBS -o '+job.outfile+'\n'
        c+='#PBS -e '+job.errfile+'\n'
        if job.user_env:
            c+='#PBS -V\n'
        #end if
        c+='''
echo $PBS_O_WORKDIR
cd $PBS_O_WORKDIR
'''
        return c
    #end def write_job_header
#end class BlueWatersXE



#Decommissioned
class Titan(Supercomputer):

    name = 'titan'
    requires_account = True
    batch_capable    = True

    def write_job_header(self,job):
        if job.queue is None:
            job.queue = 'batch'
        #end if
        c= '#!/bin/bash\n'
        c+='#PBS -A {0}\n'.format(job.account)
        c+='#PBS -q {0}\n'.format(job.queue)
        c+='#PBS -N {0}\n'.format(job.name)
        c+='#PBS -o {0}\n'.format(job.outfile)
        c+='#PBS -e {0}\n'.format(job.errfile)
        c+='#PBS -l walltime={0}\n'.format(job.pbs_walltime())
        c+='#PBS -l nodes={0}\n'.format(job.nodes)
        #c+='#PBS -l gres=widow3\n'
        c+='#PBS -l gres=atlas1\n'
        if job.user_env:
            c+='#PBS -V\n'
        #end if
        c+='''
echo $PBS_O_WORKDIR
cd $PBS_O_WORKDIR
'''
        return c
    #end def write_job_header
#end class Titan


# Active
class EOS(Supercomputer):

    name = 'eos'
    requires_account = True
    batch_capable    = True

    def post_process_job(self,job):
        if job.threads>1:
            if job.threads<=8:
                job.run_options.add(ss='-ss')
            #end if
            job.run_options.add(cc='-cc numa_node')
        #end if
    #end def post_process_job


    def write_job_header(self,job):
        if job.queue is None:
            job.queue = 'batch'
        #end if
        c= '#!/bin/bash\n'
        c+='#PBS -A {0}\n'.format(job.account)
        c+='#PBS -q {0}\n'.format(job.queue)
        c+='#PBS -N {0}\n'.format(job.name)
        c+='#PBS -o {0}\n'.format(job.outfile)
        c+='#PBS -e {0}\n'.format(job.errfile)
        c+='#PBS -l walltime={0}\n'.format(job.pbs_walltime())
        c+='#PBS -l nodes={0}\n'.format(job.nodes)
        c+='#PBS -l gres=atlas1\n'
        if job.user_env:
            c+='#PBS -V\n'
        #end if
        c+='''
echo $PBS_O_WORKDIR
cd $PBS_O_WORKDIR
'''
        return c
    #end def write_job_header
#end class EOS



class ALCF_Machine(Supercomputer):
    requires_account   = True
    batch_capable      = True
    executable_subfile = True

    prefixed_output    = True
    outfile_extension  = '.output'
    errfile_extension  = '.error'

    base_partition = None

    def post_process_job(self,job):
        job.sub_options.add(
            env  = '--env BG_SHAREDMEMSIZE=32',
            mode = '--mode script'
            )
        #if job.processes<job.nodes: # seems like a good idea, but breaks idempotency
        #    job.processes_per_node=1
        ##end if
        if job.nodes<self.base_partition:
            self.warn('!!! ATTENTION !!!\n  number of nodes on {0} should not be less than {1}\n  you requested: {2}'.format(self.name,self.base_partition,job.nodes))
        else:
            partition = log(float(job.nodes)/self.base_partition)/log(2.)
            if abs(partition-int(partition))>1e-6:
                self.warn('!!! ATTENTION !!!\n  number of nodes on {0} must be {1} times a power of two\n  you requested: {2}\n  nearby valid node count: {3}'.format(self.name,self.base_partition,job.nodes,self.base_partition*2**int(round(partition))))
            #end if
        #end if
        valid_ppn = (1,2,4,8,16,32,64)
        if job.processes_per_node is None:
            self.warn('job may not run properly\nplease specify processes_per_node in each job to be launched with runjob on {0}'.format(self.name))
        elif job.processes_per_node not in valid_ppn:
            self.warn('job may not run properly\nprocesses_per_node is not a valid value for {0}\nprocesses_per_node provided: {1}\nvalid options are: {2}'.format(self.name,job.processes_per_node,valid_ppn))
        #end if
    #end def post_process_job

    def write_job_header(self,job):
        if job.queue is None:
            job.queue = 'default'
        #end if
        c= '#!/bin/bash\n'
        c+='#COBALT -q {0}\n'.format(job.queue)
        c+='#COBALT -A {0}\n'.format(job.account)
        c+='#COBALT -n {0}\n'.format(job.nodes)
        c+='#COBALT -t {0}\n'.format(job.total_minutes())
        c+='#COBALT -O {0}\n'.format(job.identifier)
        c+='\nLOCARGS="--block $COBALT_PARTNAME ${COBALT_CORNER:+--corner} $COBALT_CORNER ${COBALT_SHAPE:+--shape} $COBALT_SHAPE"\n'
        c+='echo "Cobalt location args: $LOCARGS" >&2\n\n'
        return c
    #end def write_job_header
#end class ALCF_Machine

# Decommissioned
class Vesta(ALCF_Machine):
    name = 'vesta'
    base_partition = 32
#end class Vesta

# Decommissioned
class Cetus(ALCF_Machine):
    name = 'cetus'
    base_partition = 128
#end class Cetus

# Decommissioned
class Mira(ALCF_Machine):
    name = 'mira'
    base_partition = 512
#end class Mira


# Active
class Cooley(Supercomputer):
    name = 'cooley'
    requires_account   = True
    batch_capable      = True
    executable_subfile = True

    prefixed_output    = True
    outfile_extension  = '.output'
    errfile_extension  = '.error'

    def post_process_job(self,job):
        #if job.processes_per_node is None and job.threads!=1:
        #    self.error('threads must be 1,2,3,4,6, or 12 on Cooley\nyou provided: {0}'.format(job.threads))
        ##end if

        #job.run_options.add(
            #f   = '-f $COBALT_NODEFILE',
            #ppn = '-ppn {0}'.format(job.processes_per_node),
            #)
        return
    #end def post_process_job

    def write_job_header(self,job):
        if job.queue is None:
            job.queue = 'default'
        #end if
        c= '#!/bin/bash\n'
        c+='#COBALT -q {0}\n'.format(job.queue)
        c+='#COBALT -A {0}\n'.format(job.account)
        c+='#COBALT -n {0}\n'.format(job.nodes)
        c+='#COBALT -t {0}\n'.format(job.total_minutes())
        c+='#COBALT -O {0}\n'.format(job.identifier)
        return c
    #end def write_job_header
#end class Cooley

# Active
class Theta(Supercomputer):
    name = 'theta'
    requires_account   = True
    batch_capable      = True
    executable_subfile = True

    prefixed_output    = True
    outfile_extension  = '.output'
    errfile_extension  = '.error'

    def post_process_job(self,job):
        if job.hyperthreads is None:
            job.hyperthreads = 1
        #end if
        job.run_options.add(
            N  = '-N {0}'.format(job.processes_per_node),
            cc = '-cc depth',
            d  = '-d {0}'.format(job.threads),
            j  = '-j {0}'.format(job.hyperthreads),
            e  = '-e OMP_NUM_THREADS={0}'.format(job.threads),
            )
    #end def post_process_job

    def write_job_header(self,job):
        if job.queue is None:
            job.queue = 'default'
        #end if
        c= '#!/bin/bash\n'
        c+='#COBALT -q {0}\n'.format(job.queue)
        c+='#COBALT -A {0}\n'.format(job.account)
        c+='#COBALT -n {0}\n'.format(job.nodes)
        c+='#COBALT -t {0}\n'.format(job.total_minutes())
        c+='#COBALT -O {0}\n'.format(job.identifier)
        c+='#COBALT --attrs mcdram=cache:numa=quad\n'
        return c
    #end def write_job_header
#end class Theta


# Active
class Lonestar(Supercomputer):  # Lonestar contribution from Paul Young

    name = 'lonestar' # will be converted to lowercase anyway
    requires_account = False
    batch_capable    = True

    def write_job_header(self,job):
        if job.queue is None:
            job.queue = 'batch'
        #end if
        c= '#!/bin/bash\n'
        #c+='#$ -A {0}\n'.format(job.account)
        c+='#$ -q {0}\n'.format(job.queue)
        c+='#$ -N {0}\n'.format(job.name)
        c+='#$ -o {0}\n'.format(job.outfile)
        c+='#$ -e {0}\n'.format(job.errfile)
        c+='#$ -l h_rt={0}\n'.format(job.pbs_walltime())
        c+='#$ -pe 12way {0}\n'.format(job.nodes*12)
        c+='#$ -cwd\n'
        if job.user_env:
            c+='#$ -V\n'
        #end if
        return c
    #end def write_job_header


    def read_process_id(self,output):
        pid = None
        lines = output.splitlines()

        for line in lines:
            if 'Your job' in line:
                spid = line.split(' ')[2]
                if spid.isdigit():
                    pid = int(spid)
                #end if
            #end if
        #end for
        return pid
    #end def read_process_id
#end class Lonestar



class ICMP_Machine(Supercomputer): # ICMP and Amos contributions from Ryan McAvoy
    batch_capable      = True
    executable_subfile = True

    prefixed_output    = True
    outfile_extension  = '.output'
    errfile_extension  = '.error'

    def write_job_header(self,job):
        if job.queue is None:
            job.queue = 'defq'
        #end if
        c= '#!/bin/bash -x\n'
        c+='#SBATCH --export=ALL\n'
        c+='#SBATCH -J {0}\n'.format(job.identifier)
        c+='#SBATCH -p {0}\n'.format(job.queue)
        c+='#SBATCH -o {0}\n'.format(job.outfile)
        c+='#SBATCH -e {0}\n'.format(job.errfile)
        c+='#SBATCH --nodes {0}\n'.format(job.nodes)
        c+='#SBATCH --ntasks-per-node={0}\n'.format(job.processes_per_node)
        c+='#SBATCH --cpus-per-task={0}\n'.format(job.threads)
        c+='#SBATCH -t {0}:{1}:{2}\n'.format(str(job.hours+24*job.days).zfill(2),str(job.minutes).zfill(2),str(job.seconds).zfill(2))
        return c
    #end def write_job_header
#end class ICMP_Machine

# Unknown
class Komodo(ICMP_Machine):
    name = 'komodo'
#end class Komodo
# Unknown
class Matisse(ICMP_Machine):
    name = 'matisse'
#end class Matisse


# Active
class Amos(Supercomputer):
    name = 'amos'

    #requires_account   = True
    batch_capable      = True
    executable_subfile = True

    prefixed_output    = True
    outfile_extension  = '.output'
    errfile_extension  = '.error'

    def write_job_header(self,job):
        if job.queue is None:
            job.queue = 'debug'
        #end if
        if job.queue == 'debug':
            base_partition = 1
            max_partition = 32
            max_time =1
        elif job.queue == 'small':
            base_partition = 1
            max_partition = 64
            max_time =24
        elif job.queue == 'medium':
            base_partition = 128
            max_partition = 512
            max_time =12
        elif job.queue == 'large':
            base_partition = 1024
            max_partition = 2048
            max_time =6
        elif job.queue == 'verylarge':
            base_partition = 3072
            max_partition = 4096
            max_time =6
        #end if
        job.total_hours = job.days*24 + job.hours + job.minutes/60.0 + job.seconds/3600.0
        if job.total_hours > max_time:
            self.warn('!!! ATTENTION !!!\n  the maximum runtime on {0} should not be more than {1}\n  you requested: {2}'.format(job.queue,max_time,job.total_hours))
            job.hours   = max_time
            job.minutes =0
            job.seconds =0
        #end if
        if job.nodes<base_partition:
            self.warn('!!! ATTENTION !!!\n  number of nodes in {0} should not be less than {1}\n  you requested: {2}'.format(job.queue,base_partition,job.nodes))
        elif job.nodes>max_partition:
            self.warn('!!! ATTENTION !!!\n  number of nodes in {0} should not be more than {1}\n  you requested: {2}'.format(job.queue,max_partition,job.nodes))
        else:
            if job.queue != 'verylarge':
                partition = log(float(job.nodes)/base_partition)/log(2.)
                if abs(partition-int(partition))>1e-6:
                    self.warn('!!! ATTENTION !!!\n  number of nodes on {0} must be {1} times a power of two\n  you requested: {2}\n  nearby valid node count: {3}'.format(self.name,base_partition,job.nodes,base_partition*2**int(round(partition))))
            elif job.nodes != 3072 and job.nodes != 4096:
                self.warn('!!! ATTENTION !!!\n  number of nodes on {0} must be 3072 or 4096 you requested {1}'.format(self.name,job.nodes))
            #end if
        #end if

        c= '#!/bin/bash -x\n'
        c+='#SBATCH --export=ALL\n'
        #c+=#SBATCH -D /gpfs/sb/data/<project>/<user>/
        c+='#SBATCH -J {0}\n'.format(job.identifier)
        c+='#SBATCH -p {0}\n'.format(job.queue)
        c+='#SBATCH -o {0}\n'.format(job.outfile)
        c+='#SBATCH -e {0}\n'.format(job.errfile)
        c+='#SBATCH --nodes {0}\n'.format(job.nodes)
        c+='#SBATCH --ntasks-per-node={0}\n'.format(job.processes_per_node)
        c+='#SBATCH --cpus-per-task={0}\n'.format(job.threads)
        c+='#SBATCH -t {0}:{1}:{2}\n'.format(str(job.hours+24*job.days).zfill(2),str(job.minutes).zfill(2),str(job.seconds).zfill(2))
        # c+='#SBATCH --mail-type=ALL'
        # c+='#SBATCH --mail-user=<{0}>'

        return c
    #end def write_job_header
#end class Amos


class SnlMachine(Supercomputer):
    requires_account   = True
    batch_capable      = True
    #executable_subfile = True

    prefixed_output    = True
    outfile_extension  = '.output'
    errfile_extension  = '.error'

    #for mpiexec
    def post_process_job(self,job):
        job.run_options.add(bindto="--bind-to none",npernode="--npernode {}".format(job.processes_per_node))

    def write_job_header(self,job):
        if job.queue is None:
            job.queue='batch'
        #end if

        cpus_per_task = int(floor(float(self.cores_per_node)/job.processes_per_node))

        if job.qos == 'long':
            max_time = 96
        elif 'short' in job.queue:
            max_time = 4
        else:
            max_time = 48
        #end if

        job.total_hours = job.days*24 + job.hours + job.minutes/60.0 + job.seconds/3600.0
        if job.total_hours > max_time:   # warn if job will take more than 48 hrs.
            if job.qos == 'long':
                self.warn('!!! ATTENTION !!!\n  the maximum runtime on {0} should not be more than {1} with --qos=\'long\'\n  you requested: {2}'.format(job.queue,max_time,job.total_hours))
            elif 'short' in job.queue:
                self.warn('!!! ATTENTION !!!\n  the maximum runtime on {0} should not be more than {1} with -p short[,batch]\n  you requested: {2}'.format(job.queue,max_time,job.total_hours))
            else:
                self.warn('!!! ATTENTION !!!\n  the maximum runtime on {0} should not be more than {1}\n  you requested: {2}'.format(job.queue,max_time,job.total_hours))
            #end if
            job.hours   = max_time
            job.minutes = 0
            job.seconds = 0
        #end if


        c='#!/bin/bash\n'
        c+='#SBATCH -p '+str(job.queue)+'\n'
        c+='#SBATCH --job-name '+str(job.name)+'\n'
        c+='#SBATCH --account='+str(job.account)+'\n'
        c+='#SBATCH -N '+str(job.nodes)+'\n'
        c+='#SBATCH -t {0}:{1}:{2}\n'.format(str(job.hours+24*job.days).zfill(2),str(job.minutes).zfill(2),str(job.seconds).zfill(2))
        c+='#SBATCH -o {0}\n'.format(job.outfile)
        c+='#SBATCH -e {0}\n'.format(job.errfile)
        if job.qos:
            c+='#SBATCH --qos={}\n'.format(job.qos)
        c+='\n'
        return c
    #end def write_job_header
#end class SnlMachine
#Unknown
class Chama(SnlMachine):
    name = 'chama'
#end class Chama
#Unknown
class Skybridge(SnlMachine):
    name = 'skybridge'
#end class Skybridge
#Unknown
class Eclipse(SnlMachine):
    name = 'eclipse'
#end class Eclipse
#Unknown
class Attaway(SnlMachine):
    name = 'attaway'
#end class Attaway
#Unknown
class Manzano(SnlMachine):
    name = 'manzano'
#end class Manzano
#Unknown
class Ghost(SnlMachine):
    name = 'ghost'
#end class Ghost
#Unknown
class Amber(SnlMachine):
    name = 'amber'
#end class Amber
#Unknown
class Uno(SnlMachine):
    name = 'uno'
#end class Uno
#Unknown
class Solo(SnlMachine):
    name = 'solo'
#end class Solo


# Active
# machines at LRZ  https://www.lrz.de/english/
class SuperMUC(Supercomputer):
    name = 'supermuc'
    requires_account    = False
    batch_capable       = True
    query_with_username = False

    def write_job_header(self,job):
        if job.queue is None:
            job.queue = 'general'
        #end if
        if job.type is None:
            job.type = 'MPICH'
        else:
            job.type = job.type.lower()
            if job.type=='mpich':
                job.type=job.type.upper()
            #end if
        #end if
        ibm   = job.type=='parallel'
        intel = job.type=='MPICH'
        omp   = isinstance(job.threads,int) and job.threads>1
        if not ibm and not intel:
            self.error('the only types of MPI supported are "parallel" and "MPICH"\nreceived MPI with type: {0}'.format(job.type))
        #end if
        c ='#!/bin/bash\n'
        c+='#@ job_name         = {0}\n'.format(job.name)
        c+='#@ job_type         = {0}\n'.format(job.type)
        c+='#@ class            = {0}\n'.format(job.queue)
        c+='#@ node             = {0}\n'.format(job.nodes)
        if job.nodes<512:
            icmin = 1
            icmax = 1
        else:
            icmin = int(job.nodes//512)+1
            icmax = icmin+1
        #end if
        c+='#@ island_count     = {0},{1}\n'.format(icmin,icmax)
        if intel and omp:
            c+='#@ tasks_per_node   = {0}\n'.format(job.processes_per_node)
        else:
            c+='#@ total_tasks      = {0}\n'.format(job.processes)
        #end if
        c+='#@ wall_clock_limit = {0}\n'.format(job.ll_walltime())
        c+='#@ network.MPI      = sn_all,not_shared,us\n'
        c+='#@ initialdir       = {0}\n'.format(job.abs_dir)
        c+='#@ output           = {0}\n'.format(job.outfile)
        c+='#@ error            = {0}\n'.format(job.errfile)
        c+='#@ energy_policy_tag = my_energy_tag\n'
        c+='#@ minimize_time_to_solution = yes\n'
        if job.email is None:
            c+='#@ notification     = never\n'
        else:
            c+='#@ notification     = always\n'
            c+='#@ notify_user      = {0}\n'.format(job.email)
        #end if
        c+='#@ queue\n'
        c+='. /etc/profile\n'
        c+='. /etc/profile.d/modules.sh\n'
        if ibm and omp:
            c+='export MP_SINGLE_THREAD=no\n'
            c+='export MP_TASK_AFFINITY=core:{0}\n'.format(job.threads)
        elif intel and not omp:
            c+='module unload mpi.ibm\n'
            c+='module load mpi.intel\n'
        elif intel and omp:
            c+='module unload mpi.ibm\n'
            c+='module load mpi.intel\n'
            c+='export OMP_NUM_THREADS={0}\n'.format(job.threads)
            #c+='module load mpi_pinning/hybrid_blocked\n'
        #end if
        return c
    #end def write_job_header
#end class SuperMUC


# Active
class SuperMUC_NG(Supercomputer):
    name                = 'supermucng'
    requires_account    = True
    batch_capable       = True
    query_with_username = True

    def write_job_header(self,job):
        if job.queue is None:
            job.queue = 'general'
        #end if
        if job.hyperthreads is None:
            job.hyperthreads = job.processes_per_node//48
            if job.hyperthreads==0:
                job.hyperthreads=None
            #end if
        #end if
        if job.constraint is None:
            job.contraint = 'scratch&work'
        #end if
        c ='#!/bin/bash\n'
        c+='#SBATCH --account={}\n'.format(job.account)
        c+='#SBATCH --partition={}\n'.format(job.queue)
        c+='#SBATCH -J {}\n'.format(job.name)
        c+='#SBATCH --time={}\n'.format(job.sbatch_walltime())
        c+='#SBATCH -o ./{}\n'.format(job.outfile)
        c+='#SBATCH -e ./{}\n'.format(job.errfile)
        if job.switches is not None:
            c+='#SBATCH --switches={}\n'.format(job.switches)
        #end if
        c+='#SBATCH --nodes={}\n'.format(job.nodes)
        c+='#SBATCH --ntasks-per-node={}\n'.format(job.processes_per_node)
        if job.ntasks_per_core is not None:
            c+='#SBATCH --ntasks-per-core={}\n'.format(job.ntasks_per_core)
        elif job.hyperthreads is not None:
            c+='#SBATCH --ntasks-per-core={}\n'.format(job.hyperthreads)
        #end if
        if not job.default_cpus_per_task:
            if job.cpus_per_task is None:
                c+='#SBATCH --cpus-per-task={}\n'.format(job.threads)
            else:
                c+='#SBATCH --cpus-per-task={}\n'.format(job.cpus_per_task)
            #end if
        #end if
        c+='#SBATCH -D ./\n'
        c+='#SBATCH --no-requeue\n'
        if job.constraint is not None:
            c+='#--constraint="{}"\n'.format(job.constraint)
        #end if
        if job.email is not None:
            c+='#SBATCH --mail-type=ALL\n'
            c+='#SBATCH --mail-user={}\n'.format(job.email)
        #end if
        c+='#SBATCH --export=NONE\n'
        if job.user_env:
            c+='#SBATCH --get-user-env\n'
        #end if
        return c
    #end def write_job_header
#end class SuperMUC_NG


# Decommissioned
class Stampede2(Supercomputer):
    name = 'stampede2'

    requires_account   = True
    batch_capable      = True
    #executable_subfile = True

    prefixed_output    = True
    outfile_extension  = '.output'
    errfile_extension  = '.error'

    def write_job_header(self,job):
        if job.queue is None:
            job.queue='normal'
        #end if
        
        if job.queue == 'development':
            max_nodes = 16
            max_time = 2
        elif job.queue == 'normal':
            max_nodes = 256
            max_time = 48
        elif job.queue == 'large':
            max_nodes = 2048
            max_time = 48
        elif job.queue == 'long':
            max_nodes = 32
            max_time = 96
        elif job.queue == 'flat_quadrant':
            max_nodes = 24
            max_time = 48
        elif job.queue == 'skx-dev':
            max_nodes = 4
            max_time = 2
        elif job.queue == 'skx-normal':
            max_nodes = 128
            max_time = 48
        elif job.queue == 'skx-large':
            max_nodes = 868
            max_time = 48
        #end if
        
        if 'skx' in job.queue:
            max_processes_per_node = 48
        else:
            max_processes_per_node = 68
        #end if
        job.total_hours = job.days*24 + job.hours + job.minutes/60.0 + job.seconds/3600.0
        if job.total_hours > max_time:
            self.warn('!!! ATTENTION !!!\n  the maximum runtime on {0} should not be more than {1}\n  you requested: {2}'.format(job.queue,max_time,job.total_hours))
            job.hours   = max_time
            job.minutes =0
            job.seconds =0
        #end if
        
        if job.nodes > max_nodes:
            self.warn('!!! ATTENTION !!!\n  the maximum nodes on {0} should not be more than {1}\n  you requested: {2}'.format(job.queue,max_nodes,job.nodes))
            job.nodes = max_nodes
        #end if
        
        if job.processes_per_node > max_processes_per_node:
            self.warn('!!! ATTENTION !!!\n  the maximum number of processes per node on {0} should not be more than {1}\n  you requested: {2}'.format(job.queue,max_processes_per_node,job.processes_per_node))
            job.processes_per_node = max_processes_per_node
        #end if
        
        c='#!/bin/bash\n'
        c+='#SBATCH --job-name '+str(job.name)+'\n'
        c+='#SBATCH --account='+str(job.account)+'\n'
        c+='#SBATCH -N '+str(job.nodes)+'\n'
        c+='#SBATCH --ntasks-per-node={0}\n'.format(job.processes_per_node)
        c+='#SBATCH --cpus-per-task={0}\n'.format(job.threads)
        c+='#SBATCH -t {0}:{1}:{2}\n'.format(str(job.hours+24*job.days).zfill(2),str(job.minutes).zfill(2),str(job.seconds).zfill(2))
        c+='#SBATCH -o {0}\n'.format(job.outfile)
        c+='#SBATCH -e {0}\n'.format(job.errfile)
        c+='#SBATCH -p {0}\n'.format(job.queue)
        c+='\n'
        return c
    #end def write_job_header
#end class Stampede2


# Decommissioned
# CADES at ORNL
class CadesMoab(Supercomputer):
    name = 'cades_moab'
    requires_account = True
    batch_capable    = True

    def post_process_job(self,job):
        ppn = job.processes_per_node
        if job.threads>1 and ppn is not None and ppn>1:
            processes_per_socket = int(floor(job.processes_per_node/2))
            job.run_options.add(npersocket='--npersocket {0}'.format(processes_per_socket))
        #end if
    #end def post_process_job

    def write_job_header(self,job):
        if job.queue is None:
            job.queue = 'skylake'
        #end if
        if job.qos is None:
            job.qos = 'std'
        #end if
        if job.group_list is None:
            job.group_list = 'cades-'+job.account
        #end if
        c= '#!/bin/bash\n'
        c+='#PBS -A {0}\n'.format(job.account)
        c+='#PBS -W group_list={0}\n'.format(job.group_list)
        c+='#PBS -q {0}\n'.format(job.queue)
        c+='#PBS -N {0}\n'.format(job.name)
        c+='#PBS -o {0}\n'.format(job.outfile)
        c+='#PBS -e {0}\n'.format(job.errfile)
        c+='#PBS -l qos={0}\n'.format(job.qos) # This could be qos=burst as well, but then it can be cancelled by others
        c+='#PBS -l walltime={0}\n'.format(job.pbs_walltime())
        c+='#PBS -l nodes={0}:ppn={1}\n'.format(job.nodes, job.ppn)
        c+='''
echo $PBS_O_WORKDIR
cd $PBS_O_WORKDIR
'''
        return c
    #end def write_job_header
#end class CadesMoab


# Active
class CadesSlurm(Supercomputer):
    name = 'cades'
    requires_account = True
    batch_capable    = True

    def write_job_header(self,job):
        if job.queue is None:
            job.queue = 'skylake'
        #end if

        c  = '#!/bin/bash\n'
        c += '#SBATCH -A {}\n'.format(job.account)
        c += '#SBATCH -p {}\n'.format(job.queue)
        c += '#SBATCH -J {}\n'.format(job.name)
        c += '#SBATCH -t {}\n'.format(job.sbatch_walltime())
        c += '#SBATCH -N {}\n'.format(job.nodes)
        c += '#SBATCH --ntasks-per-node={0}\n'.format(job.processes_per_node)
        c += '#SBATCH --cpus-per-task={0}\n'.format(job.threads)
        c += '#SBATCH --mem=0\n' # required on Cades
        c += '#SBATCH -o '+job.outfile+'\n'
        c += '#SBATCH -e '+job.errfile+'\n'
        c += '#SBATCH --exclusive\n'
        if job.user_env:
            c += '#SBATCH --export=ALL\n'   # equiv to PBS -V
        else:
            c += '#SBATCH --export=NONE\n'
        #end if

        return c
    #end def write_job_header
#end class CadesSlurm


# Active
# Inti at ORNL
class Inti(Supercomputer):
    name = 'inti'
    requires_account = False
    batch_capable    = True

    def write_job_header(self,job):
        if job.queue is None:
            job.queue = 'QMCREGULAR'
        #end if

        c  = '#!/bin/bash\n'
        # c += '#SBATCH -A {}\n'.format(job.account)
        c += '#SBATCH -p {}\n'.format(job.queue)
        c += '#SBATCH -J {}\n'.format(job.name)
        c += '#SBATCH -t {}\n'.format(job.sbatch_walltime())
        c += '#SBATCH -N {}\n'.format(job.nodes)
        c += '#SBATCH --ntasks-per-node={0}\n'.format(job.processes_per_node)
        c += '#SBATCH --cpus-per-task={0}\n'.format(job.threads)
        c += '#SBATCH --mem=0\n' # required on Cades
        c += '#SBATCH -o '+job.outfile+'\n'
        c += '#SBATCH -e '+job.errfile+'\n'
        c += '#SBATCH --exclusive\n'
        if job.user_env:
            c += '#SBATCH --export=ALL\n'   # equiv to PBS -V
        else:
            c += '#SBATCH --export=NONE\n'
        #end if

        return c
    #end def write_job_header
#end class Inti


# Active
# Baseline at ORNL https://docs.cades.olcf.ornl.gov/baseline_user_guide/baseline_user_guide.html
class Baseline(Supercomputer):
    name = 'baseline'
    requires_account = True
    batch_capable    = True
    queue_configs={
        'default': 'batch_cnms',
        'batch': {
            'max_nodes': 138,
            'max_walltime': '24:00:00',
        },
        'batch_low_memory': {
            'max_nodes': 68,
            'max_walltime': '24:00:00',
        },
        'batch_high_memory': {
            'max_nodes': 70,
            'max_walltime': '24:00:00',
        },
        'batch_ccsi': {
            'max_nodes': 20,
            'max_walltime': '24:00:00',
        },
        'batch_cnms': {
            'max_nodes': 20,
            'max_walltime': '24:00:00',
        },
        'gpu_acmhs': {
            'max_nodes': 1,
            'max_walltime': '24:00:00',
        }
    }
    def write_job_header(self,job):
        self.validate_queue_config(job)

        c  = '#!/bin/bash\n'
        c += '#SBATCH -A {}\n'.format(job.account)
        c += '#SBATCH -p {}\n'.format(job.queue)
        c += '#SBATCH -J {}\n'.format(job.name)
        c += '#SBATCH -t {}\n'.format(job.sbatch_walltime())
        c += '#SBATCH -N {}\n'.format(job.nodes)
        c += '#SBATCH --ntasks-per-node={0}\n'.format(job.processes_per_node)
        c += '#SBATCH --cpus-per-task={0}\n'.format(job.threads)
        c += '#SBATCH -o '+job.outfile+'\n'
        c += '#SBATCH -e '+job.errfile+'\n'
        if job.user_env:
            c += '#SBATCH --export=ALL\n'   # equiv to PBS -V
        else:
            c += '#SBATCH --export=NONE\n'
        #end if

        return c
    #end def write_job_header
#end class Baseline

class Frontier(Supercomputer):
    name = 'frontier'
    requires_account = True
    batch_capable    = True

    queue_configs = {
        'default': 'batch',
        'batch': {
            'max_nodes': 9664,
            'max_walltime': '12:00:00',
        },
        'extended': {
            'max_nodes': 1894,
            'max_walltime': '24:00:00',
        },
        'debug': {
            'max_nodes': 1,
            'max_walltime': '01:00:00',
        }
    }

    def pre_process_job(self,job):
        # Set default queue and node type
        if job.queue is None:
            job.queue = 'batch'
        #end if
        if job.constraint is None:
            job.constraint = 'cpu'
        #end if
        # Account for dual nature of Frontier
        if 'cpu' in job.constraint:
            self.cores_per_node = 56
        elif 'gpu' in job.constraint:
            self.cores_per_node = 56
            self.gpus_per_node  = 4
        else:
            self.error('SLURM input "constraint" must contain either "cpu" or "gpu" on Frontier\nyou provided: {0}'.format(job.constraint))
        #end if
    #end def pre_process_job
    

    def post_process_job(self, job):
        if 'cpu' in job.constraint:
            job.run_options.add(
                cpu_bind='--cpu-bind=threads',
                threads_per_core='--threads-per-core={0}'.format(job.threads)
            )            
        elif 'gpu' in job.constraint:
            gpus_per_task = int(floor(float(self.gpus_per_node)/job.processes_per_node))
            job.run_options.add(
                gpu_bind='--gpu-bind=closest',
                gpus_per_task='--gpus-per-task={0}'.format(gpus_per_task)
            )
        #end if
        job.run_options.add(
            N='-N {}'.format(job.nodes),
            n='-n {}'.format(job.processes),
            c='-c {}'.format(job.threads),

        )

    def write_job_header(self, job):
        self.validate_queue_config(job)
        if job.queue is None:
            job.queue = 'batch'
        elif job.queue == 'debug':
            job.qos = 'debug'
            job.queue = 'batch'
        #end if

        c = '#!/bin/sh\n'
        c += '#SBATCH -A {account}\n'.format(account=job.account)
        c += '#SBATCH -p {queue}\n'.format(queue=job.queue)
        c += '#SBATCH -J {name}\n'.format(name=job.name)
        c += '#SBATCH -t {time}\n'.format(time=job.lsf_walltime())
        c += '#SBATCH -N {nodes}\n'.format(nodes=job.nodes)
        c += '#SBATCH -S 8\n' # Uses default low-noise mode layout(default), reduces number of cores from 64 to 56. 
        c += '#SBATCH -o {name}.out\n'.format(name=job.name)
        c += '#SBATCH -e {name}.err\n'.format(name=job.name)
        return c


#end class Frontier


# Active 
# BESMS is at ORNL 
class Besms(Supercomputer):
    name = 'besms'
    requires_account = True
    batch_capable    = True
    # Using sinfo to get the queue configs
    queue_configs={
        'default': 't92',
        't92': {
            'max_nodes': 10,
            'max_walltime': '672:00:00',
        },
        't92-burst': {
            'max_nodes': 25,
            'max_walltime': '672:00:00',
        },
        'burst': {
            'max_nodes': 166,
            'max_walltime': '48:00:00',
        },
    }
    def write_job_header(self,job):
        self.validate_queue_config(job)

        c  = '#!/bin/bash\n'
        c += '#SBATCH -A {}\n'.format(job.account)
        c += '#SBATCH -p {}\n'.format(job.queue)
        c += '#SBATCH -J {}\n'.format(job.name)
        c += '#SBATCH -t {}\n'.format(job.sbatch_walltime())
        c += '#SBATCH -N {}\n'.format(job.nodes)
        c += '#SBATCH --ntasks-per-node={0}\n'.format(job.processes_per_node)
        c += '#SBATCH --cpus-per-task={0}\n'.format(job.threads)
        c += '#SBATCH -o '+job.outfile+'\n'
        c += '#SBATCH -e '+job.errfile+'\n'
        c += '#SBATCH --mem=350G\n'
        c += '#SBATCH --exclusive\n'

        if job.user_env:
            c += '#SBATCH --export=ALL\n'   # equiv to PBS -V
        else:
            c += '#SBATCH --export=NONE\n'
        #end if

        return c
    #end def write_job_header
#end class Baseline

# Decommissioned
# Summit at ORNL
class Summit(Supercomputer):

    name = 'summit'
    requires_account = True
    batch_capable    = True

    def post_process_job(self,job):
        # add the options only if the user has not supplied options
        if len(job.run_options)==0:
            opt = obj(
                launch_dist = '-d packed',
                bind        = '-b rs',
                )
            if job.gpus is None:
                job.gpus = 6 # gpus to use per node
            #end if
            if job.alloc_flags is None:
                job.alloc_flags = 'smt1'
            #end if
            if job.gpus==0:
                if job.processes%2==0:
                    resource_sets_per_node = 2
                else:
                    resource_sets_per_node = 1
                #end if
                nrs   = job.nodes*resource_sets_per_node
                pprs  = job.processes_per_node//resource_sets_per_node
                gpurs = 0
            else:
                ppn = job.processes_per_node
                if ppn is None:
                    self.warn('job may not run properly on Summit\nat least one mpi process should be present for each node\nplease check the generated bsub file for correctness')
                    ppn = 0
                #end if
                if ppn%job.gpus!=0:
                    self.warn('job may not run properly on Summit\nprocesses per node should divide evenly into number of gpus requested\nprocesses per node requested: {0}\ngpus per node requested: {1}\nplease check the generated bsub file for correctness'.format(job.processes_per_node,job.gpus))
                #end if
                resource_sets_per_node = job.gpus
                nrs   = job.nodes*resource_sets_per_node
                pprs  = ppn//resource_sets_per_node
                gpurs = 1
            #end if
            opt.set(
                resource_sets= '-n {0}'.format(nrs),
                rs_per_node  = '-r {0}'.format(resource_sets_per_node),
                tasks_per_rs = '-a {0}'.format(pprs),
                cpus_per_rs  = '-c {0}'.format(pprs*job.threads),
                gpus_per_rs  = '-g {0}'.format(gpurs),
                )
            job.run_options.add(**opt)
        #end if
    #end def post_process_job


    def write_job_header(self,job):
        c ='#!/bin/bash\n'
        c+='#BSUB -P {0}\n'.format(job.account)
        if job.queue is not None:
            c+='#BSUB -q {0}\n'.format(job.queue)
        #end if
        c+='#BSUB -J {0}\n'.format(job.name)
        c+='#BSUB -o {0}\n'.format(job.outfile)
        c+='#BSUB -e {0}\n'.format(job.errfile)
        c+='#BSUB -W {0}\n'.format(job.lsf_walltime())
        c+='#BSUB -nnodes {0}\n'.format(job.nodes)
        if job.alloc_flags is not None:
            c+='#BSUB -alloc_flags "{0}"\n'.format(job.alloc_flags)
        #end if
        return c
    #end def write_job_header


    def read_process_id(self,output):
        pid = None
        tokens = output.split()
        for t in tokens:
            if t.startswith('<'):
                spid = t.strip('<>').strip()
                if spid.isdigit():
                    pid = int(spid)
                    break
                #end if
            #end if
        #end for
        return pid
    #end def read_process_id
#end class Summit

# Unknown
## Added 28/11/2019 by A Zen
class Rhea(Supercomputer):

    name = 'rhea'
    requires_account   = True
    batch_capable      = True
    #executable_subfile = True
    prefixed_output    = True
    outfile_extension  = '.output'
    errfile_extension  = '.error'

    def post_process_job(self,job):
        job.run_options.add(
            N='-N {}'.format(job.nodes),
            n='-n {}'.format(job.processes),
            )
        if job.threads>1:
            job.run_options.add(
                c = '-c {}'.format(job.threads),
                )
            if 'cpu_bind' not in job.run_options:
                if job.processes_per_node==self.cores_per_node:
                    cpu_bind = '--cpu-bind=threads'
                else:
                    cpu_bind = '--cpu-bind=cores'
                #end if
                job.run_options.add(
                    cpu_bind = cpu_bind
                    )
            #end if
        #end if
    #end def post_process_job

    def write_job_header(self,job):
        if job.queue is None:
            job.queue='batch'
        #end if
        base_partition = None
        max_partition = 384
        if job.nodes <= 16:
            max_time = 48
        elif job.nodes <= 64:
            max_time = 36
        else:
            max_time = 3
        job.total_hours = job.days*24 + job.hours + job.minutes/60.0 + job.seconds/3600.0
        if job.total_hours > max_time:   # warn if job will take more than 96 hrs.
            self.warn('!!! ATTENTION !!!\n  the maximum runtime on {0} should not be more than {1}\n  you requested: {2}'.format(job.queue,max_time,job.total_hours))
            job.hours   = max_time
            job.minutes =0
            job.seconds =0
        #end if

        c='#!/bin/bash\n'
        c+='#SBATCH --job-name '+str(job.name)+'\n'
        c+='#SBATCH --account='+str(job.account)+'\n'
        c+='#SBATCH -N '+str(job.nodes)+'\n'
        c+='#SBATCH -t {0}:{1}:{2}\n'.format(str(job.hours+24*job.days).zfill(2),str(job.minutes).zfill(2),str(job.seconds).zfill(2))
        c+='#SBATCH -o {0}\n'.format(job.outfile)
        c+='#SBATCH -e {0}\n'.format(job.errfile)
        if job.email is not None:
            c+='#SBATCH --mail-user {}\n'.format(job.email)
            c+='#SBATCH --mail-type ALL\n'
            #c+='#SBATCH --mail-type FAIL\n'
        #end if
        c+='\n'
        c+='cd $SLURM_SUBMIT_DIR\n'
        c+='\n'
        c+='echo JobID : $SLURM_JOBID \n'
        c+='echo Number of nodes requested: $SLURM_JOB_NUM_NODES \n'
        c+='echo List of nodes assigned to the job: $SLURM_NODELIST \n'
        c+='\n'
        return c
    #end def write_job_header
#end class Rhea

# Active
# Andes at ORNL
## Added 19/03/2021 by A Zen
class Andes(Supercomputer):

    name = 'andes'
    requires_account   = True
    batch_capable      = True
    #executable_subfile = True
    prefixed_output    = True
    outfile_extension  = '.output'
    errfile_extension  = '.error'

    def post_process_job(self,job):
        if job.threads>1:
            job.run_options.add(
                c = '-c {}'.format(job.threads),
                )
            if 'cpu_bind' not in job.run_options:
                if job.processes_per_node==self.cores_per_node:
                    cpu_bind = '--cpu-bind=threads'
                else:
                    cpu_bind = '--cpu-bind=cores'
                #end if
                job.run_options.add(
                    cpu_bind = cpu_bind
                    )
            #end if
        #end if
        job.run_options.add(
            N='-N {}'.format(job.nodes),
            n='-n {}'.format(job.processes),
            )
    #end def post_process_job

    def write_job_header(self,job):
        if job.queue is None:
            job.queue='batch'
        #end if
        base_partition = None
        max_partition = 384
        if job.nodes <= 16:
            max_time = 48
        elif job.nodes <= 64:
            max_time = 36
        else:
            max_time = 3
        job.total_hours = job.days*24 + job.hours + job.minutes/60.0 + job.seconds/3600.0
        if job.total_hours > max_time:   # warn if job will take more than 96 hrs.
            self.warn('!!! ATTENTION !!!\n  the maximum runtime on {0} should not be more than {1}\n  you requested: {2}'.format(job.queue,max_time,job.total_hours))
            job.hours   = max_time
            job.minutes =0
            job.seconds =0
        #end if

        c='#!/bin/bash\n'
        c+='#SBATCH --job-name '+str(job.name)+'\n'
        c+='#SBATCH --account='+str(job.account)+'\n'
        c+='#SBATCH -N '+str(job.nodes)+'\n'
        c+='#SBATCH -t {0}:{1}:{2}\n'.format(str(job.hours+24*job.days).zfill(2),str(job.minutes).zfill(2),str(job.seconds).zfill(2))
        c+='#SBATCH -o {0}\n'.format(job.outfile)
        c+='#SBATCH -e {0}\n'.format(job.errfile)
        if job.email is not None:
            c+='#SBATCH --mail-user {}\n'.format(job.email)
            c+='#SBATCH --mail-type ALL\n'
            #c+='#SBATCH --mail-type FAIL\n'
        #end if
        c+='\n'
        c+='cd $SLURM_SUBMIT_DIR\n'
        c+='\n'
        c+='echo JobID : $SLURM_JOBID \n'
        c+='echo Number of nodes requested: $SLURM_JOB_NUM_NODES \n'
        c+='echo List of nodes assigned to the job: $SLURM_NODELIST \n'
        c+='\n'
        return c
    #end def write_job_header
#end class Andes

# Active
## Added 05/04/2022 by A Zen
class Archer2(Supercomputer):
    # https://docs.archer2.ac.uk/user-guide/hardware/

    name = 'archer2'
    requires_account   = True
    batch_capable      = True
    #executable_subfile = True
    prefixed_output    = True
    outfile_extension  = '.output'
    errfile_extension  = '.error'

    def post_process_job(self,job):
        job.run_options.add(
            distribution='--distribution=block:block',
            hint='--hint=nomultithread',
            N='-N {}'.format(job.nodes),
            n='-n {}'.format(job.processes),
            )
        if job.threads>1:
            job.run_options.add(
                c = '-c {}'.format(job.threads),
                )
#           if 'cpu_bind' not in job.run_options:
#               if job.processes_per_node==self.cores_per_node:
#                   cpu_bind = '--cpu-bind=threads'
#               else:
#                   cpu_bind = '--cpu-bind=cores'
#               #end if
#               job.run_options.add(
#                   cpu_bind = cpu_bind
#                   )
            #end if
        #end if
    #end def post_process_job

    def write_job_header(self,job):
        if job.qos is None:
            job.qos='standard'
        #end if
        base_partition = None
        if job.qos == 'long':
            max_time = 48
            max_partition = 64
        elif 'short' in job.qos:
            max_time = 20.0/60.0
            max_partition = 32
        else:
            max_time = 24
            max_partition = 1024
        #end if
        job.total_hours = job.days*24 + job.hours + job.minutes/60.0 + job.seconds/3600.0
        if job.total_hours > max_time:   
            self.warn('!!! ATTENTION !!!\n  the maximum runtime on {0} should not be more than {1}\n  you requested: {2}'.format(job.queue,max_time,job.total_hours))
            job.hours   = max_time
            job.minutes =0
            job.seconds =0
        #end if
        if job.nodes > max_partition:   
            self.warn('!!! ATTENTION !!!\n  the maximum nodes on {0} should not be more than {1}\n  you requested: {2}'.format(job.queue,max_partition,job.nodes))
            job.nodes   = max_partition
        #end if

        c='#!/bin/bash\n'
        c+='#SBATCH --job-name '+str(job.name)+'\n'
        c+='#SBATCH --account='+str(job.account)+'\n'
        c+='#SBATCH -N '+str(job.nodes)+'\n'
        c+='#SBATCH --ntasks-per-node={0}\n'.format(job.processes_per_node)
        c+='#SBATCH --cpus-per-task={0}\n'.format(job.threads)
        c+='#SBATCH -t {0}:{1}:{2}\n'.format(str(job.hours+24*job.days).zfill(2),str(job.minutes).zfill(2),str(job.seconds).zfill(2))
        c+='#SBATCH -o {0}\n'.format(job.outfile)
        c+='#SBATCH -e {0}\n'.format(job.errfile)
        c+='#SBATCH --partition=standard\n'
        c+='#SBATCH --qos={0}\n'.format(job.qos)
        if job.email is not None:
            c+='#SBATCH --mail-user {}\n'.format(job.email)
            c+='#SBATCH --mail-type ALL\n'
            #c+='#SBATCH --mail-type FAIL\n'
        #end if
        c+='\n'
        #c+='cd $SLURM_SUBMIT_DIR\n'
        #c+='\n'
        c+='echo JobID : $SLURM_JOBID\n'
        c+='echo Number of nodes requested: $SLURM_JOB_NUM_NODES\n'
        c+='echo List of nodes assigned to the job: $SLURM_NODELIST\n'
        c+='\n'
        return c
    #end def write_job_header
#end class Archer2

# Unknown
class Tomcat3(Supercomputer):
    name             = 'tomcat3'
    requires_account = False
    batch_capable    = True
    redirect_output  = True

    def write_job_header(self,job):
        if job.queue is None:
            job.queue = 'tomcat'
        #end if
        c = '#!/bin/bash -l\n'
        c+='#SBATCH -J {}\n'.format(job.name)
        c+='#SBATCH -N {}\n'.format(job.nodes)
        c+='#SBATCH -t {}\n'.format(job.sbatch_walltime())
        c+='#SBATCH -p {}\n'.format(job.queue)
        if job.email is not None:
            c+='#SBATCH --mail-user {}\n'.format(job.email)
            c+='#SBATCH --mail-type ALL\n'
        #end if
        c+='#. /home/rcohen/.bashrc\n'
        if len(job.presub)==0:
            c+='unalias cd; source /mnt/beegfs/intel/parallel_studio_xe_2019.3.062/bin/psxevars.sh\n'
            c+='ulimit -a\n'
        #end if
        return c
    #end def write_job_header
#end class Tomcat3


# Active 
# Polaris at ANL
class Polaris(Supercomputer):
    name = 'polaris'
    requires_account = True
    batch_capable    = True
    special_bundling = True

    def post_process_job(self,job):
        if len(job.run_options)==0: 
            opt = obj(
                ppn     = '--ppn {}'.format(job.processes_per_node),
                depth   = '--depth={}'.format(job.threads),
                cpubind = '--cpu-bind depth',
                threads = '--env OMP_NUM_THREADS={}'.format(job.threads),
                )
            job.run_options.add(**opt)
        #end if
    #end def post_process_job

    def write_job_header(self,job):
        if job.queue is None:
            job.queue = 'prod'
        #end if
        if job.filesystems is None:
            job.filesystems = 'home:eagle:grand'
        #end if
        c= '#!/bin/sh\n'
        c+='#PBS -l select={}:system=polaris\n'.format(job.nodes)
        c+='#PBS -l place=scatter\n'
        c+='#PBS -l filesystems={}\n'.format(job.filesystems)
        c+='#PBS -l walltime={}\n'.format(job.pbs_walltime())
        c+='#PBS -A {}\n'.format(job.account)
        c+='#PBS -q {}\n'.format(job.queue)
        c+='#PBS -N {0}\n'.format(job.name)
        c+='#PBS -k doe\n'
        c+='#PBS -o {0}\n'.format(job.outfile)
        c+='#PBS -e {0}\n'.format(job.errfile)
        c+='\n'
        c+='cd ${PBS_O_WORKDIR}\n'

        return c
    #end def write_job_header

    def specialized_bundle_commands(self,job,launcher,serial):
        c = ''
        j0 = job.bundled_jobs[0]
        c+='split --lines={} --numeric-suffixes=1 --suffix-length=3 $PBS_NODEFILE local_hostfile.\n'.format(j0.nodes)
        c+='\n'
        lhfiles = ['local_hostfile.'+str(n+1).zfill(3) for n in range(len(job.bundled_jobs))]
        for j,lh in zip(job.bundled_jobs,lhfiles):
            c+='cp {} {}\n'.format(lh,j.abs_subdir)
        #end for
        for j,lh in zip(job.bundled_jobs,lhfiles):
            j.run_options.add(hostfile='--hostfile '+lh)
            c+='\ncd '+j.abs_subdir+'\n'
            c+=j.run_command(launcher,redirect=True,serial=serial)+'\n'
        #end for
        c+='\nwait\n'
        return c
    #end def specialized_bundle_commands
#end class Polaris

# Active 
# Aurora at ANL
class Aurora(Supercomputer):
    name = 'aurora'
    requires_account = True
    batch_capable    = True
    special_bundling = True

    def pre_process_job(self,job):
        # Set default queue and node type
        if job.queue is None:
            job.queue = 'prod'
        #end if
        if job.filesystems is None:
            job.filesystems = 'flare'
        #end if
        if job.constraint is None:
            job.constraint = 'cpu'
        #end if
        # Account for dual nature of Perlmutter
        if 'cpu' in job.constraint:
            self.procs_per_node = 2
            self.cores_per_node = 102
            self.ram_per_node   = 1024
        elif 'gpu' in job.constraint:
            self.procs_per_node = 1
            self.cores_per_node = 102
            self.ram_per_node   = 768
            self.gpus_per_node  = 6
        else:
            self.error('CPU or GPU constraint must be specified for Aurora')
        #end if
    #end def pre_process_job

    def post_process_job(self,job):
        if len(job.run_options)==0:
            if 'cpu' in job.constraint:
                threads = '--env OMP_NUM_THREADS={} --env OMP_PLACES=cores'.format(job.threads)
                cpubind = '--cpu-bind depth'
            elif 'gpu' in job.constraint:
                threads = '--env OMP_NUM_THREADS={}'.format(job.threads)
                cpubind = '--cpu-bind=list'
                ind = 1
                for _ in range(job.processes_per_node):
                    cpubind += ':{}-{}'.format(ind, ind + job.threads - 1)
                    if ind + 2*(job.threads - 1) > 53 and ind < 53:
                        ind = 53
                    else:
                        ind += job.threads
                # end for

            else:
                self.error('CPU or GPU constraint must be specified for Aurora')
            #end if

            opt = obj(
                ppn     = '--ppn {}'.format(job.processes_per_node),
                depth   = '--depth={}'.format(job.threads),
                cpubind = cpubind,
                threads = threads,
                # affinity= '/soft/tools/mpi_wrapper_utils/gpu_tile_compact.sh',
                )
            job.run_options.add(**opt)
        #end if
    #end def post_process_job

    def write_job_header(self,job):
        c= '#!/bin/sh\n'
        c+='#PBS -l select={}\n'.format(job.nodes)
        c+='#PBS -l place=scatter\n'
        c+='#PBS -l filesystems={}\n'.format(job.filesystems)
        c+='#PBS -l walltime={}\n'.format(job.pbs_walltime())
        c+='#PBS -A {}\n'.format(job.account)
        c+='#PBS -q {}\n'.format(job.queue)
        c+='#PBS -N {0}\n'.format(job.name)
        c+='#PBS -k doe\n'
        c+='#PBS -o {0}\n'.format(job.outfile)
        c+='#PBS -e {0}\n'.format(job.errfile)
        c+='\n'
        c+='cd ${PBS_O_WORKDIR}\n'
        if 'gpu' in job.constraint:
            c+='export MPICH_GPU_SUPPORT_ENABLED=1\n'
        return c
    #end def write_job_header

    def specialized_bundle_commands(self,job,launcher,serial):
        c = ''
        j0 = job.bundled_jobs[0]
        c+='split --lines={} --numeric-suffixes=1 --suffix-length=3 $PBS_NODEFILE local_hostfile.\n'.format(j0.nodes)
        c+='\n'
        lhfiles = ['local_hostfile.'+str(n+1).zfill(3) for n in range(len(job.bundled_jobs))]
        for j,lh in zip(job.bundled_jobs,lhfiles):
            c+='cp {} {}\n'.format(lh,j.abs_subdir)
        #end for
        for j,lh in zip(job.bundled_jobs,lhfiles):
            j.run_options.add(hostfile='--hostfile '+lh)
            c+='\ncd '+j.abs_subdir+'\n'
            c+=j.run_command(launcher,redirect=True,serial=serial)+'\n'
        #end for
        c+='\nwait\n'
        return c
    #end def specialized_bundle_commands
#end class Aurora

# Active 
# Improv at ANL (LCRC)
class Improv(Supercomputer):
    name = 'improv'
    requires_account = True
    batch_capable    = True

    def post_process_job(self,job):
        if len(job.run_options)==0: 
            opt = obj(
                mapby   = '--map-by ppr:{}:package'.format(job.processes_per_proc),
                bindto = '--bind-to socket',
                )
            job.run_options.add(**opt)
        #end if
        if job.threads>1 and job.env is None:
            job.set_environment(
                OMP_PLACES    = 'cores',
                OMP_PROC_BIND = 'close',
                )
        #end if

    #end def post_process_job

    def write_job_header(self,job):
        if job.queue is None:
            job.queue = 'compute'
        #end if
        c= '#!/bin/bash -l\n'
        if job.threads>1:
            c+='#PBS -l select={}:ncpus=128:mpiprocs={}:ompthreads={}\n'.format(job.nodes,job.processes_per_node,job.threads)
        else:
            c+='#PBS -l select={}:ncpus=128:mpiprocs={}\n'.format(job.nodes,job.processes)
        #end if
        c+='#PBS -l walltime={}\n'.format(job.pbs_walltime())
        c+='#PBS -A {}\n'.format(job.account)
        c+='#PBS -q {}\n'.format(job.queue)
        c+='#PBS -N {0}\n'.format(job.name)
        c+='#PBS -k doe\n'
        c+='#PBS -o {0}\n'.format(job.outfile)
        c+='#PBS -e {0}\n'.format(job.errfile)
        c+='\n'
        c+='cd ${PBS_O_WORKDIR}\n'
        return c
    #end def write_job_header
#end class Improv


# Active
# Kagayaki at JAIST
## Added 05/04/2023 by Tom Ichibha
class Kagayaki(Supercomputer):
    name = 'kagayaki'
    requires_account = False
    batch_capable    = True
    special_bundling = False

    def process_job_options(self,job):
        # job.run_options.add(nodefile='-machinefile $PBS_NODEFILE', np='-np '+str(job.processes))
        opt = obj(
            nodefile='-machinefile $PBS_NODEFILE',
            omp='-x OMP_NUM_THREADS',
            np='-np {}'.format(job.processes),
            )
        job.run_options.add(**opt)

    def write_job_header(self,job):
        ppn = 16 if job.queue in ['Default', 'SINGLE', 'LONG', 'DEFAULT'] else 128
        c=''
        c+='#!/bin/bash\n'
        if (job.queue is not None):
            c+='#PBS -q ' + job.queue + '\n'
        c+='#PBS -N ' + job.name + '\n'
        c+='#PBS -o ' + job.outfile +'\n'
        c+='#PBS -e ' + job.errfile + '\n'
        c+='#PBS -l select={0}:ncpus={1}:mpiprocs={1}\n'.format(job.nodes, ppn)  
        c+='cd $PBS_O_WORKDIR\n'
        c+='export OMP_NUM_THREADS=' + str(job.threads) + '\n'
        return c
    #end def write_job_header                                                                       
#end class Kagayaki


# Active
# Kestrel at NREL
class Kestrel(Supercomputer):
    name = 'kestrel'
    requires_account = True
    batch_capable    = True

    def write_job_header(self,job):
        if job.queue is None:
            job.queue = 'regular'
        #end if
        c='#!/bin/bash\n'
        c+='#SBATCH -A '+job.account+'\n'
        c+='#SBATCH -p '+job.queue+'\n'
        c+='#SBATCH -J '+str(job.name)+'\n'
        c+='#SBATCH -t '+job.sbatch_walltime()+'\n'
        c+='#SBATCH -N '+str(job.nodes)+'\n'
        c+='#SBATCH --ntasks-per-node={0}\n'.format(job.processes_per_node)
        c+='#SBATCH --cpus-per-task={0}\n'.format(job.threads)
        c+='#SBATCH -o '+job.outfile+'\n'
        c+='#SBATCH -e '+job.errfile+'\n'
        if job.user_env:
            c+='#SBATCH --export=ALL\n'   # equiv to PBS -V
        else:
            c+='#SBATCH --export=NONE\n'
        #end if
        c+='''
echo $SLURM_SUBMIT_DIR
cd $SLURM_SUBMIT_DIR
'''
        return c
    #end def write_job_header
#end class Kestrel


# Active
# Lassen at LLNL
class Lassen(Supercomputer):

    name = 'lassen'
    requires_account = True
    batch_capable    = True

    def post_process_job(self,job):
        # add the options only if the user has not supplied options
        if len(job.run_options)==0:
            opt = obj(
                nodes = '-N {}'.format(job.nodes),
                tasks = '-T {}'.format(job.processes_per_node),
                )
            if job.threads>1:
                opt.threads = '--threads={}'.format(job.threads)
            #end if
            if job.gpus is None:
                job.gpus = 4# gpus to use per node
            #end if
            #if job.alloc_flags is None:
            #    job.alloc_flags = 'smt1'
            ##end if
            if job.gpus==0:
                None
            else:
                opt.mgpu = '-M "-gpu"'
            #end if
            job.run_options.add(**opt)
        #end if
    #end def post_process_job

    def write_job_header(self,job):
        c ='#!/bin/bash\n'
        c+='#BSUB -G {0}\n'.format(job.account)
        if job.queue is not None:
            c+='#BSUB -q {0}\n'.format(job.queue)
        #end if
        c+='#BSUB -J {0}\n'.format(job.name)
        c+='#BSUB -o {0}\n'.format(job.outfile)
        c+='#BSUB -e {0}\n'.format(job.errfile)
        c+='#BSUB -W {0}\n'.format(job.lsf_walltime())
        c+='#BSUB -nnodes {0}\n'.format(job.nodes)
        #if job.alloc_flags is not None:
        #    c+='#BSUB -alloc_flags "{0}"\n'.format(job.alloc_flags)
        ##end if
        return c
    #end def write_job_header

    def read_process_id(self,output):
        pid = None
        tokens = output.split()
        for t in tokens:
            if t.startswith('<'):
                spid = t.strip('<>').strip()
                if spid.isdigit():
                    pid = int(spid)
                    break
                #end if
            #end if
        #end for
        return pid
    #end def read_process_id
#end class Lassen


# Ruby at LLNL
class Ruby(Supercomputer):
    name = 'ruby'
    requires_account = True
    batch_capable    = True

    def write_job_header(self,job):
        if job.queue is None:
            job.queue = 'regular'
        #end if
        c='#!/bin/bash\n'
        c+='#SBATCH -A '+job.account+'\n'
        c+='#SBATCH -p '+job.queue+'\n'
        c+='#SBATCH -J '+str(job.name)+'\n'
        c+='#SBATCH -t '+job.sbatch_walltime()+'\n'
        c+='#SBATCH -N '+str(job.nodes)+'\n'
        c+='#SBATCH --ntasks-per-node={0}\n'.format(job.processes_per_node)
        c+='#SBATCH --cpus-per-task={0}\n'.format(job.threads)
        c+='#SBATCH -o '+job.outfile+'\n'
        c+='#SBATCH -e '+job.errfile+'\n'
        if job.user_env:
            c+='#SBATCH --export=ALL\n' 
        else:
            c+='#SBATCH --export=NONE\n'
        #end if
        c+='''
echo $SLURM_SUBMIT_DIR
cd $SLURM_SUBMIT_DIR
'''
        return c
    #end def write_job_header
#end class Ruby





#Known machines
#  workstations
for cores in range(1,128+1):
    Workstation('ws'+str(cores),cores,'mpirun'),
#end for
#  supercomputers and clusters
#            nodes sockets cores ram qslots  qlaunch  qsubmit     qstatus    qdelete
Jaguar(      18688,   2,     8,   32,  100,  'aprun',     'qsub',   'qstat',    'qdel')
Kraken(       9408,   2,     6,   16,  100,  'aprun',     'qsub',   'qstat',    'qdel')
Golub(          512,  2,     6,   32, 1000, 'mpirun',     'qsub',   'qstat',    'qdel')
OIC5(           28,   2,    16,  128, 1000, 'mpirun',     'qsub',   'qstat',    'qdel')
Cori(         9688,   1,    68,   96,  100,   'srun',   'sbatch',  'squeue', 'scancel')
BlueWatersXK( 3072,   1,    16,   32,  100,  'aprun',     'qsub',   'qstat',    'qdel')
BlueWatersXE(22640,   2,    16,   64,  100,  'aprun',     'qsub',   'qstat',    'qdel')
Titan(       18688,   1,    16,   32,  100,  'aprun',     'qsub',   'qstat',    'qdel')
EOS(           744,   2,     8,   64, 1000,  'aprun',     'qsub',   'qstat',    'qdel')
Vesta(        2048,   1,    16,   16,   10, 'runjob',     'qsub',  'qstata',    'qdel')
Cetus(        1024,   1,    16,   16,   10, 'runjob',     'qsub',  'qstata',    'qdel')
Mira(        49152,   1,    16,   16,   10, 'runjob',     'qsub',  'qstata',    'qdel')
Cooley(        126,   2,     6,  384,   10, 'mpirun',     'qsub',  'qstata',    'qdel')
Theta(        4392,   1,    64,  192, 1000,  'aprun',     'qsub',  'qstata',    'qdel')
Lonestar(    22656,   2,     6,   12,  128,  'ibrun',     'qsub',   'qstat',    'qdel')
Matisse(        20,   2,     8,   64,    2, 'mpirun',   'sbatch',   'sacct', 'scancel')
Komodo(         24,   2,     6,   48,    2, 'mpirun',   'sbatch',   'sacct', 'scancel')
Amos(         5120,   1,    16,   16,  128,   'srun',   'sbatch',   'sacct', 'scancel')
Chama(        1232,   2,     8,   64, 1000,'mpiexec',   'sbatch',  'squeue', 'scancel')
Uno(           168,   2,     8,  128, 1000,'mpiexec',   'sbatch',  'squeue', 'scancel')
Eclipse(      1488,   2,    18,  128, 1000,'mpiexec',   'sbatch',  'squeue', 'scancel')
Attaway(      1488,   2,    18,  192, 1000,'mpiexec',   'sbatch',  'squeue', 'scancel')
Manzano(      1488,   2,    24,  192, 1000,'mpiexec',   'sbatch',  'squeue', 'scancel')
Ghost(         740,   2,    18,  128, 1000,'mpiexec',   'sbatch',  'squeue', 'scancel')
Amber(        1496,   2,    56,  256, 1000,'mpiexec',   'sbatch',  'squeue', 'scancel')
Skybridge(    1848,   2,     8,   64, 1000,'mpiexec',   'sbatch',  'squeue', 'scancel')
Solo(          374,   2,    18,  128, 1000,'mpiexec',   'sbatch',  'squeue', 'scancel')
SuperMUC(      512,   1,    28,  256,    8,'mpiexec', 'llsubmit',     'llq','llcancel')
Stampede2(    4200,   1,    68,   96,   50,  'ibrun',   'sbatch',  'squeue', 'scancel')
CadesMoab(     156,   2,    18,  128,  100, 'mpirun',     'qsub',   'qstat',    'qdel')
CadesSlurm(    156,   2,    18,  128,  100, 'mpirun',   'sbatch',  'squeue', 'scancel')
Summit(       4608,   2,    21,  512,  100,  'jsrun',     'bsub',   'bjobs',   'bkill')
Rhea(          512,   2,     8,  128, 1000,   'srun',   'sbatch',  'squeue', 'scancel')
Andes(         704,   2,    16,  256, 1000,   'srun',   'sbatch',  'squeue', 'scancel')
Tomcat3(         8,   1,    64,  192, 1000, 'mpirun',   'sbatch',   'sacct', 'scancel')
SuperMUC_NG(  6336,   1,    48,   96, 1000,'mpiexec',   'sbatch',   'sacct', 'scancel')
Archer2(      5860,   2,    64,  512, 1000,   'srun',   'sbatch',  'squeue', 'scancel')
Polaris(       560,   1,    32,  512,    8,'mpiexec',     'qsub',   'qstat',    'qdel')
Kagayaki(      240,   2,    64,  512,   20, 'mpirun',     'qsub',   'qstat',    'qdel')
Perlmutter(   3072,   2,   128,  512, 5000,   'srun',   'sbatch',  'squeue', 'scancel')
Improv(        825,   2,    64,  256,  100, 'mpirun',     'qsub',   'qstat',    'qdel')
Lassen(        756,   2,    21,  512,  100,   'lrun',     'bsub',   'bjobs',   'bkill')
Ruby(         1480,   2,    28,  192,  100,   'srun',   'sbatch',  'squeue', 'scancel')
Kestrel(      2144,   2,    52,  256,  100,   'srun',   'sbatch',  'squeue', 'scancel')
Inti(           13,   2,    64,  256,  100,   'srun',   'sbatch',  'squeue', 'scancel')
Baseline(      128,   2,    64,  512,  100,   'srun',   'sbatch',  'squeue', 'scancel')
Besms(         166,   1,    96,  768, 1000,   'srun',   'sbatch',  'squeue', 'scancel')
Frontier(     9856,   4,    14,   64, 1000,   'srun',   'sbatch',  'squeue', 'scancel')
Aurora(      10624,   2,   102,  512, 1000,'mpiexec',     'qsub',   'qstat',    'qdel')

#machine accessor functions
get_machine_name = Machine.get_hostname
get_machine      = Machine.get

#rename Job with lowercase
job=Job






