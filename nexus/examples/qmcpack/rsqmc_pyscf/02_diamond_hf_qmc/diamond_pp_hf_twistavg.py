#! /usr/bin/env python3

from nexus import settings,job,run_project,obj
from nexus import generate_physical_system
from nexus import generate_pyscf
from nexus import generate_convert4qmc
from nexus import generate_qmcpack

settings(
    pseudo_dir = '../../pseudopotentials',
    results    = '',
    sleep      = 3,
    machine    = 'ws16',
    )

system = generate_physical_system(
    units    = 'A',
    axes     = '''1.785   1.785   0.000
                  0.000   1.785   1.785
                  1.785   0.000   1.785''',
    elem_pos = '''
               C  0.0000  0.0000  0.0000
               C  0.8925  0.8925  0.8925
               ''',
    tiling   = (2,1,1),
    kgrid    = (3,1,1),
    kshift   = (0,0,0),
    C        = 4,
    )

scf = generate_pyscf(
    identifier = 'scf',                      # log output goes to scf.out
    path       = 'diamond_ta/scf',              # directory to run in
    job        = job(serial=True,threads=16),# pyscf must run w/o mpi
    template   = './scf_template.py',        # pyscf template file
    system     = system,
    cell       = obj(                        # used to make Cell() inputs
        basis         = 'bfd-vdz',
        ecp           = 'bfd',
        drop_exponent = 0.1,
        verbose       = 5,
        ),
    save_qmc   = True,                # save wfn data for qmcpack
    )

c4q = generate_convert4qmc(
    identifier   = 'c4q',
    path         = 'diamond_ta/scf',
    job          = job(cores=1),
    dependencies = (scf,'orbitals'),
    )

qmc = generate_qmcpack(
    driver       = 'legacy',
    block        = True,
    identifier   = 'vmc',
    path         = 'diamond_ta/vmc',
    job          = job(cores=12,threads=4,app='qmcpack_complex'),
    input_type   = 'basic',
    system       = system,
    pseudos      = ['C.BFD.xml'],
    corrections  = [],
    qmc          = 'vmc',
    dependencies = [(c4q,'determinantset'),
                    (scf,'orbitals')],
    )

run_project()
