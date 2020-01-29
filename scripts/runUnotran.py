from getGeometry import setGeometry
from getGroupStructure import getGroupBounds
import sys
sys.path.append('/homes/rlreed/workspace/unotran/src')
import os
import pydgm
import numpy as np
np.set_printoptions(precision=4, linewidth=132)
import matplotlib.pyplot as plt
import time
import resource


def setSolver(solverOption):
    pydgm.control.solver_type = solverOption.ljust(256)
    pydgm.control.allow_fission = (solverOption == 'eigen')
    pydgm.control.source_value = 0.0 if solverOption == 'eigen' else 1.0


def setGroup(G, contig, order=None, basisType='dlp'):
    groupBounds, counts, mask, bName = getGroupBounds(G, contig)

    pydgm.control.xs_name = 'XS/{}gXS.anlxs'.format(G).ljust(256)

    print 'using basis basis/{}_{}_{}g'.format(basisType, bName, G)
    pydgm.control.dgm_basis_name = 'basis/{}_{}_{}g'.format(basisType, bName, G).ljust(256)

    if order is not None:
        print basisType, G, counts.values(), order, [min(d - 1, order) for d in counts.values()]
        pydgm.control.truncation_map = [min(d - 1, order) for d in counts.values()]

    if len(groupBounds) > 0:
        pydgm.control.energy_group_map = np.array(groupBounds) + 1


def setVariables(G, geoOption, solverOption, DGM, contig, basisType='dlp', truncationOrder=None, homogOption=None):
    '''
    Set the input options to run Unotran

    homogOption:
        0 - No Homogenization
        1 - Flat Delta approximation
        2 - Linear Delta approximation
        3 - Spatial homogenization
        4 - Both linear delta and spatial homogenization
    '''
    # Set the variables

    fm, cm, mm, bounds = setGeometry(geoOption)
    setSolver(solverOption)
    setGroup(G, contig, truncationOrder, basisType)
    pydgm.control.fine_mesh = fm
    pydgm.control.coarse_mesh = cm
    pydgm.control.material_map = mm

    pydgm.control.boundary_type = bounds
    pydgm.control.angle_order = 8
    pydgm.control.angle_option = pydgm.angle.gl
    pydgm.control.recon_print = 1
    pydgm.control.eigen_print = 1
    pydgm.control.outer_print = 0
    pydgm.control.recon_tolerance = 1e-12
    pydgm.control.eigen_tolerance = 1e-12
    pydgm.control.outer_tolerance = 1e-12
    pydgm.control.max_recon_iters = 100000
    pydgm.control.max_eigen_iters = 1000
    pydgm.control.max_outer_iters = 100
    pydgm.control.lamb = 0.8 if G == 44 else 0.8 if G == 238 else 0.4
    pydgm.control.use_dgm = DGM
    pydgm.control.store_psi = True
    pydgm.control.equation_type = 'DD'
    pydgm.control.scatter_legendre_order = 0
    pydgm.control.ignore_warnings = True

    if homogOption == 0 or homogOption is None:
        pydgm.control.truncate_delta = False
    elif homogOption == 1:
        pydgm.control.truncate_delta = True
        pydgm.control.delta_legendre_order = 0
    elif homogOption == 2:
        pydgm.control.truncate_delta = True
        pydgm.control.delta_legendre_order = 1
    elif homogOption == 3:
        pydgm.control.truncate_delta = False
        assert geoOption == 'full'
        pydgm.control.homogenization_map = [i + 1 for i, f in enumerate(fm) for j in range(f)]
    elif homogOption == 4:
        pydgm.control.truncate_delta = True
        pydgm.control.delta_legendre_order = 1
        assert geoOption == 'full'
        pydgm.control.homogenization_map = [i + 1 for i, f in enumerate(fm) for j in range(f)]

def solve(phi=None, k=None):
    # Initialize and solve
    if pydgm.control.use_dgm:
        pydgm.dgmsolver.initialize_dgmsolver()

        if phi is not None:
            pydgm.state.phi = phi
            for a in range(pydgm.control.number_angles * 2):
                pydgm.state.psi[:, a, :] = phi[0]
        if k is not None:
            pydgm.state.k_eff = k
        pydgm.dgmsolver.dgmsolve()
    else:
        pydgm.solver.initialize_solver()
        pydgm.solver.solve()

    # Output the scalar flux
    phi = np.copy(pydgm.state.phi)
    k = np.copy(pydgm.state.keff)
    dens = np.copy(pydgm.state.mg_density)
    count = np.copy(pydgm.state.sweep_count)

    # Cleanup
    if pydgm.control.use_dgm:
        pydgm.dgmsolver.finalize_dgmsolver()
    else:
        pydgm.solver.finalize_solver()

    return phi, k, dens, count


def run(task):
    solver = 'eigen'
    geo, G, basisType, order, contig, homog = getInfo(task)

    path = 'data'

    print 'Task={} Group={} Order={} Basis={} Contig={} Geo={} homog={}'.format(task + 1, G, order, basisType, contig, geo, homog)

    # Solve the DGM Problem
    setVariables(G, geo, solver, True, contig, basisType, order, homog)
    t_start = time.time()
    dgm_phi, dgm_k, dgm_dens, dgm_count = solve()
    t_end = time.time()

    # Save the flux arrays
    np.save('{}/dgm_phi_{}_{}_g{}_c{}_o{}_h{}'.format(path, basisType, geo, G, contig, order, homog), dgm_phi)
    np.save('{}/dgm_k_{}_{}_g{}_c{}_o{}_h{}'.format(path, basisType, geo, G, contig, order, homog), dgm_k)
    np.save('{}/dgm_dens_{}_{}_g{}_c{}_o{}_h{}'.format(path, basisType, geo, G, contig, order, homog), dgm_dens)
    np.save('{}/dgm_count_{}_{}_g{}_c{}_o{}_h{}'.format(path, basisType, geo, G, contig, order, homog), dgm_count)
    np.save('{}/dgm_time_{}_{}_g{}_c{}_o{}_h{}'.format(path, basisType, geo, G, contig, order, homog), t_end - t_start)
    np.save('{}/dgm_mem_{}_{}_g{}_c{}_o{}_h{}'.format(path, basisType, geo, G, contig, order, homog), resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)

def getInfo(task):
    Gs = [44, 238, 1968]
    geos = ['full']
    contigs = [3]

    item = 0
    for i, G in enumerate(Gs):
        for homogOption in [0, 1, 2, 3, 4]:
            for contig in contigs:
                maxOrder = max(getGroupBounds(G, contig)[1].values())
                for geo in geos:
                    # Select the basis based on the geometry
                    if geo == 0:
                        basisTypes = ['dlp']
                    elif geo == 'full':
                        basisTypes = ['dlp', 'klt_full', 'klt_combine', 'klt_uo2', 'klt_mox', 'klt_pins_full']
                    elif geo == 'core1':
                        basisTypes = ['dlp', 'klt_assay_12', 'klt_assay_all', 'klt_core1', 'klt_pins_core1']
                    elif geo == 'core2':
                        basisTypes = ['dlp', 'klt_assay_13', 'klt_assay_all', 'klt_core2', 'klt_pins_core2']

                    for j, basisType in enumerate(basisTypes):
                        for order in range(maxOrder):
                            if item == task:
                                return geo, G, basisType, order, contig, homogOption
                            else:
                                item += 1

task = os.environ['SLURM_ARRAY_TASK_ID']
task = int(task) - 1

run(task)

