from getGeometry import setGeometry
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


def setGroup(G):
    pydgm.control.xs_name = 'XS/{}gXS.anlxs'.format(G).ljust(256)


def setVariables(G, geoOption, solverOption):
    # Set the variables

    fm, cm, mm, bounds = setGeometry(geoOption)
    setSolver(solverOption)
    setGroup(G)
    pydgm.control.fine_mesh = fm
    pydgm.control.coarse_mesh = cm
    pydgm.control.material_map = mm
    pydgm.control.boundary_type = bounds
    pydgm.control.angle_order = 8
    pydgm.control.angle_option = pydgm.angle.gl
    pydgm.control.eigen_print = 1
    pydgm.control.outer_print = 1 if G > 100 else 0
    pydgm.control.eigen_tolerance = 1e-12
    pydgm.control.outer_tolerance = 1e-12
    pydgm.control.max_eigen_iters = 100000
    pydgm.control.max_outer_iters = 1000
    pydgm.control.store_psi = False
    pydgm.control.equation_type = 'DD'
    pydgm.control.scatter_legendre_order = 0
    pydgm.control.ignore_warnings = True


def solve(phi=None):
    # Initialize and solve
    pydgm.solver.initialize_solver()
    pydgm.solver.solve()

    # Output the scalar flux
    phi = np.copy(pydgm.state.phi)
    k = np.copy(pydgm.state.keff)
    dens = np.copy(pydgm.state.mg_density)
    count = np.copy(pydgm.state.sweep_count)

    # Cleanup
    pydgm.solver.finalize_solver()

    return phi, k, dens, count


def getInfo(task):
    Gs = [44, 238, 1968]
    item = 0
    for i, G in enumerate(Gs):
        for geo in ['uo2-1', 'uo2-2', 'uo2-Gd', 'mox', 'junction', 'full', 'assay1', 'assay2', 'assay3', 'core1', 'core2']:
            if item == task:
                return G, geo
            else:
                item += 1


task = os.environ['SLURM_ARRAY_TASK_ID']
task = int(task) - 1
G, geo = getInfo(task)

solver = 'eigen'
path = 'reference'
print 'Group={}'.format(G)
# Solve the NON-DGM Problem
setVariables(G, geo, solver)
t_start = time.time()
phi, k, dens, count = solve()
t_end = time.time()
# Save the reference flux
np.save('{}/ref_phi_{}_{}'.format(path, geo, G), phi)
np.save('{}/ref_k_{}_{}'.format(path, geo, G), k)
np.save('{}/ref_dens_{}_{}'.format(path, geo, G), dens)
np.save('{}/ref_count_{}_{}'.format(path, geo, G), count)
np.save('{}/ref_time_{}_{}'.format(path, geo, G), t_end - t_start)
np.save('{}/ref_mem_{}_{}'.format(path, geo, G), resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)

