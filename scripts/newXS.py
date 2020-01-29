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


def setVariables(G, geoOption, solverOption, DGM, contig, basisType='dlp', truncationOrder=None):
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
    pydgm.control.eigen_print = 1 if not DGM else 1 if G != 1968 else 1
    pydgm.control.outer_print = 0
    pydgm.control.recon_tolerance = 1e-12
    pydgm.control.eigen_tolerance = 1e-12
    pydgm.control.outer_tolerance = 1e-12
    pydgm.control.max_recon_iters = 100000
    pydgm.control.max_eigen_iters = 1000
    pydgm.control.max_outer_iters = 100
    pydgm.control.lamb = 1.0 if G == 44 else 0.8 if G == 238 else 0.4
    pydgm.control.use_dgm = DGM
    pydgm.control.store_psi = True
    pydgm.control.equation_type = 'DD'
    pydgm.control.scatter_legendre_order = 0
    pydgm.control.ignore_warnings = True



def run():
    solver = 'eigen'
    geo = 'full'
    G = 44
    basisType = 'dlp'
    contig = 3
    order = 12

    path = 'data'

    # Solve the DGM Problem
    setVariables(G, geo, solver, True, contig, basisType, order)

    # Setup the variable containers
    pydgm.dgmsolver.initialize_dgmsolver()

    # Give some dummy phi
    pydgm.state.phi = np.arange(np.product(pydgm.state.phi.shape)).reshape(pydgm.state.phi.shape) + 1

    # Get the right shapes from the inputs
    nL, nG, nC = pydgm.state.phi.shape

    # Compute the flux momenst
    pydgm.dgmsolver.compute_flux_moments()

    # Get the flux moments and shape
    m_phi = pydgm.state.mg_phi[0].T
    nCG = m_phi.shape[1]

    # Compute the test cross sections
    pydgm.dgmsolver.compute_xs_moments()
    sig_t_test = pydgm.state.mg_sig_t

    # Get the group bounds
    groupBounds, counts, mask, bName = getGroupBounds(G, contig)
    nO = np.max(counts.values())
    bounds = [0]
    for o in range(np.max(counts.keys()) + 1):
        bounds.append(bounds[o] + counts[o])

    # Load the basis
    P = np.loadtxt('basis/{}_{}_{}g'.format(basisType, bName, G))

    # Get the fine group flux
    phi = pydgm.state.phi[0,:,:]

    nMat = np.max(pydgm.mesh.mmap)

    # Create the XS containers
    ex_sig_t = np.zeros((nMat, nCG, nO))
    m_sig_t = np.zeros((nC, nCG, nO))

    # Loop over the spatial cells
    for m in range(nMat):
        # Get the fine XS
        sig_t = pydgm.material.sig_t[:, m]

        # Loop over the energy group
        for cg in range(nCG):
            P_G = P.T[bounds[cg]: bounds[cg + 1]]
            for gp in range(bounds[cg], bounds[cg + 1]):
                for o, p in enumerate(P_G):
                    ex_sig_t[m, cg, o] += p[gp] * sig_t[gp] * P_G[0, gp]


    print ex_sig_t.shape

    print m_phi.shape

    for c in range(nC):
        pass

    exit()


    m_sig_t[j] = m_sig_t.dot(m_phi) / m_phi[0]


    print m_phi
    print m_sig_t.T
    print sig_t_test

run()

