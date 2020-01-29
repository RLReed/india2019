from getGroupStructure import getGroupBounds
from collections import OrderedDict
import sys
import os
sys.path.append('/homes/rlreed/workspace/unotran/src')
import pydgm
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font', **{'family': 'serif'})
from matplotlib import rcParams
rcParams['xtick.direction'] = 'out'
rcParams['ytick.direction'] = 'out'
rcParams['xtick.labelsize'] = 18
rcParams['ytick.labelsize'] = 18
rcParams['lines.linewidth'] = 1.85
rcParams['axes.labelsize'] = 20
rcParams.update({'figure.autolayout': True})
np.set_printoptions(linewidth=132)


def getEnergy(G):
    with open('XS/{0}gXS.anlxs'.format(G), 'r') as f:
        f.readline()
        s = f.readline()
    return np.array(s.split()).astype(float) * 1e6


def makePlot(G):
    E = getEnergy(G)
    meanE = 0.5 * (E[1:] + E[:-1])
    diffE = E[:-1] - E[1:]
    snapshots = getSnapshots(G, 'full')
    snapMean = np.mean(snapshots, axis=1) * meanE / diffE
    norm = np.linalg.norm(snapMean)
    snapMean /= norm

    EE, minsnap = barchart(E, np.min(snapshots, axis=1) / norm / diffE * meanE)
    EE, maxsnap = barchart(E, np.max(snapshots, axis=1) / norm / diffE * meanE)

    full = KLT(snapshots)[:, 0] * meanE / diffE
    uo2 = KLT(getSnapshots(G, 'uo2-1'))[:, 0] * meanE / diffE
    mox = KLT(getSnapshots(G, 'mox'))[:, 0] * meanE / diffE
    com = KLT(np.concatenate((getSnapshots(G, 'uo2-1'), getSnapshots(G, 'mox'), getSnapshots(G, 'junction')), axis=1))[:, 0] * meanE / diffE
    pin = KLT(np.concatenate((getSnapshots(G, 'uo2-1'), getSnapshots(G, 'mox')), axis=1))[:, 0] * meanE / diffE

    norm = np.mean(snapMean / full)

    plt.plot(*barchart(E, full * norm), label='full', c='r', ls='-')
    plt.plot(*barchart(E, uo2 * norm), label='UO$_2$', c='m', ls='-')
    plt.plot(*barchart(E, mox * norm), label='MOX', c='c', ls='-')
    plt.plot(*barchart(E, com * norm), label='combine', c='b', ls='-')
    plt.plot(*barchart(E, com * norm), label='pin', c='g', ls='--')
    plt.fill_between(EE, minsnap, maxsnap, alpha=0.5)
    plt.xlim([min(E), max(E)])
    plt.ylim([1e-5, 1e0])
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Energy [eV]')
    plt.ylabel('Flux per unit lethargy')
    plt.legend(loc=0)
    plt.savefig('plots/{}_snapshots.png'.format(G))
    plt.clf()


def getSnapshots(G, geoOption):
    return np.load('reference/ref_phi_{}_{}.npy'.format(geoOption, G))[0]


def barchart(x, y):
    X = np.zeros(2 * len(y))
    Y = np.zeros(2 * len(y))
    for i in range(len(y)):
        X[2 * i] = x[i]
        X[2 * i + 1] = x[i + 1]
        Y[2 * i] = y[i]
        Y[2 * i + 1] = y[i]
    return X, Y


def modGramSchmidt(A):
    m, n = A.shape
    A = A.copy()
    Q = np.zeros((m, n))
    R = np.zeros((n, n))

    for k in range(n):
        R[k, k] = np.linalg.norm(A[:, k:k + 1].reshape(-1), 2)
        Q[:, k:k + 1] = A[:, k:k + 1] / R[k, k]
        R[k:k + 1, k + 1:n + 1] = np.dot(Q[:, k:k + 1].T, A[:, k + 1:n + 1])
        A[:, k + 1:n + 1] = A[:, k + 1:n + 1] - np.dot(Q[:, k:k + 1], R[k:k + 1, k + 1:n + 1])

    return Q, R


def plotBasis(G, basisType, contig, bName):
    # load the basis from file
    basis = np.loadtxt(bName)

    vectors = np.zeros((3, G))
    for g in range(G):
        b = np.trim_zeros(basis[g], trim='f')
        if len(b) >= 3:
            b = b[:3]
        else:
            bb = np.zeros(3)
            bb[:b.shape[0]] = b
            b = bb
        vectors[:, g] = b
    plot(vectors, basisType, contig, bName)


def plot(A, basisType, contig, bName):
    colors = ['b', 'g', 'm']
    plt.clf()
    G = A.shape[1]

    groupMask, counts, mask, extraName = getGroupBounds(G, contig)

    bounds = [0]
    for i, a in enumerate(A):
        ming = 0
        for CG, nGroups in counts.items():
            maxg = ming + nGroups
            plt.plot(range(ming, maxg), a[ming:maxg], c=colors[i], label='order {}'.format(i))
            bounds.append(maxg)
            ming += nGroups

    bounds = np.array(bounds)
    plt.vlines(bounds[1:-1] - 0.5, -1, 1)
    plt.xlim([0, G - 1])
    plt.ylim([-1, 1])
    plt.xlabel('Energy group')
    plt.ylabel('Normalized basis')
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = OrderedDict(zip(labels, handles))
    legend = plt.legend(by_label.values(), by_label.keys(), loc='upper center', ncol=3, fancybox=True, framealpha=0.0, bbox_to_anchor=(0.5, 1.1))
    plt.savefig('{}.png'.format(bName), transparent=True, additional_artists=[legend])
    return


def KLT(Data, flat=False):
    A = Data.T.dot(Data)
    # Perform the SVD
    w, A = np.linalg.eig(A)
    idx = w.argsort()[::-1]  # sorts eigenvalues and vectors max->min
    w = w[idx]
    A = A[:, idx]
    A = np.real(Data.dot(A))
    # Orthogonalize the matrix
    A, r = np.linalg.qr(A, mode='complete')

    if flat:
        # ignore the smallest eigenpair to make room for the flat zeroth
        A = np.insert(A.T[:-1], 0, np.ones(len(A[0])) / np.sqrt(len(A[0])), axis=0)  # Add a flat zeroth

        # Redo the SVD if the zeroth was changed
        A, r = np.linalg.qr(A.T, mode='complete')

    if np.sum(A[:, 0]) < 0:
        A *= -1

    # Make sure that the basis is orthogonal
    np.testing.assert_array_almost_equal(A.T.dot(A), np.eye(len(A)), 12)

    return A


def makeBasis(G, basisType, contig, flat=False):
    if basisType == 'full':
        name = 'klt_full'
        data = getSnapshots(G, 'full')
    elif basisType == 'mox':
        name = 'klt_mox'
        data = getSnapshots(G, 'mox')
    elif basisType == 'uo2':
        name = 'klt_uo2'
        data = getSnapshots(G, 'uo2-1')
    elif basisType == 'pins_full':
        name = 'klt_pins_full'
        data = np.concatenate((getSnapshots(G, 'uo2-1'), getSnapshots(G, 'mox')), axis=1)
    elif basisType == 'pins_core1':
        name = 'klt_pins_core1'
        data = np.concatenate((getSnapshots(G, 'uo2-1'), getSnapshots(G, 'uo2-2'), getSnapshots(G, 'uo2-Gd')), axis=1)
    elif basisType == 'pins_core2':
        name = 'klt_pins_core2'
        data = np.concatenate((getSnapshots(G, 'uo2-1'), getSnapshots(G, 'uo2-2'), getSnapshots(G, 'mox')), axis=1)
    elif basisType == 'combine':
        name = 'klt_combine'
        data = np.concatenate((getSnapshots(G, 'uo2-1'), getSnapshots(G, 'mox'), getSnapshots(G, 'junction')), axis=1)
    elif basisType == 'core1':
        name = 'klt_core1'
        data = getSnapshots(G, 'core1')
    elif basisType == 'core2':
        name = 'klt_core2'
        data = getSnapshots(G, 'core2')
    elif basisType == 'assay_12':
        name = 'klt_assay_12'
        data = np.concatenate((getSnapshots(G, 'assay1'), getSnapshots(G, 'assay2')), axis=1)
    elif basisType == 'assay_13':
        name = 'klt_assay_13'
        data = np.concatenate((getSnapshots(G, 'assay1'), getSnapshots(G, 'assay3')), axis=1)
    elif basisType == 'assay_all':
        name = 'klt_assay_all'
        data = np.concatenate((getSnapshots(G, 'assay1'), getSnapshots(G, 'assay2'), getSnapshots(G, 'assay3')), axis=1)
    else:
        raise NotImplementedError('The type: {} has not been implemented'.format(basisType))

    # Get the coarse group bounds
    groupMask, counts, mask, extraName = getGroupBounds(G, contig)
    eName = ('_' + extraName) if extraName is not None else ''

    # Initialize the basis lists
    basis = np.zeros((G, G))

    # Compute the basis for each coarse group
    for group, order in counts.items():
        # Get the mask for the currect group
        m = groupMask == group

        # Get the DLP basis for the given order
        A = KLT(data[m], flat)

        # Slice into the basis with the current group
        basis[np.ix_(m, m)] = A

    bName = 'basis/{0}{2}_{1}g'.format(name, G, eName)
    # Save the basis to file
    np.savetxt(bName, basis)

    plotBasis(G, basisType, contig, bName)


def getInfo(task):
    Gs = [44, 238, 1968]
    item = 0
    for G in Gs:
        makeSnapPlot = True
        for contig in [3]:
            # makePlot(G)
            for basisType in ['uo2', 'mox', 'pins_full', 'combine', 'full', 'pins_core1', 'pins_core2', 'assay_12', 'assay_13', 'assay_all', 'core1', 'core2']:
                if makeSnapPlot:
                    makePlot(G)
                    makeSnapePlot = False

                if item == task:
                    x = getGroupBounds(G, contig)[1]
                    print x, np.cumsum(x.values())
                    return G, basisType, contig
                else:
                    item += 1

if __name__ == '__main__':
    task = os.environ['SLURM_ARRAY_TASK_ID']
    task = int(task) - 1

    G, basisType, contig = getInfo(task)
    print 'G={}, contig={}, type={}'.format(G, contig, basisType)
    makeBasis(G, basisType, contig, flat=True)
    print 'complete'

