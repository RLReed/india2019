from getGroupStructure import getGroupBounds
from collections import OrderedDict
import sys
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


def barchart(x, y):
    X = np.zeros(2 * len(y))
    Y = np.zeros(2 * len(y))
    for i in range(0, len(y)):
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


def plotBasis(G, bName):

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
    plot(vectors, bName)


def plot(A, bName):
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
    plt.legend(by_label.values(), by_label.keys(), loc='lower center', ncol=3)
    plt.savefig('{}.png'.format(bName))
    plt.clf()
    return


def DLP(size):
    order = size
    A = np.ones((size, order))
    if order > 1:
        for j in range(size):
            A[j, 1] = (size - 1 - (2.0 * j)) / (size - 1)
        for i in range(2, order):
            for j in range(size):
                c0 = (i - 1) * (size - 1 + i)
                c1 = (2 * i - 1) * (size - 1 - 2 * j)
                c2 = i * (size - i)
                A[j, i] = (c1 * A[j, i - 1] - c0 * A[j, i - 2]) / c2
    return modGramSchmidt(A)[0]

def testBasis(basis):
    np.testing.assert_array_almost_equal(basis.T.dot(basis), np.eye(len(basis)), 10)

def makeBasis(G, contig):
    # Get the coarse group bounds
    groupMask, counts, mask, extraName = getGroupBounds(G, contig)
    eName = ('_' + extraName) if extraName is not None else ''

    print 'building {0}{2}_{1}g'.format('dlp', G, eName)
    print counts

    # Initialize the basis lists
    basis = np.zeros((G, G))

    # Compute the basis for each coarse group
    for group, order in counts.items():
        # Get the DLP basis for the given order
        A = DLP(order)

        # Get the mask for the currect group
        m = groupMask == group

        #plt.plot(A[np.argsort(mask[m]), :3])
        # plt.show()

        # Slice into the basis with the current group
        basis[np.ix_(m, m)] = A[np.argsort(mask[m])]

    testBasis(basis)


    # Save the basis to file
    bName = 'basis/{0}{2}_{1}g'.format('dlp', G, eName)
    np.savetxt(bName, basis)

    plotBasis(G, bName)

if __name__ == '__main__':
    for G in [44, 238, 1968]:
        for contig in [3]:
            makeBasis(G, contig)
            # plotBasis(G)
