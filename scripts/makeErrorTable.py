import numpy as np
from getGroupStructure import getGroupBounds
from getGeometry import setGeometry
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
plt.ioff()
from matplotlib import rc
rc('font', **{'family': 'serif'})
from matplotlib import rcParams
rcParams['xtick.direction'] = 'out'
rcParams['ytick.direction'] = 'out'
rcParams['xtick.labelsize'] = 18
rcParams['ytick.labelsize'] = 18
rcParams['lines.linewidth'] = 1.85
rcParams['axes.labelsize'] = 20
rcParams['legend.numpoints'] = 1
rcParams.update({'figure.autolayout': True})

def readData(G):
    '''
    Reads data from the data directory that has been make via Unotran

    Returns the % relative error for both k and pin-cell fission density
    '''
    groupMask, counts, mask, extraName = getGroupBounds(G, 3)
    maxOrder = max(getGroupBounds(G, 3)[1].values())
    data = {}
    for d in ['dens', 'k']:
        # Load the reference
        try:
            ref = np.load('reference/ref_{}_{}_{}.npy'.format(d, 'full', G))
            if d == 'dens':
                # Normalize the reference
                ref /= np.linalg.norm(ref)
                # Get the density for each pin
                ref = np.sum(ref.reshape(10,-1), axis=1)
                ref /= np.mean(ref)
        except IOError:
            raise IOError('Missing reference for {} with {} groups'.format(d, G))
        data[d] = {}
        # Add the reference to the data
        data[d]['ref'] = ref
        for basis in ['dlp', 'klt_full', 'klt_combine', 'klt_uo2', 'klt_mox', 'klt_pins_full']:
            data[d][basis] = {}
            for homog in range(5):
                data[d][basis][homog] = []
                for o in range(maxOrder):
                    # Load the data
                    try:
                        sol = np.load('data/dgm_{}_{}_{}_g{}_c{}_o{}_h{}.npy'.format(d, basis, 'full', G, 3, o, homog))
                        if d == 'dens':
                            sol /= np.linalg.norm(sol)
                            sol = np.sum(sol.reshape(10, -1), axis=1)
                            sol /= np.mean(sol)
                    except IOError as e:
                        print('Missing data/dgm_{}_{}_{}_g{}_c{}_o{}_h{}.npy'.format(d, basis, 'full', G, 3, o, homog))
                        sol = np.zeros(ref.shape)
                    # Compute the error
                    err = (sol - ref) / ref * 100

                    data[d][basis][homog].append(err)

    return data

def makeTable(G, trunc=False):
    '''

    '''

    data = readData(G)
    groupMask, counts, mask, extraName = getGroupBounds(G, 3)
    maxOrder = max(getGroupBounds(G, 3)[1].values())
    order = 2

    table = np.zeros((11, 11))

    # Load the reference
    table[0,0] = data['k']['ref']
    table[1:,0] = data['dens']['ref']
    for homog in range(5):
        # Get the non-energy truncated results
        table[0, homog + 1] = data['k']['klt_combine'][homog][maxOrder-1]
        table[1:, homog + 1] = data['dens']['klt_combine'][homog][maxOrder-1]
        # Get the energy truncated results
        table[0, homog + 6] = data['k']['klt_combine'][homog][order]
        table[1:, homog + 6] = data['dens']['klt_combine'][homog][order]

    table = table[:,[0, 6, 7, 8, 9, 10]] if trunc else table[:,[0, 2, 3, 4, 5]]
    label = ['$k_{\\text{eff}}$', 'Cell {}']
    print '\nData table for Truncation={}'.format(trunc)
    for i, row in enumerate(table):
        head = label[0] if i == 0 else label[1].format(i)
        print head + '&' + ' & '.join(['{: 4.3f}'.format(r) for r in row]) + '\\\\'
        if i == 0:
            print '\\hline'

def makePlots(G, pType='phi'):
    '''

    '''

    groupMask, counts, mask, extraName = getGroupBounds(G, 3)
    maxOrder = max(getGroupBounds(G, 3)[1].values())
    lowOrder = 2
    fm, cm, mm, bounds = setGeometry('full')

    dx = []
    for i, c in enumerate(cm):
        if i == len(cm) - 1:
            break
        for x in np.linspace(c, cm[i+1], fm[i]+1)[:-1]:
            dx.append(x)
    dx.append(cm[-1])
    dx = np.array(dx)
    x = 0.5 * (dx[1:] + dx[:-1])

    labels = ['Case 1', 'Case 2', 'Case 3', 'Case 4']
    ls = [':', '-', '-.', '--']
    for order in [lowOrder, maxOrder - 1]:
        for group in (range(G) if pType == 'phi' else [0]):
            ref = np.load('reference/ref_{}_{}_{}.npy'.format(pType, 'full', G))
            if pType == 'phi':
                ref = ref[0,group]
            ref /= np.linalg.norm(ref)
            plt.plot(x, ref, label='Reference')

            for homog in range(1, 5):
                sol = np.load('data/dgm_{}_{}_{}_g{}_c{}_o{}_h{}.npy'.format(pType, 'klt_combine', 'full', G, 3, order, homog))
                if pType == 'phi':
                    sol = sol[0, group]
                sol /= np.linalg.norm(sol)
                plt.plot(x, sol, ls=ls[homog-1], label=labels[homog-1])
            plt.legend(fancybox=True, framealpha=0.0)
            plt.xlabel('length [cm]')
            plt.ylabel('normalized scalar flux' if pType == 'phi' else 'normalized fission density')
            plt.xlim([0.0, 12.8])
            plt.savefig('plots/{}_g{}_o{}.png'.format(pType, group, order), transparent=True)
            plt.clf()

if __name__ == '__main__':
    makePlots(238, 'phi')
    makePlots(238, 'dens')
    makeTable(238, False)
    makeTable(238, True)

