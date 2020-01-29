from getGroupStructure import getGroupBounds
import numpy as np
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

# Select the geometry
geo = 1
if geo == 0:  # Inf medium
    basisTypes = ['dlp']
    figureNames = {key: val for key, val in zip(basisTypes, ['DLP'])}
    styles = ['ko:']
elif geo == 1:  # 10-pin
    basisTypes = ['dlp', 'klt_uo2', 'klt_mox', 'klt_pins_full', 'klt_combine', 'klt_full']
    figureNames = {key: val for key, val in zip(basisTypes, ['DLP', 'POD_uo2', 'POD_mox', 'POD_pins', 'POD_combine', 'POD_full'])}
    styles = {key: val for key, val in zip(basisTypes, ['ko:', 'mv:', 'c*:', 'g^:', 'bd:', 'rs:'])}
    basisTypes = ['dlp', 'klt_pins_full', 'klt_combine', 'klt_full']
elif geo == 2:  # BWR core 1
    basisTypes = ['dlp', 'klt_pins_core1', 'klt_assay_12', 'klt_assay_all', 'klt_core1']
    figureNames = {key: val for key, val in zip(basisTypes, ['DLP', 'POD_pins', 'POD_assay12', 'POD_assay_all', 'POD_core'])}
    styles = ['ko:', 'm^:', 'bd:', 'gv:', 'rs:']
elif geo == 3:  # BWR core 2
    basisTypes = ['dlp', 'klt_pins_core2', 'klt_assay_12', 'klt_assay_all', 'klt_core2']
    figureNames = {key: val for key, val in zip(basisTypes, ['DLP', 'POD_pins', 'POD_assay12', 'POD_assay_all', 'POD_core'])}
    styles = ['ko:', 'm^:', 'bd:', 'gv:', 'rs:']
else:
    basisTypes = ['dlp', 'mdlp']
    figureNames = {key: val for key, val in zip(basisTypes, ['DLP', 'mDLP'])}
    styles = ['ko:', 'ms:']

assert geo in [1, 2, 3]

def getGeoMap(geoOption):
    if geoOption == 1:
        # Define the number of pins of UO2
        nPins = 5
        # Define the fine mesh
        fm = [3, 22, 3] * nPins * 2
        # Define the coarse mesh
        x = [0.0, 0.09, 1.17, 1.26]
        cm = [xx + i * x[-1] for i in range(10) for xx in x[:-1]] + [2 * nPins * 1.26]
        # Define the material map
        mm = [5, 1, 5] * nPins + [5, 4, 5] * nPins
    elif geoOption == 2:
        ass1 = [5, 1, 2, 2, 1, 5]
        ass2 = [5, 1, 3, 3, 1, 5]
        fm_ass = [9, 27, 27, 27, 27, 9]
        cm_ass = [1.1176, 3.2512, 3.2512, 3.2512, 3.2512, 1.1176]
        ass_map = [1, 2, 1, 2, 1, 2, 1]
        cm = [0.0]
        fm = []
        mm = []
        for i, ass in enumerate(ass_map):
            mm += ass1 if ass == 1 else ass2
            cm += cm_ass
            fm += fm_ass
        cm = np.cumsum(cm)
    elif geoOption == 3:
        ass1 = [5, 1, 2, 2, 1, 5]
        ass2 = [5, 1, 4, 4, 1, 5]
        fm_ass = [9, 27, 27, 27, 27, 9]
        cm_ass = [1.1176, 3.2512, 3.2512, 3.2512, 3.2512, 1.1176]
        ass_map = [1, 2, 1, 2, 1, 2, 1]
        cm = [0.0]
        fm = []
        mm = []
        for i, ass in enumerate(ass_map):
            mm += ass1 if ass == 1 else ass2
            cm += cm_ass
            fm += fm_ass
        cm = np.cumsum(cm)

    return fm, cm, mm


# Plotting options
legend_loc = 'best'


def getOrder(G, o, contig=0):
    '''

    '''
    groupMask, counts, mask, extraName = getGroupBounds(G, contig)
    return sum([min(g, o) + 1 for g in counts.values()])


def makePhiPlot(G, data, orders, homog):
    basisTypes = ['dlp', 'klt_combine', 'klt_full', 'klt_pins_full']
    linestyles = [':', '-', '-.', '--']

    savedata = []
    header = '{:23s}'.format('order')
    for i in range(len(data[basisTypes[0]][0])):
        header += '{:25s}'.format('group ' + str(i+1))

    for i, basis in enumerate(basisTypes):
        np.savetxt('plots/phi_{}_{}_data_h{}.txt'.format(basis, G, homog), np.concatenate(((np.array(orders)-1)[:,np.newaxis], np.array(data[basis])), axis=1), header=header)

        d = np.array(data[basis]) * 100
        plt.fill_between(orders, np.min(d, axis=1), np.max(d, axis=1), color=styles[basis][0], alpha=0.2, linestyle=linestyles[i], linewidth=2.5)
        #plt.semilogy(orders, d, styles[basisTypes.index(basis)], alpha=0.2, marker='')

        plt.semilogy(orders, np.mean(d, axis=1), styles[basis], label=figureNames[basis], linestyle=linestyles[i])

    plt.xlim([orders[0], orders[-1]])
    plt.ylim([1e-8, 1e2])
    plt.xlabel('Degrees of Freedom')
    plt.ylabel('$max \left( abs(\phi_g - \phi^{ref}_g) \\right)$ / $\phi^{ref}_g$ [%]')
    plt.legend(loc=legend_loc, ncol=2, fancybox=True, framealpha=0.0)
    plt.savefig('plots/fill_phi_{}_h{}.png'.format(G, homog), transparent=True)
    plt.clf()


def makeFigure(G, data, plotType, maxOrder, orders, homog):

    if plotType == 'phi':
        makePhiPlot(G, data, orders, homog)
        return

        ylabel = 'max scalar flux error'
        mul = 1
        ylim = [1e-8, 1e2]
    elif plotType == 'dens':
        ylabel = 'Max rel FD error [%]'
        mul = 100
        ylim = [1e-8, 1e2]
    elif plotType == 'k':
        ylabel = 'Rel eigenvalue error [%]'
        mul = 100
        ylim = [1e-8, 1e2]
    elif plotType == 'count':
        ylabel = 'Number of sweeps'
        mul = 1
        top = np.max([np.max(data[b]) for b in basisTypes])
        ylim = [1, 10 ** (np.ceil(np.log10(top)))]
    elif plotType == 'time':
        ylabel = 'CPU time [s]'
        mul = 1
        top = np.max([np.max(data[b]) for b in basisTypes])
        ylim = [1, 10 ** (np.ceil(np.log10(top)))]
    elif plotType == 'mem':
        ylabel = 'Memory use [MB]'
        mul = 1.0 / 1000.0
        top = np.max([mul * np.max(data[b]) for b in basisTypes])
        bot = np.max([mul * np.min(data[b]) for b in basisTypes])
        ylim = [10 ** (np.floor(np.log10(bot))), 10 ** (np.ceil(np.log10(top)))]


    if plotType in ['phi', 'dens', 'k']:
        # Make plot of maximum errors
        data = [(orders, mul * np.array(data[b]), styles[b], figureNames[b]) for b in basisTypes]
    else:
        d = [(orders, mul * np.array(data[b])[:,0], styles[b], figureNames[b]) for b in basisTypes]
        d.append((orders, mul * np.array(data[basisTypes[0]])[:,1], 'k-', 'reference'))
        data = d
    labels = ['Degrees of Freedom', ylabel]
    name = 'plots/{}_error_{}_h{}.png'.format(plotType, G, homog)

    savedata = []
    header = '{:23s}'.format('groups')

    for d in data:
        x, y, s, l = d
        savedata.append(y[:])
        header += '{:25s}'.format(l)
        plt.semilogy(x, y, s, label=l)

    np.savetxt('plots/{}_{}_data_h{}.txt'.format(plotType, G, homog), np.array([x] + savedata).T, header=header)

    plt.xlim([0, G])
    plt.ylim(ylim)
    plt.xlabel(labels[0])
    plt.ylabel(labels[1])
    plt.legend(loc=legend_loc, ncol=2, fancybox=True, framealpha=0.0)
    plt.savefig(name, transparent=True)
    plt.clf()

def readData(plotTypes, maxOrder, basisTypes, G, contig, homogs):
    '''
    Load the files and import the necessary data, then compute the error.

    Inputs:
        plotTypes
        errTypes
        maxOrder
        basisTypes
        G
    '''

    data = {}

    for m, plotType in enumerate(plotTypes):
        data[plotType] = {}
        for h in homogs:
            data[plotType][h] = {}
            for k, b in enumerate(basisTypes):
                data[plotType][h][b] = [computeError(plotType, G, o, b, contig, h) for o in range(maxOrder)]
    return data


def computeError(plotType, groups, order, basis, contig, homog):
    '''
    Given the files to compare, compute the error

    Inputs:
        path
        errType
        groups
        order
        basis

    Outputs:
        max_error
        mean_error
        min_error
    '''

    geoName = 'full' if geo == 1 else 'core1' if geo == 2 else 'core2'

    if plotType in ['dens', 'phi']:
        # Get only the fuel cells
        fm, cm, mm = getGeoMap(geo)
        mask = np.array(mm) != 5
        mask = np.array([m for i, m in enumerate(mask) for j in range(fm[i])])

    # Load the reference
    try:
        ref = np.load('{}/ref_{}_{}_{}.npy'.format(ref_path, plotType, geoName, groups))

        # Normalize if scalar flux or fission density
        if plotType in ['phi', 'dens']:
            ref /= np.linalg.norm(ref)
    except IOError:
        raise IOError('Missing reference for {} with {} groups'.format(geoName, groups))

    # Load the truncated results
    try:
        sol = np.load('{}/dgm_{}_{}_{}_g{}_c{}_o{}_h{}.npy'.format(data_path, plotType, basis, geoName, groups, contig, order, homog))
        # Normalize if scalar flux or fission density
        if plotType in ['phi', 'dens']:
            sol /= np.linalg.norm(sol)
    except IOError as e:
        print('Missing {}/dgm_{}_{}_{}_g{}_c{}_o{}_h{}.npy'.format(data_path, plotType, basis, geoName, groups, contig, order, homog))
        sol = np.zeros(ref.shape)

    if plotType == 'phi':
        # Only compute the error for the zeroth moment of phi
        ref = ref[0]
        sol = sol[0]
    elif plotType == 'dens':
        # Only compute the error for fuel cells
        ref = ref[mask]
        sol = sol[mask]

    if plotType in ['phi', 'dens', 'k']:
        # Compute the difference
        err = np.abs(sol - ref)

        # Divide by the mean in each group
        if plotType == 'phi':
            err /= np.mean(ref, axis=1)[:, np.newaxis]
            return np.max(err, axis=1)

        elif plotType in ['k', 'dens']:
            err /= ref

        return np.max(err)
    else:
        return np.max(sol), np.max(ref)


# Which groups to plot
Gs = [44, 238, 1968]
contig = 3
homogs = [0, 1, 2, 3, 4]

# Define the path to the data
data_path = 'data'
ref_path = 'reference'

# Loop throught the group structures
for i, G in enumerate(Gs):
    # Compute the maximum order for group G
    groupMask, counts, mask, extraName = getGroupBounds(G, contig)
    maxOrder = max(counts.values())

    # Create the list for all the orders in the group structure
    orders = [sum([min(o + 1, g) for g in counts.values()]) for o in range(max(counts.values()))]

    # Define which type of plots to make
    plotTypes = ['phi', 'dens', 'count', 'k', 'time', 'mem']

    if G == 1968:
        try:
            basisTypes.remove('klt_uo2')
            basisTypes.remove('klt_mox')
        except ValueError:
            pass
        plotTypes = plotTypes[:-2]

    data = readData(plotTypes, maxOrder, basisTypes, G, contig, homogs)

    for k, plotType in enumerate(plotTypes):
        print(G, plotType)
        for h in homogs:
            makeFigure(G, data[plotType][h], plotType, maxOrder - 1, orders, h)

