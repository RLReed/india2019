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
np.set_printoptions(threshold=10000)


def getMats(geoOption):
    if geoOption == 'inf_med':
        return [0]
    elif geoOption == 'full':
        return [0, 3]
    elif geoOption == 'core1':
        return [0, 1, 2]
    elif geoOption == 'core2':
        return [0, 1, 3]


def getXS(fname):
    '''
    Read the input file and extract the total cross sections

    input
        fname: name of the file to read

    output
        sig_t: numpy array containing the total cross sections for each material
    '''

    # Initialize the total cross section list
    sig_t = []
    with open(fname) as f:
        # Get the number of materials and number of groups
        nMat, nGroups = [int(d) for d in f.readline().split()[:2]]
        # Skip 2 lines to pass the energy bounds and velocity
        f.readline()
        f.readline()

        # Loop over the materials
        for m in range(nMat):
            # Initialize temp list for total XS for material m
            t = []
            # Skip the material name
            f.readline()
            # Get the number of legendre moments
            nLegendre = int(f.readline().split()[0])
            # Loop over groups and extract the total XS
            for g in range(nGroups):
                t.append(float(f.readline().split()[0]))
            # Skip the scattering moments
            for g in range(nGroups * nLegendre):
                f.readline()
            # Append the temp XS list to the main list
            sig_t.append(t)

    return np.array(sig_t)


def findContiguousBounds(sig_t, minCutoff, ratioCutoff, groupCutoff):
    '''
    Find the coarse group structure given a set of cross sections

    Inputs:
        sig_t: array of the total XS indexed by material [mat, G]
        minCutoff: XS below which the ratio cutoff is ignored
        ratioCutoff: Ratio of largest to smallest XS in a coarse group
        groupCutoff: Maximum number of fine groups in a coarse group
    Outputs:
        bounds: Starting fine group indices for each coarse group
    '''

    def reset(xs):
        '''Reset the min and max cutoffs to the input cross section'''
        minXS = np.maximum(xs, minCutoff)
        maxXS = xs
        return minXS, maxXS

    sig_t = np.array(sig_t)

    # Get the number of groups
    nGroup = len(sig_t[0])

    # Initialize the boundary at the lowest energy group
    bounds = [nGroup + 1]

    # Initialize the cutoff bounds
    minXS, maxXS = reset(sig_t[:,-1])

    # Loop through the cross sections from greatest to least
    for i, xs in enumerate(sig_t.T[::-1]):
        group = nGroup - i

        # Check if the xs is below the min or above the max
        minXS = np.minimum(xs, minXS)
        maxXS = np.maximum(xs, maxXS)
        ratio = maxXS / minXS

        # Check for a ratio that is too large
        if (max(ratio) > ratioCutoff and max(maxXS) > minCutoff) or bounds[-1] - group > groupCutoff:
            bounds.append(group + 1)
            # Reset the cutoff bounds
            minXS, maxXS = reset(xs)

    # Add the highest energy group bound
    bounds.append(1)
    # Reverse to the natural order of structures
    bounds = np.array(bounds[::-1])

    return np.array([i for i, b in enumerate(bounds[1:] - bounds[:-1]) for j in range(b)])

def findNoncontiguousBounds(sig_t, minCutoff=1.0, ratioCutoff=2.0, groupCutoff=60):
    '''
    Find the coarse group structure given a set of cross sections

    Inputs:
        sig_t: array of the total XS indexed by material [mat, G]
        minCutoff: XS below which the ratio cutoff is ignored
        ratioCutoff: Ratio of largest to smallest XS in a coarse group
        groupCutoff: Maximum number of fine groups in a coarse group
    Outputs:
        structure: mapping of fine group to coarse group
    '''

    G = len(sig_t[0])

    # Get the minimum and maximum of all materials
    minXS = np.min(sig_t, axis=0)
    minmap = np.argsort(minXS)[::-1]
    maxXS = np.max(sig_t, axis=0)
    maxmap = np.argsort(maxXS)[::-1]

    unsorted_groups = list(range(G))

    # Initialize a dictionary
    structure = np.zeros(G) - 1
    coarse_group_index = 0
    while len(unsorted_groups) > 0:
        # Get the largest unsorted cross section
        maximum = maxXS[maxmap[0]]
        # If above the minimum cross section cutoff
        if maximum > minCutoff:
            # Find the ratio for the current group
            #cutoff = max(ratioCutoff, maximum / minXS[maxmap[0]])
            cutoff = ratioCutoff
            # Select all groups with a smaller ratio than the cutoff
            cutmap = minmap[maximum / minXS[minmap] <= cutoff]
            if len(cutmap) == 0:
                cutmap = [maxmap[0]]
        else:
            # Select all remaining groups
            cutmap = minmap
        # Only allow groupCutoff groups into the current group
        cutmap = cutmap[:groupCutoff]
        # Assign selected groups to the coarse_group_index
        structure[cutmap] = coarse_group_index
        # Remove the assigned groups from minmap, maxmap, and unsorted_groups
        minmap = np.array([g for g in minmap if g not in cutmap])
        maxmap = np.array([g for g in maxmap if g not in cutmap])
        unsorted_groups = [g for g in unsorted_groups if g not in cutmap]
        # Increment the coarse_group_index
        coarse_group_index += 1

    return (structure * -1 + max(structure)).astype(int)

def computeBounds(G, fname, geoOption, contig, sig_t, minCutoff, ratioCutoff, groupCutoff):
    '''
    This function determines which fine groups should belong to each coarse group

    The boundaries will not necessarily be contiguous

    Inputs:
        fname: name of the cross section file in anlxs format
    '''

    # Contig options
    names = {0: 'contiguous_min',
             1: 'contiguous_max',
             2: 'contiguous_mean',
             3: 'continuous_difference',
             4: 'noncontiguous_min',
             5: 'noncontiguous_max',
             6: 'noncontiguous_mean',
             7: 'noncontiguous_difference'}

    mats = getMats(geoOption)

    mask = np.array([r in mats for r in range(len(sig_t))]).astype(bool)
    sig_t = sig_t[mask]

    # Number the energy groups
    groups = np.arange(G)

    # Get the minimum XS across all materials
    minXS = np.min(sig_t, axis=0)
    # Get the maximum XS across all materials
    maxXS = np.max(sig_t, axis=0)
    # Get the mean XS across all materials
    meanXS = np.mean(sig_t, axis=0)

    # Select which XS upon which to compute the bounds
    if contig == 0:  # contiguous min
        mask = findContiguousBounds([minXS], minCutoff, ratioCutoff, groupCutoff)
    elif contig == 1:  # contiguous max
        mask = findContiguousBounds([maxXS], minCutoff, ratioCutoff, groupCutoff)
    elif contig == 2:  # contiguous mean
        mask = findContiguousBounds([meanXS], minCutoff, ratioCutoff, groupCutoff)
    elif contig == 3:  # continuous difference
        mask = findContiguousBounds(sig_t, minCutoff, ratioCutoff, groupCutoff)
    elif contig == 4:  # non-contiguous min
        mask = findNoncontiguousBounds([minXS], minCutoff, ratioCutoff, groupCutoff)
    elif contig == 5:  # non-contiguous max
        mask = findNoncontiguousBounds([maxXS], minCutoff, ratioCutoff, groupCutoff)
    elif contig == 6:  # non-contiguous mean
        mask = findNoncontiguousBounds([meanXS], minCutoff, ratioCutoff, groupCutoff)
    elif contig == 7:  # non-contiguous difference
        mask = findNoncontiguousBounds(sig_t, minCutoff, ratioCutoff, groupCutoff)

    # Get the coarse group bounds given the total cross sections
    s = '{}\n{}'.format(names[contig], mask.tolist())

    with open('XS/structure{}-{}-m{}-r{}-g{}'.format(G, contig, minCutoff, ratioCutoff, groupCutoff), 'w') as f:
        f.write(s)

    # Plot the cross sections
    colors = ['b', 'g', 'r', 'c', 'm', 'orange', 'purple', 'grey']
    ls = '-' if contig < 4 else '--'

    for i, m in enumerate(mask):
        #submask = mask == m
        #plt.fill_between(groups[submask], minXS[submask], maxXS[submask], label='Group {}'.format(m), linestyle=ls, color=colors[m%len(colors)])
        plt.semilogy(groups[i], minXS[i], label='Group {}'.format(m), ls=ls, marker='o', c=colors[m%len(colors)])
        plt.semilogy(groups[i], maxXS[i], label='Group {}'.format(m), ls=ls, marker='s', c=colors[m%len(colors)])
    plt.yscale('log')
    plt.xlabel('Fine group number')
    plt.ylabel('Total cross section [cm$^{-1}$]')
    plt.grid(True)
    plt.title(names[contig])
    plt.savefig('{}/{}g_{}-m{}-r{}-g{}.pdf'.format('/'.join(fname.split('/')[:-1]), G, contig, minCutoff, ratioCutoff, groupCutoff))
    plt.clf()

    return s


if __name__ == '__main__':
    geo = 1
    Gs = [44, 238, 1968]
    contigs = range(8)
    for G in Gs:
        fname = 'XS/{0}gXS.anlxs'.format(G)

        # Get the XS for the materials in the problem geometry
        sig_t = getXS(fname)

        for contig in contigs:
            print 'Working on structure for group {} config {}'.format(G, contig)
            computeBounds(G, fname, geo, contig, sig_t)
