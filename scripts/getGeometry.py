import numpy as np

def setGeometry(geoOption):
    if geoOption == 'inf_med':
        N = 1
        fm = [N]
        cm = [0.0, float(N)]
        mm = [1]
        bounds = [1.0, 1.0]
    elif geoOption == 'uo2-1':
        fm = [3, 22, 3]
        cm = [0.0, 0.09, 1.17, 1.26]
        mm = [5, 1, 5]
        bounds = [1.0, 1.0]
    elif geoOption == 'uo2-2':
        fm = [3, 22, 3]
        cm = [0.0, 0.09, 1.17, 1.26]
        mm = [5, 2, 5]
        bounds = [1.0, 1.0]
    elif geoOption == 'uo2-Gd':
        fm = [3, 22, 3]
        cm = [0.0, 0.09, 1.17, 1.26]
        mm = [5, 3, 5]
        bounds = [1.0, 1.0]
    elif geoOption == 'mox':
        fm = [3, 22, 3]
        cm = [0.0, 0.09, 1.17, 1.26]
        mm = [5, 4, 5]
        bounds = [1.0, 1.0]
    elif geoOption == 'junction':
        fm = [3, 22, 3, 3, 22, 3]
        cm = [0.0, 0.09, 1.17, 1.26, 1.35, 2.43, 2.52]
        mm = [5, 1, 5, 5, 4, 5]
        bounds = [1.0, 1.0]
    elif geoOption == 'full':
        # Define the number of pins of UO2
        nPins = 5
        # Define the fine mesh
        fm = [3, 22, 3] * nPins * 2
        # Define the coarse mesh
        x = [0.0, 0.09, 1.17, 1.26]
        cm = [xx + i * x[-1] for i in range(10) for xx in x[:-1]] + [2 * nPins * 1.26]
        cm = cm
        # Define the material map
        mm = [5, 1, 5] * nPins + [5, 4, 5] * nPins
        bounds = [1.0, 1.0]
    elif geoOption == 'assay1':
        fm = [6, 18, 18, 18, 18, 6]
        cm = [0.0, 1.1176, 4.3688, 7.62, 10.8712, 14.1224, 15.24]
        mm = [5, 1, 2, 2, 1, 5]
        bounds = [1.0, 1.0]
    elif geoOption == 'assay2':
        fm = [6, 18, 18, 18, 18, 6]
        cm = [0.0, 1.1176, 4.3688, 7.62, 10.8712, 14.1224, 15.24]
        mm = [5, 1, 3, 3, 1, 5]
        bounds = [1.0, 1.0]
    elif geoOption == 'assay3':
        fm = [6, 18, 18, 18, 18, 6]
        cm = [0.0, 1.1176, 4.3688, 7.62, 10.8712, 14.1224, 15.24]
        mm = [5, 1, 4, 4, 1, 5]
        bounds = [1.0, 1.0]
    elif geoOption == 'core1':
        ass1 = [5, 1, 2, 2, 1, 5]
        ass2 = [5, 1, 3, 3, 1, 5]
        fm_ass = [6, 18, 18, 18, 18, 6]
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
        bounds = [0.0, 0.0]
    elif geoOption == 'core2':
        ass1 = [5, 1, 2, 2, 1, 5]
        ass2 = [5, 1, 4, 4, 1, 5]
        fm_ass = [6, 18, 18, 18, 18, 6]
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
        bounds = [0.0, 0.0]

    return fm, cm, mm, bounds

