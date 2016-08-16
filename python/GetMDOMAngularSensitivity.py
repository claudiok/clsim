#!/usr/bin/env python


import numpy as np
import os
from icecube import icetray, dataclasses
from icecube.clsim import I3CLSimFunctionFromTable
from icecube.clsim import util
from I3Tray import I3Units
from os.path import expandvars


def GetMDOMAngularSensitivity():

    angularAcceptanceFile = os.environ['I3_SRC'] + '/clsim/resources/tablemaker/mDOM/AngularAcceptancePMT_FlatDisk_mDOM_cosineCorrected_Jul16.txt'
    angularAcceptanceFromFile = np.loadtxt(angularAcceptanceFile, dtype = {'names' : ('cos', 'acceptance'),
                                                                           'formats' : ('f8', 'f8')})
    AngularAcceptanceTmpFnc = util.interpolate.interp1d(angularAcceptanceFromFile['cos'], angularAcceptanceFromFile['acceptance'])
    cos_bins = np.linspace(-1.0, 1.0, 1001)
    angularAcceptanceList = list()
    for cos_bin in cos_bins:
        if cos_bin <= 0:
            angularAcceptanceList.append(0.0)
        else:
            acceptance = AngularAcceptanceTmpFnc(cos_bin) * 1.0
            angularAcceptanceList.append(acceptance)
   
    AngularAcceptanceFnc = I3CLSimFunctionFromTable(cos_bins[0], (cos_bins[1] - cos_bins[0]) , angularAcceptanceList)
    
    onesArray = np.ones(len(angularAcceptanceList))
    onesFnc = I3CLSimFunctionFromTable(cos_bins[0], (cos_bins[1] - cos_bins[0]) , onesArray)
    #return onesFnc
    
    return AngularAcceptanceFnc


