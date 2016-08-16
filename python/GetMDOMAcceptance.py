#!/usr/bin/env python

import numpy as np
import os
from icecube import icetray, dataclasses
from icecube.clsim import I3CLSimFunctionFromTable
from I3Tray import I3Units
from os.path import expandvars

def GetMDOMAcceptance():

    ### this' the new shit with flat disc 'n' everything
    
    filename = os.environ['I3_SRC'] + '/clsim/resources/tablemaker/mDOM/Acceptance_PMT_in_mDOM_IceTray.txt'
    acceptancemDOM = np.loadtxt(filename, dtype = { 'names' : ('wlen', 'acceptance'),
                                                       'formats' : ('f8', 'f8')})    
    efficiencyFnc = I3CLSimFunctionFromTable(250.0 * I3Units.nanometer,  2.0 * I3Units.nanometer, acceptancemDOM['acceptance'])


    ### this is the old stuff

    '''
    filename = os.environ['I3_SRC'] + '/clsim/resources/tablemaker/mDOM/EffectiveAreaPMTinmDOM_IceTray.txt'
    effectiveAreamDOM = np.loadtxt(filename, dtype = { 'names' : ('wlen', 'area'),
                                                       'formats' : ('f8', 'f8')})
    R_ICDOM = 0.1651
    ICDOM_Area = np.pi * R_ICDOM**2 
    R_MDOM = 0.178 * icetray.I3Units.m
    MDOMArea = R_MDOM**2 * np.pi
    ICDOMArea = 0.0856335639674771 
    normedEfficiency = effectiveAreamDOM['area'] / MDOMArea
    efficiencyFnc = I3CLSimFunctionFromTable(250.0 * I3Units.nanometer,  2.0 * I3Units.nanometer, normedEfficiency)
    '''

    onesArray = np.ones(len(acceptancemDOM['acceptance']))
    onesFnc = I3CLSimFunctionFromTable(250.0 * I3Units.nanometer,  2.0 * I3Units.nanometer, onesArray)
    #return onesFnc
    
    return efficiencyFnc

        

