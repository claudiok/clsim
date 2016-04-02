#!/usr/bin/env python
"""
Test the Spice-Lea vector transformations against
code directly ported to python from PPC.
"""

from __future__ import print_function
import numpy
import math

from icecube import icetray, dataclasses, clsim, phys_services
from I3Tray import I3Units

# test parameters
numberOfTrials = 10000
maximumDeviation = 1e-14

# anisotropy parameters
thx=216.
logk1=0.04
logk2=-0.08


# derived values for the PPC calculation
azx = numpy.cos(thx*I3Units.deg) # x-component of the "tilt" direction
azy = numpy.sin(thx*I3Units.deg) # y-component of the "tilt" direction
k1  = numpy.exp(logk1)           # coefficient of anisotropy parallel to "tilt" direction
k2  = numpy.exp(logk2)           # coefficient of anisotropy perpendicular to "tilt" direction
kz  = 1./(k1*k2)                 # a normalizing factor for the z direction


def evaluateVectorTransformationCLSim(vector, transform):
    return numpy.array(transform.ApplyTransform(numpy.array(vector)))

def evaluateVectorTransformationPPCPre(vector, azx, azy, k1, k2, kz):
    """
    Straight from PPC: pro.cu:621
    """
    retvec = numpy.zeros(3)

    n1=( azx*vector[0]+azy*vector[1])*k1
    n2=(-azy*vector[0]+azx*vector[1])*k2
    nx=n1*azx-n2*azy
    ny=n1*azy+n2*azx
    nz=vector[2]*kz
    r=1./numpy.sqrt(nx*nx+ny*ny+nz*nz)
    retvec[0]=r*nx
    retvec[1]=r*ny
    retvec[2]=r*nz

    return retvec

def evaluateVectorTransformationPPCPost(vector, azx, azy, k1, k2, kz):
    """
    Straight from PPC: pro.cu:636
    """
    retvec = numpy.zeros(3)

    n1=( azx*vector[0]+azy*vector[1])/k1
    n2=(-azy*vector[0]+azx*vector[1])/k2
    nx=n1*azx-n2*azy
    ny=n1*azy+n2*azx
    nz=vector[2]/kz
    r=1./numpy.sqrt(nx*nx+ny*ny+nz*nz)
    retvec[0]=r*nx
    retvec[1]=r*ny
    retvec[2]=r*nz

    return retvec

dummy, preTransform, postTransform = clsim.util.GetSpiceLeaAnisotropyTransforms(
    anisotropyDirAzimuth=thx*I3Units.deg, 
    magnitudeAlongDir=logk1,
    magnitudePerpToDir=logk2)


zeniths = numpy.arccos(numpy.random.uniform(0.,1.,numberOfTrials)*2. - 1.)
azimuths = numpy.random.uniform(0.,2.*math.pi,numberOfTrials)
xVals = numpy.sin(zeniths)*numpy.cos(azimuths)
yVals = numpy.sin(zeniths)*numpy.sin(azimuths)
zVals = numpy.cos(zeniths)
vectors = numpy.array([xVals, yVals, zVals]).T

maxDeviationInPre=0.
maxDeviationInPost=0.

for i in range(len(vectors)):
    preV_clsim = evaluateVectorTransformationCLSim(vectors[i], preTransform)
    preV_ppc   = evaluateVectorTransformationPPCPre(vectors[i], azx, azy, k1, k2, kz)

    postV_clsim = evaluateVectorTransformationCLSim(vectors[i], postTransform)
    postV_ppc   = evaluateVectorTransformationPPCPost(vectors[i], azx, azy, k1, k2, kz)

    # print "clsim(pre) :", vectors[i], "->", preV_clsim
    # print "  ppc(pre) :", vectors[i], "->", preV_ppc
    # print "clsim(post):", vectors[i], "->", postV_clsim
    # print "  ppc(post):", vectors[i], "->", postV_ppc

    diff_pre  = preV_clsim  - preV_ppc
    diff_post = postV_clsim - postV_ppc

    if numpy.amax(numpy.abs(diff_pre)) > maxDeviationInPre:
        maxDeviationInPre = numpy.amax(numpy.abs(diff_pre))
    if numpy.amax(numpy.abs(diff_post)) > maxDeviationInPost:
        maxDeviationInPost = numpy.amax(numpy.abs(diff_post))

print("maximum absolute deviation in pre-rotation transformation: ", maxDeviationInPre)
print("maximum absolute deviation in post-rotation transformation:", maxDeviationInPost)

if maxDeviationInPre > maximumDeviation:
    raise RuntimeError("clsim implementation results differ from ppc reference implementation results (in pre-rotation transformation)!")

if maxDeviationInPost > maximumDeviation:
    raise RuntimeError("clsim implementation results differ from ppc reference implementation results (in post-rotation transformation)!")

print("test successful!")


