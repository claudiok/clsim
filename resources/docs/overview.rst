Overview
========

The first step in simulating the light yield from a given I3Particle at a DOM
is to convert the particle into a series of light-emitting ``steps``. Each step is
assumed to have a constant speed beta=v/c, determining the Cerenkov angle and a 
constant number of photons that should be emitted over a given length.
By default, steps are created using a full *Geant4* simulation, but alternative
parameterizations can be used to speed up the process. (*Geant4* tends to get
rather slow at higher energies (E>10TeV).) The current version of *clsim* comes
with a parameterization that is compatible to *ppc*. Parameterizations can
be used for only a sub-set of particles and energies.

Once a set of steps is generated, they are uploaded to the compute device
(i.e. the GPU). The GPU runs "kernels" in parallel that are responsible for
creating photons from the steps, propagating them through ice layers and
checking for collisions with DOMs. All photons that collided are saved with
their full information (including direction, position on the DOM, wavelength,
number of scatters, properties at its point of creation, ...). Those photons
are then sent back to the host, converted to I3FrameObjects (I3Photon) and
saved in the frame. In order to keep the GPU busy, multiple frames (==events)
are simulated in parallel. The *clsim* module takes care of correctly 
re-assembling the events afterwards.

The output of the GPU simulation step is thus a list of photons at the DOM
surface. These photons still need to be converted into hits, which is done
using a dedicated module. This module finally writes a I3MCHitSeriesMap,
compatible to all existing photon simulation output.
