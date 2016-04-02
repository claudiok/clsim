from icecube import icetray, dataclasses, phys_services
from I3Tray import I3Units

import math

def I3SimpleProtonDecayGenerator(frame,
                                 posZRange=(-327.*I3Units.m,-502.*I3Units.m),
                                 radius=50.*I3Units.m,
                                 randomService=None,
                                 centerX=0.*I3Units.m,
                                 centerY=0.*I3Units.m):
    if randomService is None: raise RuntimeError("No random number generator!")
    
    massProton = 938.271996*I3Units.MeV
    massPi0 = 134.9766*I3Units.MeV
    massPositron = 0.510998910*I3Units.MeV
    
    # 2-body decay in proton rest-frame from 4-momentum conservation (TODO: check)
    ETotPi0 = (massProton**2. + massPi0**2. - massPositron**2.)/(2.*massProton)
    ETotPositron = (massProton**2. + massPositron**2. - massPi0**2.)/(2.*massProton)
    EKinPi0 = ETotPi0-massPi0
    EKinPositron = ETotPositron-massPositron
    
    posZ = randomService.uniform(posZRange[0], posZRange[1])
    while True:
        posX = randomService.uniform(-radius, radius)
        posY = randomService.uniform(-radius, radius)
        posR = math.sqrt(posX**2. + posY**2.)
        if posR <= radius: break
    
    dirZenith = math.acos(randomService.uniform(-1.,1.))
    dirAzimuth = randomService.uniform(0., 2.*math.pi)

    #posZ = (posZRange[0]+posZRange[1])/2.
    #posX = 0.
    #posY = 0.
    #dirZenith = 90.*I3Units.deg
    

    posX += centerX
    posY += centerY
    
    # remove existing I3MCTree
    if "I3MCTree" in frame: del frame["I3MCTree"]
    
    mcTree = dataclasses.I3MCTree()
    
    # clsim interprets I3Particle.Get/SetEnergy() as the kinetic energy
    primaryProton = dataclasses.I3Particle()
    primaryProton.location_type = dataclasses.I3Particle.LocationType.Anywhere # not InIce!
    primaryProton.shape = dataclasses.I3Particle.ParticleShape.Primary
    primaryProton.type = dataclasses.I3Particle.ParticleType.PPlus
    primaryProton.speed = 0.
    primaryProton.energy = 0.*I3Units.GeV
    primaryProton.time = 10.*I3Units.ns
    primaryProton.pos = dataclasses.I3Position(posX, posY, posZ)
    primaryProton.dir = dataclasses.I3Direction(dirZenith, dirAzimuth)
    mcTree.add_primary(primaryProton)

    childPiZero = dataclasses.I3Particle()
    childPiZero.location_type = dataclasses.I3Particle.LocationType.InIce
    childPiZero.shape = dataclasses.I3Particle.ParticleShape.Cascade
    childPiZero.type = dataclasses.I3Particle.ParticleType.Pi0
    childPiZero.energy = EKinPi0 # kinetic energy
    childPiZero.time = primaryProton.time
    childPiZero.pos = dataclasses.I3Position(posX, posY, posZ)
    childPiZero.dir = primaryProton.dir
    mcTree.append_child(primaryProton, childPiZero)

    childPositron = dataclasses.I3Particle()
    childPositron.location_type = dataclasses.I3Particle.LocationType.InIce
    childPositron.shape = dataclasses.I3Particle.ParticleShape.Cascade
    childPositron.type = dataclasses.I3Particle.ParticleType.EPlus
    #childPositron.type = dataclasses.I3Particle.ParticleType.MuPlus
    childPositron.energy = EKinPositron # kinetic energy
    childPositron.time = primaryProton.time
    childPositron.pos = dataclasses.I3Position(posX, posY, posZ)
    childPositron.dir = dataclasses.I3Direction(-primaryProton.dir.x,
                                                -primaryProton.dir.y,
                                                -primaryProton.dir.z)
    mcTree.append_child(primaryProton, childPositron)
    
    frame["I3MCTree"] = mcTree

