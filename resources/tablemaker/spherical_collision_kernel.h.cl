inline bool checkForCollision(const floating4_t photonPosAndTime,
    const floating4_t photonDirAndWlen,
    floating_t inv_groupvel,
    floating_t photonTotalPathLength,
    uint photonNumScatters,
    const floating4_t photonStartPosAndTime,
    const floating4_t photonStartDirAndWlen,
    const struct I3CLSimStep *step,
#ifdef STOP_PHOTONS_ON_DETECTION
    floating_t *thisStepLength,
#else
    floating_t thisStepLength,
#endif
    __global uint* hitIndex,
    uint maxHitIndex,
    __global struct I3CLSimPhoton *outputPhotons,
    __local const unsigned short *geoLayerToOMNumIndexPerStringSetLocal
    );
