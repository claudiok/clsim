    
inline void checkForCollision_OnString(
    const unsigned short stringNum,
    const floating_t photonDirLenXYSqr,
    const floating4_t photonPosAndTime,
    const floating4_t photonDirAndWlen,
#ifdef STOP_PHOTONS_ON_DETECTION
    floating_t *thisStepLength,
    bool *hitRecorded,
    unsigned short *hitOnString,
    unsigned short *hitOnDom,
#else
    floating_t thisStepLength,
    floating_t inv_groupvel,
    floating_t photonTotalPathLength,
    uint photonNumScatters,
    const floating4_t photonStartPosAndTime,
    const floating4_t photonStartDirAndWlen,
    const struct I3CLSimStep *step,
    __global uint* hitIndex,
    uint maxHitIndex,
    __write_only __global struct I3CLSimPhoton *outputPhotons,
#endif
    __local const unsigned short *geoLayerToOMNumIndexPerStringSetLocal
    );
    
inline void checkForCollision_InCell(
    const floating_t photonDirLenXYSqr,
    const floating4_t photonPosAndTime,
    const floating4_t photonDirAndWlen,
#ifdef STOP_PHOTONS_ON_DETECTION
    floating_t *thisStepLength,
    bool *hitRecorded,
    unsigned short *hitOnString,
    unsigned short *hitOnDom,
#else
    floating_t thisStepLength,
    floating_t inv_groupvel,
    floating_t photonTotalPathLength,
    uint photonNumScatters,
    const floating4_t photonStartPosAndTime,
    const floating4_t photonStartDirAndWlen,
    const struct I3CLSimStep *step,
    __global uint* hitIndex,
    uint maxHitIndex,
    __write_only __global struct I3CLSimPhoton *outputPhotons,
#endif
    __local const unsigned short *geoLayerToOMNumIndexPerStringSetLocal,
    
    __constant unsigned short *this_geoCellIndex,
    const floating_t this_geoCellStartX,
    const floating_t this_geoCellStartY,
    const floating_t this_geoCellWidthX,
    const floating_t this_geoCellWidthY,
    const int this_geoCellNumX,
    const int this_geoCellNumY
    );
    
inline void checkForCollision_InCells(
    const floating_t photonDirLenXYSqr,
    const floating4_t photonPosAndTime,
    const floating4_t photonDirAndWlen,
#ifdef STOP_PHOTONS_ON_DETECTION
    floating_t *thisStepLength,
    bool *hitRecorded,
    unsigned short *hitOnString,
    unsigned short *hitOnDom,
#else
    floating_t thisStepLength,
    floating_t inv_groupvel,
    floating_t photonTotalPathLength,
    uint photonNumScatters,
    const floating4_t photonStartPosAndTime,
    const floating4_t photonStartDirAndWlen,
    const struct I3CLSimStep *step,
    __global uint* hitIndex,
    uint maxHitIndex,
    __write_only __global struct I3CLSimPhoton *outputPhotons,
#endif
    __local const unsigned short *geoLayerToOMNumIndexPerStringSetLocal
    );

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
    __write_only __global struct I3CLSimPhoton *outputPhotons,
    __local const unsigned short *geoLayerToOMNumIndexPerStringSetLocal
    );


