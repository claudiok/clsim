
#define R2_NBINS 5
#define R2_STEP (25./R2_NBINS)
#define COSTHETA_NBINS 6
#define COSTHETA_STEP (2./COSTHETA_NBINS)
#define AZIMUTH_NBINS 6
#define AZIMUTH_STEP (PI/COSTHETA_NBINS)
#define MIN(a, b) ((a < b) ? a : b);
#define MAX(a, b) ((a > b) ? a : b);

// #define OM_RADIUS 0.2;

__constant const floating4_t sourceOrigin = {0., 0., 0., 0.};
__constant const floating4_t sourceDir = {0., 0., 1., 0.};
__constant const floating4_t sourceDirP = {0., 1., 0., 0.};
__constant const floating4_t sourceDirPx = {1., 0., 0., 0.};

inline void checkForCollision_inCell(
    unsigned short rbin, unsigned short cbin, unsigned short abin,
    const floating4_t photonPosAndTime, const floating4_t photonDirAndWlen,
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
    __global struct I3CLSimPhoton *outputPhotons
#endif
    )
{
    floating_t sina, cosa;
    sina = sincos(AZIMUTH_STEP*(abin+0.5), &cosa);
    const floating_t cost = COSTHETA_STEP*(cbin+0.5)-1.;
    const floating_t sint = my_sqrt(1.-sqr(cost));
    const floating_t r = sqr(R2_STEP*(rbin+0.5));
    const floating4_t domPos = r*(sint*(sina*sourceDirPx + cosa*sourceDirP) + cost*sourceDir)
        + sourceOrigin + (floating4_t)(ZERO, ZERO, ZERO, photonPosAndTime.w);

    const floating4_t drvec = photonPosAndTime - domPos;
    const floating_t dr2 = dot(drvec,drvec);
    const floating_t urdot = dot(drvec, photonDirAndWlen);
    floating_t discr = sqr(urdot) - dr2 + OM_RADIUS*OM_RADIUS;
    
    if ((discr >= ZERO) && (dr2 > OM_RADIUS*OM_RADIUS)) {
        
        floating_t smin = -urdot - discr;
        if (smin < ZERO) smin = -urdot + discr;
        
        if ((smin >= ZERO) && (smin < thisStepLength)) {
                const unsigned short domish = rbin;
                const unsigned short stringish = abin*COSTHETA_NBINS + cbin;
                saveHit(photonPosAndTime,
                        photonDirAndWlen,
                        smin, // this is the limited thisStepLength
                        inv_groupvel,
                        photonTotalPathLength,
                        photonNumScatters,
                        photonStartPosAndTime,
                        photonStartDirAndWlen,
                        step,
                        stringish, // "String" number
                        domish, // "OM" number
                        hitIndex,
                        maxHitIndex,
                        outputPhotons
                        );
        }
    }
}  

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
    )
{
#ifdef DEBUG_STORE_GENERATED_PHOTONS
    saveHit(photonPosAndTime,
            photonDirAndWlen,
            ZERO,
            inv_groupvel,
            photonTotalPathLength,
            photonNumScatters,
            photonStartPosAndTime,
            photonStartDirAndWlen,
            step,
            0,
            0,
            hitIndex,
            maxHitIndex,
            outputPhotons
            );
    return true;
#else // DEBUG_STORE_GENERATED_PHOTONS

    // check for collisions
    const floating_t photonDirLenXYSqr = sqr(photonDirAndWlen.x) + sqr(photonDirAndWlen.y);
    if (photonDirLenXYSqr <= ZERO) return false;

#ifdef STOP_PHOTONS_ON_DETECTION
    bool hitRecorded=false;
    unsigned short hitOnString;
    unsigned short hitOnDom;
#endif

    floating4_t photonDisplacement = photonPosAndTime-sourceOrigin;
    photonDisplacement.w = 0.;
    
    const floating_t photonR2 = my_sqrt(length(photonDisplacement));
    photonDisplacement = normalize(photonDisplacement);
    const floating_t photonCosTheta = dot(sourceDir, photonDisplacement);
    const floating_t photonAzimuth = atan2(fabs(dot(sourceDirP, photonDisplacement)), dot(sourceDirPx, photonDisplacement));
    
    photonDisplacement = photonStartPosAndTime-sourceOrigin;
    photonDisplacement.w = 0.;
    
    const floating_t photonStartR2 = my_sqrt(length(photonDisplacement));
    photonDisplacement = normalize(photonDisplacement);
    const floating_t photonStartCosTheta = dot(sourceDir, photonDisplacement);
    const floating_t photonStartAzimuth = atan2(fabs(dot(sourceDirP, photonDisplacement)), dot(sourceDirPx, photonDisplacement));
    
    const unsigned short r2Bin = photonR2/R2_STEP;
    const unsigned short cosThetaBin = (photonCosTheta+1.)/COSTHETA_STEP;
    const unsigned short azimuthBin = (photonAzimuth)/AZIMUTH_STEP;
    
    const unsigned short r2StartBin = photonStartR2/R2_STEP;
    const unsigned short cosThetaStartBin = (photonStartCosTheta+1.)/COSTHETA_STEP;
    const unsigned short azimuthStartBin = (photonStartAzimuth)/AZIMUTH_STEP;
    
    const unsigned short r0 = MIN(r2Bin, r2StartBin); const unsigned short r1 = MAX(r2Bin, r2StartBin);
    const unsigned short c0 = MIN(cosThetaBin, cosThetaStartBin); const unsigned short c1 = MAX(cosThetaBin, cosThetaStartBin);
    const unsigned short a0 = MIN(azimuthBin, azimuthStartBin); const unsigned short a1 = MAX(azimuthBin, azimuthStartBin);
    
    for (unsigned short r=r0; r <= r1; r++) {
        for (unsigned short c=c0; c <= c1; c++) {
            for (unsigned short a=a0; a <= a1; a++)
            {
                checkForCollision_inCell(r, c, a, photonPosAndTime, photonDirAndWlen,
                    thisStepLength,
                    inv_groupvel,
                    photonTotalPathLength,
                    photonNumScatters,
                    photonStartPosAndTime,
                    photonStartDirAndWlen,
                    step,
                    hitIndex,
                    maxHitIndex,
                    outputPhotons);
            }
            
        }
    
    }


#ifdef STOP_PHOTONS_ON_DETECTION
    // In case photons are stopped on detection
    // (i.e. absorbed by the DOM), we need to record
    // them here (after all possible DOM intersections
    // have been checked). 
    //
    // Otherwise, the recording is done right after
    // the intersection detection further down in
    // checkForCollision_*().
    if (hitRecorded) {
        saveHit(photonPosAndTime,
                photonDirAndWlen,
                *thisStepLength,
                inv_groupvel,
                photonTotalPathLength,
                photonNumScatters,
                photonStartPosAndTime,
                photonStartDirAndWlen,
                step,
                hitOnString,
                hitOnDom,
                hitIndex,
                maxHitIndex,
                outputPhotons
                );
    }
    return hitRecorded;
#else // STOP_PHOTONS_ON_DETECTION
    // in case photons should *not* be absorbed when they
    // hit a DOM, this will always return false (i.e.
    // no detection.)
    return false;
#endif // STOP_PHOTONS_ON_DETECTION
#endif // DEBUG_STORE_GENERATED_PHOTONS
}
