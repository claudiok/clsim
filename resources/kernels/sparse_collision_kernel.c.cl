/**
 * Copyright (c) 2011, 2012
 * Claudio Kopper <claudio.kopper@icecube.wisc.edu>
 * and the IceCube Collaboration <http://www.icecube.wisc.edu>
 *
 * Permission to use, copy, modify, and/or distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
 * WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
 * SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION
 * OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN
 * CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 *
 *
 * $Id$
 *
 * @file sparse_collision_kernel.c.cl
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

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
    floating_t distanceTraveledInAbsorptionLengths,
    const floating4_t photonStartPosAndTime,
    const floating4_t photonStartDirAndWlen,
    const struct I3CLSimStep *step,
    __global uint* hitIndex,
    uint maxHitIndex,
    __global struct I3CLSimPhoton *outputPhotons,
#ifdef SAVE_PHOTON_HISTORY
    __global float4 *photonHistory,
    float4 *currentPhotonHistory,
#endif
#endif
    __local const unsigned short *geoLayerToOMNumIndexPerStringSetLocal
    )
{
    // find the string set for this string
    unsigned char stringSet = geoStringInStringSet[stringNum];

    { // check intersection with string cylinder
        // only use test if uhat lateral component is bigger than about 0.1 (NEED to check bigger than zero)
        const floating_t smin = my_divide(sqr(((photonPosAndTime.x - convert_floating_t(geoStringPosX[stringNum]))*photonDirAndWlen.y - (photonPosAndTime.y - convert_floating_t(geoStringPosY[stringNum]))*photonDirAndWlen.x)), photonDirLenXYSqr);
        //if (smin > sqr(convert_floating_t(geoStringRadius[stringNum]))) return;  // NOTE: smin == distance squared
        if (smin > sqr(convert_floating_t(GEO_STRING_MAX_RADIUS))) return;  // NOTE: smin == distance squared
    }

    { // check if photon is above or below the string (geoStringMaxZ and geoStringMinZ do not include the OM radius!)
        if ((photonDirAndWlen.z > ZERO) && (photonPosAndTime.z > geoStringMaxZ[stringNum]+OM_RADIUS)) return;
        if ((photonDirAndWlen.z < ZERO) && (photonPosAndTime.z < geoStringMinZ[stringNum]-OM_RADIUS)) return;
    }

    // this photon could potentially be hitting an om
    // -> check them all

    int lowLayerZ = convert_int((photonPosAndTime.z-geoLayerStartZ[stringSet])/geoLayerHeight[stringSet]);
#ifdef STOP_PHOTONS_ON_DETECTION
    int highLayerZ = convert_int((photonPosAndTime.z+photonDirAndWlen.z*(*thisStepLength)-geoLayerStartZ[stringSet])/geoLayerHeight[stringSet]);
#else
    int highLayerZ = convert_int((photonPosAndTime.z+photonDirAndWlen.z*thisStepLength-geoLayerStartZ[stringSet])/geoLayerHeight[stringSet]);
#endif
    if (highLayerZ<lowLayerZ) {int tmp=lowLayerZ; lowLayerZ=highLayerZ; highLayerZ=tmp;}
    lowLayerZ = min(max(lowLayerZ, 0), geoLayerNum[stringSet]-1);
    highLayerZ = min(max(highLayerZ, 0), geoLayerNum[stringSet]-1);

#ifndef STOP_PHOTONS_ON_DETECTION
    // the number of 64bit integers needed to store bits for all doms
    #define numComponents ((GEO_MAX_DOM_INDEX + 64 - 1)/64)
    ulong dom_bitmask[numComponents];
    for (uint i=0;i<numComponents;++i) dom_bitmask[i]=0;
    #undef numComponents
#endif

    //__constant const unsigned short *geoLayerToOMNumIndex=geoLayerToOMNumIndexPerStringSet + (convert_uint(stringSet)*GEO_LAYER_STRINGSET_MAX_NUM_LAYERS) + lowLayerZ;
    __local const unsigned short *geoLayerToOMNumIndex=geoLayerToOMNumIndexPerStringSetLocal + (convert_uint(stringSet)*GEO_LAYER_STRINGSET_MAX_NUM_LAYERS) + lowLayerZ;
    for (int layer_z=lowLayerZ;layer_z<=highLayerZ;++layer_z,++geoLayerToOMNumIndex)
    {
        const unsigned short domNum = *geoLayerToOMNumIndex;
        if (domNum==0xFFFF) continue; // empty layer for this string

#ifndef STOP_PHOTONS_ON_DETECTION
        // prevent strings from being checked twice
        if (dom_bitmask[stringNum/64] & (1 << convert_ulong(domNum%64))) continue;  // already check this string
        dom_bitmask[stringNum/64] |= (1 << convert_ulong(domNum%64));               // mark this string as checked
#endif
        
        floating_t domPosX, domPosY, domPosZ;
        geometryGetDomPosition(stringNum, domNum, &domPosX, &domPosY, &domPosZ);

        floating_t urdot, discr;
        {
            const floating4_t drvec = (const floating4_t)(domPosX - photonPosAndTime.x,
                                                          domPosY - photonPosAndTime.y,
                                                          domPosZ - photonPosAndTime.z,
                                                          ZERO);
            const floating_t dr2 = dot(drvec,drvec);

            urdot = dot(drvec, photonDirAndWlen); // this assumes drvec.w==0
            discr   = sqr(urdot) - dr2 + OM_RADIUS*OM_RADIUS;   // (discr)^2
        }
        
        if (discr < ZERO) continue; // no intersection with this DOM
        
#ifdef PANCAKE_FACTOR        
        discr = my_sqrt(discr)/PANCAKE_FACTOR;
#else
        discr = my_sqrt(discr);
#endif        

        // by construction: smin1 < smin2

        {
            // distance from current point along the track to second intersection
            const floating_t smin2 = urdot + discr;
            if (smin2 < ZERO) continue; // implies smin1 < 0, so no intersection
        }

        // distance from current point along the track to first intersection
        const floating_t smin1 = urdot - discr;

        // smin2 > 0 && smin1 < 0 means that there *is* an intersection, but we are starting inside the DOM. 
        // This allows photons starting inside a DOM to leave (necessary for flashers):
        if (smin1 < ZERO) continue;

        // if we get here, there *is* an intersection with the DOM (there are two actually, one for
        // the ray enetering the DOM and one when it leaves again). We are interested in the one where enters
        // the ray enters the DOM.
        
        
        // check if distance to intersection <= thisStepLength; if not then no detection 
#ifdef STOP_PHOTONS_ON_DETECTION
        if (smin1 < *thisStepLength)
#else
        if (smin1 < thisStepLength)
#endif
        {
#ifdef STOP_PHOTONS_ON_DETECTION
            // record a hit (for later, the actual recording is done
            // in checkForCollision().)
            *thisStepLength=smin1; // limit step length
            *hitOnString=stringNum;
            *hitOnDom=domNum;
            *hitRecorded=true;
            // continue searching, maybe we hit a closer OM..
            // (in that case, no hit will be saved for this one)
#else //STOP_PHOTONS_ON_DETECTION
            // save the hit right here
            saveHit(photonPosAndTime,
                    photonDirAndWlen,
                    smin1, // this is the limited thisStepLength
                    inv_groupvel,
                    photonTotalPathLength,
                    photonNumScatters,
                    distanceTraveledInAbsorptionLengths,
                    photonStartPosAndTime,
                    photonStartDirAndWlen,
                    step,
                    stringNum,
                    domNum,
                    hitIndex,
                    maxHitIndex,
                    outputPhotons
#ifdef SAVE_PHOTON_HISTORY
                    , photonHistory,
                    currentPhotonHistory
#endif //SAVE_PHOTON_HISTORY
                    );
#endif //STOP_PHOTONS_ON_DETECTION
        }

    }

}

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
    floating_t distanceTraveledInAbsorptionLengths,
    const floating4_t photonStartPosAndTime,
    const floating4_t photonStartDirAndWlen,
    const struct I3CLSimStep *step,
    __global uint* hitIndex,
    uint maxHitIndex,
    __global struct I3CLSimPhoton *outputPhotons,
#ifdef SAVE_PHOTON_HISTORY
    __global float4 *photonHistory,
    float4 *currentPhotonHistory,
#endif
#endif
    __local const unsigned short *geoLayerToOMNumIndexPerStringSetLocal,
    
    __constant unsigned short *this_geoCellIndex,
    const floating_t this_geoCellStartX,
    const floating_t this_geoCellStartY,
    const floating_t this_geoCellWidthX,
    const floating_t this_geoCellWidthY,
    const int this_geoCellNumX,
    const int this_geoCellNumY
    )
{
    int lowCellX = convert_int((photonPosAndTime.x-this_geoCellStartX)/this_geoCellWidthX);
    int lowCellY = convert_int((photonPosAndTime.y-this_geoCellStartY)/this_geoCellWidthY);

#ifdef STOP_PHOTONS_ON_DETECTION
    int highCellX = convert_int((photonPosAndTime.x+photonDirAndWlen.x*(*thisStepLength)-this_geoCellStartX)/this_geoCellWidthX);
    int highCellY = convert_int((photonPosAndTime.y+photonDirAndWlen.y*(*thisStepLength)-this_geoCellStartY)/this_geoCellWidthY);
#else
    int highCellX = convert_int((photonPosAndTime.x+photonDirAndWlen.x*thisStepLength-this_geoCellStartX)/this_geoCellWidthX);
    int highCellY = convert_int((photonPosAndTime.y+photonDirAndWlen.y*thisStepLength-this_geoCellStartY)/this_geoCellWidthY);
#endif

    if (highCellX<lowCellX) {int tmp=lowCellX; lowCellX=highCellX; highCellX=tmp;}
    if (highCellY<lowCellY) {int tmp=lowCellY; lowCellY=highCellY; highCellY=tmp;}

    lowCellX = min(max(lowCellX, 0), this_geoCellNumX-1);
    lowCellY = min(max(lowCellY, 0), this_geoCellNumY-1);
    highCellX = min(max(highCellX, 0), this_geoCellNumX-1);
    highCellY = min(max(highCellY, 0), this_geoCellNumY-1);

#ifndef STOP_PHOTONS_ON_DETECTION
    // the number of 64bit integers needed to store bits for all strings
    #define numComponents ((NUM_STRINGS + 64 - 1)/64)
    ulong string_bitmask[numComponents];
    for (uint i=0;i<numComponents;++i) string_bitmask[i]=0;
    #undef numComponents
#endif

    for (int cell_y=lowCellY;cell_y<=highCellY;++cell_y)
    {
        for (int cell_x=lowCellX;cell_x<=highCellX;++cell_x)
        {
            const unsigned short stringNum = this_geoCellIndex[cell_y*this_geoCellNumX+cell_x];
            if (stringNum==0xFFFF) continue; // empty cell
        
#ifndef STOP_PHOTONS_ON_DETECTION
            // prevent strings from being checked twice
            if (string_bitmask[stringNum/64] & (1 << convert_ulong(stringNum%64))) continue;    // already check this string
            string_bitmask[stringNum/64] |= (1 << convert_ulong(stringNum%64));             // mark this string as checked
#endif
        
            checkForCollision_OnString(
                stringNum,
                photonDirLenXYSqr,
                photonPosAndTime,
                photonDirAndWlen,
#ifdef STOP_PHOTONS_ON_DETECTION
                thisStepLength,
                hitRecorded,
                hitOnString,
                hitOnDom,
#else // STOP_PHOTONS_ON_DETECTION
                thisStepLength,
                inv_groupvel,
                photonTotalPathLength,
                photonNumScatters,
                distanceTraveledInAbsorptionLengths,
                photonStartPosAndTime,
                photonStartDirAndWlen,
                step,
                hitIndex,
                maxHitIndex,
                outputPhotons,
#ifdef SAVE_PHOTON_HISTORY
                photonHistory,
                currentPhotonHistory,
#endif // SAVE_PHOTON_HISTORY
#endif // STOP_PHOTONS_ON_DETECTION
                geoLayerToOMNumIndexPerStringSetLocal
                );
        }
    }
    
}

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
    floating_t distanceTraveledInAbsorptionLengths,
    const floating4_t photonStartPosAndTime,
    const floating4_t photonStartDirAndWlen,
    const struct I3CLSimStep *step,
    __global uint* hitIndex,
    uint maxHitIndex,
    __global struct I3CLSimPhoton *outputPhotons,
#ifdef SAVE_PHOTON_HISTORY
    __global float4 *photonHistory,
    float4 *currentPhotonHistory,
#endif
#endif
    __local const unsigned short *geoLayerToOMNumIndexPerStringSetLocal
    )
{
    // using macros and hard-coded names is
    // not really the best thing to do here..
    // replace with a loop sometime.
    
#ifdef STOP_PHOTONS_ON_DETECTION
#define DO_CHECK(subdetectorNum)                \
    checkForCollision_InCell(                   \
        photonDirLenXYSqr,                      \
        photonPosAndTime,                       \
        photonDirAndWlen,                       \
        thisStepLength,                         \
        hitRecorded,                            \
        hitOnString,                            \
        hitOnDom,                               \
        geoLayerToOMNumIndexPerStringSetLocal,  \
                                                \
        geoCellIndex_ ## subdetectorNum,        \
        GEO_CELL_START_X_ ## subdetectorNum,    \
        GEO_CELL_START_Y_ ## subdetectorNum,    \
        GEO_CELL_WIDTH_X_ ## subdetectorNum,    \
        GEO_CELL_WIDTH_Y_ ## subdetectorNum,    \
        GEO_CELL_NUM_X_ ## subdetectorNum,      \
        GEO_CELL_NUM_Y_ ## subdetectorNum       \
        );
#else // STOP_PHOTONS_ON_DETECTION
#ifdef SAVE_PHOTON_HISTORY
#define DO_CHECK(subdetectorNum)                \
    checkForCollision_InCell(                   \
        photonDirLenXYSqr,                      \
        photonPosAndTime,                       \
        photonDirAndWlen,                       \
        thisStepLength,                         \
        inv_groupvel,                           \
        photonTotalPathLength,                  \
        photonNumScatters,                      \
        distanceTraveledInAbsorptionLengths,    \
        photonStartPosAndTime,                  \
        photonStartDirAndWlen,                  \
        step,                                   \
        hitIndex,                               \
        maxHitIndex,                            \
        outputPhotons,                          \
        photonHistory,                          \
        currentPhotonHistory,                   \
        geoLayerToOMNumIndexPerStringSetLocal,  \
                                                \
        geoCellIndex_ ## subdetectorNum,        \
        GEO_CELL_START_X_ ## subdetectorNum,    \
        GEO_CELL_START_Y_ ## subdetectorNum,    \
        GEO_CELL_WIDTH_X_ ## subdetectorNum,    \
        GEO_CELL_WIDTH_Y_ ## subdetectorNum,    \
        GEO_CELL_NUM_X_ ## subdetectorNum,      \
        GEO_CELL_NUM_Y_ ## subdetectorNum       \
        );
#else //SAVE_PHOTON_HISTORY
#define DO_CHECK(subdetectorNum)                \
    checkForCollision_InCell(                   \
        photonDirLenXYSqr,                      \
        photonPosAndTime,                       \
        photonDirAndWlen,                       \
        thisStepLength,                         \
        inv_groupvel,                           \
        photonTotalPathLength,                  \
        photonNumScatters,                      \
        distanceTraveledInAbsorptionLengths,    \
        photonStartPosAndTime,                  \
        photonStartDirAndWlen,                  \
        step,                                   \
        hitIndex,                               \
        maxHitIndex,                            \
        outputPhotons,                          \
        geoLayerToOMNumIndexPerStringSetLocal,  \
                                                \
        geoCellIndex_ ## subdetectorNum,        \
        GEO_CELL_START_X_ ## subdetectorNum,    \
        GEO_CELL_START_Y_ ## subdetectorNum,    \
        GEO_CELL_WIDTH_X_ ## subdetectorNum,    \
        GEO_CELL_WIDTH_Y_ ## subdetectorNum,    \
        GEO_CELL_NUM_X_ ## subdetectorNum,      \
        GEO_CELL_NUM_Y_ ## subdetectorNum       \
        );
#endif //SAVE_PHOTON_HISTORY
#endif // STOP_PHOTONS_ON_DETECTION

    // argh..
#if GEO_CELL_NUM_SUBDETECTORS > 0
    DO_CHECK(0);
#endif        

#if GEO_CELL_NUM_SUBDETECTORS > 1
    DO_CHECK(1);
#endif        

#if GEO_CELL_NUM_SUBDETECTORS > 2
    DO_CHECK(2);
#endif        

#if GEO_CELL_NUM_SUBDETECTORS > 3
    DO_CHECK(3);
#endif        

#if GEO_CELL_NUM_SUBDETECTORS > 4
    DO_CHECK(4);
#endif        

#if GEO_CELL_NUM_SUBDETECTORS > 5
    DO_CHECK(5);
#endif        

#if GEO_CELL_NUM_SUBDETECTORS > 6
    DO_CHECK(6);
#endif        

#if GEO_CELL_NUM_SUBDETECTORS > 7
    DO_CHECK(7);
#endif        

#if GEO_CELL_NUM_SUBDETECTORS > 8
    DO_CHECK(8);
#endif        

#if GEO_CELL_NUM_SUBDETECTORS > 9
    #error more than 9 subdetectors are currently not supported.
#endif

#undef DO_CHECK
}

inline bool checkForCollision(const floating4_t photonPosAndTime,
    const floating4_t photonDirAndWlen,
    floating_t inv_groupvel,
    floating_t photonTotalPathLength,
    uint photonNumScatters,
    floating_t distanceTraveledInAbsorptionLengths,
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
#ifdef SAVE_PHOTON_HISTORY
    __global float4 *photonHistory,
   float4 *currentPhotonHistory,
#endif
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
            distanceTraveledInAbsorptionLengths,
            photonStartPosAndTime,
            photonStartDirAndWlen,
            step,
            0,
            0,
            hitIndex,
            maxHitIndex,
            outputPhotons
#ifdef SAVE_PHOTON_HISTORY
          , photonHistory,
            currentPhotonHistory
#endif
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

    checkForCollision_InCells(
        photonDirLenXYSqr,
        photonPosAndTime,
        photonDirAndWlen,
#ifdef STOP_PHOTONS_ON_DETECTION
        thisStepLength,
        &hitRecorded,
        &hitOnString,
        &hitOnDom,
#else // STOP_PHOTONS_ON_DETECTION
        thisStepLength,
        inv_groupvel,
        photonTotalPathLength,
        photonNumScatters,
        distanceTraveledInAbsorptionLengths,
        photonStartPosAndTime,
        photonStartDirAndWlen,
        step,
        hitIndex,
        maxHitIndex,
        outputPhotons,
#ifdef SAVE_PHOTON_HISTORY
        photonHistory,
        currentPhotonHistory,
#endif // SAVE_PHOTON_HISTORY
#endif // STOP_PHOTONS_ON_DETECTION
        geoLayerToOMNumIndexPerStringSetLocal);

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
                distanceTraveledInAbsorptionLengths,
                photonStartPosAndTime,
                photonStartDirAndWlen,
                step,
                hitOnString,
                hitOnDom,
                hitIndex,
                maxHitIndex,
                outputPhotons
#ifdef SAVE_PHOTON_HISTORY
              , photonHistory,
                currentPhotonHistory
#endif
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
