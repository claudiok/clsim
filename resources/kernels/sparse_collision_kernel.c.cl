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

        floating_t smin1 = INFINITY;
#ifndef OM_HEIGHT // check if DOM/pDOM or GEN2 OM
        {
            const floating4_t drvec = (const floating4_t)(domPosX - photonPosAndTime.x,
                                                          domPosY - photonPosAndTime.y,
                                                          domPosZ - photonPosAndTime.z,
                                                          ZERO);
            const floating_t dr2 = dot(drvec,drvec);

            const floating_t urdot = dot(drvec, photonDirAndWlen); // this assumes drvec.w==0
            floating_t discr   = sqr(urdot) - dr2 + OM_RADIUS*OM_RADIUS;   // (discr)^2  // new constant introduced 

            if (discr < ZERO) continue; // no intersection with this DOM
#ifdef PANCAKE_FACTOR        
            discr = my_sqrt(discr)/PANCAKE_FACTOR;
#else
            discr = my_sqrt(discr);
#endif   
            // I removed smin2 because if smin1 < ZERO => either smin2 < ZERO or flasher photon starting in om in both cases a "continue" follows. 
            // The only exception should be tangential photons that are very improbable 1/(floatmax) and should be suppressed by the OM angular acceptance anyway
            // However if any unexpected behavior occurs try uncommenting the following block
            /*{// by construction: smin1 < smin2
                // distance from current point along the track to second intersection
                // smin2 > 0 && smin1 < 0 means that there *is* an intersection, but we are starting inside the DOM. 
                // This allows photons starting inside a DOM to leave (necessary for flashers):
                const floating_t smin2 = urdot + discr;
                if (smin2 < ZERO) continue; // implies smin1 < 0, so no intersection
            }*/

            // distance from current point along the track to first intersection
            smin1 = urdot - discr;
        }
#else
        {
            //Calculate intersection for cylindrical modules::
            //spheres calculated as before but twice with +/- height/2 offset
                //get smin1 each
            //for the cylinder side area calculate the intersection in the horizontal plane and ignore height for starters
                //also get smin1 but only keep photons where where photonDirAndWlen.z*smin1 < height/2. && photonDirAndWlen.z*smin > -height/2 relative to OM position.
            // keep the shortest (earliest) smin1 intersection
            // FYI smin 2 is actually not needed and has been excluded for simplification and speed
            floating_t urdot_us, urdot_ls, urdot_m, discr_us, discr_ls, discr_m, h_norm;
            {//upper sphere
                const floating4_t drvec_us = (const floating4_t)(domPosX - photonPosAndTime.x,
                                                              domPosY - photonPosAndTime.y,
                                                              domPosZ - photonPosAndTime.z + OM_HEIGHT/2.,
                                                              ZERO);
                const floating_t dr2_us = dot(drvec_us, drvec_us);

                urdot_us = dot(drvec_us, photonDirAndWlen); // this assumes drvec.w==0 / actually it does not need to be 0, it is ignored anyway 
                discr_us   = sqr(urdot_us) - dr2_us + OM_RADIUS*OM_RADIUS;   // (discr)^2  // new constant introduced 
            }
            {//lower sphere
                const floating4_t drvec_ls = (const floating4_t)(domPosX - photonPosAndTime.x,
                                                              domPosY - photonPosAndTime.y,
                                                              domPosZ - photonPosAndTime.z - OM_HEIGHT/2.,
                                                              ZERO);
                const floating_t dr2_ls = dot(drvec_ls, drvec_ls);

                urdot_ls = dot(drvec_ls, photonDirAndWlen); // this assumes drvec.w==0 
                discr_ls   = sqr(urdot_ls) - dr2_ls + OM_RADIUS*OM_RADIUS;   // (discr)^2  // new constant introduced 
            }
            {//cylinder (side area) / mantle (it's a German thing ;) )
                const floating4_t drvec_m = (const floating4_t)(domPosX - photonPosAndTime.x,
                                                              domPosY - photonPosAndTime.y,
                                                              ZERO, //assume infinite long mantle for starters
                                                              ZERO);
                const floating_t dr2_m = dot(drvec_m, drvec_m);

                h_norm = my_sqrt(sqr(photonDirAndWlen.x) + sqr(photonDirAndWlen.y)); //length of the two dimensional direction vector
                urdot_m = (h_norm > ZERO )  ? dot(drvec_m, photonDirAndWlen/h_norm) : ZERO ; // avoid division by zero, h_norm==0 means straight up going photon thus no intersection with cylinder side area anyway / this assumes drvec.w==0 
                discr_m   = sqr(urdot_m) - dr2_m + OM_RADIUS*OM_RADIUS;   // (discr)^2 
            }
            if (discr_us < ZERO && discr_ls < ZERO && discr_m < ZERO) continue; // no intersection with the upper sphere and lower sphere and infinite mantle aka no intersection with this OM

            //how can I distribute the following to 3 different worker groups?
            if (discr_us > ZERO ){
                discr_us = my_sqrt(discr_us);
                floating_t smin_temp  = urdot_us - discr_us;
                smin1 = (smin_temp < smin1 )  ? smin_temp : smin1 ;
                if (smin1 < ZERO) continue; // only first intersection needed
            }            
            if (discr_ls > ZERO ){
                discr_ls = my_sqrt(discr_ls);
                floating_t smin_temp  = urdot_ls - discr_ls;
                smin1 = (smin_temp < smin1 )  ? smin_temp : smin1 ; // take the shortest distance until intersection
                if (smin1 < ZERO) continue; // only first intersection needed
            }
            if (discr_m > ZERO ){
                discr_m = my_sqrt(discr_m);
                floating_t smin_temp  = (urdot_m - discr_m)/h_norm; // calculating distance for 3 Dimensions
                smin1 = (smin_temp < smin1 
                            && (domPosZ - photonPosAndTime.z - photonDirAndWlen.z*smin_temp) < OM_HEIGHT/2. //check lower constrain
                            && (domPosZ - photonPosAndTime.z - photonDirAndWlen.z*smin_temp) > -OM_HEIGHT/2. ) //check upper constrain / height constrain for mantle intersection
                            ? smin_temp : smin1 ;
                if (smin1 < ZERO) continue; // only first intersection needed
            }
        }
#endif

        //by cosntruction flasher must be in the upper or lower sphere, if they are within the cylinder but not any of the spheres the code needs to be modified!!!
        if (smin1 < ZERO || smin1 == INFINITY) continue; // reject if no hit or inside any of the 3 geometries 2 spheres 1 cylinder

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
