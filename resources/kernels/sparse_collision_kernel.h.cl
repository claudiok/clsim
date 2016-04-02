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
 * @file sparse_collision_kernel.h.cl
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
    );

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
    );


