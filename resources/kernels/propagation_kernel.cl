#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics : enable
#pragma OPENCL EXTENSION cl_khr_local_int32_base_atomics : enable
#pragma OPENCL EXTENSION cl_khr_byte_addressable_store : enable
//#pragma OPENCL EXTENSION cl_khr_fp16 : enable

// debug mode: store all photons right after they are generated on string0/OM0
//#define DEBUG_STORE_GENERATED_PHOTONS

//#define PRINTF_ENABLED

// disable dbg_printf for GPU
//#define dbg_printf(format, ...)

#ifdef PRINTF_ENABLED
// enable printf for CPU
#pragma OPENCL EXTENSION cl_amd_printf : enable
#define dbg_printf(format, ...) printf(format, ##__VA_ARGS__)
#endif

__constant float speedOfLight = 0.299792458f; // [m/ns]
__constant float recip_speedOfLight = 3.33564095f; // [ns/m]
__constant float PI = 3.14159265359f;


#ifdef USE_NATIVE_MATH
inline float my_divide(float a, float b) {return native_divide(a,b);}
inline float my_recip(float a) {return native_recip(a);}
inline float my_powr(float a, float b) {return native_powr(a,b);}
inline float my_sqrt(float a) {return native_sqrt(a);}
inline float my_rsqrt(float a) {return native_rsqrt(a);}
inline float my_cos(float a) {return native_cos(a);}
inline float my_sin(float a) {return native_sin(a);}
inline float my_log(float a) {return native_log(a);}
inline float my_exp(float a) {return native_exp(a);}
#else
inline float my_divide(float a, float b) {return a/b;}
inline float my_recip(float a) {return 1.f/a;}
inline float my_powr(float a, float b) {return powr(a,b);}
inline float my_sqrt(float a) {return sqrt(a);}
inline float my_rsqrt(float a) {return rsqrt(a);}
inline float my_cos(float a) {return cos(a);}
inline float my_sin(float a) {return sin(a);}
inline float my_log(float a) {return log(a);}
inline float my_exp(float a) {return exp(a);}
#endif

inline float sqr(float a) {return a*a;}


struct __attribute__ ((packed)) I3CLSimStep 
{
    float4 posAndTime;   // x,y,z,time
    float4 dirAndLengthAndBeta; // theta,phi,length,beta
    uint numPhotons;
    float weight;
    uint identifier;
    uint dummy;
};

struct __attribute__ ((packed)) I3CLSimPhoton 
{
    float4 posAndTime;   // x,y,z,time
    float2 dir; // theta,phi
    float wavelength; // photon wavelength
    float cherenkovDist; // Cherenkov distance travelled 
    uint numScatters; // number of scatters
    float weight;
    uint identifier;
    short stringID;
    ushort omID;
    float4 startPosAndTime;
    float2 startDir;
    float groupVelocity;
    uint dummy;
};


inline int findLayerForGivenZPos(float posZ)
{
    return convert_int((posZ-(float)MEDIUM_LAYER_BOTTOM_POS)/(float)MEDIUM_LAYER_THICKNESS);
}

inline float mediumLayerBoundary(int layer)
{
    return (convert_float(layer)*((float)MEDIUM_LAYER_THICKNESS)) + (float)MEDIUM_LAYER_BOTTOM_POS;
}

void scatterDirectionByAngle(float cosa,
    float sina,
    float4 *direction,
    float randomNumber)
{
    //printf("direction before=(%f,%f,%f) len^2=%f  -> cos=%f, sin=%f, r=%f\n",
    //       (*direction).x, (*direction).y, (*direction).z,
    //       (*direction).x*(*direction).x + (*direction).y*(*direction).y + (*direction).z*(*direction).z,
    //       cosa, sina, randomNumber);

    // randomize direction of scattering (rotation around old direction axis)
    const float b=2.0f*PI*randomNumber;
    const float cosb=my_cos(b);
    const float sinb=my_sin(b);

    // Rotate new direction into absolute frame of reference 
    const float sinth = my_sqrt(max(0.0f, 1.0f-(*direction).z*(*direction).z));

    if(sinth>0.f){  // Current direction not vertical, so rotate 
        const float4 oldDir = *direction;

        (*direction).x=oldDir.x*cosa-my_divide((oldDir.y*cosb+oldDir.z*oldDir.x*sinb)*sina,sinth);
        (*direction).y=oldDir.y*cosa+my_divide((oldDir.x*cosb-oldDir.z*oldDir.y*sinb)*sina,sinth);
        (*direction).z=oldDir.z*cosa+sina*sinb*sinth;
    }else{         // Current direction is vertical, so this is trivial
        (*direction).x=sina*cosb;
        (*direction).y=sina*sinb;
        (*direction).z=cosa*sign((*direction).z);
    }

    {
        const float recip_length = my_recip(fast_length((float4)((*direction).x, (*direction).y, (*direction).z, 0.0f)));

        (*direction).x *= recip_length;
        (*direction).y *= recip_length;
        (*direction).z *= recip_length;
    }

    //printf("direction after=(%f,%f,%f) len^2=%f\n",
    //       (*direction).x, (*direction).y, (*direction).z,
    //       (*direction).x*(*direction).x + (*direction).y*(*direction).y + (*direction).z*(*direction).z);

}


inline void createPhotonFromTrack(struct I3CLSimStep *step,
    const float4 stepDir,
    RNG_ARGS,
    float4 *photonPosAndTime,
    float4 *photonDirAndWlen)
{
    float shiftMultiplied = step->dirAndLengthAndBeta.z*RNG_CALL_UNIFORM_CO;
    float inverseParticleSpeed = my_recip(speedOfLight*step->dirAndLengthAndBeta.w);

    // move along the step direction
    *photonPosAndTime = (float4)
        (
        step->posAndTime.x+stepDir.x*shiftMultiplied,
        step->posAndTime.y+stepDir.y*shiftMultiplied,
        step->posAndTime.z+stepDir.z*shiftMultiplied,
        step->posAndTime.w+inverseParticleSpeed*shiftMultiplied
        );

    // determine the photon layer (clamp if necessary)
    unsigned int layer = min(max(findLayerForGivenZPos( (*photonPosAndTime).z ), 0), MEDIUM_LAYERS-1);

    // our photon still needs a wavelength. create one!
    const float wavelength = generateWavelength(RNG_ARGS_TO_CALL);

    const float cosCherenkov = my_recip(step->dirAndLengthAndBeta.w*getPhaseRefIndex(layer, wavelength)); // cos theta = 1/(beta*n)
    const float sinCherenkov = my_sqrt(1.0f-cosCherenkov*cosCherenkov);

    // determine the photon direction

    // start with the track direction
    (*photonDirAndWlen).xyz = stepDir.xyz;
    (*photonDirAndWlen).w = wavelength;

    // and now rotate to cherenkov emission direction
    //printf("gen:\n");
    scatterDirectionByAngle(cosCherenkov, sinCherenkov, photonDirAndWlen, RNG_CALL_UNIFORM_CO);
    //printf("endgen.\n");
}

inline float2 sphDirFromCar(float4 carDir)
{
    // Calculate Spherical coordinates from Cartesian
    const float r_inv = my_rsqrt(carDir.x*carDir.x+carDir.y*carDir.y+carDir.z*carDir.z);

    float theta = 0.f;
    if (fabs(carDir.z*r_inv)<=1.f) {
        theta=acos(carDir.z*r_inv);
    } else {
        if (carDir.z<0.f) theta=PI;
    }
    if (theta<0.f) theta+=2.f*PI;

    float phi=atan2(carDir.y,carDir.x);
    if (phi<0.f) phi+=2.f*PI;

    return (float2)(theta, phi);
}


inline void checkForCollision_OnString(
    const unsigned short stringNum,
    const float photonDirLenXYSqr,
    const float4 photonPosAndTime,
    const float4 photonDirAndWlen,
    float *thisStepLength,
    bool *hitRecorded,
    unsigned short *hitOnString,
    unsigned short *hitOnDom,
    __local const unsigned short *geoLayerToOMNumIndexPerStringSetLocal
    )
{
    // find the string set for this string
    unsigned char stringSet = geoStringInStringSet[stringNum];

    { // check intersection with string cylinder
        // only use test if uhat lateral component is bigger than about 0.1 (NEED to check bigger than zero)
        const float smin = my_divide(sqr(((photonPosAndTime.x - convert_float(geoStringPosX[stringNum]))*photonDirAndWlen.y - (photonPosAndTime.y - convert_float(geoStringPosY[stringNum]))*photonDirAndWlen.x)), photonDirLenXYSqr);
        //if (smin > sqr(convert_float(geoStringRadius[stringNum]))) return;  // NOTE: smin == distance squared
        if (smin > sqr(convert_float(GEO_STRING_MAX_RADIUS))) return;  // NOTE: smin == distance squared
    }

    { // check if photon is above or below the string
        if ((photonDirAndWlen.z > 0.f) && (photonPosAndTime.z > geoStringMaxZ[stringNum])) return;
        if ((photonDirAndWlen.z < 0.f) && (photonPosAndTime.z < geoStringMinZ[stringNum])) return;
    }

    // this photon could potentially be hitting an om
    // -> check them all

    int lowLayerZ = convert_int((photonPosAndTime.z-geoLayerStartZ[stringSet])/geoLayerHeight[stringSet]);
    int highLayerZ = convert_int((photonPosAndTime.z+photonDirAndWlen.z*(*thisStepLength)-geoLayerStartZ[stringSet])/geoLayerHeight[stringSet]);
    if (highLayerZ<lowLayerZ) {int tmp=lowLayerZ; lowLayerZ=highLayerZ; highLayerZ=tmp;}
    lowLayerZ = min(max(lowLayerZ, 0), geoLayerNum[stringSet]-1);
    highLayerZ = min(max(highLayerZ, 0), geoLayerNum[stringSet]-1);

    //__constant const unsigned short *geoLayerToOMNumIndex=geoLayerToOMNumIndexPerStringSet + (convert_uint(stringSet)*GEO_LAYER_STRINGSET_MAX_NUM_LAYERS) + lowLayerZ;
    __local const unsigned short *geoLayerToOMNumIndex=geoLayerToOMNumIndexPerStringSetLocal + (convert_uint(stringSet)*GEO_LAYER_STRINGSET_MAX_NUM_LAYERS) + lowLayerZ;
    for (unsigned int layer_z=lowLayerZ;layer_z<=highLayerZ;++layer_z,++geoLayerToOMNumIndex)
    {
        const unsigned short domNum = *geoLayerToOMNumIndex;
        if (domNum==0xFFFF) continue; // empty layer for this string

        float domPosX, domPosY, domPosZ;
        geometryGetDomPosition(stringNum, domNum, &domPosX, &domPosY, &domPosZ);

        const float4 drvec = (const float4)(photonPosAndTime.x - domPosX,
            photonPosAndTime.y - domPosY,
            photonPosAndTime.z - domPosZ,
            0.f);

        const float dr2     = dot(drvec,drvec);
        const float urdot   = dot(drvec, photonDirAndWlen); // this assumes drvec.w==0

        float discr   = sqr(urdot) - dr2 + OM_RADIUS*OM_RADIUS;   // (discr)^2

        //if (dr2 < OM_RADIUS*OM_RADIUS) // start point inside the OM
        //{
        //    *thisStepLength=0.f;
        //
        //    // record a hit
        //    *hitOnString=stringNum;
        //    *hitOnDom=domNum;
        //
        //    *hitRecorded=true;
        //}
        //else
        if (discr >= 0.0f) 
        {
            discr = my_sqrt(discr);

            float smin = -urdot - discr;
            if (smin < 0.0f) smin = -urdot + discr;

            // check if distance to intersection <= thisStepLength; if not then no detection 
            if ((smin >= 0.0f) && (smin < *thisStepLength))
            {
                *thisStepLength=smin; // limit step length

                // record a hit
                *hitOnString=stringNum;
                *hitOnDom=domNum;

                *hitRecorded=true;
                // continue searching, maybe we hit a closer OM..
            }
        }
    }

}

inline void checkForCollision_InCell(
    const float photonDirLenXYSqr,
    const float4 photonPosAndTime,
    const float4 photonDirAndWlen,
    float *thisStepLength,
    bool *hitRecorded,
    unsigned short *hitOnString,
    unsigned short *hitOnDom,
    __local const unsigned short *geoLayerToOMNumIndexPerStringSetLocal,
    
    __constant unsigned short *this_geoCellIndex,
    const float this_geoCellStartX,
    const float this_geoCellStartY,
    const float this_geoCellWidthX,
    const float this_geoCellWidthY,
    const int this_geoCellNumX,
    const int this_geoCellNumY
    )
{
    int lowCellX = convert_int((photonPosAndTime.x-this_geoCellStartX)/this_geoCellWidthX);
    int lowCellY = convert_int((photonPosAndTime.y-this_geoCellStartY)/this_geoCellWidthY);

    int highCellX = convert_int((photonPosAndTime.x+photonDirAndWlen.x*(*thisStepLength)-this_geoCellStartX)/this_geoCellWidthX);
    int highCellY = convert_int((photonPosAndTime.y+photonDirAndWlen.y*(*thisStepLength)-this_geoCellStartY)/this_geoCellWidthY);

    if (highCellX<lowCellX) {int tmp=lowCellX; lowCellX=highCellX; highCellX=tmp;}
    if (highCellY<lowCellY) {int tmp=lowCellY; lowCellY=highCellY; highCellY=tmp;}

    lowCellX = min(max(lowCellX, 0), this_geoCellNumX-1);
    lowCellY = min(max(lowCellY, 0), this_geoCellNumY-1);
    highCellX = min(max(highCellX, 0), this_geoCellNumX-1);
    highCellY = min(max(highCellY, 0), this_geoCellNumY-1);

    for (unsigned int cell_y=lowCellY;cell_y<=highCellY;++cell_y)
    {
        for (unsigned int cell_x=lowCellX;cell_x<=highCellX;++cell_x)
        {
            const unsigned short stringNum = this_geoCellIndex[cell_y*this_geoCellNumX+cell_x];
            if (stringNum==0xFFFF) continue; // empty cell
        
            checkForCollision_OnString(
                stringNum,
                photonDirLenXYSqr,
                photonPosAndTime,
                photonDirAndWlen,
                thisStepLength,
                hitRecorded,
                hitOnString,
                hitOnDom,
                geoLayerToOMNumIndexPerStringSetLocal
                );
        }
    }
    
}

inline void checkForCollision_InCells(
    const float photonDirLenXYSqr,
    const float4 photonPosAndTime,
    const float4 photonDirAndWlen,
    float *thisStepLength,
    bool *hitRecorded,
    unsigned short *hitOnString,
    unsigned short *hitOnDom,
    __local const unsigned short *geoLayerToOMNumIndexPerStringSetLocal
    )
{
    // using macros and hard-coded names is
    // not really the best thing to do here..
    // replace with a loop sometime.
    
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
        );                                      \

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
    
inline bool checkForCollision(const float4 photonPosAndTime,
    const float4 photonDirAndWlen,
    float inv_groupvel,
    float photonTotalPathLength,
    uint photonNumScatters,
    float4 photonStartPosAndTime,
    float4 photonStartDirAndWlen,
    const struct I3CLSimStep *step,
    float *thisStepLength,
    __global uint* hitIndex,
    uint maxHitIndex,
    __write_only __global struct I3CLSimPhoton *outputPhotons,
    __local const unsigned short *geoLayerToOMNumIndexPerStringSetLocal
    )
{
    bool hitRecorded=false;
    unsigned short hitOnString;
    unsigned short hitOnDom;
    
    // check for collisions
    const float photonDirLenXYSqr = sqr(photonDirAndWlen.x) + sqr(photonDirAndWlen.y);
#ifdef DEBUG_STORE_GENERATED_PHOTONS
    hitRecorded=true;
    hitOnString=0;
    hitOnDom=0;

#else
    if (photonDirLenXYSqr <= 0.f) return false;

    checkForCollision_InCells(
        photonDirLenXYSqr,
        photonPosAndTime,
        photonDirAndWlen,
        thisStepLength,
        &hitRecorded,
        &hitOnString,
        &hitOnDom,
        geoLayerToOMNumIndexPerStringSetLocal);
#endif    
    
    if (hitRecorded)
    {
        uint myIndex = atom_inc(hitIndex);
        if (myIndex < maxHitIndex)
        {
#ifdef PRINTF_ENABLED
            dbg_printf("     -> photon record added at position %u.\n",
                myIndex);
#endif

            outputPhotons[myIndex].posAndTime = (float4)
                (
                photonPosAndTime.x+(*thisStepLength)*photonDirAndWlen.x,
                photonPosAndTime.y+(*thisStepLength)*photonDirAndWlen.y,
                photonPosAndTime.z+(*thisStepLength)*photonDirAndWlen.z,
                photonPosAndTime.w+(*thisStepLength)*inv_groupvel
                );

            outputPhotons[myIndex].dir = sphDirFromCar(photonDirAndWlen);
            outputPhotons[myIndex].wavelength = photonDirAndWlen.w;

            outputPhotons[myIndex].cherenkovDist = photonTotalPathLength+(*thisStepLength);
            outputPhotons[myIndex].numScatters = photonNumScatters;
            outputPhotons[myIndex].weight = step->weight / getWavelengthBias(photonDirAndWlen.w);
            outputPhotons[myIndex].identifier = step->identifier;

            outputPhotons[myIndex].stringID = convert_short(hitOnString);
            outputPhotons[myIndex].omID = convert_ushort(hitOnDom);

            outputPhotons[myIndex].startPosAndTime=photonStartPosAndTime;
            outputPhotons[myIndex].startDir = sphDirFromCar(photonStartDirAndWlen);

            outputPhotons[myIndex].groupVelocity = my_recip(inv_groupvel);


#ifdef PRINTF_ENABLED
            dbg_printf("     -> stored photon: p=(%f,%f,%f), d=(%f,%f), t=%f, wlen=%fnm\n",
                outputPhotons[myIndex].posAndTime.x, outputPhotons[myIndex].posAndTime.y, outputPhotons[myIndex].posAndTime.z,
                outputPhotons[myIndex].dir.x, outputPhotons[myIndex].dir.y,
                outputPhotons[myIndex].posAndTime.w, outputPhotons[myIndex].wavelength/1e-9f);
#endif

        }
    }   

    return hitRecorded;
}


__kernel void propKernel(__global uint *hitIndex,   // deviceBuffer_CurrentNumOutputPhotons
    const uint maxHitIndex2,    // maxNumOutputPhotons_
    __read_only __global unsigned short *geoLayerToOMNumIndexPerStringSet,

    __read_only __global struct I3CLSimStep *inputSteps, // deviceBuffer_InputSteps
    __write_only __global struct I3CLSimPhoton *outputPhotons, // deviceBuffer_OutputPhotons

    __global ulong* MWC_RNG_x,
    __global uint* MWC_RNG_a)
{
#ifdef PRINTF_ENABLED
    dbg_printf("Start kernel... (work item %u of %u)\n", get_global_id(0), get_global_size(0));
#endif

    __local unsigned short geoLayerToOMNumIndexPerStringSetLocal[GEO_geoLayerToOMNumIndexPerStringSet_BUFFER_SIZE];

    // copy the geo data to our local memory (this is done by a whole work group in parallel)
    event_t copyFinishedEvent =
        async_work_group_copy(geoLayerToOMNumIndexPerStringSetLocal,
        geoLayerToOMNumIndexPerStringSet, 
        (size_t)GEO_geoLayerToOMNumIndexPerStringSet_BUFFER_SIZE,
        0);
    wait_group_events(1, &copyFinishedEvent);
    //barrier(CLK_LOCAL_MEM_FENCE);

    unsigned int i = get_global_id(0);
    unsigned int global_size = get_global_size(0);

    //download MWC RNG state
    ulong real_rnd_x = MWC_RNG_x[i];
    uint real_rnd_a = MWC_RNG_a[i];
    ulong *rnd_x = &real_rnd_x;
    uint *rnd_a = &real_rnd_a;

    const uint maxHitIndex = get_global_size(0);

    // download the step
    struct I3CLSimStep step = inputSteps[i];

    float4 stepDir;
    {
        const float rho = my_sin(step.dirAndLengthAndBeta.x); // sin(theta)
        stepDir = (float4)(rho*my_cos(step.dirAndLengthAndBeta.y), // rho*cos(phi)
            rho*my_sin(step.dirAndLengthAndBeta.y), // rho*sin(phi)
            my_cos(step.dirAndLengthAndBeta.x),    // cos(phi)
            0.f);
    }

#ifdef PRINTF_ENABLED
    dbg_printf("Step at: p=(%f,%f,%f), d=(%f,%f,%f), t=%f, l=%f, N=%u\n",
        step.posAndTime.x,
        step.posAndTime.y,
        step.posAndTime.z,
        stepDir.x, stepDir.y, stepDir.z,
        step.posAndTime.w,
        step.dirAndLengthAndBeta.z,
        step.numPhotons);
#endif

    #define EPSILON 0.00001f

    uint photonsLeftToPropagate=step.numPhotons;
    float abs_lens_left=0.f;
    
    float4 photonStartPosAndTime;
    float4 photonStartDirAndWlen;
    float4 photonPosAndTime;
    float4 photonDirAndWlen;
    uint photonNumScatters=0;
    float photonTotalPathLength=0.f;
    int currentPhotonLayer=0;

#ifndef FUNCTION_getGroupVelocity_DOES_NOT_DEPEND_ON_LAYER
#error This kernel only works with a constant group velocity (constant per layer)
#endif
    float inv_groupvel=0.f;


    while (photonsLeftToPropagate > 0)
    {
        if (abs_lens_left < EPSILON)
        {
            // create a new photon
            createPhotonFromTrack(&step,
                stepDir,
                RNG_ARGS_TO_CALL,
                &photonPosAndTime,
                &photonDirAndWlen);
            
            // save the start position and time
            photonStartPosAndTime=photonPosAndTime;
            photonStartDirAndWlen=photonDirAndWlen;
            
            photonNumScatters=0;
            photonTotalPathLength=0.f;
            
#ifdef PRINTF_ENABLED
            dbg_printf("   created photon %u at: p=(%f,%f,%f), d=(%f,%f,%f), t=%f, wlen=%fnm\n",
                photonsLeftToPropagate-step.numPhotons,
                photonPosAndTime.x, photonPosAndTime.y, photonPosAndTime.z,
                photonDirAndWlen.x, photonDirAndWlen.y, photonDirAndWlen.z,
                photonPosAndTime.w, photonDirAndWlen.w/1e-9f);
#endif
            
            currentPhotonLayer = min(max(findLayerForGivenZPos(photonPosAndTime.z), 0), MEDIUM_LAYERS-1);
            //currentPhotonLayer = findLayerForGivenZPos(photonPosAndTime.z);
#ifdef PRINTF_ENABLED
            dbg_printf("   in layer %i (valid between 0 and up to including %u)\n", currentPhotonLayer, MEDIUM_LAYERS-1);
#endif
            
            inv_groupvel = my_recip(getGroupVelocity(0, photonDirAndWlen.w));
            
            // the photon needs a lifetime. determine distance to next scatter and absorption
            // (this is in units of absorption/scattering lengths)
            abs_lens_left = -my_log(RNG_CALL_UNIFORM_OC);
            
            //if ((currentPhotonLayer < 0) || (currentPhotonLayer >= MEDIUM_LAYERS)) abs_lens_left=0.f; // outside, do not track
            
#ifdef PRINTF_ENABLED
            dbg_printf("   - total track length will be %f absorption lengths\n", abs_lens_left);
#endif
        }

        // this block is along the lines of the PPC kernel
        float distancePropagated;
        {
            const float photon_dz=photonDirAndWlen.z;
            
            // the "next" medium boundary (either top or bottom, depending on step direction)
            float mediumBoundary = (photon_dz<0.f)?(mediumLayerBoundary(currentPhotonLayer)):(mediumLayerBoundary(currentPhotonLayer)+(float)MEDIUM_LAYER_THICKNESS);

            // track this thing to the next scattering point
            float sca_step_left = -my_log(RNG_CALL_UNIFORM_OC);
#ifdef PRINTF_ENABLED
            dbg_printf("   - next scatter in %f scattering lengths\n", sca_step_left);
#endif
            
            float currentScaLen = getScatteringLength(currentPhotonLayer, photonDirAndWlen.w);
            float currentAbsLen = getAbsorptionLength(currentPhotonLayer, photonDirAndWlen.w);
            
            float ais=( photon_dz*sca_step_left - my_divide((mediumBoundary-photonPosAndTime.z),currentScaLen) )*(1.f/(float)MEDIUM_LAYER_THICKNESS);
            float aia=( photon_dz*abs_lens_left - my_divide((mediumBoundary-photonPosAndTime.z),currentAbsLen) )*(1.f/(float)MEDIUM_LAYER_THICKNESS);

#ifdef PRINTF_ENABLED
            dbg_printf("   - ais=%f, aia=%f, j_initial=%i\n", ais, aia, currentPhotonLayer);
#endif
        
            // propagate through layers
            int j=currentPhotonLayer;
            if(photon_dz<0) {
                for (; (j>0) && (ais<0.f) && (aia<0.f); 
                     mediumBoundary-=(float)MEDIUM_LAYER_THICKNESS,
                     currentScaLen=getScatteringLength(j, photonDirAndWlen.w),
                     currentAbsLen=getAbsorptionLength(j, photonDirAndWlen.w),
                     ais+=my_recip(currentScaLen),
                     aia+=my_recip(currentAbsLen)) --j;
            } else {
                for (; (j<MEDIUM_LAYERS-1) && (ais>0.f) && (aia>0.f);
                     mediumBoundary+=(float)MEDIUM_LAYER_THICKNESS,
                     currentScaLen=getScatteringLength(j, photonDirAndWlen.w),
                     currentAbsLen=getAbsorptionLength(j, photonDirAndWlen.w),
                     ais-=my_recip(currentScaLen),
                     aia-=my_recip(currentAbsLen)) ++j;
            }
        
#ifdef PRINTF_ENABLED
            dbg_printf("   - j_final=%i\n", j);
#endif
        
            float tot;
            if ((currentPhotonLayer==j) || fabs(photon_dz)<EPSILON) {
                distancePropagated=sca_step_left*currentScaLen;
                tot=abs_lens_left*currentAbsLen;
            } else {
                const float recip_photon_dz = my_recip(photon_dz);
                distancePropagated=(ais*((float)MEDIUM_LAYER_THICKNESS)*currentScaLen+mediumBoundary-photonPosAndTime.z)*recip_photon_dz;
                tot=(aia*((float)MEDIUM_LAYER_THICKNESS)*currentAbsLen+mediumBoundary-photonPosAndTime.z)*recip_photon_dz;
            }
            currentPhotonLayer=j;
            
#ifdef PRINTF_ENABLED
            dbg_printf("   - distancePropagated=%f\n", distancePropagated);
#endif
        
            // get overburden for distance
            if (tot<distancePropagated) {
                distancePropagated=tot;
                abs_lens_left=0.f;
            } else {
                abs_lens_left=my_divide(tot-distancePropagated, currentAbsLen);
            }
        }

        // the photon is now either being absorbed or scattered.
        // Check for collisions in its way
#ifdef DEBUG_STORE_GENERATED_PHOTONS
        bool collided;
        if (RNG_CALL_UNIFORM_OC > 0.9)  // prescale: 10%
#else
        bool
#endif
        collided = checkForCollision(photonPosAndTime, 
            photonDirAndWlen, 
            inv_groupvel,
            photonTotalPathLength,
            photonNumScatters,
            photonStartPosAndTime,
            photonStartDirAndWlen,
            &step,
            &distancePropagated, 
            hitIndex, 
            maxHitIndex, 
            outputPhotons, 
            geoLayerToOMNumIndexPerStringSetLocal
            );
            
#ifdef DEBUG_STORE_GENERATED_PHOTONS
        collided = true;
#endif
        if (collided) {
            // get rid of the photon if we detected it
            abs_lens_left = 0.f;

#ifdef PRINTF_ENABLED
            dbg_printf("    . colission detected, step limited to thisStepLength=%f!\n", 
                distancePropagated);
#endif
        }

        // update the track to its next position
        photonPosAndTime.x += photonDirAndWlen.x*distancePropagated;
        photonPosAndTime.y += photonDirAndWlen.y*distancePropagated;
        photonPosAndTime.z += photonDirAndWlen.z*distancePropagated;
        photonPosAndTime.w += inv_groupvel*distancePropagated;
        photonTotalPathLength += distancePropagated;


        // absorb or scatter the photon
        if (abs_lens_left < EPSILON) 
        {
            // photon was absorbed.
            // a new one will be generated at the begin of the loop.
            --photonsLeftToPropagate;
        }
        else
        {
            // photon was NOT absorbed. scatter it and re-start the loop
            
            // calculate a new direction
#ifdef PRINTF_ENABLED
            dbg_printf("   - photon is not yet absorbed (abs_len_left=%f)! Scattering!\n", abs_lens_left);
#endif

#ifdef PRINTF_ENABLED
            dbg_printf("    . photon direction before: d=(%f,%f,%f), wlen=%f\n",
                photonDirAndWlen.x, photonDirAndWlen.y, photonDirAndWlen.z,
                photonDirAndWlen.w/1e-9f);
#endif

            const float cosScatAngle = makeScatteringCosAngle(RNG_ARGS_TO_CALL);
            const float sinScatAngle = my_sqrt(1.0f - sqr(cosScatAngle));

            scatterDirectionByAngle(cosScatAngle, sinScatAngle, &photonDirAndWlen, RNG_CALL_UNIFORM_CO);

#ifdef PRINTF_ENABLED
            dbg_printf("    . cos(scat_angle)=%f sin(scat_angle)=%f\n",
                cosScatAngle, sinScatAngle);
#endif

#ifdef PRINTF_ENABLED
            dbg_printf("    . photon direction after:  d=(%f,%f,%f), wlen=%f\n",
                photonDirAndWlen.x, photonDirAndWlen.y, photonDirAndWlen.z,
                photonDirAndWlen.w/1e-9f);
#endif

            ++photonNumScatters;

#ifdef PRINTF_ENABLED
            dbg_printf("    . the photon has now been scattered %u time(s).\n", photonNumScatters);
#endif
        }


    }

#ifdef PRINTF_ENABLED
    dbg_printf("Stop kernel... (work item %u of %u)\n", i, global_size);
    dbg_printf("Kernel finished.\n");
#endif

    //upload MWC RNG state
    MWC_RNG_x[i] = real_rnd_x;
    MWC_RNG_a[i] = real_rnd_a;
}
