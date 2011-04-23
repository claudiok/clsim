
#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics : enable
#pragma OPENCL EXTENSION cl_khr_local_int32_base_atomics : enable
#pragma OPENCL EXTENSION cl_khr_byte_addressable_store : enable
//#pragma OPENCL EXTENSION cl_khr_fp16 : enable

// disable dbg_printf for GPU
#define dbg_printf(format, ...)

// enable printf for CPU
//#pragma OPENCL EXTENSION cl_amd_printf : enable
//#define dbg_printf(format, ...) printf(format, ##__VA_ARGS__)

__constant float speedOfLight = 0.299792458f; // [m/ns]
__constant float recip_speedOfLight = 3.33564095f; // [ns/m]
__constant float PI = 3.14159265359f;


#ifdef USE_NATIVE_MATH
inline float my_divide(float a, float b) {return native_divide(a,b);}
inline float my_recip(float a) {return native_recip(a);}
inline float my_powr(float a, float b) {return native_powr(a,b);}
inline float my_sqrt(float a) {return native_sqrt(a);}
inline float my_cos(float a) {return native_cos(a);}
inline float my_sin(float a) {return native_sin(a);}
inline float my_log(float a) {return native_log(a);}
inline float my_exp(float a) {return native_exp(a);}
#else
inline float my_divide(float a, float b) {return a/b;}
inline float my_recip(float a) {return 1.f/a;}
inline float my_powr(float a, float b) {return powr(a,b);}
inline float my_sqrt(float a) {return sqrt(a);}
inline float my_cos(float a) {return cos(a);}
inline float my_sin(float a) {return sin(a);}
inline float my_log(float a) {return log(a);}
inline float my_exp(float a) {return exp(a);}
#endif

inline float sqr(float a) {return a*a;}


struct I3CLSimStep 
{
    float4 posAndTime;   // x,y,z,time
    float4 dirAndLengthAndBeta; // theta,phi,length,beta
    uint numPhotons;
    float weight;
    uint identifier;
    uint dummy;
};

struct I3CLSimPhoton 
{
    float4 posAndTime;   // x,y,z,time
    float2 dir; // theta,phi
    float wavelength; // photon wavelength
    float cherenkovDist; // Cherenkov distance travelled 
    uint numScatters; // number of scatters
    float weight;
    uint identifier;
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
    const float r = my_sqrt(carDir.x*carDir.x+carDir.y*carDir.y+carDir.z*carDir.z);

    float theta = 0.f;
    if (fabs(carDir.z/r)<=1.f) {
        theta=acos(carDir.z/r);
    } else {
        if (carDir.z<0.f) theta=PI;
    }
    if (theta<0.f) theta+=2.f*PI;

    float phi=atan2(carDir.x,carDir.y);
    if (phi<0.f) phi+=2.f*PI;

    return (float2)(theta, phi);
}


inline bool checkForCollision(const float4 photonPosAndTime,
                              const float4 photonDirAndWlen,
                              float inv_groupvel,
                              float photonTotalPathLength,
                              uint photonNumScatters,
                              const struct I3CLSimStep *step,
                              float *thisStepLength,
                              __global uint* hitIndex,
                              uint maxHitIndex,
                              __write_only __global struct I3CLSimPhoton *outputPhotons,
                              __local const unsigned char *geoLayerToOMNumIndexPerStringSetLocal
                              )
{
    bool hitRecorded=false;
    unsigned short hitOnString;
    unsigned char hitOnDom;

    // check for collisions
    const float photonDirLenXYSqr = sqr(photonDirAndWlen.x) + sqr(photonDirAndWlen.y);
    if (photonDirLenXYSqr <= 0.f) return false;

    int lowCellX = convert_int((photonPosAndTime.x-GEO_CELL_START_X)/GEO_CELL_WIDTH_X);
    int lowCellY = convert_int((photonPosAndTime.y-GEO_CELL_START_Y)/GEO_CELL_WIDTH_Y);
    
    int highCellX = convert_int((photonPosAndTime.x+photonDirAndWlen.x*(*thisStepLength)-GEO_CELL_START_X)/GEO_CELL_WIDTH_X);
    int highCellY = convert_int((photonPosAndTime.y+photonDirAndWlen.y*(*thisStepLength)-GEO_CELL_START_Y)/GEO_CELL_WIDTH_Y);
    
    if (highCellX<lowCellX) {int tmp=lowCellX; lowCellX=highCellX; highCellX=tmp;}
    if (highCellY<lowCellY) {int tmp=lowCellY; lowCellY=highCellY; highCellY=tmp;}

    lowCellX = min(max(lowCellX, 0), GEO_CELL_NUM_X-1);
    lowCellY = min(max(lowCellY, 0), GEO_CELL_NUM_X-1);
    highCellX = min(max(highCellX, 0), GEO_CELL_NUM_X-1);
    highCellY = min(max(highCellY, 0), GEO_CELL_NUM_Y-1);
    
    for (unsigned int cell_y=lowCellY;cell_y<=highCellY;++cell_y)
    for (unsigned int cell_x=lowCellX;cell_x<=highCellX;++cell_x)
    {
        const unsigned short stringNum = geoCellIndex[cell_y*GEO_CELL_NUM_X+cell_x];
        if (stringNum==0xFFFF) continue; // empty cell

        // find the string set for this string
        unsigned char stringSet = geoStringInStringSet[stringNum];

        { // check intersection with string cylinder
            // only use test if uhat lateral component is bigger than about 0.1 (NEED to check bigger than zero)
            const float smin = my_divide(sqr(((photonPosAndTime.x - convert_float(geoStringPosX[stringNum]))*photonDirAndWlen.y - (photonPosAndTime.y - convert_float(geoStringPosY[stringNum]))*photonDirAndWlen.x)), photonDirLenXYSqr);
            //if (smin > sqr(convert_float(geoStringRadius[stringNum]))) continue;  // NOTE: smin == distance squared
            if (smin > sqr(convert_float(GEO_STRING_MAX_RADIUS))) continue;  // NOTE: smin == distance squared
        }

        { // check if photon is above or below the string
            if ((photonDirAndWlen.z > 0.f) && (photonPosAndTime.z > geoStringMaxZ[stringNum])) continue;
            if ((photonDirAndWlen.z < 0.f) && (photonPosAndTime.z < geoStringMinZ[stringNum])) continue;
        }
        
        // this photon could potentially be hitting an om
        // -> check them all
        
        int lowLayerZ = convert_int((photonPosAndTime.z-geoLayerStartZ[stringSet])/geoLayerHeight[stringSet]);
        int highLayerZ = convert_int((photonPosAndTime.z+photonDirAndWlen.z*(*thisStepLength)-geoLayerStartZ[stringSet])/geoLayerHeight[stringSet]);
        if (highLayerZ<lowLayerZ) {int tmp=lowLayerZ; lowLayerZ=highLayerZ; highLayerZ=tmp;}
        lowLayerZ = min(max(lowLayerZ, 0), geoLayerNum[stringSet]-1);
        highLayerZ = min(max(highLayerZ, 0), geoLayerNum[stringSet]-1);

        //__constant const unsigned char *geoLayerToOMNumIndex=geoLayerToOMNumIndexPerStringSet + (stringSet*GEO_LAYER_STRINGSET_MAX_NUM_LAYERS) + lowLayerZ;
        __local const unsigned char *geoLayerToOMNumIndex=geoLayerToOMNumIndexPerStringSetLocal + (stringSet*GEO_LAYER_STRINGSET_MAX_NUM_LAYERS) + lowLayerZ;
        for (unsigned int layer_z=lowLayerZ;layer_z<=highLayerZ;++layer_z,++geoLayerToOMNumIndex)
        {
            const unsigned char domNum = *geoLayerToOMNumIndex;
            if (domNum==0xFF) continue; // empty layer for this string
            
            float domPosX, domPosY, domPosZ;
            geometryGetDomPosition(stringNum, domNum, &domPosX, &domPosY, &domPosZ);
            
            const float4 drvec = (const float4)(photonPosAndTime.x - domPosX,
                                                photonPosAndTime.y - domPosY,
                                                photonPosAndTime.z - domPosZ,
                                                0.f);
            
            const float dr2     = dot(drvec,drvec);
            const float urdot   = dot(drvec, photonDirAndWlen); // this assumes drvec.w==0
            
            float discr   = sqr(urdot) - dr2 + OM_RADIUS*OM_RADIUS;   // (discr)^2
            
            if (dr2 < OM_RADIUS*OM_RADIUS) // start point inside the OM
            {
                *thisStepLength=0.f;
            
                // record a hit
                hitOnString=stringNum;
                hitOnDom=domNum;
                
                hitRecorded=true;
            }
            else
            if (discr >= 0.0f) 
            {
                discr = my_sqrt(discr); // /5.f;
                
                float smin = -urdot - discr;
                if (smin < 0.0f) smin = -urdot + discr;
                
                // check if distance to intersection <= thisStepLength; if not then no detection 
                if ((smin >= 0.0f) && (smin < *thisStepLength))
                {
                    *thisStepLength=smin; // limit step length

                    // record a hit
                    hitOnString=stringNum;
                    hitOnDom=domNum;

                    hitRecorded=true;
                    // continue searching, maybe we hit a closer OM..
                }
            }
        }
        
    } // for (strings/cells)

    if (hitRecorded)
    {
        uint myIndex = atom_inc(hitIndex);
        if (myIndex < maxHitIndex)
        {
            dbg_printf("     -> photon record added at position %u.\n",
                       myIndex);

        
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
            
            outputPhotons[myIndex].dummy = convert_int(hitOnString)*1000+convert_int(hitOnDom);
            
            
            dbg_printf("     -> stored photon: p=(%f,%f,%f), d=(%f,%f), t=%f, wlen=%fnm\n",
                        outputPhotons[myIndex].posAndTime.x, outputPhotons[myIndex].posAndTime.y, outputPhotons[myIndex].posAndTime.z,
                        outputPhotons[myIndex].dir.x, outputPhotons[myIndex].dir.y,
                        outputPhotons[myIndex].posAndTime.w, outputPhotons[myIndex].wavelength/1e-9f);

            
        }
    }   

    return hitRecorded;
}



__kernel void propKernel(__global uint *hitIndex,   // deviceBuffer_CurrentNumOutputPhotons
                         const uint maxHitIndex2,    // maxNumOutputPhotons_
                         __read_only __global unsigned char *geoLayerToOMNumIndexPerStringSet,

                         __read_only __global struct I3CLSimStep *inputSteps, // deviceBuffer_InputSteps
                         __write_only __global struct I3CLSimPhoton *outputPhotons, // deviceBuffer_OutputPhotons

                         __global ulong* MWC_RNG_x,
                         __global uint* MWC_RNG_a)
{
    dbg_printf("Start kernel... (work item %u of %u)\n", get_global_id(0), get_global_size(0));

    __local unsigned char geoLayerToOMNumIndexPerStringSetLocal[GEO_geoLayerToOMNumIndexPerStringSet_BUFFER_SIZE];

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

    dbg_printf("Step at: p=(%f,%f,%f), d=(%f,%f,%f), t=%f, l=%f, N=%u\n",
               step.posAndTime.x,
               step.posAndTime.y,
               step.posAndTime.z,
               stepDir.x, stepDir.y, stepDir.z,
               step.posAndTime.w,
               step.dirAndLengthAndBeta.z,
               step.numPhotons);


    for (uint jj=0;jj<step.numPhotons;++jj)
    {
        float4 photonPosAndTime;
        float4 photonDirAndWlen;

        createPhotonFromTrack(&step,
                              stepDir,
                              RNG_ARGS_TO_CALL,
                              &photonPosAndTime,
                              &photonDirAndWlen);

        uint photonNumScatters=0;
        float photonTotalPathLength=0.f;

        dbg_printf("   created photon %u at: p=(%f,%f,%f), d=(%f,%f,%f), t=%f, wlen=%fnm\n",
                   jj,
                   photonPosAndTime.x, photonPosAndTime.y, photonPosAndTime.z,
                   photonDirAndWlen.x, photonDirAndWlen.y, photonDirAndWlen.z,
                   photonPosAndTime.w, photonDirAndWlen.w/1e-9f);

        int currentPhotonLayer = findLayerForGivenZPos(photonPosAndTime.z);
        dbg_printf("   in layer %i (valid between 0 and up to including %u)\n", currentPhotonLayer, MEDIUM_LAYERS-1);
        if ((currentPhotonLayer < 0) || (currentPhotonLayer >= MEDIUM_LAYERS)) continue; // outside, do not track

        // the photon needs a lifetime. determine distance to next scatter and absorption
        // (this is in units of absorption/scattering lengths)
        float abs_lens_left = -my_log(RNG_CALL_UNIFORM_OC);
        dbg_printf("   - total track length will be %f absorption lengths\n", abs_lens_left);
        
        #define EPSILON 0.00001f

        //for (;;) // main photon propagation loop
        while (abs_lens_left > EPSILON)
        {
            float sca_step_left = -my_log(RNG_CALL_UNIFORM_OC);
            dbg_printf("   - next scatter in %f scattering lengths\n", sca_step_left);
            
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef FUNCTION_getGroupVelocity_DOES_NOT_DEPEND_ON_LAYER
            // optimized version for constant group refractive index
            // (collisions are not checked for in every layer)
            const float inv_groupvel = my_recip(getGroupVelocity(0, photonDirAndWlen.w));
            
            float totalStepLength=0.f;
            float steppedToZ=photonPosAndTime.z;
            const float photon_dz=photonDirAndWlen.z;
            
            float boundaryCurrentLayerBottom = mediumLayerBoundary(currentPhotonLayer);
            float boundaryCurrentLayerTop = boundaryCurrentLayerBottom+(float)MEDIUM_LAYER_THICKNESS;
            
            while ( (sca_step_left > EPSILON) && (abs_lens_left > EPSILON) ) 
            {
                // retrieve the absorption and scattering lengths for the current layer
                const float abslen       = getAbsorptionLength(currentPhotonLayer, photonDirAndWlen.w);
                const float scatlen      = getScatteringLength(currentPhotonLayer, photonDirAndWlen.w);
                const float inv_abslen   = my_recip(abslen);
                const float inv_scatlen  = my_recip(scatlen);
                
                float thisStepLength = min(sca_step_left*scatlen, abs_lens_left*abslen);
                
                float zBefore = steppedToZ;
                steppedToZ += thisStepLength*photon_dz; // where did our step end up?
                
                // downward crossing?
                if (steppedToZ < boundaryCurrentLayerBottom)
                {
                    // limit the current step length to the layer boundary
                    steppedToZ = boundaryCurrentLayerBottom;
                    thisStepLength = my_divide((boundaryCurrentLayerBottom-zBefore),photon_dz);
                    --currentPhotonLayer;
                    
                    boundaryCurrentLayerTop-=(float)MEDIUM_LAYER_THICKNESS;
                    boundaryCurrentLayerBottom-=(float)MEDIUM_LAYER_THICKNESS;
                }
                // upward crossing?
                else if (steppedToZ > boundaryCurrentLayerTop)
                {
                    // limit the current step length to the layer boundary
                    steppedToZ = boundaryCurrentLayerTop;
                    thisStepLength = my_divide((boundaryCurrentLayerTop-zBefore),photon_dz);
                    ++currentPhotonLayer;

                    boundaryCurrentLayerTop+=(float)MEDIUM_LAYER_THICKNESS;
                    boundaryCurrentLayerBottom+=(float)MEDIUM_LAYER_THICKNESS;
                }
                // stays within the same layer: do nothing special

                // perform the step
                abs_lens_left -= thisStepLength*inv_abslen;
                sca_step_left -= thisStepLength*inv_scatlen;
                totalStepLength += thisStepLength;
                
                if ((currentPhotonLayer < 0) || (currentPhotonLayer >= MEDIUM_LAYERS))
                {
                    // we left the known world. absorb.
                    abs_lens_left = 0.f;
                    sca_step_left = 0.f;
                }
                
            }
            
            // the photon is now either being absorbed or scattered.
            // Check for collisions in its way

            bool collided = checkForCollision(photonPosAndTime, 
                                              photonDirAndWlen, 
                                              inv_groupvel,
                                              photonTotalPathLength,
                                              photonNumScatters,
                                              &step,
                                              &totalStepLength, 
                                              hitIndex, 
                                              maxHitIndex, 
                                              outputPhotons, 
                                              geoLayerToOMNumIndexPerStringSetLocal
                                              );
            if (collided) {
                // get rid of the photon if we detected it
                abs_lens_left = 0.f;
                sca_step_left = 0.f;
                
                steppedToZ=photonPosAndTime.z+photonDirAndWlen.z*totalStepLength; // this needs to be updated, the old value has probably changed
                
                dbg_printf("    . colission detected, step limited to thisStepLength=%f, steppedToZ=%f!\n", 
                           totalStepLength, steppedToZ);
            }
            
            // update the track to its next position
            photonPosAndTime.x += photonDirAndWlen.x*totalStepLength;
            photonPosAndTime.y += photonDirAndWlen.y*totalStepLength;
            photonPosAndTime.z  = steppedToZ; // we already calculated that..
            photonPosAndTime.w += inv_groupvel*totalStepLength;
            photonTotalPathLength += totalStepLength;
            
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
#else
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
            while ( (sca_step_left > EPSILON) && (abs_lens_left > EPSILON) ) 
            {
                dbg_printf("   - stepping...\n");
                
                // retrieve the absorption and scattering lengths for the current layer
                const float abslen       = getAbsorptionLength(currentPhotonLayer, photonDirAndWlen.w);
                const float scatlen      = getScatteringLength(currentPhotonLayer, photonDirAndWlen.w);
                const float inv_groupvel = my_recip(getGroupVelocity(currentPhotonLayer, photonDirAndWlen.w));
                const float inv_abslen   = my_recip(abslen);
                
                dbg_printf("    . in this layer (%u): abslen=%f, scatlen=%f, 1/c_gr=%f, 1/abslen=%f\n",
                           currentPhotonLayer, abslen, scatlen, inv_groupvel, inv_abslen);
                       
                // determine the current step length in meters
                float thisStepLength = min(sca_step_left*scatlen, abs_lens_left*abslen);
                
                float steppedToZ = photonPosAndTime.z+thisStepLength*photonDirAndWlen.z; // where did our step end up?
                dbg_printf("    . trying a step of %fm (to z=%f)\n", thisStepLength, steppedToZ);
                
                {
                    const float boundaryCurrentLayerBottom = mediumLayerBoundary(currentPhotonLayer);
                    const float boundaryCurrentLayerTop = boundaryCurrentLayerBottom+(float)MEDIUM_LAYER_THICKNESS;
                    
                    // downward crossing?
                    if (steppedToZ < boundaryCurrentLayerBottom)
                    {
                        // limit the current step length to the layer boundary
                        steppedToZ = boundaryCurrentLayerBottom;
                        thisStepLength = my_divide((boundaryCurrentLayerBottom-photonPosAndTime.z),photonDirAndWlen.z);
                        
                        // perform the step
                        abs_lens_left -= thisStepLength*inv_abslen;
                        sca_step_left -= thisStepLength*my_recip(scatlen);
                        
                        --currentPhotonLayer;
                        
                        dbg_printf("      -> just crossed a layer boundary (downwards): now in layer %u\n", currentPhotonLayer);
                        dbg_printf("         step limited: thisStepLength=%f, steppedToZ=%f, abs_lens_left=%f, sca_step_left=%f\n",
                                   thisStepLength, steppedToZ, abs_lens_left, sca_step_left);
                    }
                    // upward crossing?
                    else if (steppedToZ > boundaryCurrentLayerTop)
                    {
                        // limit the current step length to the layer boundary
                        steppedToZ = boundaryCurrentLayerTop;
                        thisStepLength = my_divide((boundaryCurrentLayerTop-photonPosAndTime.z),photonDirAndWlen.z);
                        
                        // perform the step
                        abs_lens_left -= thisStepLength*inv_abslen;
                        sca_step_left -= thisStepLength*my_recip(scatlen);
                        
                        ++currentPhotonLayer;
                        
                        dbg_printf("      -> just crossed a layer boundary (upwards): now in layer %u\n", currentPhotonLayer);
                        dbg_printf("         step limited: thisStepLength=%f, steppedToZ=%f, abs_lens_left=%f, sca_step_left=%f\n",
                                   thisStepLength, steppedToZ, abs_lens_left, sca_step_left);
                    }
                    // stays within the same layer
                    else
                    {
                        // perform the step
                        abs_lens_left -= thisStepLength*inv_abslen;
                        sca_step_left  = 0.f;
                        
                        dbg_printf("      -> we ended up in the same layer! The photon will either be scattered or absorbed now.\n");
                        dbg_printf("         abs_len_left -= %f  =>  abs_len_left=%f\n", thisStepLength*inv_abslen, abs_lens_left);
                        dbg_printf("          step done: thisStepLength=%f, steppedToZ=%f, abs_lens_left=%f, sca_step_left=%f\n",
                                   thisStepLength, steppedToZ, abs_lens_left, sca_step_left);
                    }
                }
                
                bool collided = checkForCollision(photonPosAndTime, 
                                                  photonDirAndWlen, 
                                                  inv_groupvel,
                                                  photonTotalPathLength,
                                                  photonNumScatters,
                                                  &step,
                                                  &thisStepLength, 
                                                  hitIndex, 
                                                  maxHitIndex, 
                                                  outputPhotons, 
                                                  geoLayerToOMNumIndexPerStringSetLocal
                                                  );
                if (collided) {
                    // get rid of the photon if we detected it
                    abs_lens_left = 0.f;
                    sca_step_left = 0.f;
                    
                    steppedToZ=photonPosAndTime.z+photonDirAndWlen.z*thisStepLength; // this needs to be updated, the old value has probably changed
                    
                    dbg_printf("    . colission detected, step limited to thisStepLength=%f, steppedToZ=%f!\n", 
                               thisStepLength, steppedToZ);
                }
                
                // update the track to its next position
                photonPosAndTime.x += photonDirAndWlen.x*thisStepLength;
                photonPosAndTime.y += photonDirAndWlen.y*thisStepLength;
                photonPosAndTime.z  = steppedToZ; // we already calculated that..
                photonPosAndTime.w += inv_groupvel*thisStepLength;
                photonTotalPathLength += thisStepLength;
                
                dbg_printf("    . photon position updated: p=(%f,%f,%f), d=(%f,%f,%f), t=%f, wlen=%fnm\n",
                           photonPosAndTime.x, photonPosAndTime.y, photonPosAndTime.z,
                           photonDirAndWlen.x, photonDirAndWlen.y, photonDirAndWlen.z,
                           photonPosAndTime.w, photonDirAndWlen.w/1e-9f);
                
                if ((currentPhotonLayer < 0) || (currentPhotonLayer >= MEDIUM_LAYERS))
                {
                    // we left the known world. absorb.
                    abs_lens_left = 0.f;
                    sca_step_left = 0.f;
                    
                    dbg_printf("    . photon left the world (upper or lower layer boundary). Killing it!\n");
                }
                
            } // while()
#endif
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

            dbg_printf("   - step performed! abs_lens_left=%f, sca_step_left=%f\n", abs_lens_left, sca_step_left);
            
            // if we got here, the photon was either absorbed or it needs to be scattered.
            
            if (abs_lens_left > EPSILON)
            {
                // it was not absorbed. calculate a new direction
                dbg_printf("   - photon is not yet absorbed (abs_len_left=%f)! Scattering!\n", abs_lens_left);
                
                dbg_printf("    . photon direction before: d=(%f,%f,%f), wlen=%f\n",
                       photonDirAndWlen.x, photonDirAndWlen.y, photonDirAndWlen.z,
                       photonDirAndWlen.w/1e-9f);
                
                const float cosScatAngle = makeScatteringCosAngle(RNG_ARGS_TO_CALL);
                const float sinScatAngle = my_sqrt(1.0f - sqr(cosScatAngle));
                
                scatterDirectionByAngle(cosScatAngle, sinScatAngle, &photonDirAndWlen, RNG_CALL_UNIFORM_CO);

                dbg_printf("    . cos(scat_angle)=%f sin(scat_angle)=%f\n",
                       cosScatAngle, sinScatAngle);
                
                dbg_printf("    . photon direction after:  d=(%f,%f,%f), wlen=%f\n",
                       photonDirAndWlen.x, photonDirAndWlen.y, photonDirAndWlen.z,
                       photonDirAndWlen.w/1e-9f);
                
                ++photonNumScatters;
                
                dbg_printf("    . the photon has now been scattered %u time(s).\n", photonNumScatters);
            }
        
        } // while()
        
        dbg_printf(" * photon #%u finished.\n", jj);
        
    }
    
    
    dbg_printf("Stop kernel... (work item %u of %u)\n", i, global_size);
    dbg_printf("Kernel finished.\n");
    
    //upload MWC RNG state
    MWC_RNG_x[i] = real_rnd_x;
    MWC_RNG_a[i] = real_rnd_a;
}
