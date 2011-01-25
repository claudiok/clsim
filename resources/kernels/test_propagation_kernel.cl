

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
__constant float av_cos = 0.1f;

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



// fast
float makeRayleighScatteringCosAngle(float randomNumber)
{
	const float b = 0.835f;
	//const float p = 1.f/b;
	const float p = 1.f/0.835f;

	const float q = (b+3.f)*randomNumber-0.5f/b;
	const float d = q*q + p*p*p;

	const float u1 = -q+my_sqrt(d);
	const float u = my_powr((fabs(u1)),(1.f/3.f)) * sign(u1);
	//if (u1<0.f) u = -u;

	const float v1 = -q-my_sqrt(d);
	const float v = my_powr((fabs(v1)),(1.f/3.f)) * sign(v1);
	//if (v1 < 0.f) v = -v;

	return clamp(u+v, -1.f, 1.f);
	//return max(-1.f, min(1.f, u+v));
}

/*
// slow (whiles and ifs..)
float makeParticScatteringCosAngle(float randomNumber)
{
	unsigned int k=0;
	while (randomNumber > scatAngleDistParticIntValue[k]) {k++;}

	float angulo;

	if (k==0){
		angulo = randomNumber * scatAngleDistParticAngles[0] / scatAngleDistParticIntValue[0];
	} else {
		angulo = scatAngleDistParticAngles[k-1] + (randomNumber-scatAngleDistParticIntValue[k-1])*
		(scatAngleDistParticAngles[k] - scatAngleDistParticAngles[k-1] )/
		(scatAngleDistParticIntValue[k] - scatAngleDistParticIntValue[k-1]);
	}

	return my_cos(angulo);
}
*/

// Henyey-Greenstein with <cos theta>==0.924
inline float makeParticScatteringCosAngle(float random)
{
    const float g = 0.92400000000000004f;
    const float g2 = 0.85377600000000009f;
    
    // a random number [-1;+1]
    const float s = 2.f*(random)-1.f;
    
    const float ii = ((1.f - g2)/(1.f + g*s));
    return clamp((1.f + g2 - ii*ii) / (2.f*g), -1.f, 1.f);
}

//inline void swap(float *a, float *b) {float tmp=*a; *a=*b; *b=tmp;}
inline float sqr(float a) {return a*a;}


void scatterDirectionByAngle(float cosa,
									float sina,
									float4 *direction,
									float randomNumber)
{
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
}

float gammaDistributedNumber(float shape, 
							 ulong* rnd_x,
							 uint* rnd_a)
{
	float x;
	if(shape<1.f){  // Weibull algorithm
		//float c=1.f/shape;
		float c=my_recip(shape);
		float d=(1.f-shape)*my_powr(shape, my_divide(shape,(1.f-shape)));
		float z, e;
		do
		{
			z=-my_log(rand_MWC_oc(rnd_x, rnd_a));
			e=-my_log(rand_MWC_oc(rnd_x, rnd_a));
			x=my_powr(z, c);
		} while(z+e<d+x); // or here
	}
	else  // Cheng's algorithm
	{
		float b=shape-my_log(4.0f);
		float l=my_sqrt(2.f*shape-1.0f);
		const float cheng=1.0f+my_log(4.5f);

		//float u, v;
		float y, z, r;
		do
		{
			const float rx = rand_MWC_oc(rnd_x, rnd_a);
			const float ry = rand_MWC_oc(rnd_x, rnd_a);

			y=my_divide(my_log(my_divide(ry,(1.f-ry))),l);
			x=shape*my_exp(y);
			z=rx*ry*ry;
			r=b+(shape+l)*y-x;
		} while(r<4.5f*z-cheng && r<my_log(z));
	}
	
	return x;
}


inline float2 getGammaConstAndShape(const int type, // 0==em. shower, 1==hadronic shower
									const float energy)
{
	float2 constAndShape;

	//const float meanNumberOfPhotonsAt100GeV = 17443625.237f;
	//const float longitudinalShapeConstant = 1.91241153282223f;
	//const float longitudinalShapeFactor = 0.619279453443994f;
	//const float longitudinalScaleConstant = 0.54f;

	// from PPC
	//const float ppm=2450.08f;     // photons per meter
	const float m0=0.105658389f;  // muon rest mass [GeV]
	const float rho=0.9216f;      // density of ice [mwe]
	const float Lrad=0.358f/rho;

	// determine shower shape from energy
	const float logE=my_log(max(m0, energy));

	//float longitudinalShape, longitudinalConstant;
	if(type>0){  // hardonic shower
		constAndShape.y=1.49f+0.359f*logE;  // "a"
		constAndShape.x=Lrad/0.772f;     // "b"
	}else{  // em shower
		constAndShape.y=2.03f+0.604f*logE;  // "a"
		constAndShape.x = Lrad/0.633f;   // "b"
	}

	return constAndShape;
}


inline int findLayerForGivenZPos(float posZ)
{
	return convert_int((posZ-(float)MEDIUM_LAYER_BOTTOM_POS)/(float)MEDIUM_LAYER_THICKNESS);
}

inline float mediumLayerBoundary(int layer)
{
    return (convert_float(layer)*((float)MEDIUM_LAYER_THICKNESS)) + (float)MEDIUM_LAYER_BOTTOM_POS;
}


inline void createPhotonShower(const float4 sourcePosAndTime,
						 const float4 sourceDirAndEnergy,
						 //const float cosCherenkov,
						 //const float sinCherenkov,
						 const float shiftAlongTrack,
						 const float4 directionRandoms,
						 float4 *photonPosAndTime,
						 float4 *photonDirAndWlen)
{
	// move along the shower direction
	float sourcePosX = sourcePosAndTime.x+sourceDirAndEnergy.x*shiftAlongTrack;
	float sourcePosY = sourcePosAndTime.y+sourceDirAndEnergy.y*shiftAlongTrack;
	float sourcePosZ = sourcePosAndTime.z+sourceDirAndEnergy.z*shiftAlongTrack;
	float sourceTime = sourcePosAndTime.w+shiftAlongTrack*recip_speedOfLight;

	// determine the photon layer (clamp if necessary)
	unsigned int layer = min(max(findLayerForGivenZPos(sourcePosZ), 0), MEDIUM_LAYERS-1);
	
	// our photon still needs a wavelength. create one!
	(*photonDirAndWlen).w = my_recip(MEDIUM_MIN_PHOTON_ENERGY + directionRandoms.w * (MEDIUM_MAX_PHOTON_ENERGY-MEDIUM_MIN_PHOTON_ENERGY));
	
	const float cosCherenkov = getRecipPhaseRefIndex(layer, (*photonDirAndWlen).w);
	const float sinCherenkov = my_sqrt(1.0f-sqr(cosCherenkov));
	
	// determine the photon direction

	// start with the shower direction
	//float photonDirX = sourceDirAndEnergy.x;
	//float photonDirY = sourceDirAndEnergy.y;
	//float photonDirZ = sourceDirAndEnergy.z;
	(*photonDirAndWlen).xyz = sourceDirAndEnergy.xyz;
	
	// rotate to shower particle direction
	const float a=0.39f, b=2.61f;
	const float I=1.0f-my_exp(-b*my_powr(2.f, a));
	float cs=max(1.0f-my_powr(-my_log(1.f-directionRandoms.x*I)/b, 1.f/a), -1.0f);
	float si=my_sqrt(1.0f-cs*cs);
	scatterDirectionByAngle(cs, si, photonDirAndWlen, directionRandoms.y);

	// and now rotate to cherenkov emission direction
	scatterDirectionByAngle(cosCherenkov, sinCherenkov, photonDirAndWlen, directionRandoms.z);



	(*photonPosAndTime).x = sourcePosX;
	(*photonPosAndTime).y = sourcePosY;
	(*photonPosAndTime).z = sourcePosZ;
	(*photonPosAndTime).w = sourceTime;

}

inline void createPhotonTrack(const float4 sourcePosAndTime,
						 const float4 sourceDirAndLength,
                         const float sourceTimeFinal,
                         ulong *rnd_x,
                         uint *rnd_a,
						 float4 *photonPosAndTime,
						 float4 *photonDirAndWlen)
{
    float shiftAlongTrack = rand_MWC_co(rnd_x,rnd_a);
    float shiftMultiplied = sourceDirAndLength.w*shiftAlongTrack;

	// move along the shower direction
	float sourcePosX = sourcePosAndTime.x+sourceDirAndLength.x*shiftMultiplied;
	float sourcePosY = sourcePosAndTime.y+sourceDirAndLength.y*shiftMultiplied;
	float sourcePosZ = sourcePosAndTime.z+sourceDirAndLength.z*shiftMultiplied;
	float sourceTime = sourcePosAndTime.w+shiftAlongTrack*(sourceTimeFinal-sourcePosAndTime.w);

	// determine the photon layer (clamp if necessary)
	unsigned int layer = min(max(findLayerForGivenZPos(sourcePosZ), 0), MEDIUM_LAYERS-1);
	
	// our photon still needs a wavelength. create one!
	(*photonDirAndWlen).w = my_recip(MEDIUM_MIN_PHOTON_ENERGY + rand_MWC_oc(rnd_x, rnd_a) * (MEDIUM_MAX_PHOTON_ENERGY-MEDIUM_MIN_PHOTON_ENERGY));
	
	const float cosCherenkov = getRecipPhaseRefIndex(layer, (*photonDirAndWlen).w);
	const float sinCherenkov = my_sqrt(1.0f-sqr(cosCherenkov));
	
	// determine the photon direction

	// start with the track direction
	(*photonDirAndWlen).xyz = sourceDirAndLength.xyz;
	
	// and now rotate to cherenkov emission direction
	scatterDirectionByAngle(cosCherenkov, sinCherenkov, photonDirAndWlen, rand_MWC_co(rnd_x, rnd_a));

	(*photonPosAndTime).x = sourcePosX;
	(*photonPosAndTime).y = sourcePosY;
	(*photonPosAndTime).z = sourcePosZ;
	(*photonPosAndTime).w = sourceTime;
}

inline bool checkForCollision(const float4 photonPosAndTime,
							   const float4 photonDirAndWlen,
							   float inv_groupvel, 
							   float *thisStepLength,
							   __global uint* hitIndex,
							   uint maxHitIndex,
							   __global float4* outputRandoms,
							   //__read_only __global float* geoDomPosLocal,
							   //__read_only __global int* geoDomIDsLocal,
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
		
			const unsigned int domIndex = stringNum*MAX_NUM_DOMS_PER_STRINGS+domNum;
			const float4 drvec = (const float4)(photonPosAndTime.x - convert_float(geoDomPosX[domIndex])*(GEO_DOM_POS_MAX_ABS_X/32767.f),
												photonPosAndTime.y - convert_float(geoDomPosY[domIndex])*(GEO_DOM_POS_MAX_ABS_Y/32767.f),
												photonPosAndTime.z - convert_float(geoDomPosZ[domIndex])*(GEO_DOM_POS_MAX_ABS_Z/32767.f),
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
			else if (discr >= 0.0f) 
			{
				discr = my_sqrt(discr);
				
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
			outputRandoms[myIndex*2].x = photonPosAndTime.x+(*thisStepLength)*photonDirAndWlen.x;
			outputRandoms[myIndex*2].y = photonPosAndTime.y+(*thisStepLength)*photonDirAndWlen.y;
			outputRandoms[myIndex*2].z = photonPosAndTime.z+(*thisStepLength)*photonDirAndWlen.z;
			outputRandoms[myIndex*2].w = photonPosAndTime.w+(*thisStepLength)*inv_groupvel;
			outputRandoms[myIndex*2+1].x = photonDirAndWlen.x;
			outputRandoms[myIndex*2+1].y = photonDirAndWlen.y;
			outputRandoms[myIndex*2+1].z = photonDirAndWlen.z;
			outputRandoms[myIndex*2+1].w = convert_int(hitOnString)*1000+convert_int(hitOnDom);
		}
	}	

	return hitRecorded;
}



__kernel void testKernel(const float4 posAndTime,
						 const float4 dirAndEnergy,
						 __global uint* hitIndex,
						 const uint maxHitIndex,
                         __read_only __global float* photonData,
                         //__read_only __global uint* photonMultiplicity,
						 __read_only __global int* domIDs,
						 __read_only __global float* domPositions,
						 __read_only __global unsigned char* geoLayerToOMNumIndexPerStringSet,
						 
						 __write_only __global float4* outputRandoms,
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
	barrier(CLK_LOCAL_MEM_FENCE);

	unsigned int i = get_global_id(0);
	unsigned int global_size = get_global_size(0);


	//download MWC RNG state
	ulong rnd_x = MWC_RNG_x[i];
	uint rnd_a = MWC_RNG_a[i];

    //float4 showerPosAndTime = photonData[2*i];
    //float4 showerDirAndEnergy = photonData[2*i + 1];

    // decode individual photon data
    float4 showerPosAndTime   = (float4)(photonData[10*i+0], photonData[10*i+1], photonData[10*i+2], photonData[10*i+3]);
    float4 showerDirAndEnergy = (float4)(photonData[10*i+4], photonData[10*i+5], photonData[10*i+6], photonData[10*i+7]);
    const float showerTimeFinal = photonData[10*i+8];
    //const float showerWeight = photonData[10*i+9];
    const uint showerPhotons = ((__read_only __global uint*)photonData)[10*i+9];
    
    //const uint showerPhotons = photonMultiplicity[i];

	//float4 showerPosAndTime=posAndTime;
	//float4 showerDirAndEnergy=dirAndEnergy;

	dbg_printf("Vertex at: p=(%f,%f,%f), d=(%f,%f,%f), t=%f, E=%f\n",
		   showerPosAndTime.x, showerPosAndTime.y, showerPosAndTime.z,
		   showerDirAndEnergy.x, showerDirAndEnergy.y, showerDirAndEnergy.z,
		   showerPosAndTime.w, showerDirAndEnergy.w);

	//barrier(CLK_GLOBAL_MEM_FENCE);

//	const float2 gammaConstAndShape = getGammaConstAndShape(0, // 0==em. shower, 1==hadronic shower
//													  showerDirAndEnergy.w);
//	const float shiftAlongTrack = gammaConstAndShape.x*gammaDistributedNumber(gammaConstAndShape.y, &rnd_x, &rnd_a);
	//const float shiftAlongTrack = 0.f;

    //if (showerPhotons != 51)
    //    printf("Shower photons: %u\n", showerPhotons);

	//for (unsigned int jj=0;jj<NUM_PHOTONS_PER_KERNEL;++jj)
	for (unsigned int jj=0;jj<showerPhotons;++jj)
	{
		dbg_printf(" * photon #%u:\n", jj);

		float4 photonPosAndTime;
		float4 photonDirAndWlen;
		unsigned int photonNumScatters;
		float photonTotalPathLength;

		//barrier(CLK_LOCAL_MEM_FENCE);
		//float4 directionRandoms;
		//directionRandoms.x = rand_MWC_co(&rnd_x, &rnd_a);
		//directionRandoms.y = rand_MWC_co(&rnd_x, &rnd_a);
		//directionRandoms.z = rand_MWC_co(&rnd_x, &rnd_a);
		//directionRandoms.w = rand_MWC_oc(&rnd_x, &rnd_a);
		//createPhotonShower(showerPosAndTime, showerDirAndEnergy, shiftAlongTrack, directionRandoms, &photonPosAndTime, &photonDirAndWlen);
		
		//barrier(CLK_LOCAL_MEM_FENCE);

		// create a new photon!
		createPhotonTrack(showerPosAndTime, showerDirAndEnergy, showerTimeFinal, &rnd_x, &rnd_a, &photonPosAndTime, &photonDirAndWlen);
		photonNumScatters=0;
		photonTotalPathLength=0.f;

		dbg_printf("   created photon at: p=(%f,%f,%f), d=(%f,%f,%f), t=%f, wlen=%f\n",
			   photonPosAndTime.x, photonPosAndTime.y, photonPosAndTime.z,
			   photonDirAndWlen.x, photonDirAndWlen.y, photonDirAndWlen.z,
			   photonPosAndTime.w, photonDirAndWlen.w);

		int currentPhotonLayer = findLayerForGivenZPos(photonPosAndTime.z);
		dbg_printf("   in layer %i (valid between 0 and up to including %u)\n", currentPhotonLayer, MEDIUM_LAYERS-1);
		if ((currentPhotonLayer < 0) || (currentPhotonLayer >= MEDIUM_LAYERS)) continue; // outside, do not track

		// the photon needs a lifetime. determine distance to next scatter and absorption
		// (this is in units of absorption/scattering lengths)
		float abs_lens_left = -my_log(rand_MWC_oc(&rnd_x,&rnd_a));
		dbg_printf("   - total track length will be %f absorption lengths\n", abs_lens_left);
		
#define EPSILON 0.00001f

		//for (;;) // main photon propagation loop
		while (abs_lens_left > EPSILON)
		{
			float sca_step_left = -my_log(rand_MWC_oc(&rnd_x,&rnd_a));
			dbg_printf("   - next scatter in %f scattering lengths\n", sca_step_left);
					
			while ( (sca_step_left > EPSILON) && (abs_lens_left > EPSILON) ) 
			{
				dbg_printf("   - stepping...\n");

				// retrieve the absorption and scattering lengths for the current layer
				const float abslen       = getAbsorptionLengthForLayer(currentPhotonLayer, photonDirAndWlen.w);
				const float scatlen      = getScatteringLengthForLayer(currentPhotonLayer, photonDirAndWlen.w);
				const float inv_groupvel = getRecipGroupVelocityForLayer(currentPhotonLayer, photonDirAndWlen.w);
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
												  &thisStepLength, 
												  hitIndex, 
												  maxHitIndex, 
												  outputRandoms, 
												  //domPositions, 
												  //domIDs,
												  //geoLayerToOMNumIndexPerStringSet
												  geoLayerToOMNumIndexPerStringSetLocal
												  );
				if (collided) {
					// get rid of the photon if we detected it
					abs_lens_left = 0.f;
					sca_step_left = 0.f;
					
					steppedToZ=photonDirAndWlen.z*thisStepLength; // this needs to be updated, the old value has probably changed

					dbg_printf("    . colission detected, step limited to thisStepLength=%f, steppedToZ=%f!\n", 
					       thisStepLength, steppedToZ);
				}

				// update the track to its next position
				photonPosAndTime.x += photonDirAndWlen.x*thisStepLength;
				photonPosAndTime.y += photonDirAndWlen.y*thisStepLength;
				photonPosAndTime.z  = steppedToZ; // we already calculated that..
				photonPosAndTime.w += inv_groupvel*thisStepLength;
				photonTotalPathLength += thisStepLength;

				dbg_printf("    . photon position updated: p=(%f,%f,%f), d=(%f,%f,%f), t=%f, wlen=%f\n",
					   photonPosAndTime.x, photonPosAndTime.y, photonPosAndTime.z,
					   photonDirAndWlen.x, photonDirAndWlen.y, photonDirAndWlen.z,
					   photonPosAndTime.w, photonDirAndWlen.w);

				if ((currentPhotonLayer < 0) || (currentPhotonLayer >= MEDIUM_LAYERS))
				{
					// we left the known world. absorb.
					abs_lens_left = 0.f;
					sca_step_left = 0.f;

					dbg_printf("    . photon left the world (upper or lower layer boundary). Killing it!\n");
				}

			} // while()
			
			dbg_printf("   - step performed! abs_lens_left=%f, sca_step_left=%f\n", abs_lens_left, sca_step_left);
			
			// if we got here, the photon was either absorbed or it needs to be scattered.
			
			if (abs_lens_left > EPSILON)
			{
				// it was not absorbed. calculate a new direction
				dbg_printf("   - photon is not yet absorbed (abs_len_left=%f)! Scattering!\n", abs_lens_left);
				
				// TODO: get the scattering model right

				dbg_printf("    . photon direction before: d=(%f,%f,%f), wlen=%f\n",
					   photonDirAndWlen.x, photonDirAndWlen.y, photonDirAndWlen.z,
					   photonDirAndWlen.w);
				
                const float rr = rand_MWC_co(&rnd_x, &rnd_a);
				//const bool doRayleigh = (rand_MWC_co(&rnd_x, &rnd_a) < 0.17f);
				//const bool doRayleigh = true;
                // re-use the perfectly good random number (within its range, scaled back to [0,1])
				float cosScatAngle = (rr<0.17f)?makeRayleighScatteringCosAngle(rr/0.17f):makeParticScatteringCosAngle((1.f-rr)/(1.f-0.17f));
				float sinScatAngle = my_sqrt(1.0f - sqr(cosScatAngle));

				scatterDirectionByAngle(cosScatAngle, sinScatAngle, &photonDirAndWlen, rand_MWC_co(&rnd_x, &rnd_a));

				//if (doRayleigh)
				//	dbg_printf("    . doing Rayleigh scattering.\n");
				//else
				//	dbg_printf("    . doing particle scattering.\n");
				
				dbg_printf("    . photon direction after:  d=(%f,%f,%f), wlen=%f\n",
					   photonDirAndWlen.x, photonDirAndWlen.y, photonDirAndWlen.z,
					   photonDirAndWlen.w);
				
				++photonNumScatters;

				dbg_printf("    . the photon has now been scattered %u time(s).\n", photonNumScatters);
			}

		} // while()
		
		dbg_printf(" * photon #%u finished.\n", jj);

	}

	dbg_printf("Kernel finished.\n");

	//barrier(CLK_GLOBAL_MEM_FENCE);


	//upload MWC RNG state
	MWC_RNG_x[i] = rnd_x;
	MWC_RNG_a[i] = rnd_a;
}

