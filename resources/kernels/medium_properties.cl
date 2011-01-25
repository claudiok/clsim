// "hard-code" the ice/water properties.
// The nice thing about OpenCL is that we
// compile this stuff at program startup, so the
// host code will be able to change this on the fly :-)

#define MEDIUM_LAYERS 1

#define MEDIUM_WLEN_NUM 17
#define MEDIUM_WLEN_START 290
#define MEDIUM_WLEN_STEP 20

#define MEDIUM_MIN_PHOTON_ENERGY 1.f/(float)(MEDIUM_WLEN_START+MEDIUM_WLEN_NUM*MEDIUM_WLEN_STEP)
#define MEDIUM_MAX_PHOTON_ENERGY 1.f/(float)MEDIUM_WLEN_START

// ANTARES/KM3NeT:
// one layer from -1200m to 1200m
#define MEDIUM_LAYER_BOTTOM_POS -1200.0f
#define MEDIUM_LAYER_THICKNESS   2400.0f

//// IceCube:
//// 123 layers from -600m to 630m
//#define MEDIUM_LAYER_BOTTOM_POS -600.0f
//#define MEDIUM_LAYER_THICKNESS    10.0f

// in meter
__constant float mediumLayerBoundary[MEDIUM_LAYERS+1] = {
	-1200.f, 
	 1200.f
};

//// Maarten:
//const float refind_a0 = 1.3201f;   // offset
//const float refind_a1 =  1.4e-5f;  // dn/dP
//const float refind_a2 = 16.2566f;  // d^1n/(dx)^1
//const float refind_a3 = -4383.0f;  // d^2n/(dx)^2
//const float refind_a4 = 1.1455e6f; // d^3n/(dx)^3
//const float refind_P  = 240.0f;    // ambient pressure [atm]


// Quan&Fry (taken from W. Schuster's thesis):
// salinity in ppt:
#define refind_S   38.44f    
// temperature in degC:
#define refind_T   13.1f     
// ambient pressure [atm]:
#define refind_P   240.0f    
// offset:
#define refind_n0  1.31405f  
#define refind_n1  1.45e-5f
#define refind_n2  1.779e-4f
#define refind_n3  1.05e-6f
#define refind_n4  1.6e-8f
#define refind_n5  2.02e-6f
#define refind_n6  15.868f
#define refind_n7  0.01155f
#define refind_n8  0.00423f
#define refind_n9  4382.f
#define refind_n10 1.1455e6f

// these get used in the calculation:
__constant float refind_a0 = refind_n0+(refind_n2-refind_n3*refind_T+refind_n4*refind_T*refind_T)*refind_S-refind_n5*refind_T*refind_T;
__constant float refind_a1 = refind_n1;
__constant float refind_a2 = refind_n6+refind_n7*refind_S-refind_n8*refind_T;
__constant float refind_a3 = -refind_n9;
__constant float refind_a4 = refind_n10;

inline float getRecipPhaseRefIndex(unsigned int layer, float wavelength)
{
	
#ifdef USE_NATIVE_MATH
	const float x = native_recip(wavelength);
	return native_recip(refind_a0  +  refind_a1*refind_P  +  x*(refind_a2 + x*(refind_a3 + x*refind_a4)));
#else
	const float x = 1.f/wavelength;
	return 1.f/(refind_a0  +  refind_a1*refind_P  +  x*(refind_a2 + x*(refind_a3 + x*refind_a4)));
#endif
}

inline float getDispersionPhase(const float wavelength)
{

#ifdef USE_NATIVE_MATH
	const float x = native_recip(wavelength);
#else
	const float x = 1.f/(wavelength);
#endif
	return -x*x*(refind_a2 + x*(2.0f*refind_a3 + x*3.0f*refind_a4));
}

inline float getRecipGroupVelocityForLayer(unsigned int layer, float wavelength)
{
	const float c_light = 0.299792458f;

	const float n_inv = getRecipPhaseRefIndex(layer, wavelength);
	const float y = getDispersionPhase(wavelength);
		
	return c_light * (1.0f + y*wavelength*n_inv) * n_inv;
}


inline float getScatteringLengthForLayer(unsigned int layer, float wavelength)
{
	const float volumeConcentrationSmallParticles=0.0075f;   // in ppm
	const float volumeConcentrationLargeParticles=0.0075f;   // in ppm
	const float refWlen = 550.f;
	
#ifdef USE_NATIVE_MATH
	const float x = native_divide(refWlen,wavelength);
	const float scatCoeff = 0.0017f * native_powr(x,4.3f) + 
							1.34f  * volumeConcentrationSmallParticles * native_powr(x,1.7f) + 
							0.312f * volumeConcentrationLargeParticles * native_powr(x,0.3f);
	return native_recip(scatCoeff);
#else
	const float x = refWlen/wavelength;
	const float scatCoeff = 0.0017f * powr(x,4.3f) + 
							1.34f  * volumeConcentrationSmallParticles * powr(x,1.7f) + 
							0.312f * volumeConcentrationLargeParticles * powr(x,0.3f);
	return 1.f/scatCoeff;
#endif
}


//#define INTERP_NON_RECIP_ABSSLEN  // interpolate the absorption length instead of its reciprocal
// 1/(absorption length) in 1/meter
__constant float mediumRecipAbsorptionLength[MEDIUM_LAYERS*MEDIUM_WLEN_NUM] = {
// layer 0 (-1200m to 1200m)
//  290nm           310nm           330nm           350nm           370nm           390nm           410nm           430nm           450nm           470nm           490nm           510nm           530nm           550nm           570nm           590nm           610nm
	1.f/4.750413286f, 1.f/7.004812306f, 1.f/9.259259259f, 1.f/14.92537313f, 1.f/20.00000000f, 1.f/26.31578947f, 1.f/34.48275862f, 1.f/43.47826087f, 1.f/50.00000000f, 1.f/62.50000000f, 1.f/58.82352941f, 1.f/50.00000000f, 1.f/29.41176471f, 1.f/17.85714286f, 1.f/16.12903226f, 1.f/8.849557522f, 1.f/4.504504505f,
// layer 1 (1200m to ...)
    // does not exist
};


inline void getInterpolationBinAndFraction(float wavelength, unsigned int layer, 
										   int *bin, float *fraction)
{
	float fbin;
	*fraction = modf((wavelength-(float)MEDIUM_WLEN_START)/(float)MEDIUM_WLEN_STEP, &fbin);

	int ibin=(int)fbin;

	if (ibin<0) {
		ibin=0;
		*fraction=0.f;
	} else if (ibin>=MEDIUM_WLEN_NUM-1) {
		ibin=MEDIUM_WLEN_NUM-2;
		*fraction=1.f;
	}
	
	*bin = layer*MEDIUM_WLEN_NUM+ibin;
}


inline float getAbsorptionLengthForLayer(unsigned int layer, float wavelength)
{
	int bin; float fraction;
	getInterpolationBinAndFraction(wavelength, layer, &bin, &fraction);
	
	return mix(1.f/mediumRecipAbsorptionLength[bin], 
			   1.f/mediumRecipAbsorptionLength[bin+1],
			   fraction);
}

