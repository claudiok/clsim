
#ifndef CLSIM_I3CLSIMTABULATOR_H_INCLUDED
#define CLSIM_I3CLSIMTABULATOR_H_INCLUDED

#include <phys-services/I3RandomService.h>

#include <clsim/I3Photon.h>
#include <clsim/function/I3CLSimFunction.h>

class I3Position;
class I3Particle;

class I3Tabulator {
public:
	I3Tabulator(const std::vector<std::vector<double> > &binEdges);
	
	virtual ~I3Tabulator();
	
	off_t GetBinIndex(const I3Particle &source, const I3Position &pos, double time) const;
	double GetBinVolume(off_t idx) const;
	
	void RecordPhoton(const I3Particle &source, const I3Photon &photon);
	
private:
	void Normalize();
	
	std::vector<std::vector<double> > binEdges_;
	std::vector<off_t> strides_;
	std::vector<size_t> dims_;
	
	float *values_;
	float *weights_;
	
	double stepLength_;
	double domArea_;
	I3CLSimFunctionConstPtr angularAcceptance_;
	I3CLSimFunctionConstPtr wavelengthAcceptance_;
	I3RandomServicePtr rng_;

};

#endif /* CLSIM_I3CLSIMTABULATOR_H_INCLUDED */
