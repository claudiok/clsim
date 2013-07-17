
#include <icetray/I3ServiceBase.h>
#include <photonics-service/I3PhotonicsService.h>

#include "clsim/function/I3CLSimFunction.h"
#include "clsim/I3CLSimSimpleGeometryFromI3Geometry.h"
#include "clsim/I3CLSimStepToPhotonConverterOpenCL.h"
#include "clsim/I3CLSimLightSourceToStepConverterGeant4.h"

class I3CLShim : public I3PhotonicsService, public I3ServiceBase {
public:
	I3CLShim(const I3Context &);
	virtual ~I3CLShim();
	void Configure();
	
	/// This is the only interface we support. All others will throw.
	virtual void GetMeanAmplitudes(const std::vector<std::pair<PhotonicsSource, double> > &, std::vector<I3PhotonicsService::Receiver> &);

	virtual void SelectSource(double &meanPEs,
	    double &emissionPointDistance, double &geoTime,
	    PhotonicsSource const &source) { log_fatal("Unimplemented!"); }
	virtual void GetTimeDelay(double random, double &timeDelay) { log_fatal("Unimplemented!"); }
	virtual void GetProbabilityDensity(double &density, double timeDelay) { log_fatal("Unimplemented!"); }
	virtual bool GetProbabilityQuantiles(double *time_edges, double t_0,
	    double *amplitudes, size_t n_bins) { log_fatal("Unimplemented!"); return false; }
	virtual bool GetMeanAmplitudeGradient(double gradient[6]) { log_fatal("Unimplemented!"); return false; }
	virtual bool GetMeanAmplitudeHessian(double gradient[6], double
	    hessian[6][6]) { log_fatal("Unimplemented!"); return false; }
	virtual bool GetProbabilityQuantileGradients(double *time_edges,
	    double t_0, double gradients[][6], size_t n_bins) { log_fatal("Unimplemented!"); return false; };
	virtual bool GetProbabilityQuantileHessians(double *time_edges,
	    double t_0, double values[], double gradients[][6],
	    double hessians[][6][6], size_t n_bins) { log_fatal("Unimplemented!"); return false; };
	
	virtual const double* const GetTimeRange() const { log_fatal("Unimplemented!"); return NULL; };
	
	// I3ServiceBase should have exposed things like this in the first place
	void SetParameter(const std::string &key, const boost::python::object &value)
	{
		this->configuration_->Set(key, value);
	}

private:
	void CreateGeometry(const std::vector<I3PhotonicsService::Receiver> &);

	boost::python::object geometryFactory_;
	I3CLSimLightSourceToStepConverterPtr stepConverter_;
	I3CLSimStepToPhotonConverterOpenCLPtr photonConverter_;
	/// Parameter: Wavelength acceptance of the (D)OM as a I3WlenDependedValue object.
	I3CLSimFunctionConstPtr wavelengthAcceptance_;

	/// Parameter: Angular acceptance of the (D)OM as a I3WlenDependedValue object.
	I3CLSimFunctionConstPtr angularAcceptance_;
	
	int oversampleFactor_;
	size_t cacheDepth_;
	
	typedef std::pair<std::vector<std::pair<PhotonicsSource, double> >, std::vector<I3PhotonicsService::Receiver> > CacheEntry;
	std::list<CacheEntry> resultCache_;
	uint32_t granularity_;
	std::map<std::pair<int16_t, uint16_t>, size_t> receiverMap_;
	
	SET_LOGGER("I3CLShim");
};