
#include <icetray/I3ServiceBase.h>
#include <photonics-service/I3PhotonicsService.h>

#include "clsim/function/I3CLSimFunction.h"
#include "clsim/I3CLSimSimpleGeometryFromI3Geometry.h"
#include "clsim/I3CLSimStepToPhotonConverterOpenCL.h"
#include "clsim/I3CLSimLightSourceToStepConverterGeant4.h"

class I3CLShim : public I3BlockPhotonicsService, public I3ServiceBase {
public:
	I3CLShim(const I3Context &);
	virtual ~I3CLShim();
	void Configure();
	
	/// This is the only interface we support. All others will throw.
	virtual void GetMeanAmplitudes(std::vector<LightSource> &sources, const std::vector<Receiver> &receivers);
	
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