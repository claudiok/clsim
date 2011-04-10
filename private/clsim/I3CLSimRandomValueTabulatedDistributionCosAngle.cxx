#include <icetray/serialization.h>
#include <icetray/I3Logging.h>
#include <clsim/I3CLSimRandomValueTabulatedDistributionCosAngle.h>

#include <boost/lexical_cast.hpp>

#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include <limits>

I3CLSimRandomValueTabulatedDistributionCosAngle::
I3CLSimRandomValueTabulatedDistributionCosAngle(const std::vector<double> &angles,
                                                const std::vector<double> &values,
                                                double powerLawIndexBeforeFirstBin)
:
angles_(angles),
values_(values),
powerLawIndexBeforeFirstBin_(powerLawIndexBeforeFirstBin)
{ 
    if (angles_.size() != values_.size())
        log_fatal("The \"angles\" and \"values\" vectors must have the same size!");
    
    if (angles_.size() <= 1)
        log_fatal("At least two entries have to be specified in the vectors passed to I3CLSimRandomValueTabulatedDistributionCosAngle().");
}

I3CLSimRandomValueTabulatedDistributionCosAngle::~I3CLSimRandomValueTabulatedDistributionCosAngle() 
{ 

}

I3CLSimRandomValueTabulatedDistributionCosAngle::I3CLSimRandomValueTabulatedDistributionCosAngle() {;}

std::string I3CLSimRandomValueTabulatedDistributionCosAngle::WriteScatteringTableCode(const std::string &prefix) const
{
	typedef std::vector<double>::size_type sizeType;
	
	// sanity checks
	if (angles_.size()!=values_.size()) log_fatal("Internal error: angles_.size()!=values_.size()");
	sizeType numEntries=angles_.size();
	if (numEntries<=1) log_fatal("Internal error: insufficient number of entries.");
	
	// intermediate data containers
	std::vector<double> partic_acu(numEntries);
	std::vector<double> partic_ang(numEntries);
	std::vector<double> partic_beta(numEntries);
    
	
	// integrate scattering angle distribution (or "phase function")
	
	// normalize part_beta according to 3.8 (Light and Water, Mobley)
	for (sizeType j=0;j<numEntries;++j){
		partic_ang[j] = angles_[j];
		partic_beta[j] = 2. * M_PI * values_[j] * std::sin(partic_ang[j]);
	}
	
	// Compute first value (1e-9, practically zero) assuming
	// beta = theta**m (m = powerLawIndexBeforeFirstBin_)
	double partic_ang_m1 = 1e-9;
	double partic_beta_m1 = partic_beta[0] * std::pow((partic_ang_m1/partic_ang[0]), powerLawIndexBeforeFirstBin_);
	partic_beta_m1 = 2. * M_PI * partic_beta_m1 * std::sin(partic_ang_m1);
	
	// Integrate angular distribution (trapezoidal rule)
	// int = h*(f0+f1)/2
	partic_acu[0] = (partic_ang[0]-partic_ang_m1)*(partic_beta_m1+partic_beta[0])/2.;
	for (sizeType j=1;j<numEntries;++j){
		partic_acu[j] = partic_acu[j-1] + (partic_ang[j]-partic_ang[j-1]) * (partic_beta[j]+partic_beta[j-1]) / 2.;
	}
	
	// Normalize             
	for (sizeType j=0;j<numEntries;++j){     
		partic_acu[j] = partic_acu[j]/partic_acu[numEntries-1];
	}
	
	
	// prepare the output buffer
	std::ostringstream output(std::ostringstream::out);
	
	// write the output buffer
	output << "// this is auto-generated code written by I3CLSimRandomValueTabulatedDistributionCosAngle::WriteScatteringTableCode()" << std::endl;
	output << std::endl;
	
	
	output << "#define " << prefix << "NUM_SCAT_ANGLE_DIST_ENTRIES " << numEntries << std::endl;
	output << std::endl;
    
	output.setf(std::ios::scientific,std::ios::floatfield);
	output.precision(std::numeric_limits<float>::digits10+4); // maximum precision for a float
    
	output << "__constant float " << prefix << "scatAngleDistParticAngles[" << prefix << "NUM_SCAT_ANGLE_DIST_ENTRIES] = {" << std::endl;
	for (sizeType j=0;j<numEntries;++j){     
		output << "  " << partic_ang[j] << "f, " << std::endl;
	}
	output << "};" << std::endl;
	output << std::endl;
    
	output << "__constant float " << prefix << "scatAngleDistParticIntValue[" << prefix << "NUM_SCAT_ANGLE_DIST_ENTRIES] = {" << std::endl;
	for (sizeType j=0;j<numEntries;++j){     
		output << "  " << partic_acu[j] << "f, " << std::endl;
	}
	output << "};" << std::endl;
	output << std::endl;
	
	// return the code we just wrote to the caller
    return output.str();
}


std::string I3CLSimRandomValueTabulatedDistributionCosAngle::GetOpenCLFunction
(const std::string &functionName,
 const std::string &functionArgs,
 const std::string &functionArgsToCall,
 const std::string &uniformRandomCall_co,
 const std::string &uniformRandomCall_oc
 ) const
{
    const std::string functionDecl = std::string("inline float ") + functionName + "(" + functionArgs + ")";

    const std::string tableDecl = WriteScatteringTableCode(std::string("_") + functionName);
    const std::string scatAngleDistParticAnglesName = std::string("_") + functionName + "scatAngleDistParticAngles";
    const std::string scatAngleDistParticIntValueName = std::string("_") + functionName + "scatAngleDistParticIntValue";
    
    return tableDecl + "\n\n" + functionDecl + "\n"
    "{\n"
    "    const float randomNumber = " + uniformRandomCall_co + ";\n"
    "    \n"
    "    unsigned int k=0;\n"
    "    while (randomNumber > " + scatAngleDistParticIntValueName + "[k]) {k++;}\n"
    "    \n"
    "    float angulo;\n"
    "    \n"
    "    if (k==0){\n"
    "        angulo = randomNumber * " + scatAngleDistParticAnglesName + "[0] / " + scatAngleDistParticIntValueName + "[0];\n"
    "    } else {\n"
    "        angulo = " + scatAngleDistParticAnglesName + "[k-1] + (randomNumber-" + scatAngleDistParticIntValueName + "[k-1])*\n"
    "        (" + scatAngleDistParticAnglesName + "[k] - " + scatAngleDistParticAnglesName + "[k-1] )/\n"
    "        (" + scatAngleDistParticIntValueName + "[k] - " + scatAngleDistParticIntValueName + "[k-1]);\n"
    "    }\n"
    "    \n"
    "#ifdef USE_NATIVE_MATH\n"
    "    return native_cos(angulo);\n"
    "#else\n"
    "    return cos(angulo);\n"
    "#endif\n"
    "}\n"
    ;
}

bool I3CLSimRandomValueTabulatedDistributionCosAngle::CompareTo(const I3CLSimRandomValue &other) const
{
    try
    {
        const I3CLSimRandomValueTabulatedDistributionCosAngle &other_ = dynamic_cast<const I3CLSimRandomValueTabulatedDistributionCosAngle &>(other);

        if (other_.powerLawIndexBeforeFirstBin_ != powerLawIndexBeforeFirstBin_) return false;
        if (other_.angles_.size() != angles_.size()) return false;
        if (other_.values_.size() != values_.size()) return false;
        
        for (std::size_t i = 0; i < angles_.size(); ++i)
        {
            if (angles_[i] != other_.angles_[i]) return false;
            if (values_[i] != other_.values_[i]) return false;
        }
        
        return true;
    }
    catch (const std::bad_cast& e)
    {
        // not of the same type, treat it as non-equal
        return false;
    }
    
}

template <class Archive>
void I3CLSimRandomValueTabulatedDistributionCosAngle::serialize(Archive &ar, unsigned version)
{
    if (version>i3clsimrandomvaluetabulateddistributioncosangle_version_)
        log_fatal("Attempting to read version %u from file but running version %u of I3CLSimRandomValueTabulatedDistributionCosAngle class.",version,i3clsimrandomvaluetabulateddistributioncosangle_version_);

    ar & make_nvp("I3CLSimRandomValue", base_object<I3CLSimRandomValue>(*this));
    ar & make_nvp("angles", angles_);
    ar & make_nvp("values", values_);
    ar & make_nvp("powerLawIndexBeforeFirstBin", powerLawIndexBeforeFirstBin_);
}     


I3_SERIALIZABLE(I3CLSimRandomValueTabulatedDistributionCosAngle);
