#include <icetray/serialization.h>
#include <icetray/I3Logging.h>
#include <clsim/random_value/I3CLSimRandomValueInterpolatedDistribution.h>

#include <boost/lexical_cast.hpp>

#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include <limits>

I3CLSimRandomValueInterpolatedDistribution::
I3CLSimRandomValueInterpolatedDistribution(const std::vector<double> &x,
                                           const std::vector<double> &y)
:
x_(x),
y_(y),
constantXSpacing_(NAN),
firstX_(NAN)
{
    if (x_.size() != y_.size())
        log_fatal("The \"x\" and \"y\" vectors must have the same size!");
    
    if (x_.size() <= 1)
        log_fatal("At least two entries have to be specified in the vectors passed to I3CLSimRandomValueInterpolatedDistribution().");
}

I3CLSimRandomValueInterpolatedDistribution::
I3CLSimRandomValueInterpolatedDistribution(double xFirst, double xSpacing,
                                           const std::vector<double> &y)
:
y_(y),
constantXSpacing_(xSpacing),
firstX_(xFirst)
{
    if (isnan(firstX_)) log_fatal("\"xFirst\" must not be NaN!");
    if (isnan(constantXSpacing_)) log_fatal("\"xSpacing\" must not be NaN!");
    if (constantXSpacing_<=0.) log_fatal("\"xSpacing\" must not be <= 0!");
    
    if (y_.size() <= 1)
        log_fatal("At least two entries have to be specified in the vector passed to I3CLSimRandomValueInterpolatedDistribution().");
}


I3CLSimRandomValueInterpolatedDistribution::~I3CLSimRandomValueInterpolatedDistribution() 
{ 
}

I3CLSimRandomValueInterpolatedDistribution::I3CLSimRandomValueInterpolatedDistribution() {;}

std::string I3CLSimRandomValueInterpolatedDistribution::WriteTableCode(const std::string &prefix) const
{
    typedef std::vector<double>::size_type sizeType;
    
    // sanity checks
    if (isnan(constantXSpacing_)) {
        if (x_.size()!=y_.size()) log_fatal("Internal error: angles_.size()!=values_.size()");
    }
    
    sizeType numEntries=y_.size();
    if (numEntries<=1) log_fatal("Internal error: insufficient number of entries.");
    
    
    // intermediate data containers
    std::vector<double> data_acu(numEntries);
    std::vector<double> data_beta(numEntries);

    // simple trapezoidal integration (which is not an approximation in this case)
    data_acu[0]=0.;
    if (isnan(constantXSpacing_))
    {
        for (std::size_t j=1;j<numEntries;++j)
        {
            data_acu[j] = data_acu[j-1] + (x_[j]-x_[j-1]) * (y_[j]+y_[j-1]) / 2.;
        }
    } 
    else
    {
        for (std::size_t j=1;j<numEntries;++j)
        {
            data_acu[j] = data_acu[j-1] + (constantXSpacing_) * (y_[j]+y_[j-1]) / 2.;
        }
    }
        
    // normalize
    for (std::size_t j=0;j<numEntries;++j)
    {
        data_beta[j] = y_[j]/data_acu[numEntries-1];
        data_acu[j] = data_acu[j]/data_acu[numEntries-1];
    }
    
    
    
    // prepare the output buffer
    std::ostringstream output(std::ostringstream::out);
    
    // write the output buffer
    output << "// this is auto-generated code written by I3CLSimRandomValueInterpolatedDistribution::WriteTableCode()" << std::endl;
    output << std::endl;
    
    
    output << "#define " << prefix << "NUM_DIST_ENTRIES " << numEntries << std::endl;
    output << std::endl;
    
    output.setf(std::ios::scientific,std::ios::floatfield);
    output.precision(std::numeric_limits<float>::digits10+4); // maximum precision for a float
    
    if (isnan(constantXSpacing_))
    {
        output << "__constant float " << prefix << "distXValues[" << prefix << "NUM_DIST_ENTRIES] = {" << std::endl;
        for (sizeType j=0;j<numEntries;++j){     
            output << "  " << x_[j] << "f, " << std::endl;
        }
        output << "};" << std::endl;
        output << std::endl;
    }
    
    output << "__constant float " << prefix << "distYValues[" << prefix << "NUM_DIST_ENTRIES] = {" << std::endl;
    for (sizeType j=0;j<numEntries;++j){     
        output << "  " << data_beta[j] << "f, " << std::endl;
    }
    output << "};" << std::endl;
    output << std::endl;

    output << "__constant float " << prefix << "distYCumulativeValues[" << prefix << "NUM_DIST_ENTRIES] = {" << std::endl;
    for (sizeType j=0;j<numEntries;++j){     
        output << "  " << data_acu[j] << "f, " << std::endl;
    }
    output << "};" << std::endl;
    output << std::endl;
    
    // return the code we just generated to the caller
    return output.str();
}


std::string I3CLSimRandomValueInterpolatedDistribution::GetOpenCLFunction
(const std::string &functionName,
 const std::string &functionArgs,
 const std::string &functionArgsToCall,
 const std::string &uniformRandomCall_co,
 const std::string &uniformRandomCall_oc
 ) const
{
    const std::string functionDecl = std::string("inline float ") + functionName + "(" + functionArgs + ")";

    const std::string tableDecl = WriteTableCode(std::string("_") + functionName);
    const std::string distXValuesName = std::string("_") + functionName + "distXValues";
    const std::string distYValuesName = std::string("_") + functionName + "distYValues";
    const std::string distYCumulativeValuesName = std::string("_") + functionName + "distYCumulativeValues";

    std::string constantXSpacingValueString;
    std::string firstXValueString;
    
    if (!isnan(constantXSpacing_)) {
        std::ostringstream output(std::ostringstream::out);
        output.setf(std::ios::scientific,std::ios::floatfield);
        output.precision(std::numeric_limits<float>::digits10+4); // maximum precision for a float
        
        output << constantXSpacing_ << "f";
        constantXSpacingValueString = output.str();
    }

    if (!isnan(firstX_)) {
        std::ostringstream output(std::ostringstream::out);
        output.setf(std::ios::scientific,std::ios::floatfield);
        output.precision(std::numeric_limits<float>::digits10+4); // maximum precision for a float
        
        output << firstX_ << "f";
        firstXValueString = output.str();
    }

    std::string retString = 
    tableDecl + "\n\n" + functionDecl + ";\n\n" + functionDecl + "\n"
    "{\n"
    "    const float randomNumber = " + uniformRandomCall_oc + ";\n"
    "    \n"
    "    unsigned int k=0;\n"
    "    //float this_acu = " + distYCumulativeValuesName + "[0];\n"
    "    float this_acu = 0.f; //this is 0 by definition!\n"
    "    for (;;)\n"
    "    {\n"
    "        float next_acu = " + distYCumulativeValuesName + "[k+1];\n"
    "        if (next_acu >= randomNumber) break;\n"
    "        this_acu = next_acu;\n"
    "        ++k;\n"
    "    }\n"
    "    \n"
    "    // look between bins k and k+1\n"
    "    \n"
    "    const float b = " + distYValuesName + "[k];\n";
    
    if (isnan(constantXSpacing_)) {
        retString = retString + 
        "    const float x0 = " + distXValuesName + "[k];\n"
        "    \n"
        "    const float slope = (" + distYValuesName + "[k+1]-b)/(" + distXValuesName + "[k+1]-x0);\n";
    } else {
        retString = retString + 
        "    const float x0 = convert_float_rtz(k)*(" + constantXSpacingValueString + ") + (" + firstXValueString + "); // _rtz==\"round to zero\"\n"
        "    \n"
        "    const float slope = (" + distYValuesName + "[k+1]-b)/(" + constantXSpacingValueString + ");\n";
    }

    retString = retString + 
    "    const float dy = randomNumber-this_acu;\n"
    "    \n"
    "    \n"
    "    if ((b==0.f) && (slope==0.f))\n"
    "    {\n"
    "        return x0;\n"
    "    }\n"
    "    else if (b==0.f)\n"
    "    {\n"
    "#ifdef USE_NATIVE_MATH\n"
    "        return x0 + native_sqrt(2.f*dy/slope);\n"
    "#else\n"
    "        return x0 + sqrt(2.f*dy/slope);\n"
    "#endif\n"
    "    }\n"
    "    else if (slope==0.f)\n"
    "    {\n"
    "        return x0 + dy/b;\n"
    "    }\n"
    "    else\n"
    "    {\n"
    "#ifdef USE_NATIVE_MATH\n"
    "        return x0 + (native_sqrt(dy * (2.f*slope)/pown(b,2) + 1.f)-1.f)*b/slope;\n"
    "#else\n"
    "        return x0 + (sqrt(dy * (2.f*slope)/pown(b,2) + 1.f)-1.f)*b/slope;\n"
    "#endif\n"
    "    }\n"
    "}\n"
    ;
    
    return retString;
}

bool I3CLSimRandomValueInterpolatedDistribution::CompareTo(const I3CLSimRandomValue &other) const
{
    try
    {
        const I3CLSimRandomValueInterpolatedDistribution &other_ = dynamic_cast<const I3CLSimRandomValueInterpolatedDistribution &>(other);

        if (isnan(other_.constantXSpacing_) && (!isnan(constantXSpacing_))) return false;
        if ((!isnan(other_.constantXSpacing_)) && isnan(constantXSpacing_)) return false;

        if (isnan(other_.firstX_) && (!isnan(firstX_))) return false;
        if ((!isnan(other_.firstX_)) && isnan(firstX_)) return false;

        if ( (!isnan(constantXSpacing_)) && (other_.constantXSpacing_ != constantXSpacing_)) return false;
        if ( (!isnan(firstX_)) && (other_.firstX_ != firstX_)) return false;

        if (other_.x_.size() != x_.size()) return false;
        if (other_.y_.size() != y_.size()) return false;
        
        for (std::size_t i = 0; i < x_.size(); ++i)
        {
            if (x_[i] != other_.x_[i]) return false;
            if (y_[i] != other_.y_[i]) return false;
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
void I3CLSimRandomValueInterpolatedDistribution::serialize(Archive &ar, unsigned version)
{
    if (version>i3clsimrandomvalueinterpolateddistribution_version_)
        log_fatal("Attempting to read version %u from file but running version %u of I3CLSimRandomValueInterpolatedDistribution class.",
                  version,
                  i3clsimrandomvalueinterpolateddistribution_version_);

    ar & make_nvp("I3CLSimRandomValue", base_object<I3CLSimRandomValue>(*this));
    ar & make_nvp("x", x_);
    ar & make_nvp("y", y_);
    ar & make_nvp("constantXSpacing", constantXSpacing_);
    ar & make_nvp("firstX", firstX_);
}


I3_SERIALIZABLE(I3CLSimRandomValueInterpolatedDistribution);
