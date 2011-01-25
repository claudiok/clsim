#include <icetray/serialization.h>
#include <icetray/I3Units.h>
#include <clsim/I3CLSimMediumProperties.h>

#include <limits>

#define __STDC_FORMAT_MACROS
#include <inttypes.h>

#include <boost/shared_ptr.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/foreach.hpp>

const uint32_t I3CLSimMediumProperties::default_layersNum=1;
const double I3CLSimMediumProperties::default_layersZStart=-5000.*I3Units::m;
const double I3CLSimMediumProperties::default_layersHeight=10000.*I3Units::m;
const double I3CLSimMediumProperties::default_mediumDensity=1.*I3Units::g/I3Units::cm3;
const double I3CLSimMediumProperties::default_rockZCoordinate=I3CLSimMediumProperties::default_layersZStart;
const double I3CLSimMediumProperties::default_airZCoordinate=-I3CLSimMediumProperties::default_layersZStart;

I3CLSimMediumProperties::I3CLSimMediumProperties(double mediumDensity,
                                                 uint32_t layersNum, 
                                                 double layersZStart, 
                                                 double layersHeight,
                                                 double rockZCoordinate,
                                                 double airZCoordinate)
:
mediumDensity_(mediumDensity),
layersNum_(layersNum),
layersZStart_(layersZStart),
layersHeight_(layersHeight),
rockZCoordinate_(rockZCoordinate),
airZCoordinate_(airZCoordinate),
forcedMinWlen_(-std::numeric_limits<double>::infinity()),
forcedMaxWlen_(std::numeric_limits<double>::infinity())
{
    if (layersNum_ <= 0) log_fatal("layersNum must be at least 1 during construction of I3CLSimMediumProperties.");
    
    if (layersZStart < rockZCoordinate_)
        log_fatal("layers start inside rock");

    if (layersZStart + static_cast<double>(layersNum)*layersHeight > airZCoordinate_)
        log_fatal("layers end inside air");

    absorptionLength_.assign(layersNum_, I3CLSimWlenDependentValueConstPtr());
    scatteringLength_.assign(layersNum_, I3CLSimWlenDependentValueConstPtr());
    phaseRefractiveIndex_.assign(layersNum_, I3CLSimWlenDependentValueConstPtr());
}

I3CLSimMediumProperties::~I3CLSimMediumProperties() 
{ 

}

double I3CLSimMediumProperties::GetMinWavelength() const
{
    if (!IsReady()) return NAN; // no valid answer if not ready

    double mini=forcedMinWlen_;
    
    BOOST_FOREACH(const I3CLSimWlenDependentValueConstPtr &ptr, absorptionLength_)
    {
        if (!ptr) return NAN;
        const double wlen=ptr->GetMinWlen();
        if (wlen > mini) mini=wlen;
        
    }
    BOOST_FOREACH(const I3CLSimWlenDependentValueConstPtr &ptr, scatteringLength_)
    {
        if (!ptr) return NAN;
        const double wlen=ptr->GetMinWlen();
        if (wlen > mini) mini=wlen;
    }
    BOOST_FOREACH(const I3CLSimWlenDependentValueConstPtr &ptr, phaseRefractiveIndex_)
    {
        if (!ptr) return NAN;
        const double wlen=ptr->GetMinWlen();
        if (wlen > mini) mini=wlen;
    }
    
    return mini;
}

double I3CLSimMediumProperties::GetMaxWavelength() const
{
    if (!IsReady()) return NAN; // no valid answer if not ready

    double maxi=forcedMaxWlen_;

    BOOST_FOREACH(const I3CLSimWlenDependentValueConstPtr &ptr, absorptionLength_)
    {
        if (!ptr) return NAN;
        const double wlen=ptr->GetMaxWlen();
        if (wlen < maxi) maxi=wlen;
        
    }
    BOOST_FOREACH(const I3CLSimWlenDependentValueConstPtr &ptr, scatteringLength_)
    {
        if (!ptr) return NAN;
        const double wlen=ptr->GetMaxWlen();
        if (wlen < maxi) maxi=wlen;
    }
    BOOST_FOREACH(const I3CLSimWlenDependentValueConstPtr &ptr, phaseRefractiveIndex_)
    {
        if (!ptr) return NAN;
        const double wlen=ptr->GetMaxWlen();
        if (wlen < maxi) maxi=wlen;
    }
    
    return maxi;
}

bool I3CLSimMediumProperties::IsReady() const
{
    BOOST_FOREACH(const I3CLSimWlenDependentValueConstPtr &ptr, absorptionLength_)
    {
        if (!ptr) return false;
    }
    BOOST_FOREACH(const I3CLSimWlenDependentValueConstPtr &ptr, scatteringLength_)
    {
        if (!ptr) return false;
    }
    BOOST_FOREACH(const I3CLSimWlenDependentValueConstPtr &ptr, phaseRefractiveIndex_)
    {
        if (!ptr) return false;
    }
    if (!scatteringCosAngleDist_) return false;
    
    return true;
}

const std::vector<I3CLSimWlenDependentValueConstPtr> &I3CLSimMediumProperties::GetAbsorptionLengths() const
{
    return absorptionLength_;
}

const std::vector<I3CLSimWlenDependentValueConstPtr> &I3CLSimMediumProperties::GetScatteringLengths() const
{
    return scatteringLength_;
}

const std::vector<I3CLSimWlenDependentValueConstPtr> &I3CLSimMediumProperties::GetPhaseRefractiveIndices() const
{
    return phaseRefractiveIndex_;
}


I3CLSimWlenDependentValueConstPtr I3CLSimMediumProperties::GetAbsorptionLength(uint32_t layer) const
{
    if (layer >= layersNum_) log_fatal("Invalid layer num: %" PRIu32, layer);
    return absorptionLength_[layer];
}

I3CLSimWlenDependentValueConstPtr I3CLSimMediumProperties::GetScatteringLength(uint32_t layer) const
{
    if (layer >= layersNum_) log_fatal("Invalid layer num: %" PRIu32, layer);
    return scatteringLength_[layer];
}

I3CLSimWlenDependentValueConstPtr I3CLSimMediumProperties::GetPhaseRefractiveIndex(uint32_t layer) const
{
    if (layer >= layersNum_) log_fatal("Invalid layer num: %" PRIu32, layer);
    return phaseRefractiveIndex_[layer];
}

I3CLSimRandomValueConstPtr I3CLSimMediumProperties::GetScatteringCosAngleDistribution() const
{
    return scatteringCosAngleDist_;
}

void I3CLSimMediumProperties::SetAbsorptionLength(uint32_t layer, I3CLSimWlenDependentValueConstPtr ptr)
{
    if (layer >= layersNum_) log_fatal("Invalid layer num: %" PRIu32, layer);
    absorptionLength_[layer]=ptr;
}

void I3CLSimMediumProperties::SetScatteringLength(uint32_t layer, I3CLSimWlenDependentValueConstPtr ptr)
{
    if (layer >= layersNum_) log_fatal("Invalid layer num: %" PRIu32, layer);
    scatteringLength_[layer]=ptr;
}

void I3CLSimMediumProperties::SetPhaseRefractiveIndex(uint32_t layer, I3CLSimWlenDependentValueConstPtr ptr)
{
    if (layer >= layersNum_) log_fatal("Invalid layer num: %" PRIu32, layer);
    phaseRefractiveIndex_[layer]=ptr;
}

void I3CLSimMediumProperties::SetScatteringCosAngleDistribution(I3CLSimRandomValueConstPtr ptr)
{
    scatteringCosAngleDist_=ptr;
}


namespace {
    template <typename T, class Archive>
    void LoadFromArchiveIntoVectorConstPtr(Archive &ar, const std::string &name, std::vector<boost::shared_ptr<const T> > &arg)
    {
        std::vector<shared_ptr<T> > tmp;
        ar >> make_nvp(name.c_str(), tmp);
        arg.assign(tmp.begin(), tmp.end());
    }

    template <typename T, class Archive>
    void LoadFromArchiveIntoConstPtr(Archive &ar, const std::string &name, boost::shared_ptr<const T> &arg)
    {
        shared_ptr<T> tmp;
        ar >> make_nvp(name.c_str(), tmp);
        arg = tmp;
    }
}

template <class Archive>
void I3CLSimMediumProperties::load(Archive &ar, unsigned version)
{
    if (version>i3clsimmediumproperties_version_)
        log_fatal("Attempting to read version %u from file but running version %u of I3CLSimMediumProperties class.",version,i3clsimmediumproperties_version_);
    
    ar >> make_nvp("mediumDensity", mediumDensity_);
    ar >> make_nvp("layersNum", layersNum_);
    ar >> make_nvp("layersZStart", layersZStart_);
    ar >> make_nvp("layersHeight", layersHeight_);

    ar >> make_nvp("rockZCoordinate", rockZCoordinate_);
    ar >> make_nvp("airZCoordinate", airZCoordinate_);
    ar >> make_nvp("forcedMinWlen", forcedMinWlen_);
    ar >> make_nvp("forcedMaxWlen", forcedMaxWlen_);
    
    LoadFromArchiveIntoVectorConstPtr(ar, "absorptionLength", absorptionLength_);
    LoadFromArchiveIntoVectorConstPtr(ar, "scatteringLength", scatteringLength_);
    LoadFromArchiveIntoVectorConstPtr(ar, "phaseRefractiveIndex", phaseRefractiveIndex_);
    LoadFromArchiveIntoConstPtr(ar, "scatteringCosAngleDist", scatteringCosAngleDist_);
}     


template <class Archive>
void I3CLSimMediumProperties::save(Archive &ar, unsigned version) const
{
    ar << make_nvp("mediumDensity", mediumDensity_);
    ar << make_nvp("layersNum", layersNum_);
    ar << make_nvp("layersZStart", layersZStart_);
    ar << make_nvp("layersHeight", layersHeight_);
    
    ar << make_nvp("rockZCoordinate", rockZCoordinate_);
    ar << make_nvp("airZCoordinate", airZCoordinate_);
    ar << make_nvp("forcedMinWlen", forcedMinWlen_);
    ar << make_nvp("forcedMaxWlen", forcedMaxWlen_);
    
    ar << make_nvp("absorptionLength", absorptionLength_);
    ar << make_nvp("scatteringLength", scatteringLength_);
    ar << make_nvp("phaseRefractiveIndex", phaseRefractiveIndex_);
    ar << make_nvp("scatteringCosAngleDist", scatteringCosAngleDist_);
}     


I3_SERIALIZABLE(I3CLSimMediumProperties);
