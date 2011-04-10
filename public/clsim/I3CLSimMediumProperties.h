#ifndef I3CLSIMMEDIUMPROPERTIES_H_INCLUDED
#define I3CLSIMMEDIUMPROPERTIES_H_INCLUDED

#include "icetray/serialization.h"
#include "icetray/I3TrayHeaders.h"

#include <vector>

#include "clsim/I3CLSimWlenDependentValue.h"
#include "clsim/I3CLSimRandomValue.h"

/**
 * @brief Describes the simulation medium (i.e. water or
 * ice) with all its properties like refractive index,
 * absorption length, scattering length, ..
 */
static const unsigned i3clsimmediumproperties_version_ = 0;

struct I3CLSimMediumProperties 
{
    
public:
    static const double default_mediumDensity;
    static const uint32_t default_layersNum;
    static const double default_layersZStart;
    static const double default_layersHeight;
    static const double default_rockZCoordinate;
    static const double default_airZCoordinate;
    
    I3CLSimMediumProperties(double mediumDensity=default_mediumDensity,
                            uint32_t layersNum=default_layersNum, 
                            double layersZStart=default_layersZStart, 
                            double layersHeight=default_layersHeight,
                            double rockZCoordinate=default_rockZCoordinate,
                            double airZCoordinate=default_airZCoordinate);
    ~I3CLSimMediumProperties();

    /**
     * Returns true if media definitions for all layers are set.
     */
    bool IsReady() const;
    
    const std::vector<I3CLSimWlenDependentValueConstPtr> &GetAbsorptionLengths() const;
    const std::vector<I3CLSimWlenDependentValueConstPtr> &GetScatteringLengths() const;
    const std::vector<I3CLSimWlenDependentValueConstPtr> &GetPhaseRefractiveIndices() const;
    const std::vector<I3CLSimWlenDependentValueConstPtr> &GetGroupRefractiveIndicesOverride() const;
    
    I3CLSimWlenDependentValueConstPtr GetAbsorptionLength(uint32_t layer) const;
    I3CLSimWlenDependentValueConstPtr GetScatteringLength(uint32_t layer) const;
    I3CLSimWlenDependentValueConstPtr GetPhaseRefractiveIndex(uint32_t layer) const;
    I3CLSimWlenDependentValueConstPtr GetGroupRefractiveIndexOverride(uint32_t layer) const;
    I3CLSimRandomValueConstPtr GetScatteringCosAngleDistribution() const;

    void SetAbsorptionLength(uint32_t layer, I3CLSimWlenDependentValueConstPtr ptr);
    void SetScatteringLength(uint32_t layer, I3CLSimWlenDependentValueConstPtr ptr);
    void SetPhaseRefractiveIndex(uint32_t layer, I3CLSimWlenDependentValueConstPtr ptr);
    void SetGroupRefractiveIndexOverride(uint32_t layer, I3CLSimWlenDependentValueConstPtr ptr);
    void SetScatteringCosAngleDistribution(I3CLSimRandomValueConstPtr ptr);

    double GetMinWavelength() const;
    double GetMaxWavelength() const;
    
    // getters for constructor arguments
    inline double GetMediumDensity() const {return mediumDensity_;}
    inline uint32_t GetLayersNum() const {return layersNum_;}
    inline double GetLayersZStart() const {return layersZStart_;}
    inline double GetLayersHeight() const {return layersHeight_;}
    inline double GetRockZCoord() const {return rockZCoordinate_;}
    inline double GetAirZCoord() const {return airZCoordinate_;}
    
    // The user can specify a minimum and maximum wavelength.
    // This is necessary if none of the I3CLSimWlenDependentValues
    // specifies any.
    inline double GetForcedMinWlen() const {return forcedMinWlen_;}
    inline double GetForcedMaxWlen() const {return forcedMaxWlen_;}
    inline void SetForcedMinWlen(double val) {forcedMinWlen_=val;}
    inline void SetForcedMaxWlen(double val) {forcedMaxWlen_=val;}
    
    
private:
    double mediumDensity_;
    uint32_t layersNum_;
    double layersZStart_;
    double layersHeight_;
    
    double rockZCoordinate_;
    double airZCoordinate_;
    
    double forcedMinWlen_;
    double forcedMaxWlen_;
    
    std::vector<I3CLSimWlenDependentValueConstPtr> absorptionLength_;
    std::vector<I3CLSimWlenDependentValueConstPtr> scatteringLength_;
    std::vector<I3CLSimWlenDependentValueConstPtr> phaseRefractiveIndex_;
    std::vector<I3CLSimWlenDependentValueConstPtr> groupRefractiveIndexOverride_;
    I3CLSimRandomValueConstPtr scatteringCosAngleDist_;
    
private:
    friend class boost::serialization::access;
    template <class Archive> void load(Archive & ar, unsigned version);
    template <class Archive> void save(Archive & ar, unsigned version) const;
    BOOST_SERIALIZATION_SPLIT_MEMBER();
};

BOOST_CLASS_VERSION(I3CLSimMediumProperties, i3clsimmediumproperties_version_);

I3_POINTER_TYPEDEFS(I3CLSimMediumProperties);

#endif //I3CLSIMMEDIUMPROPERTIES_H_INCLUDED
