#ifndef I3CLSIMSIMPLEGEOMETRYFROMI3GEOMETRY_H_INCLUDED
#define I3CLSIMSIMPLEGEOMETRYFROMI3GEOMETRY_H_INCLUDED

#include "clsim/I3CLSimSimpleGeometry.h"
#include "dataclasses/geometry/I3Geometry.h"

#include <string>
#include <set>

/**
 * @brief Describes a detector geometry.
 *
 * Reads from a simple text file.
 */
class I3CLSimSimpleGeometryFromI3Geometry : public I3CLSimSimpleGeometry
{
public:
    static const std::set<int> default_ignoreStrings;
    static const std::set<unsigned int> default_ignoreDomIDs;
    static const std::set<std::string> default_ignoreSubdetectors;
    static const int32_t default_ignoreStringIDsSmallerThan;
    static const int32_t default_ignoreStringIDsLargerThan;
    static const uint32_t default_ignoreDomIDsSmallerThan;
    static const uint32_t default_ignoreDomIDsLargerThan;
    static const bool default_splitIntoPartsAccordingToPosition;
    static const bool default_useHardcodedDeepCoreSubdetector;

    I3CLSimSimpleGeometryFromI3Geometry(double OMRadius, 
                                        const I3GeometryConstPtr &geometry,
                                        const std::set<int> &ignoreStrings=default_ignoreStrings,
                                        const std::set<unsigned int> &ignoreDomIDs=default_ignoreDomIDs,
                                        const std::set<std::string> &ignoreSubdetectors=default_ignoreSubdetectors,
                                        int32_t ignoreStringIDsSmallerThan=default_ignoreStringIDsSmallerThan,
                                        int32_t ignoreStringIDsLargerThan=default_ignoreStringIDsLargerThan,
                                        uint32_t ignoreDomIDsSmallerThan=default_ignoreDomIDsSmallerThan,
                                        uint32_t ignoreDomIDsLargerThan=default_ignoreDomIDsLargerThan,
                                        bool splitIntoPartsAccordingToPosition=default_splitIntoPartsAccordingToPosition,
                                        bool useHardcodedDeepCoreSubdetector=default_useHardcodedDeepCoreSubdetector);
    virtual ~I3CLSimSimpleGeometryFromI3Geometry();

    virtual std::size_t size() const {return numOMs_;}

    virtual double GetOMRadius() const {return OMRadius_;}
    
    virtual const std::vector<int32_t> &GetStringIDVector() const {return stringIDs_;}
    virtual const std::vector<uint32_t> &GetDomIDVector() const {return domIDs_;}
    virtual const std::vector<double> &GetPosXVector() const {return posX_;}
    virtual const std::vector<double> &GetPosYVector() const {return posY_;}
    virtual const std::vector<double> &GetPosZVector() const {return posZ_;}
    virtual const std::vector<std::string> &GetSubdetectorVector() const {return subdetectors_;}

    virtual int32_t GetStringID(std::size_t pos) const {return stringIDs_.at(pos);}
    virtual uint32_t GetDomID(std::size_t pos) const {return domIDs_.at(pos);}
    virtual double GetPosX(std::size_t pos) const {return posX_.at(pos);}
    virtual double GetPosY(std::size_t pos) const {return posY_.at(pos);}
    virtual double GetPosZ(std::size_t pos) const {return posZ_.at(pos);}
    virtual std::string GetSubdetector(std::size_t pos) const {return subdetectors_.at(pos);}
    
private:
    double OMRadius_;
    bool splitIntoPartsAccordingToPosition_;
    std::size_t numOMs_;
    
    std::vector<int32_t> stringIDs_;
    std::vector<uint32_t> domIDs_;
    std::vector<double> posX_;
    std::vector<double> posY_;
    std::vector<double> posZ_;
    std::vector<std::string> subdetectors_;
    
    std::set<int> ignoreStrings_;
    std::set<unsigned int> ignoreDomIDs_;
    std::set<std::string> ignoreSubdetectors_;

    int32_t ignoreStringIDsSmallerThan_;
    int32_t ignoreStringIDsLargerThan_;
    uint32_t ignoreDomIDsSmallerThan_;
    uint32_t ignoreDomIDsLargerThan_;

    bool useHardcodedDeepCoreSubdetector_;
    
};

I3_POINTER_TYPEDEFS(I3CLSimSimpleGeometryFromI3Geometry);

#endif //I3CLSIMSIMPLEGEOMETRYFROMI3GEOMETRY_H_INCLUDED
