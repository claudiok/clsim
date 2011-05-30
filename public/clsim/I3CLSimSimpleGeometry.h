#ifndef I3CLSIMSIMPLEGEOMETRY_H_INCLUDED
#define I3CLSIMSIMPLEGEOMETRY_H_INCLUDED

#include "icetray/I3TrayHeaders.h"

#include <vector>
#include <string>

/**
 * @brief Describes a detector geometry (abstract base class)
 */

class I3CLSimSimpleGeometry 
{
    
public:
    virtual std::size_t size() const = 0;
    virtual double GetOMRadius() const = 0;

    virtual const std::vector<int32_t> &GetStringIDVector() const = 0;
    virtual const std::vector<uint32_t> &GetDomIDVector() const = 0;
    virtual const std::vector<double> &GetPosXVector() const = 0;
    virtual const std::vector<double> &GetPosYVector() const = 0;
    virtual const std::vector<double> &GetPosZVector() const = 0;
    virtual const std::vector<std::string> &GetSubdetectorVector() const = 0;
    
    virtual int32_t GetStringID(std::size_t pos) const = 0;
    virtual uint32_t GetDomID(std::size_t pos) const = 0;
    virtual double GetPosX(std::size_t pos) const = 0;
    virtual double GetPosY(std::size_t pos) const = 0;
    virtual double GetPosZ(std::size_t pos) const = 0;
    virtual std::string GetSubdetector(std::size_t pos) const = 0;
    
    
};

I3_POINTER_TYPEDEFS(I3CLSimSimpleGeometry);

#endif //I3CLSIMSIMPLEGEOMETRY_H_INCLUDED
