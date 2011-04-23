/**
 * copyright (C) 2011
 * Claudio Kopper <claudio.kopper@nikhef.nl>
 * $Id$
 *
 * @file I3CLSimPhoton.h
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 *
 *
 *  This file is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>
 *  
 */

#ifndef I3CLSIMPHOTON_H_INCLUDED
#define I3CLSIMPHOTON_H_INCLUDED

#ifndef I3CLSIM_WITHOUT_OPENCL

#ifdef __APPLE__
#include "OpenCL/opencl.h"                                                                                                                                                                       
#else
#include "CL/opencl.h"                                                                                                                                                                       
#endif

#else //I3CLSIM_WITHOUT_OPENCL

// I3CLSimPhoton is being built without OpenCL. Using our own copy
// of the cl_platform header
#include "clsim/fake_cl_platform.h"

#endif //I3CLSIM_WITHOUT_OPENCL


#include <cstring>

#include "dataclasses/I3Vector.h"
#include "dataclasses/I3Map.h"
#include "dataclasses/I3Direction.h"
#include "dataclasses/I3Position.h"



/**
 * @brief A single Cherenkov photon, either before or
 * after propagation to a target (DOM)
 */
static const unsigned i3clsimphoton_version_ = 0;

struct I3CLSimPhoton 
{
    
public:
    
    I3CLSimPhoton() {;}
    
    ~I3CLSimPhoton();

    
    inline float GetPosX() const {return ((const cl_float *)&posAndTime)[0];}
    inline float GetPosY() const {return ((const cl_float *)&posAndTime)[1];}
    inline float GetPosZ() const {return ((const cl_float *)&posAndTime)[2];}
    inline float GetTime() const {return ((const cl_float *)&posAndTime)[3];}
    inline float GetDirTheta() const {return ((const cl_float *)&dir)[0];}
    inline float GetDirPhi() const {return ((const cl_float *)&dir)[1];}

    inline float GetStartPosX() const {return ((const cl_float *)&startPosAndTime)[0];}
    inline float GetStartPosY() const {return ((const cl_float *)&startPosAndTime)[1];}
    inline float GetStartPosZ() const {return ((const cl_float *)&startPosAndTime)[2];}
    inline float GetStartTime() const {return ((const cl_float *)&startPosAndTime)[3];}
    inline float GetStartDirTheta() const {return ((const cl_float *)&startDir)[0];}
    inline float GetStartDirPhi() const {return ((const cl_float *)&startDir)[1];}

    inline float GetWavelength() const {return wavelength;}
    inline float GetCherenkovDist() const {return cherenkovDist;}
    inline uint32_t GetNumScatters() const {return numScatters;}
    inline float GetWeight() const {return weight;}
    inline uint32_t GetID() const {return identifier;}
    inline int16_t GetStringID() const {return stringID;}
    inline uint16_t GetOMID() const {return omID;}
    inline float GetGroupVelocity() const {return groupVelocity;}

    inline I3PositionPtr GetPos() const {return I3PositionPtr(new I3Position(((const cl_float *)&posAndTime)[0], ((const cl_float *)&posAndTime)[1], ((const cl_float *)&posAndTime)[2]));}
    inline I3PositionPtr GetStartPos() const
    {
        return I3PositionPtr(new I3Position(((const cl_float *)&startPosAndTime)[0],
                                            ((const cl_float *)&startPosAndTime)[1],
                                            ((const cl_float *)&startPosAndTime)[2]));
    }

    inline I3DirectionPtr GetDir() const 
    {
        I3DirectionPtr retval(new I3Direction());
        retval->SetThetaPhi(((const cl_float *)&dir)[0], ((const cl_float *)&dir)[1]);
        return retval;
    }

    inline I3DirectionPtr GetStartDir() const 
    {
        I3DirectionPtr retval(new I3Direction());
        retval->SetThetaPhi(((const cl_float *)&startDir)[0], ((const cl_float *)&startDir)[1]);
        return retval;
    }

    
    
    inline void SetPosX(const float &val) {((cl_float *)&posAndTime)[0]=val;}
    inline void SetPosY(const float &val) {((cl_float *)&posAndTime)[1]=val;}
    inline void SetPosZ(const float &val) {((cl_float *)&posAndTime)[2]=val;}
    inline void SetTime(const float &val) {((cl_float *)&posAndTime)[3]=val;}
    inline void SetDirTheta(const float &val) {((cl_float *)&dir)[0]=val;}
    inline void SetDirPhi(const float &val) {((cl_float *)&dir)[1]=val;}

    inline void SetStartPosX(const float &val) {((cl_float *)&startPosAndTime)[0]=val;}
    inline void SetStartPosY(const float &val) {((cl_float *)&startPosAndTime)[1]=val;}
    inline void SetStartPosZ(const float &val) {((cl_float *)&startPosAndTime)[2]=val;}
    inline void SetStartTime(const float &val) {((cl_float *)&startPosAndTime)[3]=val;}
    inline void SetStartDirTheta(const float &val) {((cl_float *)&startDir)[0]=val;}
    inline void SetStartDirPhi(const float &val) {((cl_float *)&startDir)[1]=val;}

    inline void SetWavelength(const float &val) {wavelength=val;}
    inline void SetCherenkovDist(const float &val) {cherenkovDist=val;}
    inline void SetNumScatters(const uint32_t &val) {numScatters=val;}
    inline void SetWeight(const float &val) {weight=val;}
    inline void SetID(const uint32_t &val) {identifier=val;}
    inline void SetStringID(const int16_t &val) {stringID=val;}
    inline void SetOMID(const uint16_t &val) {omID=val;}
    inline void SetGroupVelocity(const float &val) {groupVelocity=val;}

    inline void SetPos(const I3Position &pos)
    {
        ((cl_float *)&posAndTime)[0]=pos.GetX();
        ((cl_float *)&posAndTime)[1]=pos.GetY();
        ((cl_float *)&posAndTime)[2]=pos.GetZ();
    }
    inline void SetStartPos(const I3Position &pos)
    {
        ((cl_float *)&startPosAndTime)[0]=pos.GetX();
        ((cl_float *)&startPosAndTime)[1]=pos.GetY();
        ((cl_float *)&startPosAndTime)[2]=pos.GetZ();
    }

    inline void SetDir(const I3Direction &dir) 
    {
        ((cl_float *)&dir)[0]=dir.CalcTheta();
        ((cl_float *)&dir)[1]=dir.CalcPhi();
    }
    inline void SetDir(const double &x, const double &y, const double &z) 
    {
        const I3Direction dir(x,y,z);
        ((cl_float *)&dir)[0]=dir.CalcTheta();
        ((cl_float *)&dir)[1]=dir.CalcPhi();
    }

    inline void SetStartDir(const I3Direction &dir) 
    {
        ((cl_float *)&startDir)[0]=dir.CalcTheta();
        ((cl_float *)&startDir)[1]=dir.CalcPhi();
    }
    inline void SetStartDir(const double &x, const double &y, const double &z) 
    {
        const I3Direction dir(x,y,z);
        ((cl_float *)&startDir)[0]=dir.CalcTheta();
        ((cl_float *)&startDir)[1]=dir.CalcPhi();
    }

    
    
    // cl_float4 is a struct consisting of 4 floats named .x, .y, .z, .w
    // .. well.. On MacOS OpenCL 1.0, it's just an array[4].. so access it like that
    
    cl_float4 posAndTime;   // x,y,z,time
    cl_float2 dir; // theta,phi
    cl_float wavelength; // photon wavelength
    cl_float cherenkovDist; // Cherenkov distance travelled 
    cl_uint numScatters; // number of scatters
    cl_float weight;
    cl_uint identifier;
    cl_short stringID;
    cl_ushort omID;
    cl_float4 startPosAndTime;
    cl_float2 startDir;
    cl_float groupVelocity;
    cl_uint dummy;
    
private:
    friend class boost::serialization::access;
    template <class Archive> void load(Archive & ar, unsigned version);
    template <class Archive> void save(Archive & ar, unsigned version) const;
    BOOST_SERIALIZATION_SPLIT_MEMBER();
} __attribute__ ((packed)) ;


inline bool operator==(const I3CLSimPhoton &a, const I3CLSimPhoton &b)
{
    // compare all fields (binary) [except the last field, which is a dummy)
    return (std::memcmp(&a, &b, sizeof(I3CLSimPhoton)-sizeof(cl_uint))==0);
}

BOOST_CLASS_VERSION(I3CLSimPhoton, i3clsimphoton_version_);

typedef I3Vector<I3CLSimPhoton> I3CLSimPhotonSeries;
typedef I3Map<OMKey, I3CLSimPhotonSeries> I3CLSimPhotonSeriesMap;

I3_POINTER_TYPEDEFS(I3CLSimPhoton);
I3_POINTER_TYPEDEFS(I3CLSimPhotonSeries);
I3_POINTER_TYPEDEFS(I3CLSimPhotonSeriesMap);

#endif //I3CLSIMPHOTON_H_INCLUDED
