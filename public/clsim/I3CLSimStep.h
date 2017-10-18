/**
 * Copyright (c) 2011, 2012
 * Claudio Kopper <claudio.kopper@icecube.wisc.edu>
 * and the IceCube Collaboration <http://www.icecube.wisc.edu>
 *
 * Permission to use, copy, modify, and/or distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
 * WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
 * SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION
 * OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN
 * CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 *
 *
 * $Id$
 *
 * @file I3CLSimStep.h
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#ifndef I3CLSIMSTEP_H_INCLUDED
#define I3CLSIMSTEP_H_INCLUDED

#ifndef I3CLSIM_WITHOUT_OPENCL

#ifdef __APPLE__
#include "OpenCL/opencl.h"                                                                                                                                                                       
#else
#include "CL/opencl.h"                                                                                                                                                                       
#endif

#else //I3CLSIM_WITHOUT_OPENCL

// I3CLSimStep is being built without OpenCL. Using our own copy
// of the cl_platform header
#include "clsim/fake_cl_platform.h"

#endif //I3CLSIM_WITHOUT_OPENCL


#include <cstring>

#include "dataclasses/I3Vector.h"
#include "dataclasses/I3Map.h"
#include "dataclasses/I3Direction.h"
#include "dataclasses/I3Position.h"

 namespace boost { namespace archive {
     class portable_binary_iarchive;
     class portable_binary_oarchive;
 } }


/**
 * @brief A single step of a particle producing
 * a certain amount of Cherenkov photons. This is
 * typically produced by a full particle tracker
 * like Geant4.
 */
static const unsigned i3clsimstep_version_ = 0;

struct I3CLSimStep 
{
    
public:
    
    I3CLSimStep() {;}
    
    ~I3CLSimStep();

    
    inline float GetPosX() const {return ((const cl_float *)&posAndTime)[0];}
    inline float GetPosY() const {return ((const cl_float *)&posAndTime)[1];}
    inline float GetPosZ() const {return ((const cl_float *)&posAndTime)[2];}
    inline float GetTime() const {return ((const cl_float *)&posAndTime)[3];}
    inline float GetDirTheta() const {return ((const cl_float *)&dirAndLengthAndBeta)[0];}
    inline float GetDirPhi() const {return ((const cl_float *)&dirAndLengthAndBeta)[1];}
    inline float GetLength() const {return ((const cl_float *)&dirAndLengthAndBeta)[2];}
    inline float GetBeta() const {return ((const cl_float *)&dirAndLengthAndBeta)[3];}
    inline uint32_t GetNumPhotons() const {return numPhotons;}
    inline float GetWeight() const {return weight;}
    inline uint32_t GetID() const {return identifier;}
    inline uint8_t GetSourceType() const {return sourceType;}
    inline uint8_t GetDummy1() const {return dummy1;}
    inline uint16_t GetDummy2() const {return dummy2;}

    inline I3PositionPtr GetPos() const {return I3PositionPtr(new I3Position(((const cl_float *)&posAndTime)[0], ((const cl_float *)&posAndTime)[1], ((const cl_float *)&posAndTime)[2]));}

    inline I3DirectionPtr GetDir() const 
    {
        I3DirectionPtr retval(new I3Direction());
        retval->SetThetaPhi(((const cl_float *)&dirAndLengthAndBeta)[0], ((const cl_float *)&dirAndLengthAndBeta)[1]);
        return retval;
    }
    
    
    
    inline void SetPosX(const float &val) {((cl_float *)&posAndTime)[0]=val;}
    inline void SetPosY(const float &val) {((cl_float *)&posAndTime)[1]=val;}
    inline void SetPosZ(const float &val) {((cl_float *)&posAndTime)[2]=val;}
    inline void SetTime(const float &val) {((cl_float *)&posAndTime)[3]=val;}
    inline void SetDirTheta(const float &val) {((cl_float *)&dirAndLengthAndBeta)[0]=val;}
    inline void SetDirPhi(const float &val) {((cl_float *)&dirAndLengthAndBeta)[1]=val;}
    inline void SetLength(const float &val) {((cl_float *)&dirAndLengthAndBeta)[2]=val;}
    inline void SetBeta(const float &val) {((cl_float *)&dirAndLengthAndBeta)[3]=val;}
    inline void SetNumPhotons(const uint32_t &val) {numPhotons=val;}
    inline void SetWeight(const float &val) {weight=val;}
    inline void SetID(const uint32_t &val) {identifier=val;}
    inline void SetSourceType(const uint8_t &val) {sourceType=val;}
    inline void SetDummy1(const uint8_t &val) {dummy1=val;}
    inline void SetDummy2(const uint16_t &val) {dummy2=val;}
    
    inline void SetPos(const I3Position &pos) {((cl_float *)&posAndTime)[0]=pos.GetX(); ((cl_float *)&posAndTime)[1]=pos.GetY(); ((cl_float *)&posAndTime)[2]=pos.GetZ();}

    inline void SetDir(const I3Direction &dir) 
    {
        ((cl_float *)&dirAndLengthAndBeta)[0]=dir.CalcTheta();
        ((cl_float *)&dirAndLengthAndBeta)[1]=dir.CalcPhi();
    }
    inline void SetDir(const double &x, const double &y, const double &z) 
    {
        const I3Direction dir(x,y,z);
        ((cl_float *)&dirAndLengthAndBeta)[0]=dir.CalcTheta();
        ((cl_float *)&dirAndLengthAndBeta)[1]=dir.CalcPhi();
    }
    
    
    
    // cl_float4 is a struct consisting of 4 floats named .x, .y, .z, .w
    // .. well.. On MacOS OpenCL 1.0, it's just an array[4].. so access it like that
    
    cl_float4 posAndTime;   // x,y,z,time
    cl_float4 dirAndLengthAndBeta; // theta,phi,length,beta
    cl_uint numPhotons;
    cl_float weight;
    cl_uint identifier;
    cl_uchar sourceType;
    cl_uchar dummy1;
    cl_ushort dummy2;

private:
    friend class icecube::serialization::access;
    template <class Archive> void load(Archive & ar, unsigned version);
    template <class Archive> void save(Archive & ar, unsigned version) const;
    I3_SERIALIZATION_SPLIT_MEMBER();
} __attribute__ ((packed)) ;

template<> void I3CLSimStep::save(icecube::archive::portable_binary_oarchive &ar, unsigned version) const;
template<> void I3CLSimStep::load(icecube::archive::portable_binary_iarchive &ar, unsigned version);

inline bool operator==(const I3CLSimStep &a, const I3CLSimStep &b)
{
    // compare all fields except for the last one (which is a dummy)
    return (std::memcmp(&a, &b, sizeof(I3CLSimStep)-sizeof(cl_int))==0);
}
std::ostream& operator<<(std::ostream&, const I3CLSimStep&);

I3_CLASS_VERSION(I3CLSimStep, i3clsimstep_version_);

typedef I3Vector<I3CLSimStep> I3CLSimStepSeries;

I3_POINTER_TYPEDEFS(I3CLSimStep);
I3_POINTER_TYPEDEFS(I3CLSimStepSeries);

template<> template<> void I3Vector<I3CLSimStep>::serialize(icecube::archive::portable_binary_iarchive &ar, unsigned version);
template<> template<> void I3Vector<I3CLSimStep>::serialize(icecube::archive::portable_binary_oarchive &ar, unsigned version);

#endif //I3CLSIMSTEP_H_INCLUDED
