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
 * @file I3CLSimPhoton.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#include <icetray/serialization.h>
#include <clsim/I3CLSimPhoton.h>

#include <icetray/portable_binary_archive.hpp>

#include <boost/static_assert.hpp>
#include <boost/serialization/binary_object.hpp>

using namespace boost::archive;

namespace {
    const std::size_t blobSizeV0 = 80; // size of our structure in bytes
}

I3CLSimPhoton::~I3CLSimPhoton() { }

template <class Archive>
void I3CLSimPhoton::save(Archive &ar, unsigned version) const
{
    ar << make_nvp("x", ((const cl_float *)&posAndTime)[0]);
    ar << make_nvp("y", ((const cl_float *)&posAndTime)[1]);
    ar << make_nvp("z", ((const cl_float *)&posAndTime)[2]);
    ar << make_nvp("time", ((const cl_float *)&posAndTime)[3]);
    
    ar << make_nvp("theta", ((const cl_float *)&dir)[0]);
    ar << make_nvp("phi", ((const cl_float *)&dir)[1]);
    ar << make_nvp("wavelength", wavelength);
    ar << make_nvp("cherenkovDist", cherenkovDist);
    
    ar << make_nvp("numScatters", numScatters);
    ar << make_nvp("weight", weight);
    ar << make_nvp("id", identifier);
    ar << make_nvp("stringID", stringID);
    ar << make_nvp("omID", omID);

    ar << make_nvp("startX", ((const cl_float *)&startPosAndTime)[0]);
    ar << make_nvp("startY", ((const cl_float *)&startPosAndTime)[1]);
    ar << make_nvp("startZ", ((const cl_float *)&startPosAndTime)[2]);
    ar << make_nvp("startTime", ((const cl_float *)&startPosAndTime)[3]);
    
    ar << make_nvp("startTheta", ((const cl_float *)&startDir)[0]);
    ar << make_nvp("startPhi", ((const cl_float *)&startDir)[1]);

    ar << make_nvp("groupVelocity", groupVelocity);
    ar << make_nvp("distInAbsLens", distInAbsLens);
}     


template <class Archive>
void I3CLSimPhoton::load(Archive &ar, unsigned version)
{
    if (version > i3clsimphoton_version_)
        log_fatal("Attempting to read version %u from file but running version %u of I3CLSimPhoton class.",version,i3clsimphoton_version_);
    
    float temp; uint32_t temp_uint; int16_t temp_short; uint16_t temp_ushort;
    ar >> make_nvp("x", temp); ((cl_float *)&posAndTime)[0]=temp;
    ar >> make_nvp("y", temp); ((cl_float *)&posAndTime)[1]=temp;
    ar >> make_nvp("z", temp); ((cl_float *)&posAndTime)[2]=temp;
    ar >> make_nvp("time", temp); ((cl_float *)&posAndTime)[3]=temp;

    ar >> make_nvp("theta", temp); ((cl_float *)&dir)[0]=temp;
    ar >> make_nvp("phi", temp); ((cl_float *)&dir)[1]=temp;
    ar >> make_nvp("wavelength", temp); wavelength=temp;
    ar >> make_nvp("cherenkovDist", temp); cherenkovDist=temp;
    
    ar >> make_nvp("numScatters", temp_uint); numScatters=temp_uint;
    ar >> make_nvp("weight", temp); weight=temp;
    ar >> make_nvp("id", temp_uint); identifier=temp_uint;
    ar >> make_nvp("stringID", temp_short); stringID=temp_short;
    ar >> make_nvp("omID", temp_ushort); omID=temp_ushort;

    ar >> make_nvp("startX", temp); ((cl_float *)&startPosAndTime)[0]=temp;
    ar >> make_nvp("startY", temp); ((cl_float *)&startPosAndTime)[1]=temp;
    ar >> make_nvp("startZ", temp); ((cl_float *)&startPosAndTime)[2]=temp;
    ar >> make_nvp("startTime", temp); ((cl_float *)&startPosAndTime)[3]=temp;
    
    ar >> make_nvp("startTheta", temp); ((cl_float *)&startDir)[0]=temp;
    ar >> make_nvp("startPhi", temp); ((cl_float *)&startDir)[1]=temp;

    ar >> make_nvp("groupVelocity", temp); groupVelocity=temp;
    ar >> make_nvp("distInAbsLens", temp); distInAbsLens=temp;

}     



// just save the binary blob for binary archives (internal storage is little-endian)

template <>
void I3CLSimPhoton::save(portable_binary_oarchive &ar, unsigned version) const
{
    // check an assumption we will make throughout the code
    BOOST_STATIC_ASSERT((sizeof(I3CLSimPhoton) == blobSizeV0));
    
    ar << make_nvp("blob", boost::serialization::make_binary_object((void *)this, blobSizeV0));
}     

template <>
void I3CLSimPhoton::load(portable_binary_iarchive &ar, unsigned version)
{
    if (version > i3clsimphoton_version_)
        log_fatal("Attempting to read version %u from file but running version %u of I3CLSimPhoton class.",version,i3clsimphoton_version_);
    
    
    // check an assumption we will make throughout the code
    BOOST_STATIC_ASSERT((sizeof(I3CLSimPhoton) == blobSizeV0));
    
    ar >> make_nvp("blob", boost::serialization::make_binary_object(this, blobSizeV0));
}     

// this serialization is endian-dependent (i.e. if you serializae on a big-endian system, you will
// not be able to read this on a little-endian system). This should not matter since these objects
// are not meant to be trnamsitted over a network. TODO: To be extremely safe, an endianess flag
// should be included in the serialized output.

template<>
template<>
void I3Vector<I3CLSimPhoton>::serialize(portable_binary_iarchive &ar, unsigned version)
{
    ar >> make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));
    unsigned I3CLSimPhoton_version;
    ar >> make_nvp("i3clsimphoton_version", I3CLSimPhoton_version);
    if (I3CLSimPhoton_version != i3clsimphoton_version_)
        log_fatal("This reader can only read I3Vector<I3CLSimPhoton> version %u, but %u was provided.",i3clsimphoton_version_,I3CLSimPhoton_version);
    uint64_t size;
    ar >> make_nvp("num", size);
    
    this->resize(size);

    // read the binary blob in one go..
    ar >> make_nvp("blob", boost::serialization::make_binary_object( &((*this)[0]), blobSizeV0*size));
}

template<>
template<>
void I3Vector<I3CLSimPhoton>::serialize(portable_binary_oarchive &ar, unsigned version)
{
    ar << make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));
    ar << make_nvp("i3clsimphoton_version", i3clsimphoton_version_);
    uint64_t size = this->size();
    ar << make_nvp("num", size);
    ar << make_nvp("blob", boost::serialization::make_binary_object( &((*this)[0]), blobSizeV0*size ));
}



I3_SERIALIZABLE(I3CLSimPhoton);
I3_SERIALIZABLE(I3CLSimPhotonSeries);
I3_SERIALIZABLE(I3CLSimPhotonSeriesMap);
