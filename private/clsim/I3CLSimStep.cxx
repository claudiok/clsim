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
 * @file I3CLSimStep.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#include <icetray/serialization.h>
#include <clsim/I3CLSimStep.h>

#include <icetray/portable_binary_archive.hpp>

#include <boost/static_assert.hpp>
#include <boost/serialization/binary_object.hpp>

using namespace boost::archive;

namespace {
    const std::size_t blobSizeV0 = 48; // size of our structure in bytes
}

I3CLSimStep::~I3CLSimStep() { }

template <class Archive>
void I3CLSimStep::save(Archive &ar, unsigned version) const
{
    ar << make_nvp("x", ((const cl_float *)&posAndTime)[0]);
    ar << make_nvp("y", ((const cl_float *)&posAndTime)[1]);
    ar << make_nvp("z", ((const cl_float *)&posAndTime)[2]);
    ar << make_nvp("time", ((const cl_float *)&posAndTime)[3]);
    
    ar << make_nvp("theta", ((const cl_float *)&dirAndLengthAndBeta)[0]);
    ar << make_nvp("phi", ((const cl_float *)&dirAndLengthAndBeta)[1]);
    ar << make_nvp("length", ((const cl_float *)&dirAndLengthAndBeta)[2]);
    ar << make_nvp("beta", ((const cl_float *)&dirAndLengthAndBeta)[3]);
    
    ar << make_nvp("num", numPhotons);
    ar << make_nvp("weight", weight);
    ar << make_nvp("id", identifier);

    ar << make_nvp("sourceType", sourceType);
    ar << make_nvp("dummy1", dummy1);
    ar << make_nvp("dummy2", dummy2);
}     


template <class Archive>
void I3CLSimStep::load(Archive &ar, unsigned version)
{
    if (version > i3clsimstep_version_)
        log_fatal("Attempting to read version %u from file but running version %u of I3CLSimStep class.",version,i3clsimstep_version_);

    float temp; uint32_t temp_uint;
    uint16_t temp_uint16; uint8_t temp_uint8;
    ar >> make_nvp("x", temp); ((cl_float *)&posAndTime)[0]=temp;
    ar >> make_nvp("y", temp); ((cl_float *)&posAndTime)[1]=temp;
    ar >> make_nvp("z", temp); ((cl_float *)&posAndTime)[2]=temp;
    ar >> make_nvp("time", temp); ((cl_float *)&posAndTime)[3]=temp;

    ar >> make_nvp("theta", temp); ((cl_float *)&dirAndLengthAndBeta)[0]=temp;
    ar >> make_nvp("phi", temp); ((cl_float *)&dirAndLengthAndBeta)[1]=temp;
    ar >> make_nvp("length", temp); ((cl_float *)&dirAndLengthAndBeta)[2]=temp;
    ar >> make_nvp("beta", temp); ((cl_float *)&dirAndLengthAndBeta)[3]=temp;

    ar >> make_nvp("num", temp_uint); numPhotons=temp_uint;
    ar >> make_nvp("weight", temp); weight=temp;
    ar >> make_nvp("id", temp_uint); identifier=temp_uint;
    ar >> make_nvp("sourceType", temp_uint8); sourceType=temp_uint8;
    ar >> make_nvp("dummy1", temp_uint8); dummy1=temp_uint8;
    ar >> make_nvp("dummy2", temp_uint16); dummy2=temp_uint16;

}



// just save the binary blob for binary archives (internal storage is little-endian)

template <>
void I3CLSimStep::save(portable_binary_oarchive &ar, unsigned version) const
{
    // check an assumption we will make throughout the code
    BOOST_STATIC_ASSERT((sizeof(I3CLSimStep) == blobSizeV0));

    ar << make_nvp("blob", boost::serialization::make_binary_object((void *)this, blobSizeV0));
}

template <>
void I3CLSimStep::load(portable_binary_iarchive &ar, unsigned version)
{
    if (version > i3clsimstep_version_)
        log_fatal("Attempting to read version %u from file but running version %u of I3CLSimStep class.",version,i3clsimstep_version_);

    // check an assumption we will make throughout the code
    BOOST_STATIC_ASSERT((sizeof(I3CLSimStep) == blobSizeV0));

    ar >> make_nvp("blob", boost::serialization::make_binary_object(this, blobSizeV0));
}

template<>
template<>
void I3Vector<I3CLSimStep>::serialize(portable_binary_iarchive &ar, unsigned version)
{
    ar >> make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));
    unsigned I3CLSimStep_version;
    ar >> make_nvp("I3CLSimStep_version", I3CLSimStep_version);
    if (I3CLSimStep_version != i3clsimstep_version_)
        log_fatal("This reader can only read I3Vector<I3CLSimStep> version %u, but %u was provided.",i3clsimstep_version_,I3CLSimStep_version);
    uint64_t size;
    ar >> make_nvp("num", size);

    this->resize(size);

    // read the binary blob in one go..
    ar >> make_nvp("blob", boost::serialization::make_binary_object( &((*this)[0]), blobSizeV0*size));
}

template<>
template<>
void I3Vector<I3CLSimStep>::serialize(portable_binary_oarchive &ar, unsigned version)
{
    ar << make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));
    ar << make_nvp("I3CLSimStep_version", i3clsimstep_version_);
    uint64_t size = this->size();
    ar << make_nvp("num", size);
    ar << make_nvp("blob", boost::serialization::make_binary_object( &((*this)[0]), blobSizeV0*size ));
}



I3_SERIALIZABLE(I3CLSimStep);
I3_SERIALIZABLE(I3CLSimStepSeries);
