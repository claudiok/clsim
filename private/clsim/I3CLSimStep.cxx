/**
 * copyright (C) 2011
 * Claudio Kopper <claudio.kopper@nikhef.nl>
 * $Id$
 *
 * @file I3CLSimStep.cxx
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

#include <icetray/serialization.h>
#include <clsim/I3CLSimStep.h>

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

    ar << make_nvp("dummy", dummy);
}     


template <class Archive>
void I3CLSimStep::load(Archive &ar, unsigned version)
{
	if (version > i3clsimstep_version_)
		log_fatal("Attempting to read version %u from file but running version %u of I3CLSimStep class.",version,i3clsimstep_version_);
    
    float temp; uint32_t temp_uint;
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
    ar >> make_nvp("dummy", temp_uint); dummy=temp_uint;

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
