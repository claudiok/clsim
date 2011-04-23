/**
 * copyright (C) 2011
 * Claudio Kopper <claudio.kopper@nikhef.nl>
 * $Id$
 *
 * @file I3CLSimPhoton.cxx
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
#include <clsim/I3CLSimPhoton.h>

#include <boost/static_assert.hpp>
#include <boost/serialization/binary_object.hpp>

#include <boost/archive/portable_binary_iarchive.hpp>
#include <boost/archive/portable_binary_oarchive.hpp>

using namespace boost::archive;

namespace {
    const std::size_t blobSizeV0 = 48; // size of our structure in bytes
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

// TODO: make it work on big-endian systems!

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
