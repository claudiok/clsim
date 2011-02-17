//
//   Copyright (c) 2011  Claudio Kopper <claudio.kopper@nikhef.nl>
//   
#ifndef BOOST_SERIALIZABLE_PICKLE_SUITE_HPP_INCLUDED
#define BOOST_SERIALIZABLE_PICKLE_SUITE_HPP_INCLUDED

#include <vector>

#include "icetray/serialization.h"
#include "icetray/I3Frame.h"

#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/device/array.hpp>
#include <boost/iostreams/device/back_inserter.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filtering_stream.hpp>

#include <boost/python/tuple.hpp>
#include <boost/python/str.hpp>
#include <boost/python/dict.hpp>
#include <boost/python/extract.hpp>


namespace I3 {namespace detail {
    
template<class T>
inline void to_vector(const T &element, std::vector<char> &blobBuffer)
{
    typedef boost::iostreams::stream<boost::iostreams::back_insert_device<std::vector<char> > > vecstream_t;
    vecstream_t blobBufStream(blobBuffer);
    {
        boost::archive::portable_binary_oarchive blobBufArchive(blobBufStream);
        blobBufArchive << make_nvp("T", element);
    }
    blobBufStream.flush();
}            

template<class T>
inline void from_buffer(T &element, const char *buffer, std::size_t size)
{
    boost::iostreams::array_source src(buffer, size);
    boost::iostreams::filtering_istream fis(src);
    boost::archive::portable_binary_iarchive pia(fis);
    
    pia >> element;
}            

// specialization for I3Frame
template<>
inline void to_vector(const I3Frame &element, std::vector<char> &blobBuffer)
{
    boost::iostreams::filtering_ostream blobBufStream(boost::iostreams::back_inserter(blobBuffer));
    {
        element.save(blobBufStream);
    }
    blobBufStream.flush();
}            

// specialization for I3Frame
template<>
inline void from_buffer(I3Frame &element, const char *buffer, std::size_t size)
{
    boost::iostreams::array_source src(buffer, size);
    boost::iostreams::filtering_istream fis(src);
    
    element.load(fis);
}            

}};

template<class T>
struct boost_serializable_pickle_suite : boost::python::pickle_suite
{
    static boost::python::tuple getstate(boost::python::object element_obj)
    {
        const T &element = boost::python::extract<const T &>(element_obj)();
        
        // serialize the object into a buffer
        std::vector<char> blobBuffer;
        I3::detail::to_vector(element, blobBuffer);
        
        return boost::python::make_tuple
        (
         element_obj.attr("__dict__"),
         boost::python::str( &(blobBuffer[0]), blobBuffer.size() )
        );
    }
    
    static void setstate(boost::python::object element_obj, boost::python::tuple state)
    {
        T &element = boost::python::extract<T &>(element_obj)();
        
        if (boost::python::len(state) != 2)
        {
            PyErr_SetObject(PyExc_ValueError,
                            ("expected 2-item tuple in call to __setstate__; got %s"
                             % state).ptr()
                            );
            boost::python::throw_error_already_set();
        }
        
        // restore the object's __dict__
        boost::python::dict d = boost::python::extract<boost::python::dict>(element_obj.attr("__dict__"))();
        d.update(state[0]);
        
        // restore the C++ object
        boost::python::str buffer = boost::python::extract<boost::python::str>(state[1]);

        I3::detail::from_buffer(element, boost::python::extract<char const*>(buffer), boost::python::len(buffer));
    }
    
    static bool getstate_manages_dict() { return true; }
};

#endif //BOOST_SERIALIZABLE_PICKLE_SUITE_HPP_INCLUDED
