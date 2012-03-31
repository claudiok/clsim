#ifndef ICETRAY_PYTHON_STREAM_TO_STRING_HPP_INCLUDED
#define ICETRAY_PYTHON_STREAM_TO_STRING_HPP_INCLUDED

#include <sstream>

template <class T>
std::string stream_to_string(T t){
  std::stringstream ss;
  ss << t;
  return ss.str();
}

#endif // ICETRAY_PYTHON_STREAM_TO_STRING_HPP_INCLUDED
