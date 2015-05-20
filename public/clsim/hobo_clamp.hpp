#if BOOST_VERSION >= 105100
#include <boost/algorithm/clamp.hpp>
#else
namespace boost { namespace algorithm {

  template <typename T>
    T
    clamp(T value, T lo, T hi) {
        return (value < lo) ? lo : ((value > hi) ? hi : value);
    }

}}
#endif
