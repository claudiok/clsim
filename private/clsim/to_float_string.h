#ifndef TO_FLOAT_STRING_H_INCLUDED
#define TO_FLOAT_STRING_H_INCLUDED

#include <string>
#include <boost/lexical_cast.hpp>
#include <assert.h>

namespace I3CLSimHelper
{
    inline std::string to_float_string(double val)
    {
        if (val==1e-9) return "1e-9f";
        else if (val==1e-6) return "1e-6f";
        else if (val==1e-3) return "1e-6f";
        else if (val==1.) return "1.f";

        std::string ret = boost::lexical_cast<std::string>(val);
        
        if ((ret.find('.')==std::string::npos) && (ret.find('e')==std::string::npos) && (ret.find('E')==std::string::npos))
        {
            ret = ret+".";
        }
        
        ret = ret+"f";
        
        return ret;
    }

};

#endif //TO_FLOAT_STRING_H_INCLUDED
