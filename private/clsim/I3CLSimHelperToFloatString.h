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
 * @file I3CLSimHelperToFloatString.h
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#ifndef TO_FLOAT_STRING_H_INCLUDED
#define TO_FLOAT_STRING_H_INCLUDED

#include <string>
#include <boost/lexical_cast.hpp>
#include <assert.h>
#include <sstream>

namespace I3CLSimHelper
{
    inline std::string ToFloatString(double val)
    {
        if (val==1e-9) return "1e-9f";
        else if (val==1e-6) return "1e-6f";
        else if (val==1e-3) return "1e-3f";
        else if (val==1.) return "1.f";
        else if (val==0.) return "0.f";

        std::ostringstream output(std::ostringstream::out);
        output.setf(std::ios::scientific,std::ios::floatfield);
        output.precision(std::numeric_limits<float>::digits10+4); // maximum precision for a float

        output << val;
        std::string ret = output.str();

        if ((ret.find('.')==std::string::npos) && (ret.find('e')==std::string::npos) && (ret.find('E')==std::string::npos))
        {
            ret = ret+".";
        }

        ret = ret+"f";

        return ret;
    }

    inline std::string ToDoubleString(double val)
    {
        if (val==1e-9) return "1e-9";
        else if (val==1e-6) return "1e-6";
        else if (val==1e-3) return "1e-3";
        else if (val==1.) return "1.";
        else if (val==0.) return "0.";

        std::ostringstream output(std::ostringstream::out);
        output.setf(std::ios::scientific,std::ios::floatfield);
        output.precision(std::numeric_limits<double>::digits10+4); // maximum precision for a float

        output << val;
        std::string ret = output.str();

        if ((ret.find('.')==std::string::npos) && (ret.find('e')==std::string::npos) && (ret.find('E')==std::string::npos))
        {
            ret = ret+".";
        }

        return ret;
    }

};

#endif //TO_FLOAT_STRING_H_INCLUDED
