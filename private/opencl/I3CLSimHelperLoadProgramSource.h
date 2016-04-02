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
 * @file I3CLSimHelperLoadProgramSource.h
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#ifndef I3CLSIMHELPERLOADPROGRAMSOURCE_H_INCLUDED
#define I3CLSIMHELPERLOADPROGRAMSOURCE_H_INCLUDED

#include <string>
#include <stdexcept>
#include <cstdio>

namespace I3CLSimHelper
{
    inline std::string LoadProgramSource(const std::string &filename)
    {
        // open the OpenCL source code file
        FILE* pFileStream = fopen(filename.c_str(), "rb");
        if (pFileStream == 0) 
        {       
            throw std::runtime_error("Could not open source file.");
        }
        
        // get the length of the source code
        fseek(pFileStream, 0, SEEK_END); 
        std::size_t szSourceLength = ftell(pFileStream);
        fseek(pFileStream, 0, SEEK_SET); 
        
        // allocate a buffer for the source code string and read it in
        char* cSourceString = (char *)malloc(szSourceLength + 1); 
        if (fread((cSourceString), szSourceLength, 1, pFileStream) != 1)
        {
            fclose(pFileStream);
            free(cSourceString);
            throw std::runtime_error("Could not load source file.");
        }
        
        // close the file and return the total length of the combined (preamble + source) string
        fclose(pFileStream);

        cSourceString[szSourceLength] = '\0';

        std::string retString(cSourceString, szSourceLength);
        free(cSourceString);

        return retString;
    }
    


};

#endif //I3CLSIMHELPERLOADPROGRAMSOURCE_H_INCLUDED
