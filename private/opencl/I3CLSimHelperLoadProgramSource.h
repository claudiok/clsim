#ifndef I3CLSIMHELPERLOADPROGRAMSOURCE_H_INCLUDED
#define I3CLSIMHELPERLOADPROGRAMSOURCE_H_INCLUDED

#include <string>
#include <stdexcept>

namespace I3CLSimHelper
{
    std::string LoadProgramSource(const std::string &filename)
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
