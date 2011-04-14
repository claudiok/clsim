#include "clsim/I3CLSimSimpleGeometryTextFile.h"

#include <stdexcept>
#include <fstream>
#include <limits>

#include <boost/lexical_cast.hpp>

// TODO: these defaults are IceCube-specific!
const int32_t I3CLSimSimpleGeometryTextFile::default_ignoreStringIDsSmallerThan = 1;
const int32_t I3CLSimSimpleGeometryTextFile::default_ignoreStringIDsLargerThan = std::numeric_limits<int32_t>::max();
const uint32_t I3CLSimSimpleGeometryTextFile::default_ignoreDomIDsSmallerThan = 1;
const uint32_t I3CLSimSimpleGeometryTextFile::default_ignoreDomIDsLargerThan = 60;

I3CLSimSimpleGeometryTextFile::
I3CLSimSimpleGeometryTextFile(double OMRadius,
                              const std::string &filename,
                              int32_t ignoreStringIDsSmallerThan,
                              int32_t ignoreStringIDsLargerThan,
                              uint32_t ignoreDomIDsSmallerThan,
                              uint32_t ignoreDomIDsLargerThan
                              )
:
OMRadius_(OMRadius),
ignoreStringIDsSmallerThan_(ignoreStringIDsSmallerThan),
ignoreStringIDsLargerThan_(ignoreStringIDsLargerThan),
ignoreDomIDsSmallerThan_(ignoreDomIDsSmallerThan),
ignoreDomIDsLargerThan_(ignoreDomIDsLargerThan)
{
	std::ifstream inFile;
	inFile.open(filename.c_str(), std::ifstream::in);
    
    if (inFile.fail()) throw std::runtime_error("Could not open input file");
    
    int64_t readString; // 1 - 86
    int64_t readDom;    // 1 - 60
    double readx, ready, readz;  // dom x, y, z read from file
    
    numOMs_=0;
    while (inFile >> readString >> readDom >> readx >> ready >> readz)
    {
        int32_t string;
        uint32_t dom;
        double x;
        double y;
        double z;
        
        try
        {
            string = boost::lexical_cast<int32_t>(readString);
            dom = boost::lexical_cast<uint32_t>(readDom);
            x = boost::lexical_cast<double>(readx);
            y = boost::lexical_cast<double>(ready);
            z = boost::lexical_cast<double>(readz);
        }
        catch(boost::bad_lexical_cast &)
        {
            throw std::runtime_error("Read error (numeric conversion)!");
        }
        
        if ((string < ignoreStringIDsSmallerThan_) ||
            (string > ignoreStringIDsLargerThan_) ||
            (dom < ignoreDomIDsSmallerThan_) ||
            (dom > ignoreDomIDsLargerThan_))
            continue;
        
        stringIDs_.push_back(string);
        domIDs_.push_back(dom);
        posX_.push_back(x);
        posY_.push_back(y);
        posZ_.push_back(z);
        ++numOMs_;
    }

    inFile.close();
}

I3CLSimSimpleGeometryTextFile::
~I3CLSimSimpleGeometryTextFile()
{
    
}
