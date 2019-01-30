
#ifndef CLSIM_I3CLSimPhotonToMCPEConverter_h_INCLUDED
#define CLSIM_I3CLSimPhotonToMCPEConverter_h_INCLUDED

#include <icetray/OMKey.h>
#include <simclasses/I3MCPE.h>

class ModuleKey;
class I3CompressedPhoton;

class I3CLSimPhotonToMCPEConverter {
public:
    virtual ~I3CLSimPhotonToMCPEConverter();
    virtual std::tuple<OMKey,I3MCPE,bool> Convert(const ModuleKey&, const I3CompressedPhoton &) const = 0;
};

#endif
