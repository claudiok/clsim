#include <clsim/I3CLSimSpectrumTable.h>

I3CLSimSpectrumTable::I3CLSimSpectrumTable()
:
spectra_(1, I3CLSimFunctionConstPtr()) // add a NULL pointer at index #0 (the Cherenkov spectrum)
{ 
    
}

I3CLSimSpectrumTable::~I3CLSimSpectrumTable() 
{ 

}

std::size_t I3CLSimSpectrumTable::append(I3CLSimFunctionConstPtr newSpectrum)
{
    if (!newSpectrum)
        log_fatal("Cannot add a NULL spectrum.");
    
    // see if the spectrum already exists
    for (std::size_t i=1;i<spectra_.size();++i)
    {
        if (*(spectra_[i]) == *newSpectrum) // compare by value
        {
            return i;
        }
    }
    
    // not found, insert a reference to the spectrum
    spectra_.push_back(newSpectrum);
    
    return spectra_.size()-1;
}
