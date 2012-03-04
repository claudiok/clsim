//
//   Copyright (c) 2012  Claudio Kopper
//   
//   $Id$
//
//   This file is part of the IceTray module clsim.
//
//   clsim is free software; you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation; either version 3 of the License, or
//   (at your option) any later version.
//
//   IceTray is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#ifndef I3CLSIMSPECTRUMTABLE_H_INCLUDED
#define I3CLSIMSPECTRUMTABLE_H_INCLUDED

#include "icetray/I3TrayHeaders.h"

#include "clsim/function/I3CLSimFunction.h"

#include <vector>
#include <stdexcept>

/**
 * @brief Stores all spectra that may be used during photon generation.
 * If a spectrum that is being added should already exist, it is not 
 * added a second time. The previous version is re-used.
 *
 * Spectrum number "0" is always the Cherenkov spectrum and does
 * not need to be added by the user.
 */
struct I3CLSimSpectrumTable 
{
public:
    I3CLSimSpectrumTable();
    ~I3CLSimSpectrumTable();

    /**
     * Returns the vector of spectra. The first spectrum is the Cherenkov
     * spectrum and will be set to NULL.
     */
    const std::vector<I3CLSimFunctionConstPtr> &GetSpectra() const
    {
        return spectra_;
    }
    
    /**
     * Add a new spectrum. The spectrum index will be returned.
     * This can also be used to retrieve the index of an existing spectrum.
     */
    std::size_t append(I3CLSimFunctionConstPtr newSpectrum);

    I3CLSimFunctionConstPtr operator[](std::size_t index) const
    {
        if (index>=spectra_.size())
            throw std::out_of_range("invalid index");
        return spectra_[index];
    }

    
    std::size_t size() const {return spectra_.size();}
    
private:
    std::vector<I3CLSimFunctionConstPtr> spectra_;
    
};

I3_POINTER_TYPEDEFS(I3CLSimSpectrumTable);

#endif //I3CLSIMSPECTRUMTABLE_H_INCLUDED
