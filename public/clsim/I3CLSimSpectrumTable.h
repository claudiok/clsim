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
 * @file I3CLSimSpectrumTable.h
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

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
