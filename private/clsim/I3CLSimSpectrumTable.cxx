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
 * @file I3CLSimSpectrumTable.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

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
