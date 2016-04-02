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
 * @file I3CLSimLightSourceToStepConverter.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#include <icetray/serialization.h>
#include <clsim/I3CLSimLightSourceToStepConverter.h>

I3CLSimLightSourceToStepConverter::I3CLSimLightSourceToStepConverter() {;}
I3CLSimLightSourceToStepConverter::~I3CLSimLightSourceToStepConverter() {;}

void I3CLSimLightSourceToStepConverter::SetLightSourceParameterizationSeries
(const I3CLSimLightSourceParameterizationSeries &parameterizationSeries_)
{
    if (IsInitialized())
        throw I3CLSimLightSourceToStepConverter_exception("SetLightSourceParameterizationSeries() called after Initialize().");
        
    parameterizationSeries = parameterizationSeries_;
}

const I3CLSimLightSourceParameterizationSeries &
I3CLSimLightSourceToStepConverter::GetLightSourceParameterizationSeries() const
{ 
    return parameterizationSeries;
}

I3CLSimStepSeriesConstPtr I3CLSimLightSourceToStepConverter::GetConversionResult(double timeout)
{
    bool dummy;
    return GetConversionResultWithBarrierInfo(dummy, timeout);
}
