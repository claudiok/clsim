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
 * @file I3CLSimEventStatistics.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#include <icetray/serialization.h>
#include <clsim/I3CLSimEventStatistics.h>

I3CLSimEventStatistics::~I3CLSimEventStatistics() { }

I3CLSimEventStatistics::I3CLSimEventStatistics()
{
    Reset();
}

template <class Archive>
void I3CLSimEventStatistics::serialize (Archive &ar, unsigned version)
{
    if (version > i3clsimeventstatistics_version_)
        log_fatal("Attempting to read version %u from file but running version %u of I3CLSimEventStatistics class.",version,i3clsimeventstatistics_version_);

    ar & make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));

    ar & make_nvp("numberOfPhotonsGeneratedPerParticle",numberOfPhotonsGeneratedPerParticle_);
    ar & make_nvp("sumOfWeightsPhotonsGeneratedPerParticle",sumOfWeightsPhotonsGeneratedPerParticle_);
    ar & make_nvp("totalNumberOfPhotonsGenerated",totalNumberOfPhotonsGenerated_);
    ar & make_nvp("totalSumOfWeightsPhotonsGenerated",totalSumOfWeightsPhotonsGenerated_);

    ar & make_nvp("numberOfPhotonsAtDOMsPerParticle",numberOfPhotonsAtDOMsPerParticle_);
    ar & make_nvp("sumOfWeightsPhotonsAtDOMsPerParticle",sumOfWeightsPhotonsAtDOMsPerParticle_);
    ar & make_nvp("totalNumberOfPhotonsAtDOMs",totalNumberOfPhotonsAtDOMs_);
    ar & make_nvp("totalSumOfWeightsPhotonsAtDOMs",totalSumOfWeightsPhotonsAtDOMs_);
}

I3_SERIALIZABLE(I3CLSimEventStatistics);
