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
 * @file I3MCHitConverterWithIDs.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#include "I3MCHitConverterWithIDs.h"

void
I3MCHitConverterWithIDs::AddFields
(I3TableRowDescriptionPtr desc, const booked_type &)
{
    desc->AddField<double>  ("time", "ns", "time");
    desc->AddField<double>  ("weight", "PE", "The number of photoelectrons the hit represents.");
    desc->AddField<double>  ("cherenkov_distance", "m", "FIXME: document");
    MAKE_ENUM_VECTOR(hit_source, ::I3MCHit, HitSource, I3MCHIT_H_I3MCHit_HitSource);
    desc->AddEnumField<I3MCHit::HitSource>("source",hit_source,"","");
    desc->AddField<int>     ("id", "", "hit id");
    desc->AddField<uint64_t>("partmajorid", "", "the particle major ID");
    desc->AddField<int>     ("partminorid", "", "the particle minor ID");
}

void
I3MCHitConverterWithIDs::FillSingleRow
(const I3MCHitConverterWithIDs::booked_type &hit, I3TableRowPtr row)
{
    row->Set<double>  ("time", hit.GetTime());
    row->Set<uint64_t>("weight", hit.GetNPE());
    row->Set<double>  ("cherenkov_distance", hit.GetCherenkovDistance());
    row->Set<I3MCHit::HitSource>("source", hit.GetHitSource());
    row->Set<int>     ("id", hit.GetHitID());
    row->Set<uint64_t>("partmajorid", hit.GetParticleMajorID());
    row->Set<int>     ("partminorid", hit.GetParticleMinorID());
}

//I3_CONVERTER(I3MCHitSeriesMapConverterWithIDs, I3MCHitSeriesMapWithIDs); 




