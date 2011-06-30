/**
 * copyright  (C) 2011
 *
 * $Id$
 *
 * @version $Revision$
 * @date $LastChangedDate$
 * @author Claudio Kopper <claudio.kopper@nikhef.nl>
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
    row->Set<double>  ("weight", hit.GetWeight());
    row->Set<double>  ("cherenkov_distance", hit.GetCherenkovDistance());
    row->Set<I3MCHit::HitSource>("source", hit.GetHitSource());
    row->Set<int>     ("id", hit.GetHitID());
    row->Set<uint64_t>("partmajorid", hit.GetParticleMajorID());
    row->Set<int>     ("partminorid", hit.GetParticleMinorID());
}

I3_CONVERTER(I3MCHitSeriesMapConverterWithIDs, I3MCHitSeriesMapWithIDs); 




