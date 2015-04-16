/**
 * Copyright (c) 2013
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
 * @file I3MCPEConverterWithIDs.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#include "I3MCPEConverterWithIDs.h"

void
I3MCPEConverterWithIDs::AddFields
(I3TableRowDescriptionPtr desc, const booked_type &)
{
    desc->AddField<double>  ("time", "ns", "time");
    desc->AddField<double>  ("npe", "PE", "The number of photoelectrons the hit represents.");
    desc->AddField<uint64_t>("partmajorid", "", "the particle major ID");
    desc->AddField<int>     ("partminorid", "", "the particle minor ID");
}

void
I3MCPEConverterWithIDs::FillSingleRow
(const I3MCPEConverterWithIDs::booked_type &hit, I3TableRowPtr row)
{
    row->Set<double>  ("time", hit.time);
    row->Set<double>  ("npe", hit.npe);
    row->Set<uint64_t>("partmajorid", hit.ID.majorID);
    row->Set<int>     ("partminorid", hit.ID.minorID);
}

//I3_CONVERTER(I3MCPESeriesMapConverterWithIDs, I3MCPESeriesMapWithIDs); 




