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
 * @file I3MCHitConverterWithIDs.h
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

/// The default tableio I3MCHit converter does not store te hit
/// id. This is a similar converter that also stores hit ids
/// and particle ids.

#include "tableio/I3Converter.h"
#include "tableio/converter/I3MapConverter.h"

#include "dataclasses/physics/I3MCHit.h"

struct I3MCHitConverterWithIDs 
{
    typedef I3MCHit booked_type;
    typedef booked_type value_type;

    void AddFields(I3TableRowDescriptionPtr desc, const booked_type& = booked_type());
    void FillSingleRow(const booked_type& dl, I3TableRowPtr row);
};

typedef I3MapOMKeyVectorConverter<I3MCHitConverterWithIDs> I3MCHitSeriesMapConverterWithIDs;
