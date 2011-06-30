/**
 * copyright  (C) 2011
 *
 * $Id$
 *
 * @version $Revision$
 * @date $LastChangedDate$
 * @author Claudio Kopper <claudio.kopper@nikhef.nl>
 */

/// The default tableio I3MCHit converter does not store te hit
/// id. This is a similar converter that also stores hit ids
/// and particle ids.

#include "tableio/I3ConverterFactory.h"
#include "tableio/converter/I3MapConverter.h"

#include "dataclasses/physics/I3MCHit.h"

struct I3MCHitConverterWithIDs 
{
    typedef I3MCHit booked_type;

    void AddFields(I3TableRowDescriptionPtr desc, const booked_type& = booked_type());
    void FillSingleRow(const booked_type& dl, I3TableRowPtr row);
};

typedef I3MapOMKeyVectorConverter<I3MCHitConverterWithIDs> I3MCHitSeriesMapConverterWithIDs;
