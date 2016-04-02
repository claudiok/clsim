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
 * @file I3CLSimSimpleGeometryUserConfigurable.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#include "clsim/I3CLSimSimpleGeometryUserConfigurable.h"


I3CLSimSimpleGeometryUserConfigurable::
I3CLSimSimpleGeometryUserConfigurable(double OMRadius,
                                      std::size_t numOMs)
:
OMRadius_(OMRadius),
numOMs_(numOMs),
stringIDs_(numOMs, 0),
domIDs_(numOMs, 0),
posX_(numOMs, NAN),
posY_(numOMs, NAN),
posZ_(numOMs, NAN),
subdetectors_(numOMs, "")
{
    
}

I3CLSimSimpleGeometryUserConfigurable::
~I3CLSimSimpleGeometryUserConfigurable()
{
    
}
