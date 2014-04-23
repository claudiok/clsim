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
 * @file I3PhotonConverter.h
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#include "tableio/I3Converter.h"
#include "tableio/converter/I3MapConverter.h"

#include "clsim/I3Photon.h"

class I3PhotonConverter : public I3ConverterImplementation<I3Photon>
{
public:
    typedef booked_type value_type;

    void AddFields(I3TableRowDescriptionPtr, const value_type&  = value_type());
    void FillSingleRow(const value_type&, I3TableRowPtr);

private:
    I3TableRowDescriptionPtr CreateDescription(const I3Photon &photon); 
    std::size_t FillRows(const I3Photon &photon, I3TableRowPtr rows);
};

typedef I3MapModuleKeyVectorConverter<I3PhotonConverter, I3PhotonSeriesMap> I3PhotonSeriesMapConverter;
