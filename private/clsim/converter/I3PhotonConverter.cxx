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
 * @file I3PhotonConverter.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#include "I3PhotonConverter.h"

I3TableRowDescriptionPtr
I3PhotonConverter::CreateDescription
(const I3Photon &photon)
{
    I3TableRowDescriptionPtr desc(new I3TableRowDescription());
    
    desc->AddField<int32_t> ("id",                 "",       "photon id");
    desc->AddField<double>  ("weight",             "",       "photon weight");
    desc->AddField<uint64_t>("partmajorid",        "",       "negative log likelihood");
    desc->AddField<int32_t> ("partminorid",        "",       "negative log likelihood");

    desc->AddField<double>  ("cherenkov_distance", "m",      "full photon track length from emission to detection");
    desc->AddField<double>  ("cherenkov_time",     "ns",     "time difference between emission to detection");
    desc->AddField<double>  ("wavelength",         "nm",     "the photon wavelength");
    desc->AddField<double>  ("group_velocity",     "m/ns",   "the photon's group velocity");

    desc->AddField<double>  ("time",               "ns",     "arrival time of the photon on the DOM surface");
    desc->AddField<double>  ("hit_x",              "m",      "photon position on the DOM surface (x coordinate)");
    desc->AddField<double>  ("hit_y",              "m",      "photon position on the DOM surface (y coordinate)");
    desc->AddField<double>  ("hit_z",              "m",      "photon position on the DOM surface (z coordinate)");

    desc->AddField<double>  ("zen",                "radian", "zenith angle of the photon direction vector");
    desc->AddField<double>  ("azi",                "radian", "azimuthal angle of the photon direction vector");

    desc->AddField<double>  ("start_time",         "ns",     "time of photon emission");
    desc->AddField<double>  ("start_x",            "m",      "photon emission position (x coordinate)");
    desc->AddField<double>  ("start_y",            "m",      "photon emission position (y coordinate)");
    desc->AddField<double>  ("start_z",            "m",      "photon emission position (z coordinate)");

    desc->AddField<double>  ("start_zen",          "radian", "zenith angle of the photon direction vector at emission");
    desc->AddField<double>  ("start_azi",          "radian", "azimuthal angle of the photon direction vector at emission");

    desc->AddField<uint32_t>("num_scattered",      "",       "number of times this photon has been scattered");

    return desc;
}

std::size_t
I3PhotonConverter::FillRows
(const I3Photon &photon, I3TableRowPtr rows)
{
    rows->Set<int32_t> ("id",                 photon.GetID());
    rows->Set<double>  ("weight",             photon.GetWeight());
    rows->Set<uint64_t>("partmajorid",        photon.GetParticleMajorID());
    rows->Set<int32_t> ("partminorid",        photon.GetParticleMinorID());
    
    rows->Set<double>  ("cherenkov_distance", photon.GetCherenkovDist()/I3Units::m);
    rows->Set<double>  ("cherenkov_time",     photon.GetCherenkovTime()/I3Units::ns);
    rows->Set<double>  ("wavelength",         photon.GetWavelength()/I3Units::nanometer);
    rows->Set<double>  ("group_velocity",     photon.GetGroupVelocity()/(I3Units::m/I3Units::ns));
    
    rows->Set<double>  ("time",               photon.GetTime()/I3Units::ns);
    rows->Set<double>  ("hit_x",              photon.GetPos().GetX()/I3Units::m);
    rows->Set<double>  ("hit_y",              photon.GetPos().GetY()/I3Units::m);
    rows->Set<double>  ("hit_z",              photon.GetPos().GetZ()/I3Units::m);
    
    rows->Set<double>  ("zen",                photon.GetDir().GetZenith()/I3Units::rad);
    rows->Set<double>  ("azi",                photon.GetDir().GetAzimuth()/I3Units::rad);
    
    rows->Set<double>  ("start_time",         photon.GetStartTime()/I3Units::ns);
    rows->Set<double>  ("start_x",            photon.GetStartPos().GetX()/I3Units::m);
    rows->Set<double>  ("start_y",            photon.GetStartPos().GetY()/I3Units::m);
    rows->Set<double>  ("start_z",            photon.GetStartPos().GetZ()/I3Units::m);
    
    rows->Set<double>  ("start_zen",          photon.GetStartDir().GetZenith()/I3Units::rad);
    rows->Set<double>  ("start_azi",          photon.GetStartDir().GetAzimuth()/I3Units::rad);
    
    rows->Set<uint32_t>("num_scattered",      photon.GetNumScattered());
    
    return 1;
}

void I3PhotonConverter::AddFields(I3TableRowDescriptionPtr desc, const value_type& val)
{
    I3TableRowDescriptionPtr this_desc = CreateDescription(val);
    *desc << *this_desc;
}

void I3PhotonConverter::FillSingleRow(const value_type& val, I3TableRowPtr row)
{
    FillRows(val, row);
}

//I3_CONVERTER(I3PhotonConverter, I3Photon);
//I3_CONVERTER(I3PhotonSeriesMapConverter, I3PhotonSeriesMap); 




