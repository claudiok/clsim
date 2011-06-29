/**
 * copyright  (C) 2011
 *
 * $Id$
 *
 * @version $Revision$
 * @date $LastChangedDate$
 * @author Claudio Kopper <claudio.kopper@nikhef.nl>
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
    rows->Set<double>  ("start_y",            photon.GetStartPos().GetX()/I3Units::m);
    rows->Set<double>  ("start_z",            photon.GetStartPos().GetX()/I3Units::m);
    
    rows->Set<double>  ("start_zen",          photon.GetStartDir().GetZenith()/I3Units::rad);
    rows->Set<double>  ("start_azi",          photon.GetStartDir().GetAzimuth()/I3Units::rad);
    
    rows->Set<uint32_t>("num_scattered",      photon.GetNumScattered());
    
    return 1;
}

I3_CONVERTER(I3PhotonConverter, I3Photon);
I3_CONVERTER(I3PhotonSeriesMapConverter, I3PhotonSeriesMap); 




