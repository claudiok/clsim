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
 * @file I3CLSimMediumProperties.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS
#endif
#include <inttypes.h>

#include <icetray/serialization.h>
#include <icetray/I3Units.h>
#include <clsim/I3CLSimMediumProperties.h>

#include <limits>


#include <boost/shared_ptr.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/foreach.hpp>

const uint32_t I3CLSimMediumProperties::default_layersNum=1;
const double I3CLSimMediumProperties::default_layersZStart=-5000.*I3Units::m;
const double I3CLSimMediumProperties::default_layersHeight=10000.*I3Units::m;
const double I3CLSimMediumProperties::default_mediumDensity=1.*I3Units::g/I3Units::cm3;
const double I3CLSimMediumProperties::default_rockZCoordinate=I3CLSimMediumProperties::default_layersZStart;
const double I3CLSimMediumProperties::default_airZCoordinate=-I3CLSimMediumProperties::default_layersZStart;

I3CLSimMediumProperties::I3CLSimMediumProperties(double mediumDensity,
                                                 uint32_t layersNum, 
                                                 double layersZStart, 
                                                 double layersHeight,
                                                 double rockZCoordinate,
                                                 double airZCoordinate)
:
mediumDensity_(mediumDensity),
layersNum_(layersNum),
layersZStart_(layersZStart),
layersHeight_(layersHeight),
rockZCoordinate_(rockZCoordinate),
airZCoordinate_(airZCoordinate),
forcedMinWlen_(-std::numeric_limits<double>::infinity()),
forcedMaxWlen_(std::numeric_limits<double>::infinity())
{
    if (layersNum_ <= 0) log_fatal("layersNum must be at least 1 during construction of I3CLSimMediumProperties.");
    
    if (layersZStart < rockZCoordinate_)
        log_fatal("layers start inside rock");

    if (layersZStart + static_cast<double>(layersNum)*layersHeight > airZCoordinate_)
        log_fatal("layers end inside air");

    absorptionLength_.assign(layersNum_, I3CLSimFunctionConstPtr());
    scatteringLength_.assign(layersNum_, I3CLSimFunctionConstPtr());
    phaseRefractiveIndex_.assign(layersNum_, I3CLSimFunctionConstPtr());
    groupRefractiveIndexOverride_.assign(layersNum_, I3CLSimFunctionConstPtr());
}

I3CLSimMediumProperties::~I3CLSimMediumProperties() 
{ 

}

double I3CLSimMediumProperties::GetMinWavelength() const
{
    if (!IsReady()) return NAN; // no valid answer if not ready

    double mini=forcedMinWlen_;
    
    BOOST_FOREACH(const I3CLSimFunctionConstPtr &ptr, absorptionLength_)
    {
        if (!ptr) return NAN;
        const double wlen=ptr->GetMinWlen();
        if (wlen > mini) mini=wlen;
        
    }
    BOOST_FOREACH(const I3CLSimFunctionConstPtr &ptr, scatteringLength_)
    {
        if (!ptr) return NAN;
        const double wlen=ptr->GetMinWlen();
        if (wlen > mini) mini=wlen;
    }
    BOOST_FOREACH(const I3CLSimFunctionConstPtr &ptr, phaseRefractiveIndex_)
    {
        if (!ptr) return NAN;
        const double wlen=ptr->GetMinWlen();
        if (wlen > mini) mini=wlen;
    }
    BOOST_FOREACH(const I3CLSimFunctionConstPtr &ptr, groupRefractiveIndexOverride_)
    {
        if (!ptr) continue; // these are optional
        const double wlen=ptr->GetMinWlen();
        if (wlen > mini) mini=wlen;
    }
    
    return mini;
}

double I3CLSimMediumProperties::GetMaxWavelength() const
{
    if (!IsReady()) return NAN; // no valid answer if not ready

    double maxi=forcedMaxWlen_;

    BOOST_FOREACH(const I3CLSimFunctionConstPtr &ptr, absorptionLength_)
    {
        if (!ptr) return NAN;
        const double wlen=ptr->GetMaxWlen();
        if (wlen < maxi) maxi=wlen;
        
    }
    BOOST_FOREACH(const I3CLSimFunctionConstPtr &ptr, scatteringLength_)
    {
        if (!ptr) return NAN;
        const double wlen=ptr->GetMaxWlen();
        if (wlen < maxi) maxi=wlen;
    }
    BOOST_FOREACH(const I3CLSimFunctionConstPtr &ptr, phaseRefractiveIndex_)
    {
        if (!ptr) return NAN;
        const double wlen=ptr->GetMaxWlen();
        if (wlen < maxi) maxi=wlen;
    }
    BOOST_FOREACH(const I3CLSimFunctionConstPtr &ptr, groupRefractiveIndexOverride_)
    {
        if (!ptr) continue; // these are optional
        const double wlen=ptr->GetMaxWlen();
        if (wlen < maxi) maxi=wlen;
    }
    
    return maxi;
}

bool I3CLSimMediumProperties::IsReady() const
{
    BOOST_FOREACH(const I3CLSimFunctionConstPtr &ptr, absorptionLength_)
    {
        if (!ptr) return false;
    }
    BOOST_FOREACH(const I3CLSimFunctionConstPtr &ptr, scatteringLength_)
    {
        if (!ptr) return false;
    }
    BOOST_FOREACH(const I3CLSimFunctionConstPtr &ptr, phaseRefractiveIndex_)
    {
        if (!ptr) return false;
    }

    // either all NULL or all configured
    if (!groupRefractiveIndexOverride_.empty())
    {
        if (!groupRefractiveIndexOverride_[0])
        {
            // all entries have to be NULL
            BOOST_FOREACH(const I3CLSimFunctionConstPtr &ptr, groupRefractiveIndexOverride_)
            {
                if (ptr) return false;
            }
        }
        else
        {
            // none of the entries must be NULL
            BOOST_FOREACH(const I3CLSimFunctionConstPtr &ptr, groupRefractiveIndexOverride_)
            {
                if (!ptr) return false;
            }
        }
        
    }
    
    if (!scatteringCosAngleDist_) return false;
    if (!directionalAbsorptionLengthCorrection_) return false;
    if (!preScatterDirectionTransform_) return false;
    if (!postScatterDirectionTransform_) return false;
    if (!iceTiltZShift_) return false;
    
    return true;
}

const std::vector<I3CLSimFunctionConstPtr> &I3CLSimMediumProperties::GetAbsorptionLengths() const
{
    return absorptionLength_;
}

const std::vector<I3CLSimFunctionConstPtr> &I3CLSimMediumProperties::GetScatteringLengths() const
{
    return scatteringLength_;
}

const std::vector<I3CLSimFunctionConstPtr> &I3CLSimMediumProperties::GetPhaseRefractiveIndices() const
{
    return phaseRefractiveIndex_;
}

const std::vector<I3CLSimFunctionConstPtr> &I3CLSimMediumProperties::GetGroupRefractiveIndicesOverride() const
{
    return groupRefractiveIndexOverride_;
}


I3CLSimFunctionConstPtr I3CLSimMediumProperties::GetAbsorptionLength(uint32_t layer) const
{
    if (layer >= layersNum_) log_fatal("Invalid layer num: %" PRIu32, layer);
    return absorptionLength_[layer];
}

I3CLSimFunctionConstPtr I3CLSimMediumProperties::GetScatteringLength(uint32_t layer) const
{
    if (layer >= layersNum_) log_fatal("Invalid layer num: %" PRIu32, layer);
    return scatteringLength_[layer];
}

I3CLSimFunctionConstPtr I3CLSimMediumProperties::GetPhaseRefractiveIndex(uint32_t layer) const
{
    if (layer >= layersNum_) log_fatal("Invalid layer num: %" PRIu32, layer);
    return phaseRefractiveIndex_[layer];
}

I3CLSimFunctionConstPtr I3CLSimMediumProperties::GetGroupRefractiveIndexOverride(uint32_t layer) const
{
    if (layer >= layersNum_) log_fatal("Invalid layer num: %" PRIu32, layer);
    return groupRefractiveIndexOverride_[layer];
}

I3CLSimRandomValueConstPtr I3CLSimMediumProperties::GetScatteringCosAngleDistribution() const
{
    return scatteringCosAngleDist_;
}

I3CLSimScalarFieldConstPtr I3CLSimMediumProperties::GetDirectionalAbsorptionLengthCorrection() const
{
    return directionalAbsorptionLengthCorrection_;
}

I3CLSimVectorTransformConstPtr I3CLSimMediumProperties::GetPreScatterDirectionTransform() const
{
    return preScatterDirectionTransform_;
}

I3CLSimVectorTransformConstPtr I3CLSimMediumProperties::GetPostScatterDirectionTransform() const
{
    return postScatterDirectionTransform_;
}

I3CLSimScalarFieldConstPtr I3CLSimMediumProperties::GetIceTiltZShift() const
{
    return iceTiltZShift_;
}

void I3CLSimMediumProperties::SetAbsorptionLength(uint32_t layer, I3CLSimFunctionConstPtr ptr)
{
    if (layer >= layersNum_) log_fatal("Invalid layer num: %" PRIu32, layer);
    absorptionLength_[layer]=ptr;
}

void I3CLSimMediumProperties::SetScatteringLength(uint32_t layer, I3CLSimFunctionConstPtr ptr)
{
    if (layer >= layersNum_) log_fatal("Invalid layer num: %" PRIu32, layer);
    scatteringLength_[layer]=ptr;
}

void I3CLSimMediumProperties::SetPhaseRefractiveIndex(uint32_t layer, I3CLSimFunctionConstPtr ptr)
{
    if (layer >= layersNum_) log_fatal("Invalid layer num: %" PRIu32, layer);
    phaseRefractiveIndex_[layer]=ptr;
}

void I3CLSimMediumProperties::SetGroupRefractiveIndexOverride(uint32_t layer, I3CLSimFunctionConstPtr ptr)
{
    if (layer >= layersNum_) log_fatal("Invalid layer num: %" PRIu32, layer);
    groupRefractiveIndexOverride_[layer]=ptr;
}

void I3CLSimMediumProperties::SetScatteringCosAngleDistribution(I3CLSimRandomValueConstPtr ptr)
{
    scatteringCosAngleDist_=ptr;
}

void I3CLSimMediumProperties::SetDirectionalAbsorptionLengthCorrection(I3CLSimScalarFieldConstPtr ptr)
{
    directionalAbsorptionLengthCorrection_=ptr;
}

void I3CLSimMediumProperties::SetPreScatterDirectionTransform(I3CLSimVectorTransformConstPtr ptr)
{
    preScatterDirectionTransform_=ptr;
}

void I3CLSimMediumProperties::SetPostScatterDirectionTransform(I3CLSimVectorTransformConstPtr ptr)
{
    postScatterDirectionTransform_=ptr;
}

void I3CLSimMediumProperties::SetIceTiltZShift(I3CLSimScalarFieldConstPtr ptr)
{
    iceTiltZShift_=ptr;
}


namespace {
    template <typename T, class Archive>
    void LoadFromArchiveIntoVectorConstPtr(Archive &ar, const std::string &name, std::vector<boost::shared_ptr<const T> > &arg)
    {
        std::vector<shared_ptr<T> > tmp;
        ar >> make_nvp(name.c_str(), tmp);
        arg.assign(tmp.begin(), tmp.end());
    }

    template <typename T, class Archive>
    void LoadFromArchiveIntoConstPtr(Archive &ar, const std::string &name, boost::shared_ptr<const T> &arg)
    {
        shared_ptr<T> tmp;
        ar >> make_nvp(name.c_str(), tmp);
        arg = tmp;
    }
}

template <class Archive>
void I3CLSimMediumProperties::load(Archive &ar, unsigned version)
{
    if (version>i3clsimmediumproperties_version_)
        log_fatal("Attempting to read version %u from file but running version %u of I3CLSimMediumProperties class.",version,i3clsimmediumproperties_version_);
    
    if (version>=2)
        ar >> make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));

    ar >> make_nvp("mediumDensity", mediumDensity_);
    ar >> make_nvp("layersNum", layersNum_);
    ar >> make_nvp("layersZStart", layersZStart_);
    ar >> make_nvp("layersHeight", layersHeight_);

    ar >> make_nvp("rockZCoordinate", rockZCoordinate_);
    ar >> make_nvp("airZCoordinate", airZCoordinate_);
    ar >> make_nvp("forcedMinWlen", forcedMinWlen_);
    ar >> make_nvp("forcedMaxWlen", forcedMaxWlen_);
    
    LoadFromArchiveIntoVectorConstPtr(ar, "absorptionLength", absorptionLength_);
    LoadFromArchiveIntoVectorConstPtr(ar, "scatteringLength", scatteringLength_);
    LoadFromArchiveIntoVectorConstPtr(ar, "phaseRefractiveIndex", phaseRefractiveIndex_);
    LoadFromArchiveIntoVectorConstPtr(ar, "groupRefractiveIndexOverride", groupRefractiveIndexOverride_);
    LoadFromArchiveIntoConstPtr(ar, "scatteringCosAngleDist", scatteringCosAngleDist_);

    if (version>=1) {
        LoadFromArchiveIntoConstPtr(ar, "directionalAbsorptionLengthCorrection", directionalAbsorptionLengthCorrection_);
        LoadFromArchiveIntoConstPtr(ar, "preScatterDirectionTransform", preScatterDirectionTransform_);
        LoadFromArchiveIntoConstPtr(ar, "postScatterDirectionTransform", postScatterDirectionTransform_);
        LoadFromArchiveIntoConstPtr(ar, "iceTiltZShift", iceTiltZShift_);
    } else {
        log_fatal("Cannot load version 0 of I3CLSimMediumProperties at this time.");
    }
}     


template <class Archive>
void I3CLSimMediumProperties::save(Archive &ar, unsigned version) const
{
    // version 2:
    ar << make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));

    // version 0:
    ar << make_nvp("mediumDensity", mediumDensity_);
    ar << make_nvp("layersNum", layersNum_);
    ar << make_nvp("layersZStart", layersZStart_);
    ar << make_nvp("layersHeight", layersHeight_);
    
    ar << make_nvp("rockZCoordinate", rockZCoordinate_);
    ar << make_nvp("airZCoordinate", airZCoordinate_);
    ar << make_nvp("forcedMinWlen", forcedMinWlen_);
    ar << make_nvp("forcedMaxWlen", forcedMaxWlen_);
    
    ar << make_nvp("absorptionLength", absorptionLength_);
    ar << make_nvp("scatteringLength", scatteringLength_);
    ar << make_nvp("phaseRefractiveIndex", phaseRefractiveIndex_);
    ar << make_nvp("groupRefractiveIndexOverride", groupRefractiveIndexOverride_);
    ar << make_nvp("scatteringCosAngleDist", scatteringCosAngleDist_);

    // version 1:
    ar << make_nvp("directionalAbsorptionLengthCorrection", directionalAbsorptionLengthCorrection_);
    ar << make_nvp("preScatterDirectionTransform", preScatterDirectionTransform_);
    ar << make_nvp("postScatterDirectionTransform", postScatterDirectionTransform_);
    ar << make_nvp("iceTiltZShift", iceTiltZShift_);

}     


I3_SPLIT_SERIALIZABLE(I3CLSimMediumProperties);
