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
 * @file I3CLSimMediumProperties.h
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#ifndef I3CLSIMMEDIUMPROPERTIES_H_INCLUDED
#define I3CLSIMMEDIUMPROPERTIES_H_INCLUDED

#include "icetray/serialization.h"
#include "icetray/I3TrayHeaders.h"
#include "icetray/I3FrameObject.h"

#include <vector>

#include "clsim/function/I3CLSimFunction.h"
#include "clsim/function/I3CLSimScalarField.h"
#include "clsim/function/I3CLSimVectorTransform.h"
#include "clsim/random_value/I3CLSimRandomValue.h"

/**
 * @brief Describes the simulation medium (i.e. water or
 * ice) with all its properties like refractive index,
 * absorption length, scattering length, ..
 */
static const unsigned i3clsimmediumproperties_version_ = 2;

class I3CLSimMediumProperties : public I3FrameObject
{
public:
    static const double default_mediumDensity;
    static const uint32_t default_layersNum;
    static const double default_layersZStart;
    static const double default_layersHeight;
    static const double default_rockZCoordinate;
    static const double default_airZCoordinate;
    
    I3CLSimMediumProperties(double mediumDensity=default_mediumDensity,
                            uint32_t layersNum=default_layersNum, 
                            double layersZStart=default_layersZStart, 
                            double layersHeight=default_layersHeight,
                            double rockZCoordinate=default_rockZCoordinate,
                            double airZCoordinate=default_airZCoordinate);
    ~I3CLSimMediumProperties();

    /**
     * Returns true if media definitions for all layers are set.
     */
    bool IsReady() const;
    
    const std::vector<I3CLSimFunctionConstPtr> &GetAbsorptionLengths() const;
    const std::vector<I3CLSimFunctionConstPtr> &GetScatteringLengths() const;
    const std::vector<I3CLSimFunctionConstPtr> &GetPhaseRefractiveIndices() const;
    const std::vector<I3CLSimFunctionConstPtr> &GetGroupRefractiveIndicesOverride() const;
    
    I3CLSimFunctionConstPtr GetAbsorptionLength(uint32_t layer) const;
    I3CLSimFunctionConstPtr GetScatteringLength(uint32_t layer) const;
    I3CLSimFunctionConstPtr GetPhaseRefractiveIndex(uint32_t layer) const;
    I3CLSimFunctionConstPtr GetGroupRefractiveIndexOverride(uint32_t layer) const;
    I3CLSimRandomValueConstPtr GetScatteringCosAngleDistribution() const;
    I3CLSimScalarFieldConstPtr GetDirectionalAbsorptionLengthCorrection() const;
    I3CLSimVectorTransformConstPtr GetPreScatterDirectionTransform() const;
    I3CLSimVectorTransformConstPtr GetPostScatterDirectionTransform() const;
    I3CLSimScalarFieldConstPtr GetIceTiltZShift() const;

    void SetAbsorptionLength(uint32_t layer, I3CLSimFunctionConstPtr ptr);
    void SetScatteringLength(uint32_t layer, I3CLSimFunctionConstPtr ptr);
    void SetPhaseRefractiveIndex(uint32_t layer, I3CLSimFunctionConstPtr ptr);
    void SetGroupRefractiveIndexOverride(uint32_t layer, I3CLSimFunctionConstPtr ptr);
    void SetScatteringCosAngleDistribution(I3CLSimRandomValueConstPtr ptr);
    void SetDirectionalAbsorptionLengthCorrection(I3CLSimScalarFieldConstPtr ptr);
    void SetPreScatterDirectionTransform(I3CLSimVectorTransformConstPtr ptr);
    void SetPostScatterDirectionTransform(I3CLSimVectorTransformConstPtr ptr);
    void SetIceTiltZShift(I3CLSimScalarFieldConstPtr ptr);

    double GetMinWavelength() const;
    double GetMaxWavelength() const;
    
    // getters for constructor arguments
    inline double GetMediumDensity() const {return mediumDensity_;}
    inline uint32_t GetLayersNum() const {return layersNum_;}
    inline double GetLayersZStart() const {return layersZStart_;}
    inline double GetLayersHeight() const {return layersHeight_;}
    inline double GetRockZCoord() const {return rockZCoordinate_;}
    inline double GetAirZCoord() const {return airZCoordinate_;}
    
    // The user can specify a minimum and maximum wavelength.
    // This is necessary if none of the I3CLSimFunctions
    // specifies any.
    inline double GetForcedMinWlen() const {return forcedMinWlen_;}
    inline double GetForcedMaxWlen() const {return forcedMaxWlen_;}
    inline void SetForcedMinWlen(double val) {forcedMinWlen_=val;}
    inline void SetForcedMaxWlen(double val) {forcedMaxWlen_=val;}
    
    
private:
    double mediumDensity_;
    uint32_t layersNum_;
    double layersZStart_;
    double layersHeight_;
    
    double rockZCoordinate_;
    double airZCoordinate_;
    
    double forcedMinWlen_;
    double forcedMaxWlen_;
    
    std::vector<I3CLSimFunctionConstPtr> absorptionLength_;
    std::vector<I3CLSimFunctionConstPtr> scatteringLength_;
    std::vector<I3CLSimFunctionConstPtr> phaseRefractiveIndex_;
    std::vector<I3CLSimFunctionConstPtr> groupRefractiveIndexOverride_;
    I3CLSimRandomValueConstPtr scatteringCosAngleDist_;
    I3CLSimScalarFieldConstPtr directionalAbsorptionLengthCorrection_;
    I3CLSimVectorTransformConstPtr preScatterDirectionTransform_;
    I3CLSimVectorTransformConstPtr postScatterDirectionTransform_;
    I3CLSimScalarFieldConstPtr iceTiltZShift_;

private:
    friend class icecube::serialization::access;
    template <class Archive> void load(Archive & ar, unsigned version);
    template <class Archive> void save(Archive & ar, unsigned version) const;
    I3_SERIALIZATION_SPLIT_MEMBER();
};

I3_CLASS_VERSION(I3CLSimMediumProperties, i3clsimmediumproperties_version_);

I3_POINTER_TYPEDEFS(I3CLSimMediumProperties);

#endif //I3CLSIMMEDIUMPROPERTIES_H_INCLUDED
