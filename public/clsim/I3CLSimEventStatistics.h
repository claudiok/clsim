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
 * @file I3CLSimEventStatistics.h
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#ifndef I3CLSIMEVENTSTATISTICS_H_INCLUDED
#define I3CLSIMEVENTSTATISTICS_H_INCLUDED

#include "icetray/I3FrameObject.h"

#include <map>
#include <vector>

#include "dataclasses/physics/I3Particle.h"

static const unsigned i3clsimeventstatistics_version_ = 0;

/**
 * @brief This class collects statistics/information
 * on simulated events (for example the total number
 * of generated photons).
 */
class I3CLSimEventStatistics : public I3FrameObject
{
public:
    I3CLSimEventStatistics();
    virtual ~I3CLSimEventStatistics();

    
    //// NUMBER OF GENERATED PHOTONS
    inline uint64_t GetNumberOfPhotonsGeneratedForParticle(uint64_t majorID, int minorID) const
    {
        std::map<std::pair<uint64_t, int>, uint64_t>::const_iterator it = 
        numberOfPhotonsGeneratedPerParticle_.find(std::pair<uint64_t, int>(majorID, minorID));
        
        if (it==numberOfPhotonsGeneratedPerParticle_.end()) {
            log_error("Particle not found");
            return 0;
        }
        
        return it->second;
    }
    inline uint64_t GetNumberOfPhotonsGeneratedForParticle(const I3Particle &particle) const
    {
        return GetNumberOfPhotonsGeneratedForParticle(particle.GetMajorID(), particle.GetMinorID());
    }
    inline uint64_t GetTotalNumberOfPhotonsGenerated() const
    {
        return totalNumberOfPhotonsGenerated_;
    }

    
    //// WEIGHT SUM OF GENERATED PHOTONS
    inline double GetSumOfWeightsPhotonsGeneratedForParticle(uint64_t majorID, int minorID) const
    {
        std::map<std::pair<uint64_t, int>, double>::const_iterator it = 
        sumOfWeightsPhotonsGeneratedPerParticle_.find(std::pair<uint64_t, int>(majorID, minorID));
        
        if (it==sumOfWeightsPhotonsGeneratedPerParticle_.end()) {
            log_error("Particle not found");
            return 0;
        }
        
        return it->second;
    }
    inline double GetSumOfWeightsPhotonsGeneratedForParticle(const I3Particle &particle) const
    {
        return GetSumOfWeightsPhotonsGeneratedForParticle(particle.GetMajorID(), particle.GetMinorID());
    }
    inline double GetTotalSumOfWeightsPhotonsGenerated() const
    {
        return totalSumOfWeightsPhotonsGenerated_;
    }
    
    
    
    
    //// NUMBER OF PHOTONS AT DOMs
    inline uint64_t GetNumberOfPhotonsAtDOMsForParticle(uint64_t majorID, int minorID) const
    {
        std::map<std::pair<uint64_t, int>, uint64_t>::const_iterator it = 
        numberOfPhotonsAtDOMsPerParticle_.find(std::pair<uint64_t, int>(majorID, minorID));
        
        if (it==numberOfPhotonsAtDOMsPerParticle_.end()) {
            log_error("Particle not found");
            return 0;
        }
        
        return it->second;
    }
    inline uint64_t GetNumberOfPhotonsAtDOMsForParticle(const I3Particle &particle) const
    {
        return GetNumberOfPhotonsAtDOMsForParticle(particle.GetMajorID(), particle.GetMinorID());
    }
    inline uint64_t GetTotalNumberOfPhotonsAtDOMs() const
    {
        return totalNumberOfPhotonsAtDOMs_;
    }
    
    
    //// WEIGHT SUM OF PHOTONS AT DOMs
    inline double GetSumOfWeightsPhotonsAtDOMsForParticle(uint64_t majorID, int minorID) const
    {
        std::map<std::pair<uint64_t, int>, double>::const_iterator it = 
        sumOfWeightsPhotonsAtDOMsPerParticle_.find(std::pair<uint64_t, int>(majorID, minorID));
        
        if (it==sumOfWeightsPhotonsAtDOMsPerParticle_.end()) {
            log_error("Particle not found");
            return 0;
        }
        
        return it->second;
    }
    inline double GetSumOfWeightsPhotonsAtDOMsForParticle(const I3Particle &particle) const
    {
        return GetSumOfWeightsPhotonsAtDOMsForParticle(particle.GetMajorID(), particle.GetMinorID());
    }
    inline double GetTotalSumOfWeightsPhotonsAtDOMs() const
    {
        return totalSumOfWeightsPhotonsAtDOMs_;
    }
    
    
    //// ADD GENERATED PHOTONS
    inline void AddNumPhotonsGeneratedWithWeights(uint64_t numPhotons, double weightsForPhotons,
                                                  uint64_t majorID, int minorID)
    {
        std::pair<uint64_t, int> index(majorID, minorID);
        
        (numberOfPhotonsGeneratedPerParticle_.insert(std::make_pair(index, 0)).first->second)+=numPhotons;
        (sumOfWeightsPhotonsGeneratedPerParticle_.insert(std::make_pair(index, 0.)).first->second)+=weightsForPhotons;
        
        totalNumberOfPhotonsGenerated_+=numPhotons;
        totalSumOfWeightsPhotonsGenerated_+=weightsForPhotons;
    }

    inline void AddNumPhotonsGeneratedWithWeights(uint64_t numParticles, double weightsForParticles,
                                         const I3Particle &particle)
    {
        AddNumPhotonsGeneratedWithWeights(numParticles, weightsForParticles,
                                          particle.GetMajorID(), particle.GetMinorID());
    }

    //// ADD PHOTONS AT DOMs
    inline void AddNumPhotonsAtDOMsWithWeights(uint64_t numPhotons, double weightsForPhotons,
                                               uint64_t majorID, int minorID)
    {
        std::pair<uint64_t, int> index(majorID, minorID);
        
        (numberOfPhotonsAtDOMsPerParticle_.insert(std::make_pair(index, 0)).first->second)+=numPhotons;
        (sumOfWeightsPhotonsAtDOMsPerParticle_.insert(std::make_pair(index, 0.)).first->second)+=weightsForPhotons;
        
        totalNumberOfPhotonsAtDOMs_+=numPhotons;
        totalSumOfWeightsPhotonsAtDOMs_+=weightsForPhotons;
    }
    
    inline void AddNumPhotonsAtDOMsWithWeights(uint64_t numParticles, double weightsForParticles,
                                         const I3Particle &particle)
    {
        AddNumPhotonsAtDOMsWithWeights(numParticles, weightsForParticles,
                                       particle.GetMajorID(), particle.GetMinorID());
    }

    
    inline void Reset()
    {
        numberOfPhotonsGeneratedPerParticle_.clear();
        sumOfWeightsPhotonsGeneratedPerParticle_.clear();
        totalNumberOfPhotonsGenerated_=0;
        totalSumOfWeightsPhotonsGenerated_=0.;

        numberOfPhotonsAtDOMsPerParticle_.clear();
        sumOfWeightsPhotonsAtDOMsPerParticle_.clear();
        totalNumberOfPhotonsAtDOMs_=0;
        totalSumOfWeightsPhotonsAtDOMs_=0.;
    }
    
    
private:
    std::map<std::pair<uint64_t, int>, uint64_t> numberOfPhotonsGeneratedPerParticle_;
    std::map<std::pair<uint64_t, int>, double> sumOfWeightsPhotonsGeneratedPerParticle_;
    uint64_t totalNumberOfPhotonsGenerated_;    
    double totalSumOfWeightsPhotonsGenerated_;

    std::map<std::pair<uint64_t, int>, uint64_t> numberOfPhotonsAtDOMsPerParticle_;
    std::map<std::pair<uint64_t, int>, double> sumOfWeightsPhotonsAtDOMsPerParticle_;
    uint64_t totalNumberOfPhotonsAtDOMs_;    
    double totalSumOfWeightsPhotonsAtDOMs_;

    friend class icecube::serialization::access;
    template <class Archive> void serialize(Archive & ar, unsigned version);
};

I3_CLASS_VERSION(I3CLSimEventStatistics, i3clsimeventstatistics_version_);

I3_POINTER_TYPEDEFS(I3CLSimEventStatistics);

#endif //I3CLSIMEVENTSTATISTICS_H_INCLUDED
