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
 * @file I3CLSimLightSourceToStepConverterGeant4.h
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#ifndef I3CLSIMLIGHTSOURCEPROPAGATORGEANT4_H_INCLUDED
#define I3CLSIMLIGHTSOURCEPROPAGATORGEANT4_H_INCLUDED

#include "icetray/I3Logging.h"
#include "clsim/I3CLSimLightSourcePropagator.h"
#include "clsim/I3CLSimQueue.h"
#include "clsim/I3CLSimLightSource.h"

#include <map>
#include <string>
#include <atomic>

/**
 * @brief A particle-to-step converter using Geant4
 * for tracking.
 */

class G4RunManager;

class I3CLSimLightSourcePropagatorGeant4 : public I3CLSimLightSourcePropagator
{
private:
    static std::atomic<bool> thereCanBeOnlyOneGeant4;
    
public:
    static const std::string default_physicsListName;
    static const double default_maxBetaChangePerStep;
    static const uint32_t default_maxNumPhotonsPerStep;
    
    I3CLSimLightSourcePropagatorGeant4(std::string physicsListName=default_physicsListName,
                                         double maxBetaChangePerStep=default_maxBetaChangePerStep,
                                         uint32_t maxNumPhotonsPerStep=default_maxNumPhotonsPerStep
                                         );
    virtual ~I3CLSimLightSourcePropagatorGeant4();

    // inherited:

    virtual void SetRandomService(I3RandomServicePtr random);

    /**
     * Sets the wavelength bias. Set this to a constant value
     * of 1 if you do not need biased photon generation.
     * The Cherenkov spectrum will be multiplied by this
     * value at each wavelength.
     * This will influence the number of photons produced.
     * Will throw if used after the call to Initialize().
     */
    virtual void SetWlenBias(I3CLSimFunctionConstPtr wlenBias);

    /**
     * Sets the medium properties.
     * Will throw if used after the call to Initialize().
     */
    virtual void SetMediumProperties(I3CLSimMediumPropertiesConstPtr mediumProperties);
    
    /**
     * Initializes the simulation.
     * Will throw if already initialized.
     */
    virtual void Initialize();

    /**
     * Returns true if initialized.
     * Never throws.
     */
    virtual bool IsInitialized() const;
    
    virtual bool IsValidForLightSource(const I3CLSimLightSource &source) { return source.GetType() == I3CLSimLightSource::Particle; };
    virtual void Convert(I3CLSimLightSourceConstPtr &, uint32_t, secondary_callback, step_callback);
    
private:
    void LogGeant4Messages(bool allAsWarn=false) const;

    mutable boost::shared_ptr<I3CLSimQueue<boost::shared_ptr<std::pair<const std::string, bool> > > > queueFromGeant4Messages_;
    
    I3RandomServicePtr randomService_;
    uint32_t randomSeed_;
    std::string physicsListName_;
    double maxBetaChangePerStep_;
    uint32_t maxNumPhotonsPerStep_;
    
    bool initialized_;
    
    std::unique_ptr<G4RunManager> runManager_;
    
    I3CLSimFunctionConstPtr wlenBias_;
    I3CLSimMediumPropertiesConstPtr mediumProperties_;
    
    SET_LOGGER("I3CLSimLightSourcePropagatorGeant4");
};

I3_POINTER_TYPEDEFS(I3CLSimLightSourcePropagatorGeant4);

#endif //I3CLSIMLIGHTSOURCEPROPAGATORGEANT4_H_INCLUDED
