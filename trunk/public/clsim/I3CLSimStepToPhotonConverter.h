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
 * @file I3CLSimStepToPhotonConverter.h
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#ifndef I3CLSIMSTEPTOPHOTONCONVERTER_H_INCLUDED
#define I3CLSIMSTEPTOPHOTONCONVERTER_H_INCLUDED

#include "icetray/I3TrayHeaders.h"

#include "clsim/I3CLSimSimpleGeometry.h"

#include "clsim/I3CLSimStep.h"
#include "clsim/I3CLSimPhoton.h"
#include "clsim/I3CLSimPhotonHistory.h"
#include "clsim/I3CLSimMediumProperties.h"
#include "clsim/random_value/I3CLSimRandomValue.h"
#include "clsim/function/I3CLSimFunction.h"

#include <boost/noncopyable.hpp>

#include <vector>
#include <map>
#include <string>
#include <stdexcept>

/**
 * @brief Base class for objects that receive a vector of
 * I3CLSimSteps and convert them into I3CLSimPhotons.
 *
 * Depending on the implementation, these photons may
 * either be stored directly at their emission point or
 * after propagation to a target (e.g. a DOM)
 */

class I3CLSimStepToPhotonConverter_exception : public std::runtime_error
{
public:
    virtual ~I3CLSimStepToPhotonConverter_exception() throw() {;}

    I3CLSimStepToPhotonConverter_exception(const std::string &msg)
    :std::runtime_error(msg)
    {;}
};

struct I3CLSimStepToPhotonConverter : private boost::noncopyable
{
public:
    // simple struct that holds conversion results
    struct ConversionResult_t 
    {
        ConversionResult_t() : identifier(0) {;}
        ConversionResult_t(uint32_t identifier_,
                           I3CLSimPhotonSeriesPtr photons_=I3CLSimPhotonSeriesPtr(),
                           I3CLSimPhotonHistorySeriesPtr photonHistories_=I3CLSimPhotonHistorySeriesPtr())
        :
        identifier(identifier_),
        photons(photons_),
        photonHistories(photonHistories_)
        {;}
        
        uint32_t identifier;
        I3CLSimPhotonSeriesPtr photons;
        I3CLSimPhotonHistorySeriesPtr photonHistories;
    };
    
    //virtual ~I3CLSimStepToPhotonConverter();

    /**
     * Sets the wavelength generators. 
     * The first generator (index 0) is assumed to return a Cherenkov
     * spectrum that may have a bias applied to it. This bias factor
     * needs to be set using SetWlenBias().
     * All other generator indices are assumed to be for flasher/laser
     * light generation. During generation, no Cherenkov angle
     * rotation will be applied to those photons with indices >= 1.
     * Will throw if used after the call to Initialize().
     */
    virtual void SetWlenGenerators(const std::vector<I3CLSimRandomValueConstPtr> &wlenGenerators) = 0;

    /**
     * Sets the wavelength weights. Set this to a constant value
     * of 1 if you do not need biased photon generation.
     * The wavelength spectrum set with SetWlenGenerator()
     * is assumed to have a biassing factor already applied to it.
     * This call sets this factor in order to be able to assign
     * correct weights.
     * Will throw if used after the call to Initialize().
     */
    virtual void SetWlenBias(I3CLSimFunctionConstPtr wlenBias) = 0;

    /**
     * Sets the medium properties.
     * Will throw if used after the call to Initialize().
     */
    virtual void SetMediumProperties(I3CLSimMediumPropertiesConstPtr mediumProperties) = 0;

    /**
     * Sets the geometry.
     * Will throw if used after the call to Initialize().
     */
    virtual void SetGeometry(I3CLSimSimpleGeometryConstPtr geometry) = 0;
    
    /**
     * Initializes the simulation.
     * Will throw if already initialized.
     */
    virtual void Initialize() = 0;

    /**
     * Returns true if initialized.
     * Never throws.
     */
    virtual bool IsInitialized() const = 0;
    
    /**
     * Adds a new I3CLSimStepSeries to the queue.
     * The resulting I3CLSimPhotonSeries can be retrieved from the
     * I3CLSimStepToPhotonConverter after some processing time.
     *
     * Enqueuing a vector after calling EnqueueBarrier 
     * will throw if not all photons have been retrieved.
     *
     * Will throw if not initialized.
     */
    virtual void EnqueueSteps(I3CLSimStepSeriesConstPtr steps, uint32_t identifier) = 0;

    /**
     * Reports the current queue size. The queue works asynchronously,
     * so this value will probably have changed once you use it.
     *
     * Will throw if not initialized.
     */
    virtual std::size_t QueueSize() const = 0; 
    
    /**
     * Returns true if more photons are available.
     * If the return value is false, the current simulation is finished
     * and a new step vector may be set.
     * 
     * Will throw if not initialized.
     */
    virtual bool MorePhotonsAvailable() const = 0;

    /**
     * Returns a bunch of photons stored in a vector<I3CLSimPhoton>.
     *
     * The return value is a pair<uint, vector<I3CLSimPhoton> >.
     * The integer is the same identifier as specified in the call
     * to EnqueueSteps().
     *
     * Might block if no photons are available.
     * 
     * Will throw if not initialized.
     */
    virtual ConversionResult_t GetConversionResult() = 0;
    
protected:
};

I3_POINTER_TYPEDEFS(I3CLSimStepToPhotonConverter);

#endif //I3CLSIMSTEPTOPHOTONCONVERTER_H_INCLUDED
