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
 * @file I3CLSimFlasherPulse.h
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#ifndef I3CLSIMFLASHERPULSE_H_INCLUDED
#define I3CLSIMFLASHERPULSE_H_INCLUDED

#include "icetray/I3TrayHeaders.h"
#include "dataclasses/I3Vector.h"
#include "dataclasses/I3Direction.h"
#include "dataclasses/I3Position.h"

#include <string>
#include <vector>
#include <boost/shared_ptr.hpp>

#define I3CLSIMFLASHERPULSE_H_I3CLSimFlasherPulse_FlasherPulseType \
    (Unknown)(LED340nm)(LED370nm)(LED405nm)(LED450nm)(LED505nm)    \
    (SC1)(SC2)

static const unsigned i3clsimflasherpulse_version_ = 0;

/**
 * @brief Describes a single flasher LED (or SC) pulse.
 * I3FlasherInfo is not really suited for this, so we keep
 * our own structure here. One I3FlasherInfo object can
 * describe multiple LEDs at the same time (using a bitmask),
 * which is something we don't want here. (These LEDs could
 * have different wavelengths, ...)
 */
class I3CLSimFlasherPulse
{
public:
    ~I3CLSimFlasherPulse();
    I3CLSimFlasherPulse();

    enum FlasherPulseType {
        Unknown  = 0,
        LED340nm = 1,
        LED370nm = 2,
        LED405nm = 3,
        LED450nm = 4,
        LED505nm = 5,
        SC1      = 6,
        SC2      = 7
    };

    inline FlasherPulseType GetType() const {return flasherPulseType_;}
    inline void SetType(FlasherPulseType newType) {flasherPulseType_ = newType;}
    
    
    inline const I3Position &GetPos() const {return pos_;}
    inline const I3Direction &GetDir() const {return dir_;}

    inline void SetPos(const I3Position &newPos) {pos_=newPos;}
    inline void SetDir(const I3Direction &newDir) {dir_=newDir;}


    inline double GetTime() const {return time_;}
    inline void SetTime(double newTime) {time_=newTime;}

    inline double GetNumberOfPhotonsNoBias() const {return numberOfPhotonsNoBias_;}
    inline void SetNumberOfPhotonsNoBias(double newNumber) {numberOfPhotonsNoBias_=newNumber;}

    inline double GetPulseWidth() const {return pulseWidth_;}
    inline void SetPulseWidth(double newWidth) {pulseWidth_=newWidth;}

    inline double GetAngularEmissionSigmaPolar() const {return angularEmissionSigmaPolar_;}
    inline void SetAngularEmissionSigmaPolar(double newSigma) {angularEmissionSigmaPolar_=newSigma;}

    inline double GetAngularEmissionSigmaAzimuthal() const {return angularEmissionSigmaAzimuthal_;}
    inline void SetAngularEmissionSigmaAzimuthal(double newSigma) {angularEmissionSigmaAzimuthal_=newSigma;}

    std::ostream& Print(std::ostream&) const;
private:
    FlasherPulseType flasherPulseType_;

    I3Direction dir_;
    I3Position pos_;
    double time_;
    
    double numberOfPhotonsNoBias_;
    double pulseWidth_; // FWHM
    
    double angularEmissionSigmaPolar_;
    double angularEmissionSigmaAzimuthal_;
    

private: // static stuff
    friend class icecube::serialization::access;
    template <class Archive> void serialize(Archive & ar, unsigned version);

    friend bool operator==(const I3CLSimFlasherPulse &, const I3CLSimFlasherPulse &);
};
bool operator==(const I3CLSimFlasherPulse &a, const I3CLSimFlasherPulse &b);
std::ostream& operator<<(std::ostream&, const I3CLSimFlasherPulse&);

I3_CLASS_VERSION(I3CLSimFlasherPulse, i3clsimflasherpulse_version_);

typedef I3Vector<I3CLSimFlasherPulse> I3CLSimFlasherPulseSeries;

I3_POINTER_TYPEDEFS(I3CLSimFlasherPulse);
I3_POINTER_TYPEDEFS(I3CLSimFlasherPulseSeries);

#endif //I3CLSIMFLASHERPULSE_H_INCLUDED
