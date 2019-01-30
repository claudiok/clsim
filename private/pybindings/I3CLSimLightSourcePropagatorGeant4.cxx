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
 * @file I3CLSimLightSourceToStepConverter.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#include <geant4/I3CLSimLightSourcePropagatorGeant4.h>

namespace bp = boost::python;

void register_I3CLSimLightSourcePropagatorGeant4()
{
    bp::class_<I3CLSimLightSourcePropagatorGeant4,
               boost::shared_ptr<I3CLSimLightSourcePropagatorGeant4>,
               bp::bases<I3CLSimLightSourcePropagator>,
               boost::noncopyable>(
                 "I3CLSimLightSourcePropagatorGeant4",
                 bp::init<std::string,double,uint32_t>((
                   bp::arg("physicsListName")=I3CLSimLightSourcePropagatorGeant4::default_physicsListName,
                   bp::arg("maxBetaChangePerStep")=I3CLSimLightSourcePropagatorGeant4::default_maxBetaChangePerStep,
                   bp::arg("maxNumPhotonsPerStep")=I3CLSimLightSourcePropagatorGeant4::default_maxNumPhotonsPerStep
                ))
    );
}
