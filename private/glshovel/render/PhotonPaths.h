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
 * @file PhotonPaths.h
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#ifndef GLSHOVEL_RENDER_PHOTONPATHS_H_INCLUDED
#define GLSHOVEL_RENDER_PHOTONPATHS_H_INCLUDED

// Qt defines signal as "protected" using the preprocessor.
// This breaks boost::signal, so we need to include it first here
// in order to have everything defined before including Qt.
#include <boost/signals.hpp>

#include "clsim/I3Photon.h"

#include <QWidget>
#ifdef signals
#undef signals
#endif

#include <glshovel/gl.h>
#include <glshovel/util.h>
#include <glshovel/render_dispatch.h>
#include <glshovel/Configurable.h>


namespace render 
{
    class PhotonPaths : public render_dispatch<I3PhotonSeriesMap, PhotonPaths>
    {
        //ConfigurableDouble diameter_;
        //ConfigurableDouble shininess_;
        //ConfigurableColor color_;
        ConfigurableBool colorByTime_;
        
        double mintime_, maxtime_;
        unsigned long maxScatters_;
        //GLuint afterlist_;
        
        //HitList hit_list_;
        
    public:
        
        typedef render_dispatch<I3PhotonSeriesMap, PhotonPaths> base_t;
        
        void 
        render_one_photon(const I3Photon& p, double fromtime, double totime);
        
        PhotonPaths();
        void scan();
        void clear();
        void render(GLuint list, 
                    double fromtime,
                    double totime);
        
        typedef boost::false_type want_viewport_changes;
        
        template <class Archive>
        void
        serialize(Archive& ar, unsigned version);
        
    };
    
}
#endif

