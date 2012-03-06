/**
 *    $Id$
 *
 *    Copyright (C) 2012   Claudio Kopper <claudio.kopper@icecube.wisc.edu>
 *    and the IceCube Collaboration <http://www.icecube.wisc.edu>
 *    
 *    This file is free software; you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation; either version 3 of the License, or
 *    (at your option) any later version.
 *    
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *    
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef GLSHOVEL_RENDER_PHOTONPATHS_H_INCLUDED
#define GLSHOVEL_RENDER_PHOTONPATHS_H_INCLUDED

#include <QWidget>
#include <glshovel/gl.h>
#include <glshovel/util.h>
#include <glshovel/render_dispatch.h>
#include <dataclasses/geometry/I3Geometry.h>
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

