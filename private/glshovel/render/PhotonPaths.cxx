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
#include <limits>
#include <deque>

#include <glshovel/Renderer.h>
#include <glshovel/ConfigurableColor.h>
#include <clsim/I3Photon.h>
#include <glshovel/render/PhotonPaths.h>
#include <boost/foreach.hpp>
#include <boost/format.hpp>

namespace render
{
    
    PhotonPaths::PhotonPaths()
    { 
        colorByTime_.init(true);
        Renderer::add_param("Assign photon colors by time, not by wavelength", colorByTime_);
    }
    
    template <class Archive>
    void 
    PhotonPaths::serialize(Archive& ar, unsigned version)
    {
        RENDERER_SAVE((colorByTime_));
        ar & make_nvp("base", base_object<base_t>(*this));
    }
    
    void
    PhotonPaths::clear()
    {
    }
    
    void
    PhotonPaths::scan()
    {
        mintime_ = NAN;
        maxtime_ = NAN;
        
        for (I3PhotonSeriesMap::const_iterator it=obj_->begin();
             it!=obj_->end();++it)
        {
            const I3PhotonSeries &photons = it->second;
            
            BOOST_FOREACH(const I3Photon& photon, photons)
            {
                if (isnan(mintime_)) {
                    mintime_ = photon.GetStartTime(); 
                } else {
                    if (photon.GetStartTime() < mintime_) 
                        mintime_=photon.GetStartTime();
                }

                if (isnan(maxtime_)) {
                    maxtime_ = photon.GetTime(); 
                    maxScatters_ = photon.GetNumScattered();
                } else {
                    if (photon.GetTime() > maxtime_) 
                        maxtime_=photon.GetTime();
                }
                
                
                if (photon.GetNumScattered() > maxScatters_)
                    maxScatters_=photon.GetNumScattered();

            }
        }
    }

    namespace {
        struct PhotonPathsStep {
            PhotonPathsStep() {;}
            double fromX, fromY, fromZ, fromT;
            double toX, toY, toZ, toT;
            bool skipped;
        };
    }
    
    void 
    PhotonPaths::render_one_photon(const I3Photon& p, 
                                     double fromtime, double totime)
    {
        const std::size_t numEntries = p.GetNumPositionListEntries();
        if (numEntries < 2) return;
        
        std::deque<PhotonPathsStep> steps;
        
        const double speed = p.GetGroupVelocity();
        
        I3PositionConstPtr to = p.GetPositionListEntry(numEntries-1);
        double to_time = p.GetTime();
        if (!to) return; // error
        
        I3PositionConstPtr from;
        double from_time;
        
        bool reachedGap=false;
        
        for (uint32_t i=numEntries-2;;--i)
        {
            I3PositionConstPtr from = p.GetPositionListEntry(i);
            if (!from) {
                reachedGap=true;
                continue; // no step
            }
            
            if ((reachedGap) && (i>0)) continue;
            
            if (i==0) {
                from_time = p.GetStartTime(); 
            } else {
                double distance = std::sqrt((to->GetX()-from->GetX())*(to->GetX()-from->GetX()) +
                                            (to->GetY()-from->GetY())*(to->GetY()-from->GetY()) +
                                            (to->GetZ()-from->GetZ())*(to->GetZ()-from->GetZ()) );
                double duration = distance/speed;
                from_time = to_time-duration;
            }
            
            steps.push_front(PhotonPathsStep());
            PhotonPathsStep &path = steps.front();
            path.fromX=from->GetX();
            path.fromY=from->GetY();
            path.fromZ=from->GetZ();
            path.fromT=from_time;
            path.toX=to->GetX();
            path.toY=to->GetY();
            path.toZ=to->GetZ();
            path.toT=to_time;
            path.skipped = reachedGap;
            
            
            to=from;
            to_time=from_time;
            
            if (i==0) break;
        }
        

        double wavelength = p.GetWavelength();
        if (isnan(wavelength)) wavelength=0.; // no nans here..
        
        const double minWavelength = 265.*I3Units::nanometer;
        const double maxWavelength = 675.*I3Units::nanometer;
        
        BOOST_FOREACH(const PhotonPathsStep &step, steps)
        {
            if (step.toT <= fromtime) continue;
            if (step.fromT >= totime) continue;
            
            double startTime = std::max(fromtime,step.fromT);
            double endTime   = std::min(totime,step.toT);
            
            double hue;
            if (colorByTime_) {
                hue = ( (step.toT+step.fromT)/2. - mintime_) / (maxtime_ - mintime_);
                if (hue > 1.) hue = 1.;
                if (hue < 0.) hue = 0.;
            } else {
                hue = 1. - (( wavelength - minWavelength) / (maxWavelength - minWavelength));
                if (hue > 1.) hue = 1.;
                if (hue < 0.) hue = 0.;
                
                // some gamma correction to make things more colorful
                hue = std::pow(hue, 2.);
            }
            
            
            const double d_x = step.toX-step.fromX;
            const double d_y = step.toY-step.fromY;
            const double d_z = step.toZ-step.fromZ;
            const double d_t_from = (startTime-step.fromT)/(step.toT-step.fromT);
            const double d_t_to   = (endTime  -step.fromT)/(step.toT-step.fromT);
            
            const double dispFromX = step.fromX + d_t_from*d_x;
            const double dispFromY = step.fromY + d_t_from*d_y;
            const double dispFromZ = step.fromZ + d_t_from*d_z;

            const double dispToX = step.fromX + d_t_to*d_x;
            const double dispToY = step.fromY + d_t_to*d_y;
            const double dispToZ = step.fromZ + d_t_to*d_z;

            
            QColor qc = QColor::fromHsvF(0.66666*hue, 1, 1);
            
            set_material(GL_EMISSION, qc);
            set_material(GL_SPECULAR, qc);
            set_material(GL_AMBIENT_AND_DIFFUSE,  qc);
            
            if (step.skipped) {
                glLineStipple(1, 0xAAAA);
                glEnable(GL_LINE_STIPPLE);
            }
            glBegin(GL_LINES);
            glVertex3d(dispFromX, dispFromY, dispFromZ);
            glVertex3d(dispToX, dispToY, dispToZ);
            glEnd();

            if (step.skipped) {
                glDisable(GL_LINE_STIPPLE);
            }
            
        }
        
        
    }
    
    void
    PhotonPaths::render(GLuint list,
                        double fromtime, 
                        double totime)
    {
        log_debug("starting PhotonPaths render of %s...", frameobject_name().c_str());
        
        glNewList(list, GL_COMPILE);
        
        const unsigned brightness = 192;
        set_material(GL_EMISSION, QColor(0, 0, 0, 0));
        set_material(GL_AMBIENT_AND_DIFFUSE, QColor(0, brightness, brightness, 255));
        set_material(GL_EMISSION, QColor(0, brightness, brightness, 255));
        
        for (I3PhotonSeriesMap::const_iterator it=obj_->begin();
             it!=obj_->end();++it)
        {
            const I3PhotonSeries &photons = it->second;
            
            BOOST_FOREACH(const I3Photon& photon, photons)
            {
                PhotonPaths::render_one_photon(photon, fromtime, totime);
            }
            
            
        }
        
        glEndList();
        
        log_debug("done.");
    }
    
    RENDERER(PhotonPaths);
    
}

RENDERER_SERIALIZABLE(render::PhotonPaths);
BOOST_CLASS_EXPORT(render::PhotonPaths);


