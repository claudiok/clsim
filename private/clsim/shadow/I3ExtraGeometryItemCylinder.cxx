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
 * @file I3ExtraGeometryItemCylinder.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#include <limits>

#include <icetray/serialization.h>
#include <icetray/I3Units.h>
#include <clsim/shadow/I3ExtraGeometryItemCylinder.h>

#include <boost/foreach.hpp>

I3ExtraGeometryItemCylinder::I3ExtraGeometryItemCylinder()
:
radius_(NAN),
boundingBoxCalculated_(false)
{;}

I3ExtraGeometryItemCylinder::
I3ExtraGeometryItemCylinder(const I3Position &from, const I3Position &to, double radius)
:
from_(from),
to_(to),
radius_(radius),
boundingBoxCalculated_(false)
{
}

I3ExtraGeometryItemCylinder::~I3ExtraGeometryItemCylinder() { }

void I3ExtraGeometryItemCylinder::CalculateBoundingBox() const
{
    if (boundingBoxCalculated_) return;
    
    
    boundingBoxLower_=I3Position();
    boundingBoxUpper_=I3Position();
    
    boundingBoxCalculated_=true;
}

unsigned int I3ExtraGeometryItemCylinder::FindIntersections(const I3Position &lineStart,
                                                            const I3Position &lineEnd,
                                                            double &p0, double &p1,
                                                            double &line_segment_length
                                                            ) const
{
    // This function is based on code from http://www.geometrictools.com/
    // Geometric Tools, LLC
    // Copyright (c) 1998-2011
    // Distributed under the Boost Software License, Version 1.0.
    // http://www.boost.org/LICENSE_1_0.txt
    
    double line_dir_x = lineEnd.GetX()-lineStart.GetX();
    double line_dir_y = lineEnd.GetY()-lineStart.GetY();
    double line_dir_z = lineEnd.GetZ()-lineStart.GetZ();
    
    line_segment_length = std::sqrt(line_dir_x*line_dir_x + line_dir_y*line_dir_y + line_dir_z*line_dir_z);
    line_dir_x /= line_segment_length;
    line_dir_y /= line_segment_length;
    line_dir_z /= line_segment_length;
    
    
    // the cylinder direction
    double Wx = to_.GetX() - from_.GetX();
    double Wy = to_.GetY() - from_.GetY();
    double Wz = to_.GetZ() - from_.GetZ();
    double halfHeight;
    {
        const double W_len = std::sqrt(Wx*Wx + Wy*Wy + Wz*Wz);
        Wx /= W_len; Wy /= W_len; Wz /= W_len;
        halfHeight = W_len/2.; // the cylinder height
    }
    const double rSqr = radius_*radius_;
    
    // two orthogonal directions
    double Ux, Uy, Uz, Vx, Vy, Vz;
    {
        const double x = Wx < 0.0 ? -Wx : Wx;
        const double y = Wy < 0.0 ? -Wy : Wy;
        const double z = Wz < 0.0 ? -Wz : Wz;
        if (x < y) {
            if (x < z) {
                Ux = 0.; Uy = Wz; Uz=-Wy;
            } else {
                Ux = Wy; Uy = -Wx; Uz=0.;
            }
        } else {
            if (y < z) {
                Ux = -Wz; Uy = 0.; Uz=Wx;
            } else {
                Ux = Wy; Uy = -Wx; Uz=0.;
            }
        }
        
        const double U_len = std::sqrt(Ux*Ux + Uy*Uy + Uz*Uz);
        Ux /= U_len; Uy /= U_len; Uz /= U_len;
        
        // V = W x U  =>  U x V = W
        Vx = Wy*Uz - Wz*Uy;
        Vy = Wz*Ux - Wx*Uz;
        Vz = Wx*Uy - Wy*Ux;
    }
    
    // diff from line origin to cylinder center
    const double cylOrigin_x = from_.GetX() + Wx;
    const double cylOrigin_y = from_.GetY() + Wy;
    const double cylOrigin_z = from_.GetZ() + Wz;
    
    const double diff_x = lineStart.GetX() - cylOrigin_x;
    const double diff_y = lineStart.GetY() - cylOrigin_y;
    const double diff_z = lineStart.GetZ() - cylOrigin_z;
    
    // convert to local cylinder coordinates
    const double Px = Ux*diff_x + Uy*diff_y + Uz*diff_z;
    const double Py = Vx*diff_x + Vy*diff_y + Vz*diff_z;
    const double Pz = Wx*diff_x + Wy*diff_y + Wz*diff_z;
    
    // Get the z-value, in cylinder coordinates, of the incoming line's
    // unit-length direction.
    const double dz = Wx*line_dir_x + Wy*line_dir_y + Wz*line_dir_z;
    
    const double EPSILON=1e-8;
    
    if (std::fabs(dz) >= 1. - EPSILON)
    {
        // The line is parallel to the cylinder axis.  Determine if the line
        // intersects the cylinder end disks.
        double radialSqrDist = rSqr - Px*Px - Py*Py;
        if (radialSqrDist < 0.)
        {
            // Line outside the cylinder, no intersection.
            return 0;
        }
        
        // Line intersects the cylinder end disks.
        if (dz > 0.)
        {
            p0 = -Pz - halfHeight;
            p1 = -Pz + halfHeight;
        }
        else
        {
            p0 = Pz - halfHeight;
            p1 = Pz + halfHeight;
        }
        return 2;
    }

    // convert incoming line unit-length direction to cylinder coordinates
    const double Dx = Ux*line_dir_x + Uy*line_dir_y + Uz*line_dir_z;
    const double Dy = Vx*line_dir_x + Vy*line_dir_y + Vz*line_dir_z;
    const double Dz = dz;
    
    double a0, a1, a2, discr, root, inv, tValue;
    
    if (std::fabs(Dz) <= EPSILON)
    {
        // The line is perpendicular to the cylinder axis.
        if (std::fabs(Pz) > halfHeight)
        {
            // Line is outside the planes of the cylinder end disks.
            return 0;
        }
        
        // Test intersection of line P+t*D with infinite cylinder
        // x^2+y^2 = r^2.  This reduces to computing the roots of a
        // quadratic equation.  If P = (px,py,pz) and D = (dx,dy,dz),
        // then the quadratic equation is
        //   (dx^2+dy^2)*t^2 + 2*(px*dx+py*dy)*t + (px^2+py^2-r^2) = 0
        a0 = Px*Px + Py*Py - rSqr;
        a1 = Px*Dx + Py*Dy;
        a2 = Dx*Dx + Dy*Dy;
        discr = a1*a1 - a0*a2;
        if (discr < 0.)
        {
            // Line does not intersect cylinder.
            return 0;
        }
        else if (discr > EPSILON)
        {
            // Line intersects cylinder in two places.
            root = std::sqrt(discr);
            inv = 1./a2;
            p0 = (-a1 - root)*inv;
            p1 = (-a1 + root)*inv;
            return 2;
        }
        else
        {
            // Line is tangent to the cylinder.
            p0 = -a1/a2;
            p1 = NAN;
            return 1;
        }
    }
    
    // Test plane intersections first.
    int quantity = 0;
    inv = 1./Dz;
    
    double t0 = (-halfHeight - Pz)*inv;
    double xTmp = Px + t0*Dx;
    double yTmp = Py + t0*Dy;
    if (xTmp*xTmp + yTmp*yTmp <= rSqr)
    {
        // Planar intersection inside the top cylinder end disk.
        p0 = t0;
        quantity++;
    }
    
    double t1 = (+halfHeight - Pz)*inv;
    xTmp = Px + t1*Dx;
    yTmp = Py + t1*Dy;
    if (xTmp*xTmp + yTmp*yTmp <= rSqr)
    {
        // Planar intersection inside the bottom cylinder end disk.
        if (quantity==0) {
            p0 = t1;
        } else {
            p1 = t1;
        }
        quantity++;
    }
    
    if (quantity == 2)
    {
        // Line intersects both top and bottom cylinder end disks.
        if (p0 > p1)
        {
            double save = p0;
            p0 = p1;
            p1 = save;
        }
        return 2;
    }

    // If quantity == 1, then the line must intersect cylinder wall in a
    // single point somewhere between the end disks.  This case is detected
    // in the following code that tests for intersection between line and
    // cylinder wall.
    a0 = Px*Px + Py*Py - rSqr;
    a1 = Px*Dx + Py*Dy;
    a2 = Dx*Dx + Dy*Dy;
    discr = a1*a1 - a0*a2;
    if (discr < 0.)
    {
        // Line does not intersect cylinder wall.
        if (quantity>0) {
            log_fatal("Internal logic error: found wall intersection when there should have been none.");
        }
        return 0;
    }
    else if (discr > EPSILON)
    {
        root = std::sqrt(discr);
        inv = 1./a2;
        tValue = (-a1 - root)*inv;
        if (t0 <= t1)
        {
            if ((t0 <= tValue) && (tValue <= t1))
            {
                if (quantity==0) {
                    p0 = tValue;
                } else {
                    p1 = tValue;
                }
                quantity++;
            }
        }
        else
        {
            if ((t1 <= tValue) && (tValue <= t0))
            {
                if (quantity==0) {
                    p0 = tValue;
                } else {
                    p1 = tValue;
                }
                quantity++;
            }
        }
        
        if (quantity == 2)
        {
            // Line intersects one of the cylinder end disks and once on the
            // cylinder wall.
            if (p0 > p1)
            {
                double save = p0;
                p0 = p1;
                p1 = save;
            }
            return 2;
        }
        
        tValue = (-a1 + root)*inv;
        if (t0 <= t1)
        {
            if (t0 <= tValue && tValue <= t1)
            {
                if (quantity==0) {
                    p0 = tValue;
                } else {
                    p1 = tValue;
                }
                quantity++;
            }
        }
        else
        {
            if (t1 <= tValue && tValue <= t0)
            {
                if (quantity==0) {
                    p0 = tValue;
                } else {
                    p1 = tValue;
                }
                quantity++;
            }
        }
    }
    else
    {
        tValue = -a1/a2;
        if (t0 <= t1)
        {
            if ((t0 <= tValue) && (tValue <= t1))
            {
                if (quantity==0) {
                    p0 = tValue;
                } else {
                    p1 = tValue;
                }
                quantity++;
            }
        }
        else
        {
            if ((t1 <= tValue) && (tValue <= t0))
            {
                if (quantity==0) {
                    p0 = tValue;
                } else {
                    p1 = tValue;
                }
                quantity++;
            }
        }
    }
    
    if (quantity == 2)
    {
        if (p0 > p1)
        {
            double save = p0;
            p0 = p1;
            p1 = save;
        }
    }
    
    return quantity;
    
}


bool
I3ExtraGeometryItemCylinder::DoesLineIntersect
(const I3Position &lineStart,
 const I3Position &lineEnd) const
{
    if (isnan(radius_)) return false;
    
    double result[2];
    double line_segment_length;
    const unsigned int numResults = FindIntersections(lineStart,
                                                      lineEnd,
                                                      result[0],
                                                      result[1],
                                                      line_segment_length);


    for (unsigned int i=0; i<numResults; ++i)
    {
        if ((result[i] >= 0.) && (result[i] <= line_segment_length)) return true;
    }

    // no intersection found
    return false;
}

std::pair<I3Position, I3Position>
I3ExtraGeometryItemCylinder::GetBoundingBox
() const
{
    CalculateBoundingBox();
    return std::make_pair(boundingBoxLower_, boundingBoxUpper_);
}



template <class Archive>
void I3ExtraGeometryItemCylinder::load(Archive &ar, unsigned version)
{
    if (version > i3extrageometryitemcylinder_version_)
        log_fatal("Attempting to read version %u from file but running version %u of I3ExtraGeometryItem class.",version,i3extrageometryitemcylinder_version_);
    
    ar >> make_nvp("I3ExtraGeometryItem", base_object<I3ExtraGeometryItem>(*this));

    ar >> make_nvp("from", from_);
    ar >> make_nvp("to", to_);
    ar >> make_nvp("radius", radius_);
}


template <class Archive>
void I3ExtraGeometryItemCylinder::save(Archive &ar, unsigned version) const
{
    ar << make_nvp("I3ExtraGeometryItem", base_object<I3ExtraGeometryItem>(*this));

    ar << make_nvp("from", from_);
    ar << make_nvp("to", to_);
    ar << make_nvp("radius", radius_);
}


std::ostream& I3ExtraGeometryItemCylinder::operator<<(std::ostream& oss) const
{
    oss << "[ I3ExtraGeometryItemCylinder :" << std::endl;
    oss << "                         from : (" << from_.GetX()/I3Units::m << "," << from_.GetY()/I3Units::m << "," << from_.GetZ()/I3Units::m << ")m" << std::endl;
    oss << "                           to : (" << to_.GetX()/I3Units::m << "," << to_.GetY()/I3Units::m << "," << to_.GetZ()/I3Units::m << ")m" << std::endl;
    oss << "                       radius : " << radius_/I3Units::m << "m" << std::endl;
    oss << "]" << std::endl;
    return oss;
}

// the following line needs to go away once a more recent version of icetray makes it into a release:
#include "split_serializable_backport.h"

I3_SPLIT_SERIALIZABLE(I3ExtraGeometryItemCylinder);
