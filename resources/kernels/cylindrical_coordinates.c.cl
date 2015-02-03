
typedef floating4_t coordinate_t;

inline floating_t magnitude(floating4_t vec)
{
    return my_sqrt(vec.x*vec.x + vec.y*vec.y + vec.z*vec.z);
}

// Coordinate transform for infinite muon tracks. Here the source depth is
// degenerate with time, so we have 1 fewer dimension
inline coordinate_t
getCoordinates(const floating4_t absPos, const struct I3CLSimReferenceParticle *source)
{
    coordinate_t coords;
    
    // NB: the reference vectors sourcePos, sourceDir, and perpDir are
    //     defined as static variables at compile time
    floating4_t pos = absPos - source->posAndTime;
    floating_t l = dot(pos, source->dir);
    floating4_t rho = pos - l*source->dir;
    
    // perpendicular distance
    coords.s0 = magnitude(rho);
    // azimuth (in radians, because that's how the first implementation worked)
    coords.s1 = (coords.s0 > 0) ?
        acos(-dot(rho,source->perpDir)/coords.s0) : 0;
    // depth of closest approach
    // FIXME: doing this properly requires sampling different vertex depths in
    //        the same run, which in turn requires passing the source vertex
    //        and direction to the kernel rather than compiling it in
    coords.s2 = source->posAndTime.z + l*source->dir.z;
    // delay time
    coords.s3 = pos.w - (l + coords.s0*tan_thetaC)*recip_speedOfLight;
    
    return coords;
}
