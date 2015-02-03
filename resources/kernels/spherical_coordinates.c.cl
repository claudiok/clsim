
typedef floating4_t coordinate_t;

inline floating_t magnitude(floating4_t vec)
{
    return my_sqrt(vec.x*vec.x + vec.y*vec.y + vec.z*vec.z);
}

inline coordinate_t
getCoordinates(const floating4_t absPos)
{
    coordinate_t coords;
    
    // NB: the reference vectors sourcePos, sourceDir, and perpDir are
    //     defined as static variables at compile time
    floating4_t pos = absPos - sourcePos;
    floating_t l = dot(pos, sourceDir);
    floating4_t rho = pos - l*sourceDir;
    floating_t n_rho = magnitude(rho);
    
    // radius
    coords.s0 = magnitude(pos);
    // azimuth
    coords.s1 = (n_rho > 0) ?
        acos(-dot(rho,perpDir)/n_rho)/(PI/180) : 0;
    // cos(polar angle)
    coords.s2 = (coords.s0 > 0) ? my_divide(l, coords.s0) : 0;
    // delay time
    coords.s3 = pos.w - coords.s0*min_invPhaseVel;
    
    dbg_printf("     %4.1f %4.1f %4.2f %6.2f\n", coords.s0, coords.s1, coords.s2, coords.s3);
    
    return coords;
}
