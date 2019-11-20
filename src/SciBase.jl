#files = [
#		]

#for f in files
#    include(f * ".jl")
#end

module constants
    using Unitful
    gEarth = 9.80665u"m/s^2"
    GEarth = (6.67408*10^(-11))u"N*m^2/kg^2"
    k = (8.9875517887*10^9)u"N*m^2/C^2"
    c = 299792458u"m/s"
    e = (1.60217662*10^(-19))u"C"
    h = (6.62607004*10^(-34))u"-34"
    eps0 = (8.85418782*10^(-12))u"m^(-3)*kg^(-1)*s^(4)*A^(2)"
    mu0 = (1.25663706*10^(-6))u"m*kg*s^-2*A^-2"
    mEarth = (5.972*10^24)u"kg"
    rEarth = 6378.1u"km"
    mSun = (1.989*10^30)u"kg"
    export
        gEarth,
        GEarth,
        k,
        c,
        e,
        h,
        eps0,
        mu0,
        mEarth,
        rEarth,
        mSun
end

module physics
    module kinematics
        using ....constants
        # Displacement
        disp(x0,x) = x-x0
        dispConstAcc1(v0,a,t) = v0*t+1/2*a*t^2
        dispConstAcc2(v0,v,t) = ((v+v0)/2)*t
        # Velocity
        avgVel(x0,x,t0,t) = (x-x0)/(t-t0)
        velConstAcc1(v0,a,t) = v0+a*t
        velConstAcc2(x0,x,v0,a) = sqrt(v0^2+2*a*(x-x0))
        # Acceleration
        accConst(x0,x,v0,v) = (v^2-v0^2)/(2*(x-x0))
        range(v0,th0) = v0^2/gEarth*sin(2*th0)
        maxHeight(v0,th0) = v0^2/(2*gEarth)*sin(th0)^2
    end
end