export
    getBody,
    getMu,
    loadBodies!,
    loadPlanet!,
    loadMoon!,
    loadSun!,
    loadPoints!,
    LCoords,
    lagrangePoint,
    L45Stable,
    hillRadius,
    sun,
    obliquity,
    moon,
    moonConsts,
    planCoords_hci,
    meanMoonCoords_ecl,
    nuations_Earth,
    meanObliqHigh_Earth,
    meanObliqLow_Earth,
    trueObliq_Earth,
    planetEphemerides,
    planetRV


function loadBodies!(sys)
    nonEarthBodies = ["SUN", "MERCURY", "VENUS", "MARS", "JUPITER", "SATURN", "NEPTUNE", "URANUS"]
    for key in keys(sys)
        if in(key, sys["bodies"]) && !haskey(sys[key], "coords")
            if key == "MOON"
                loadMoon!(sys)
            elseif key == "SUN"
                loadSun!(sys)
            elseif key == "EARTH"
                sys[key]["coords"] = zeros(sys["tSteps"], 3)
                sys[key]["frame"] = "eci"
            else
                loadPlanet!(sys, key)
            end
        end
    end
end


function loadPlanet!(sys, planet::String)
    if !haskey(sys[planet], "coords")
        X⃗₀ = planetRV(planet, sys["t"][1])
        c = zeros(sys["tSteps"], 3)
        for i = 1:sys["tSteps"]
            Δt = (sys["t"][i] - sys["t"][1]) * solSecsEarth # [sec]
            c[i,:] = kepler(X⃗₀, Δt, μ=μ☉)[1]
            # c[i,:] = planCoords_hci(sys["t"][i], planet)
        end
        sys[planet]["coords"] = c
        sys[planet]["frame"] = "hci"
    end
end


function loadMoon!(sys)
    if !haskey(sys["MOON"], "coords")
        m = zeros(sys["tSteps"], 3)
        for i = 1:sys["tSteps"]
            r      = moon(sys["t"][i]) # [rad rad km]
            m[i,:] = αδr2ijk(r)        # [km km km]
            # m[i,:] = meanMoonCoords_ecl(sys["t"][i])
        end
        sys["MOON"]["coords"] = m
        sys["MOON"]["frame"]  = "ecl"
    end
end


function loadSun!(sys)
    if !haskey(sys["SUN"], "coords")
        # sys["SUN"]["coords"] = zeros(sys["tSteps"], 3)
        # sys["SUN"]["frame"] = "hci"
        m = zeros(sys["tSteps"], 3)
        for i = 1:sys["tSteps"]
            r      = sun(sys["t"][i])    # [rad rad km]
            m[i,:] = αδr2ijk(r)          # [km km km]
        end
        sys["SUN"]["coords"] = m
        sys["SUN"]["frame"]  = "ecl"
    end
end


function loadPoints!(sys)
    pointDict = groupPoints(sys["points"])
    loadPointBodies!(sys, keys(pointDict))
    for key in keys(pointDict)
        loadSystemPoints!(sys, key, sort(vec(pointDict[key])))
    end
end


function groupPoints(points)
    pointDict = Dict()
    for i = 1:length(points)
        bA = points[i][1]
        bB = points[i][2]
        n  = parse(Int, points[i][4])
        key = sort!([getBody(bA), getBody(bB)])
        pointDict[key] = haskey(pointDict, key) ? [pointDict[key] n] : [n]
    end
    return pointDict
end


function loadPointBodies!(sys, keys)
    for key in keys
        for body in key
            if !haskey(sys, body)
                sys[body] = Dict{String,Any}("type" => "body", "plot" => false)
                push!(sys["bodies"], body)
            end
        end
    end
    loadBodies!(sys)
end


function getBody(ch::Char)
    (ch == 'E') && return "EARTH"
    (ch == 'M') && return "MOON"
    (ch == 'S') && return "SUN"
    error("ERROR: Could not define body ", ch, ".")
end


function loadSystemPoints!(sys, bodyArr, numArr)
    L     = zeros(sys["tSteps"], 3, length(numArr))
    bodyA = bodyArr[1]
    bodyB = bodyArr[2]
    μA    = getMu(bodyA)
    μB    = getMu(bodyB)

    if (sys[bodyA]["frame"] ≠ sys[bodyB]["frame"])
        sys[bodyB]["coords"] = transformItem(sys["t"], sys[bodyB]["coords"], sys[bodyB]["frame"], sys[bodyA]["frame"])
        sys[bodyB]["frame"] = sys[bodyA]["frame"]
    end

    for i = 1:sys["tSteps"]
        R = norm(sys[bodyA]["coords"][i,:] - sys[bodyB]["coords"][i,:])
        l = LCoords(μA, μB, R)
        for j = 1:length(numArr)
            L[i,:,j] = l[numArr[j]]
        end
    end

    for i = 1:length(numArr)
        key = string(bodyA[1], bodyB[1], 'L', numArr[i])
        sys[key]["coords"] = L[:,:,i]
        sys[key]["frame"] = join(sort(collect(string(bodyA[1], bodyB[1])))) * "syn"
    end
end


getMu(body::String) = eval(Meta.parse("mu" * uppercasefirst(lowercase(body))))


# Return locations of Lagrange Points in 2D synodic frame. Source: Wikipedia
# function LCoords(μₐ, μB, R; order=1)
#     μ₁     = max(μₐ, μB) # Larger body
#     μ₂     = min(μₐ, μB) # Smaller body

#     o_R₁   = R*(μ₂ / (μ₁ + μ₂))
#     o_R₂   = R*(μ₁ / (μ₁ + μ₂))
#     ω      = √((μ₁+μ₂)/R^3)

#     init   = R*∛(μ₂/(3.0*μ₁)) # initial guess for L1/L2 dist.
#     fₗ₁(x)  =  (μ₂/(x^2)) + (o_R₂ - x)*(ω^2) - μ₁/((R-x)^2)
#     fₗ₂(x)  = -(μ₂/(x^2)) + (o_R₂ + x)*(ω^2) - μ₁/((R+x)^2)
#     b_rₗ₁   = [fzero(fₗ₁, -init, order=order) 0.0 0.0]
#     b_rₗ₂   = [fzero(fₗ₂, +init, order=order) 0.0 0.0]

#     init   = R * 7.0*μ₂ / (12.0*μ₁) # initial guess for L3 dist.
#     fₗ₃(x)  = -(μ₁/((R-x)^2)) - (μ₂/((2.0*R-x)^2)) + (ω^2)*(o_R₁+R-x)
#     a_rₗ₃   = [-R + fzero(fₗ₃, init, order=order) 0.0 0.0]

#     a_rₗ₄   = R*[sin(π/6.0)  cos(π/6.0) 0.0]
#     a_rₗ₅   = R*[sin(π/6.0) -cos(π/6.0) 0.0]

#     r₁      = [-o_R₁ 0.0 0.0]
#     r₂      = [ o_R₂ 0.0 0.0]
#     rₗ₁    =  r₂ + b_rₗ₁
#     rₗ₂    =  r₂ + b_rₗ₂
#     rₗ₃    = -r₁ + a_rₗ₃
#     rₗ₄    = -r₁ + a_rₗ₄
#     rₗ₅    = -r₁ + a_rₗ₅

#     rₗ₁, rₗ₂, rₗ₃, rₗ₄, rₗ₅, r₁, r₂
# end

function LCoords(μₐ::T, μB::T, R::T; tol::T=1e-8) where T <: Float64
    μ = min(μₐ,μB)/(μₐ+μB)

    rₗ₁ = R*lagrangePoint(μ, 1, tol=tol)
    rₗ₂ = R*lagrangePoint(μ, 2, tol=tol)
    rₗ₃ = R*lagrangePoint(μ, 3, tol=tol)
    rₗ₄ = R*lagrangePoint(μ, 4         )
    rₗ₅ = R*lagrangePoint(μ, 5         )
    r₁ = R*[-μ                 0 0]
    r₂ = R*[max(μₐ,μB)/(μₐ+μB) 0 0]

    return rₗ₁, rₗ₂, rₗ₃, rₗ₄, rₗ₅, r₁, r₂
end

# Return locations of Lagrange Point N (1-5) in normalized 2D synodic frame.
# Source: http://descanso.jpl.nasa.gov/monograph/series12/LunarTraj--Overall.pdf p. 357-8
# After Szebehely (1967)
function lagrangePoint(μ::Real, N::Int64; tol::Float64=1e-12, max_i::Int64=20)
    N == 1 && (f = (x) -> ∛(   μ *(x-1)^2/(3-2μ-x*(3-μ-x))); g = (y) -> 1-μ-y)
    N == 2 && (f = (x) -> ∛(   μ *(x+1)^2/(3-2μ+x*(3-μ+x))); g = (y) -> 1-μ+y)
    N == 3 && (f = (x) -> ∛((1-μ)*(x+1)^2/(1+2μ+x*(2+μ+x))); g = (y) ->  -μ-y)
    N == 4 && return [-μ+0.5  √3*0.5 0.0]
    N == 5 && return [-μ+0.5 -√3*0.5 0.0]

    γ₀, γ, i = ∛(μ*(1-μ)/3), 0.0, 1
    while true
        γ  = f(γ₀)
        (abs(γ-γ₀) ≤ tol || i > max_i) && break
        γ₀, i = γ, i+1
    end
    return [g(γ) 0.0 0.0]
end

#Source: http://map.gsfc.nasa.gov/ContentMedia/lagrange.pdf
L45Stable(μₐ::Real, μB::Real) = max(μₐ, μB)/min(μₐ, μB) ≥ 12.5(1 + √0.9936)

# Returns radius of Hill sphere for a smaller body orbiting a larger body with semi-major axis a and eccentricity e.
hillRadius(μₐ::Real, μB::Real, a::Real; e::Float64=0.0) = a*(1-e)*∛(min(μₐ, μB)/(3max(μₐ, μB)))


# Returns geocentric ecliptic positon vector of Sun
# Source: Vallado (2001) p. 266-7. Sun algorithm.
function sun(JDUT1::Float64)
    TUT1 = JD2JDC(JDUT1)

    λM☉  = 280.4606184 + 36000.77005361TUT1
    λM☉  = rem(λM☉, 360.0)                                            # [deg]

    TTDB = TUT1
    M☉   = 357.5277233 + 35999.05034TTDB                              # [deg]
    M☉   = mod2pi(deg2rad(M☉))                                        # [rad]

    λecl = λM☉ + 1.914666471sin(M☉) + 0.019994643sin(2M☉)             # [deg]
    λecl = mod2pi(deg2rad(λecl))                                      # [rad]

    r☉   = 1.000140612 - 0.016708617cos(M☉) - 0.000139589cos(2M☉)     # [AU]

    return [λecl 0.0 r☉*km_in_AU] # [rad rad km]
end

# Returns ε̄. T_TDB is the Julian centuries rom J2000.
# Source: Vallado (2001) p.209
function obliquity(T_TDB::Float64)
    ε̄ = 23.439291 - 1.30042e-2T_TDB - 1.64e-7T_TDB^2 + 5.04e-7T_TDB^3
    return deg2rad(ε̄)
end


# ------------------------------------------------------------------------------
#
#                           function moon
#
#  Calculates moon's geocentric position vector at TDB Julian date.
#
#  author        : David Vallado
#
#  revisions
#    Spicer      - Julia implementation         30 mar 2016
#
#  inputs          description                    range / [units]
#    JDTDB       - TDB julian date                 [days]
#
#  outputs       :
#    r⃗☾          - moon's ijk position vector       [AU]
#
#  references    :
#    vallado       2001, p.274, Algorithm 31
#
# r⃗☾ = moon(JDTDB)
# ------------------------------------------------------------------------------

function moon(JDTDB::Float64)
    TTDB = JD2JDC(JDTDB)

    λecl = 218.32 + 481267.8813TTDB +
                      6.29 *sind(134.9 + 477198.85TTDB) -
                      1.27 *sind(259.2 - 413335.38TTDB) +
                      0.66 *sind(235.7 + 890534.23TTDB) +
                      0.21 *sind(269.9 + 954397.70TTDB) -
                      0.19 *sind(357.5 +  35999.05TTDB) -
                      0.11 *sind(186.6 + 966404.05TTDB)       # [deg]

    ϕecl =            5.13 *sind( 93.3 + 483202.03TTDB) +
                      0.28 *sind(228.2 + 960400.87TTDB) -
                      0.28 *sind(318.3 +   6003.18TTDB) -
                      0.17 *sind(217.6 - 407332.20TTDB)       # [deg]

    P    =   0.9508 + 0.0518cosd(134.9 + 477198.85TTDB) +
                      0.0095cosd(259.2 - 413335.38TTDB) +
                      0.0078cosd(235.7 + 890534.23TTDB) +
                      0.0028cosd(269.9 + 954397.70TTDB)       # [deg]

    λecl = mod2pi(deg2rad(λecl))                           # [rad]
    ϕecl = mod2pi(deg2rad(ϕecl))                           # [rad]
    P    = mod2pi(deg2rad(P   ))                           # [rad]

    r☾   = R⨁*csc(P)

    return [λecl ϕecl r☾] # [rad rad km]
end




# ------------------------------------------------------------------------------
#
#                           function planetEphemerides
#
#  Returns planet's orbital elements at Julian century TTDB. Originally from
#  Meeus (1991, 202-204). Values for Pluto from Astronomical Almanac 1995.
#
#  author        : David Vallado
#
#  revisions
#    Spicer      - Julia implementation         10 Apr 2016
#
#  inputs          description                    range / [units]
#    planet      - string, e.g. "jupiter"
#    TTDB        - TTDB Julian century             [centuries]
#
#  outputs       :
#    a           - Planet's semimajor axis         [AU]
#    e           - Eccentricity                    [-]
#    ι           - Inclination                     [rad]
#    Ω           - Longitude of ascending node     [rad]
#    ω̃           - Longitude of periapsis          [rad]
#    λM          - Mean Longitude                  [rad]
#
#  references    :
#    vallado       2001, p.911-913, Appendix D.4
#
# a, e, ι, Ω, ω̃, λM = planetEphemerides(planet, TTDB)
# ------------------------------------------------------------------------------
function planetEphemerides(planet::String, TTDB::Float64)
    if planet == "mercury"
        a  =   0.387098310
        e  =   0.20563175  +      0.000020406*TTDB - 0.0000000284 *TTDB^2 - 0.00000000017*TTDB^3
        ι  =   7.004986    -      0.0059516  *TTDB + 0.00000081   *TTDB^2 + 0.000000041  *TTDB^3
        Ω  =  48.330893    -      0.1254229  *TTDB - 0.00008833   *TTDB^2 - 0.000000196  *TTDB^3
        ω̃  =  77.456119    +      0.1588643  *TTDB - 0.00001343   *TTDB^2 + 0.000000039  *TTDB^3
        λM = 252.250906    + 149472.6746358  *TTDB - 0.00000535   *TTDB^2 + 0.000000002  *TTDB^3
    elseif planet == "venus"
        a  =   0.723329820
        e  =   0.00677188  -     0.000047766 *TTDB + 0.0000000975 *TTDB^2 + 0.00000000044*TTDB^3
        ι  =   3.394662    -     0.0008568   *TTDB - 0.00003244   *TTDB^2 + 0.000000010  *TTDB^3
        Ω  =  76.679920    -     0.2780080   *TTDB - 0.00014256   *TTDB^2 - 0.000000198  *TTDB^3
        ω̃  = 131.563707    +     0.0048646   *TTDB - 0.00138232   *TTDB^2 - 0.000005332  *TTDB^3
        λM = 181.979801    + 58517.8156760   *TTDB + 0.00000165   *TTDB^2 - 0.000000002  *TTDB^3
    elseif planet == "earth" 
        a  =   1.000001018 
        e  =   0.01670862  -     0.000042037 *TTDB - 0.0000001236 *TTDB^2 + 0.00000000004*TTDB^3
        ι  =   0.0000000   +     0.0130546   *TTDB - 0.00000931   *TTDB^2 - 0.000000034  *TTDB^3
        Ω  = 174.873174    -     0.2410908   *TTDB + 0.00004067   *TTDB^2 - 0.000001327  *TTDB^3
        ω̃  = 102.937348    +     0.3225557   *TTDB + 0.00015026   *TTDB^2 + 0.000000478  *TTDB^3
        λM = 100.466449    + 35999.3728519   *TTDB - 0.00000568   *TTDB^2 + 0.000000000  *TTDB^3
    elseif planet == "mars" 
        a  =   1.523679342 
        e  =   0.09340062  +     0.000090483 *TTDB - 0.0000000806 *TTDB^2 - 0.00000000035*TTDB^3
        ι  =   1.849726    -     0.0081479   *TTDB - 0.00002255   *TTDB^2 - 0.000000027  *TTDB^3
        Ω  =  49.558093    -     0.2949846   *TTDB - 0.00063993   *TTDB^2 - 0.000002143  *TTDB^3
        ω̃  = 336.060234    +     0.4438898   *TTDB - 0.00017321   *TTDB^2 + 0.000000300  *TTDB^3
        λM = 355.433275    + 19140.2993313   *TTDB + 0.00000261   *TTDB^2 - 0.000000003  *TTDB^3
    elseif planet == "jupiter" 
        a  =   5.202603191 +    0.0000001913 *TTDB
        e  =   0.04849485  +    0.000163244  *TTDB - 0.0000004719 *TTDB^2 - 0.00000000197*TTDB^3
        ι  =   1.303270    -    0.0019872    *TTDB + 0.00003318   *TTDB^2 + 0.000000092  *TTDB^3
        Ω  = 100.464441    +    0.1766828    *TTDB + 0.00090387   *TTDB^2 - 0.000007032  *TTDB^3
        ω̃  =  14.331309    +    0.2155525    *TTDB + 0.00072252   *TTDB^2 - 0.000004590  *TTDB^3
        λM =  34.351484    + 3034.9056746    *TTDB - 0.00008501   *TTDB^2 + 0.000000004  *TTDB^3
    elseif planet == "saturn" 
        a  =   9.554909596 -    0.0000021389 *TTDB
        e  =   0.05550862  -    0.000346818  *TTDB - 0.0000006456 *TTDB^2 + 0.00000000338*TTDB^3
        ι  =   2.488878    +    0.0025515    *TTDB - 0.00004903   *TTDB^2 + 0.000000018  *TTDB^3
        Ω  = 113.665524    -    0.2566649    *TTDB - 0.00018345   *TTDB^2 + 0.000000357  *TTDB^3
        ω̃  =  93.056787    +    0.5665496    *TTDB + 0.00052809   *TTDB^2 + 0.000004882  *TTDB^3
        λM =  50.077471    + 1222.1137943    *TTDB + 0.00021004   *TTDB^2 - 0.000000019  *TTDB^3
    elseif planet == "uranus" 
        a  =  19.218446062 -    0.0000000372 *TTDB + 0.00000000098*TTDB^2
        e  =   0.04629590  -    0.000027337  *TTDB + 0.0000000790 *TTDB^2 + 0.00000000025*TTDB^3
        ι  =   0.773196    -    0.0016869    *TTDB + 0.00000349   *TTDB^2 + 0.000000016  *TTDB^3
        Ω  =  74.005947    +    0.0741461    *TTDB + 0.00040540   *TTDB^2 + 0.000000104  *TTDB^3
        ω̃  = 173.005159    +    0.0893206    *TTDB - 0.00009470   *TTDB^2 + 0.000000413  *TTDB^3
        λM = 314.055005    +  428.4669983    *TTDB - 0.00000486   *TTDB^2 + 0.000000006  *TTDB^3
    elseif planet == "neptune"
        a  =  30.110386869 -    0.0000001663 *TTDB + 0.00000000069*TTDB^2
        e  =   0.00898809  +    0.000006408  *TTDB - 0.0000000008 *TTDB^2
        ι  =   1.769952    +    0.0002257    *TTDB + 0.00000023   *TTDB^2 - 0.000000000  *TTDB^3
        Ω  = 131.784057    -    0.0061651    *TTDB - 0.00000219   *TTDB^2 - 0.000000078  *TTDB^3
        ω̃  =  48.123691    +    0.0291587    *TTDB + 0.00007051   *TTDB^2 - 0.000000000  *TTDB^3
        λM = 304.348665    +  218.4862002    *TTDB + 0.00000059   *TTDB^2 - 0.000000002  *TTDB^3
    elseif planet == "pluto"
        a  =  39.53758
        e  =   0.250877
        ι  =  17.13233
        Ω  = 110.4065
        ω̃  = 224.6148
        λM = 218.88735
    else
        error("Invalid entry for planet.")
    end

    return a, e, mod2pi(deg2rad(ι)), mod2pi(deg2rad(Ω)), mod2pi(deg2rad(ω̃)), mod2pi(deg2rad(λM))
end

# Finds X⃗ [km, km/s (hci, XYZ)] of a planet at specific JD.
# Source: Vallado (2001), Algorithm 33, p.283
function planetRV(planet::String, JD::Float64)
    TTDB = JD2JDC(JD)

    a, e, ι, Ω, ω̃, λM = planetEphemerides(planet, TTDB)

    M = λM - ω̃
    ω =  ω̃ - Ω

    E = kepEqtnE(M, e)
    ν = A2ν(E, e)

    p = a*(1-e^2)*km_in_AU # [km]

    X⃗_XYZ = RANDV(p, e, ι, Ω, ω, ν, μ=μ☉)
    X⃗_XYZ /= day_in_TU☉
    return X⃗_XYZ # [km, km/sec] (hci)
end


















# Moon tables from Meeus 1998
moonTableA = [
    (0,  0,  1,  0, 6288774, -20905355),
    (2,  0, -1,  0, 1274027,  -3699111),
    (2,  0,  0,  0,  658314,  -2955968),
    (0,  0,  2,  0,  213618,   -569925),
    (0,  1,  0,  0, -185116,     48888),
    (0,  0,  0,  2, -114332,     -3149),
    (2,  0, -2,  0,   58793,    246158),
    (2, -1, -1,  0,   57066,   -152138),
    (2,  0,  1,  0,   53322,   -170733),
    (2, -1,  0,  0,   45758,   -204586),
    (0,  1, -1,  0,  -40923,   -129620),
    (1,  0,  0,  0,  -34720,    108743),
    (0,  1,  1,  0,  -30383,    104755),
    (2,  0,  0, -2,   15327,     10321),
    (0,  0,  1,  2,  -12528,         0),
    (0,  0,  1, -2,   10980,     79661),
    (4,  0, -1,  0,   10675,    -34782),
    (0,  0,  3,  0,   10034,    -23210),
    (4,  0, -2,  0,    8548,    -21636),
    (2,  1, -1,  0,   -7888,     24208),
    (2,  1,  0,  0,   -6766,     30824),
    (1,  0, -1,  0,   -5163,     -8379),
    (1,  1,  0,  0,    4987,    -16675),
    (2, -1,  1,  0,    4036,    -12831),
    (2,  0,  2,  0,    3994,    -10445),
    (4,  0,  0,  0,    3861,    -11650),
    (2,  0, -3,  0,    3665,     14403),
    (0,  1, -2,  0,   -2689,     -7003),
    (2,  0, -1,  2,   -2602,         0),
    (2, -1, -2,  0,    2390,     10056),
    (1,  0,  1,  0,   -2348,      6322),
    (2, -2,  0,  0,    2236,     -9884),
    (0,  1,  2,  0,   -2120,      5751),
    (0,  2,  0,  0,   -2069,         0),
    (2, -2, -1,  0,    2048,     -4950),
    (2,  0,  1, -2,   -1773,      4130),
    (2,  0,  0,  2,   -1595,         0),
    (4, -1, -1,  0,    1215,     -3958),
    (0,  0,  2,  2,   -1110,         0),
    (3,  0, -1,  0,    -892,      3258),
    (2,  1,  1,  0,    -810,      2616),
    (4, -1, -2,  0,     759,     -1897),
    (0,  2, -1,  0,    -713,     -2117),
    (2,  2, -1,  0,    -700,      2354),
    (2,  1, -2,  0,     691,         0),
    (2, -1,  0, -2,     596,         0),
    (4,  0,  1,  0,     549,     -1423),
    (0,  0,  4,  0,     537,     -1117),
    (4, -1,  0,  0,     520,     -1571),
    (1,  0, -2,  0,    -487,     -1739),
    (2,  1,  0, -2,    -399,         0),
    (0,  0,  2, -2,    -381,     -4421),
    (1,  1,  1,  0,     351,         0),
    (3,  0, -2,  0,    -340,         0),
    (4,  0, -3,  0,     330,         0),
    (2, -1,  2,  0,     327,         0),
    (0,  2,  1,  0,    -323,      1165),
    (1,  1, -1,  0,     299,         0),
    (2,  0,  3,  0,     294,         0),
    (2,  0, -1, -2,       0,      8752)]


moonTableB = [
    (0,  0,  0,  1, 5128122),
    (0,  0,  1,  1,  280602),
    (0,  0,  1, -1,  277693),
    (2,  0,  0, -1,  173237),
    (2,  0, -1,  1,   55413),
    (2,  0, -1, -1,   46271),
    (2,  0,  0,  1,   32573),
    (0,  0,  2,  1,   17198),
    (2,  0,  1, -1,    9266),
    (0,  0,  2, -1,    8822),
    (2, -1,  0, -1,    8216),
    (2,  0, -2, -1,    4324),
    (2,  0,  1,  1,    4200),
    (2,  1,  0, -1,   -3359),
    (2, -1, -1,  1,    2463),
    (2, -1,  0,  1,    2211),
    (2, -1, -1, -1,    2065),
    (0,  1, -1, -1,   -1870),
    (4,  0, -1, -1,    1828),
    (0,  1,  0,  1,   -1794),
    (0,  0,  0,  3,   -1749),
    (0,  1, -1,  1,   -1565),
    (1,  0,  0,  1,   -1491),
    (0,  1,  1,  1,   -1475),
    (0,  1,  1, -1,   -1410),
    (0,  1,  0, -1,   -1344),
    (1,  0,  0, -1,   -1335),
    (0,  0,  3,  1,    1107),
    (4,  0,  0, -1,    1021),
    (4,  0, -1,  1,     833),
    (0,  0,  1, -3,     777),
    (4,  0, -2,  1,     671),
    (2,  0,  0, -3,     607),
    (2,  0,  2, -1,     596),
    (2, -1,  1, -1,     491),
    (2,  0, -2,  1,    -451),
    (0,  0,  3, -1,     439),
    (2,  0,  2,  1,     422),
    (2,  0, -3, -1,     421),
    (2,  1, -1,  1,    -366),
    (2,  1,  0,  1,    -351),
    (4,  0,  0,  1,     331),
    (2, -1,  1,  1,     315),
    (2, -2,  0, -1,     302),
    (0,  0,  1,  3,    -283),
    (2,  1,  1, -1,    -229),
    (1,  1,  0, -1,     223),
    (1,  1,  0,  1,     223),
    (0,  1, -2, -1,    -220),
    (2,  1, -1, -1,    -220),
    (1,  0,  1,  1,    -185),
    (2, -1, -2, -1,     181),
    (0,  1,  2,  1,    -177),
    (4,  0, -2, -1,     176),
    (4, -1, -1, -1,     166),
    (1,  0,  1, -1,    -164),
    (4,  0,  1, -1,     132),
    (1,  0, -1, -1,    -119),
    (4, -1,  0, -1,     115),
    (2, -2,  0,  1,     107)]


# Calculates constants required by other Moon functions
function moonConsts(T)
    _kL1    = Poly([deg2rad(218.3164477), deg2rad(481267.88123421), deg2rad(-0.0015786), deg2rad( 1.0/00538841), deg2rad(-1.0/065194000)])
    _kD     = Poly([deg2rad(297.8501921), deg2rad(445267.1114034),  deg2rad(-0.0018819), deg2rad( 1.0/00545868), deg2rad(-1.0/113065000)])
    _kM     = Poly([deg2rad(357.5291092), deg2rad(035999.0502909),  deg2rad(-0.0001536), deg2rad( 1.0/24490000)])
    _kM1    = Poly([deg2rad(134.9633964), deg2rad(477198.8675055),  deg2rad( 0.0087414), deg2rad( 1.0/00069699), deg2rad(-1.0/014712000)])
    _kF     = Poly([deg2rad(093.2720950), deg2rad(483202.0175233),  deg2rad(-0.0036539), deg2rad(-1.0/03526000), deg2rad( 1.0/863310000)])

    _kA1    = Poly([deg2rad(119.75), deg2rad(000131.849)])
    _kA2    = Poly([deg2rad(053.09), deg2rad(479264.290)])
    _kA3    = Poly([deg2rad(313.45), deg2rad(481266.484)])

    L1      = mod2pi(polyval(_kL1, T))
    D       = mod2pi(polyval(_kD,  T))
    M       = mod2pi(polyval(_kM,  T))
    M1      = mod2pi(polyval(_kM1, T))
    F       = mod2pi(polyval(_kF,  T))

    A1      = mod2pi(polyval(_kA1, T))
    A2      = mod2pi(polyval(_kA2, T))
    A3      = mod2pi(polyval(_kA3, T))

    E       = polyval(Poly([1.0, -0.002516, -0.0000074]), T)

    return (L1, D, M, M1, F, A1, A2, A3, E, E^2)
end


function meanMoonCoords_ecl(JD)
    T = JD2JDC(JD)
    L1, D, M, M1, F, A1, A2, A3, E, E2 = moonConsts(T)

    # Longitude and Radius
    suml, sumr = 0.0, 0.0
    for i in moonTableA
        (tD, tM, tM1, tF, tl, tr) = i
        arg = (tD*D) + (tM*M) + (tM1*M1) + (tF*F)
        if abs(tM) == 1
            tl *= E
            tr *= E
        elseif abs(tM) == 2
            tl *= E2
            tr *= E2
        end
        suml += tl*sin(arg)
        sumr += tr*cos(arg)
    end

    # Latitude
    sumb = 0.0
    for i in moonTableB
        (tD, tM, tM1, tF, tb) = i
        arg = (tD*D) + (tM*M) + (tM1*M1) + (tF*F)
        if abs(tM) == 1
            tb *= E
        elseif abs(tM) == 2
            tb *= E2
        end
        sumb += tb*sin(arg)
    end

    suml = suml + 3958 * sin(A1     ) +
                  1962 * sin(L1 - F ) +
                  0318 * sin(A2     )

    sumb = sumb - 2235 * sin(L1     ) +
                  0382 * sin(A3     ) +
                  0175 * sin(A1 - F ) +
                  0175 * sin(A1 + F ) +
                  0127 * sin(L1 - M1) -
                  0115 * sin(L1 + M1)

    λ = L1 + deg2rad(suml/1e6) # [rad]
    β = deg2rad(sumb/1e6)      # [rad]
    ∆ = (385000560 + sumr)/1e3     # [km]

    return ∆*[cos(β)*cos(λ) cos(β)*sin(λ) sin(β)]
end


# Converts moon coordinates into ecl xyz coords.
#=
function moonCoords_ecl(JD)
    lambda, beta, delta = meanMoonCoords(JD)
    O_r_M_ecl = delta*[cos(beta)*cos(lambda) cos(beta)*sin(lambda) sin(beta)]
end
=#

#=
function moonCoords_eci(JD)
    lambda, beta, Delta = meanMoonCoords(JD);
    lambda = lambda + nuations(JD)[1),
    eps = trueObliq(JD);

    alpha = atan(sin(lambda)*cos(eps) - tan(beta)*sin(eps), cos(lambda))
    delta = asin(sin(beta)*cos(eps) + cos(beta)*sin(eps)*sin(lambda))
    O_r_M_eci = Delta*[cos(delta)*cos(alpha) cos(delta)*sin(alpha) sin(delta)]
end
=#


#=
Return the nutation in right ascension (also called equation of the equinoxes.)
    Meeus-1998: page 88.
    Parameters:
        jd : Julian Day in dynamical time
    Returns:
        nutation, in radians

function nut_in_ra(jd)
    global days_per_second
    deltapsi = rad2deg(nut_in_lon(jd)) * 3600     # deltapsi in seconds
    epsilon  = true_obliquity(jd)                 # Epsilon kept in radians
    c = deltapsi * cos(epsilon)/15                # result in seconds...
    return (c * (pi * 2) * days_per_second)       # converted radians
end
=#


#=
Geocentric solar position and radius, low precision.
=#

#=
Return geocentric ecliptic longitude, latitude and radius.

        Parameters:
            JD : Julian Day in dynamical time
        Returns:
            longitude in radians
            latitude in radians
            radius in au

=#

function planCoords_hci(JD, planet)
    λ, β, R = vsop87d_dimension(JD, lowercase(planet))
    R *= km_in_AU
    return R*[cos(β)*cos(λ) cos(β)*sin(λ) sin(β)]
end

#=
_kL0 = Poly([deg2rad(280.46646),  deg2rad(36000.76983),  deg2rad( 0.0003032)])
_kM  = Poly([deg2rad(357.52911),  deg2rad(35999.05029),  deg2rad(-0.0001537)])
_kC  = Poly([deg2rad(  1.914602), deg2rad(   -0.004817), deg2rad(-0.0000140)])
_ke  = Poly([0.016708634, -0.000042037, -0.0000001267])

_ck3 = deg2rad( 0.019993)
_ck4 = deg2rad(-0.000101)
_ck5 = deg2rad( 0.000289)

# Return geometric longitude and radius vector.
# Low precision. The longitude is accurate to 0.01 degree.
# The latitude should be presumed to be 0.0. [Meeus-1998: equations 25.2 through
# 25.5, pages 165-166
#
#    Parameters:
#        JD : Julian Day in dynamical time
#    Returns:
#        longitude in radians
#        radius in au

function longitude_radius_low(JD)
    T  = JD2JDC(JD) # julian centuries from J2000.0
    L0 = polyval(_kL0, T) # geometric mean longitude of sun for mean equinox
    M  = polyval(_kM, T) # mean anomaly of sun
    e  = polyval(_ke, T) # eccentricity of Earth's orbit
    C  = (polyval(_kC, T)*sin(M))
         + ((_ck3 - _ck4*T)*sin(2*M))
         + (_ck5*sin(3*M)) # equation of the center
    L  = mod2pi(L0 + C) # true longitude
    v  = M + C # true anomaly
    R  = 1.000001018*(1 - e*e) / (1 + e*cos(v)) # radius vector
    return (L, R)
end

_lk0 = deg2rad(125.04)
_lk1 = deg2rad(1934.136)
_lk2 = deg2rad(0.00569)
_lk3 = deg2rad(0.00478)

# Correct the geometric longitude for nutation and aberration.
#    Low precision. [Meeus-1998: pg 164]
#    Parameters:
#        JD : Julian Day in dynamical time
#        L : longitude in radians
#    Returns:
#        corrected longitude in radians
#
function apparent_longitude_low(JD, L)
    T = JD2JDC(JD)
    omega = _lk0 - _lk1*T
    return mod2pi(L - _lk2 - _lk3*sin(omega))
end

_lk4 = dms2rad(0, 0, 20.4898)

# Correct for aberration; low precision, but good enough for most uses.
#
#    [Meeus-1998: pg 164]
#
#    Parameters:
#        R : radius in au
#    Returns:
#        correction in radians
function aberration_low(R)
    return -_lk4 / R
end


    # Returns the equation of time at JD
    # parameters:
    #   JD: julian day in dynamic time
    # Returns:
    #   Equation of time in minutes, can be positive or negative
function equation_time(JD)
    tau = JD2JDC(JD)/10
    _p = Poly([280.4664567, 360007.6982779, 0.03032028, 1/49931, -1/15300, -1/2000000])
    L0 = mod2pi(deg2rad(polyval(_p, tau)))
    (L, B, R)  = meanSunCoords(JD)
    (L, B)     = vsop2fk5(JD, L, B)
    deltaPsi   = nut_in_lon(JD)
    epsilon    = meanObliqHigh(JD)
    asc,decl   = ecl_to_equ(L, B, epsilon)
    eqt = L0 - deg2rad(0.0057183) - asc + deltaPsi*cos(epsilon)
    if eqt > pi/2
        eqt = -((2*pi) - eqt)
    end
    return eqt/(2*pi)*1440
end
=#








nutationTable_Earth = [
    ( 0,  0,  0,  0,  1, -171996, -1742, 92025,  89),
    (-2,  0,  0,  2,  2,  -13187,   -16,  5736, -31),
    ( 0,  0,  0,  2,  2,   -2274,    -2,   977,  -5),
    ( 0,  0,  0,  0,  2,    2062,     2,  -895,   5),
    ( 0,  1,  0,  0,  0,    1426,   -34,    54,  -1),
    ( 0,  0,  1,  0,  0,     712,     1,    -7,   0),
    (-2,  1,  0,  2,  2,    -517,    12,   224,  -6),
    ( 0,  0,  0,  2,  1,    -386,    -4,   200,   0),
    ( 0,  0,  1,  2,  2,    -301,     0,   129,  -1),
    (-2, -1,  0,  2,  2,     217,    -5,   -95,   3),
    (-2,  0,  1,  0,  0,    -158,     0,     0,   0),
    (-2,  0,  0,  2,  1,     129,     1,   -70,   0),
    ( 0,  0, -1,  2,  2,     123,     0,   -53,   0),
    ( 2,  0,  0,  0,  0,      63,     0,     0,   0),
    ( 0,  0,  1,  0,  1,      63,     1,   -33,   0),
    ( 2,  0, -1,  2,  2,     -59,     0,    26,   0),
    ( 0,  0, -1,  0,  1,     -58,    -1,    32,   0),
    ( 0,  0,  1,  2,  1,     -51,     0,    27,   0),
    (-2,  0,  2,  0,  0,      48,     0,     0,   0),
    ( 0,  0, -2,  2,  1,      46,     0,   -24,   0),
    ( 2,  0,  0,  2,  2,     -38,     0,    16,   0),
    ( 0,  0,  2,  2,  2,     -31,     0,    13,   0),
    ( 0,  0,  2,  0,  0,      29,     0,     0,   0),
    (-2,  0,  1,  2,  2,      29,     0,   -12,   0),
    ( 0,  0,  0,  2,  0,      26,     0,     0,   0),
    (-2,  0,  0,  2,  0,     -22,     0,     0,   0),
    ( 0,  0, -1,  2,  1,      21,     0,   -10,   0),
    ( 0,  2,  0,  0,  0,      17,    -1,     0,   0),
    ( 2,  0, -1,  0,  1,      16,     0,    -8,   0),
    (-2,  2,  0,  2,  2,     -16,     1,     7,   0),
    ( 0,  1,  0,  0,  1,     -15,     0,     9,   0),
    (-2,  0,  1,  0,  1,     -13,     0,     7,   0),
    ( 0, -1,  0,  0,  1,     -12,     0,     6,   0),
    ( 0,  0,  2, -2,  0,      11,     0,     0,   0),
    ( 2,  0, -1,  2,  1,     -10,     0,     5,   0),
    ( 2,  0,  1,  2,  2,      -8,     0,     3,   0),
    ( 0,  1,  0,  2,  2,       7,     0,    -3,   0),
    (-2,  1,  1,  0,  0,      -7,     0,     0,   0),
    ( 0, -1,  0,  2,  2,      -7,     0,     3,   0),
    ( 2,  0,  0,  2,  1,      -7,     0,     3,   0),
    ( 2,  0,  1,  0,  0,       6,     0,     0,   0),
    (-2,  0,  2,  2,  2,       6,     0,    -3,   0),
    (-2,  0,  1,  2,  1,       6,     0,    -3,   0),
    ( 2,  0, -2,  0,  1,      -6,     0,     3,   0),
    ( 2,  0,  0,  0,  1,      -6,     0,     3,   0),
    ( 0, -1,  1,  0,  0,       5,     0,     0,   0),
    (-2, -1,  0,  2,  1,      -5,     0,     3,   0),
    (-2,  0,  0,  0,  1,      -5,     0,     3,   0),
    ( 0,  0,  2,  2,  1,      -5,     0,     3,   0),
    (-2,  0,  2,  0,  1,       4,     0,     0,   0),
    (-2,  1,  0,  2,  1,       4,     0,     0,   0),
    ( 0,  0,  1, -2,  0,       4,     0,     0,   0),
    (-1,  0,  1,  0,  0,      -4,     0,     0,   0),
    (-2,  1,  0,  0,  0,      -4,     0,     0,   0),
    ( 1,  0,  0,  0,  0,      -4,     0,     0,   0),
    ( 0,  0,  1,  2,  0,       3,     0,     0,   0),
    ( 0,  0, -2,  2,  2,      -3,     0,     0,   0),
    (-1, -1,  1,  0,  0,      -3,     0,     0,   0),
    ( 0,  1,  1,  0,  0,      -3,     0,     0,   0),
    ( 0, -1,  1,  2,  2,      -3,     0,     0,   0),
    ( 2, -1, -1,  2,  2,      -3,     0,     0,   0),
    ( 0,  0,  3,  2,  2,      -3,     0,     0,   0),
    ( 2, -1,  0,  2,  2,      -3,     0,     0,   0)]

function nutationConsts_Earth(T)
    _kD  = Poly([deg2rad(297.85036), deg2rad(445267.111480), deg2rad(-0.0019142), deg2rad( 1.0/189474)])
    _kM  = Poly([deg2rad(357.52772), deg2rad(035999.050340), deg2rad(-0.0001603), deg2rad(-1.0/300000)])
    _kM1 = Poly([deg2rad(134.96298), deg2rad(477198.867398), deg2rad( 0.0086972), deg2rad( 1.0/ 56250)])
    _kF  = Poly([deg2rad(093.27191), deg2rad(483202.017538), deg2rad(-0.0036825), deg2rad( 1.0/327270)])
    _ko  = Poly([deg2rad(125.04452), deg2rad(-01934.136261), deg2rad( 0.0020708), deg2rad( 1.0/450000)])

    D  = mod2pi(polyval(_kD,  T))
    M  = mod2pi(polyval(_kM,  T))
    M1 = mod2pi(polyval(_kM1, T))
    F  = mod2pi(polyval(_kF,  T))
    ω  = mod2pi(polyval(_ko,  T))

    D, M, M1, F, ω
end


function nuations_Earth(JD)
    # TODO nut_in_lon() factor the /1e5 and /1e6 adjustments into the table.
    T = JD2JDC(JD)
    D, M, M1, F, omega = nutationConsts_Earth(T)
    dPsi, dEps = 0, 0
    for v in nutationTable_Earth
        (tD, tM, tM1, tF, tOmega, tPsiK, tPsiT, tEpsK, tEpsT) = v
        arg  = D*tD + M*tM + M1*tM1 + F*tF + omega*tOmega
        dPsi = dPsi + sin(arg)*(tPsiK/10000 + tPsiT*T/100000)
        dEps = dEps + cos(arg)*(tEpsK/10000 + tEpsT*T/100000)
    end
    deg2rad(dPsi/3600), deg2rad(dEps/3600) # [rad]
end


function meanObliqLow_Earth(JD)
    T = JD2JDC(JD)
    _el0 = Poly([dms2rad(23, 26,  21.448),
                 dms2rad( 0,  0, -46.8150),
                 dms2rad( 0,  0,  -0.00059),
                 dms2rad( 0,  0,   0.001813)])
    eps0 = polyval(_el0, T) # [rad]
end


function meanObliqHigh_Earth(JD)
    U = JD2JDC(JD)/100
    _el1 = Poly([   dms2rad(23, 26,    21.448),
                    dms2rad( 0,  0, -4680.93),
                    dms2rad( 0,  0,    -1.55),
                    dms2rad( 0,  0,  1999.25),
                    dms2rad( 0,  0,   -51.38),
                    dms2rad( 0,  0,  -249.67),
                    dms2rad( 0,  0,   -39.05),
                    dms2rad( 0,  0,     7.12),
                    dms2rad( 0,  0,    27.87),
                    dms2rad( 0,  0,     5.79),
                    dms2rad( 0,  0,     2.45)])
    eps0 = polyval(_el1, U) # [rad]
end

function trueObliq_Earth(JD::Float64; precision::String="l")
    if precision == "l"
        return deg2rad(23.43474888)
    elseif precision == "m"
        return meanObliqLow_Earth(JD) + nuations_Earth(JD)[2]
    else
        return meanObliqHigh_Earth(JD) + nuations_Earth(JD)[2] # [rad]
    end
end