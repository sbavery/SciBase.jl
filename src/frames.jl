export
    getPath!,
    calcPath,
    transformFrames!,
    transformItem,
    ecef2gdc,
    ecef2aer,
    gdc2gcc,
    ecl2eci,
    hci2ecl,
    eci2ecef,
    ecl2EMsyn,
    hci2ESsyn,
    ecef2sez,
    ecef2sezmat,
    eci2sez,
    coe2peri,
    peri2bci,
    rot1,
    rot2,
    rot3,
    rot1mat,
    rot2mat,
    rot3mat,
    aer2sez,
    aer2sez2,
    sez2aer,
    αδr2ijk,
    eci2topo,
    eci2geo,
    λϕr2αδr,
    RANDV,
    rsw2ijkmat,
    ntw2ijkmat,
    sph2cart,
    cart2sph

global inputs = [("eci",   "ecl",  false),
                 ("ecl",   "hci",  false),
                 ("ecef",  "eci",  false),
                 ("EMsyn", "ecl",  false),
                 ("ESsyn", "hci",  false),
                 ("gdc",   "ecef", false), # Geodetic (ground track)
                 ("gcc",   "gdc",  false), # Geocentric
                 ("aer",   "ecef", false), # Azimuth, elevation, & range
                 ("mcl",   "hci",  false), # Mars-centered ecliptic
                 ("mci",   "mcl",  false)] # Mars-centered inertial
for i in eachindex(inputs)
    global inputs = [inputs; (inputs[i][2], inputs[i][1], true)]
end


function transformFrames!(sys, toFrame::String)
    if sys["tUnits"] ≠ "JD"
        sys["t"] = datetime2julian(sys["t"])
        sys["tUnits"] = "JD"
    end
    for key in keys(sys)
        if isa(sys[key], Dict) && haskey(sys[key], "frame")
            (toFrame == "aer" && key in sys["stations"]) && continue
            (in(toFrame, ["aer", "gdc", "gcc"]) && key==getCentralBody(sys)) && continue
            (key == "EARTH" && in(sys[key]["frame"], ["eci", "ecl", "ecef"]) && in(toFrame, ["eci", "ecl", "ecef"])) && continue

            if sys[key]["frame"] ≠ toFrame
                sys[key]["coords"] = transformItem(sys["t"], sys[key]["coords"], sys[key]["frame"], toFrame)
                sys[key]["frame"] = toFrame
            end
        end
    end
end


function transformItem(T, R::Array{Float64}, f1::String, f2::String)
    (f1 == f2) && return R
    global inputs
    path, inv = String[], Bool[]
    getPath!(f1, f2, path, inv, inputs)
    str1 = join(reverse([path[i+Int(inv[i])] * "2" * path[i+Int(!inv[i])] * "(" for i in eachindex(inv)]))
    str2 = ", JD=JD, inv=" * join(inv, "), JD=JD, inv=") * ")"
    @eval f = (r, JD) -> $(Meta.parse(eval(str1 * "r" * str2)))
    c = zeros(length(T), 3)
    for i in eachindex(T)
        c[i,:] = Base.invokelatest(f, R[i:i,:], T[i])
    end
    return c
end


function getPath!(from::String, to::String, path::Array{String,1}, inv::Array{Bool,1}, inputs::Array{Tuple{String,String,Bool},1})
    if from == to
        push!(path, from)
        return true
    end
    for i in eachindex(inputs)
        if (inputs[i][1] == from && !in(inputs[i][2], path))
            push!(path, inputs[i][1]), push!(inv, !inputs[i][3])
            if !getPath!(inputs[i][2], to, path, inv, inputs)
                pop!(path), pop!(inv)
            else
                return true
            end
        end
    end
    return false
end


function ecl2eci(r::Array{Float64}; JD::Float64=2.45e6, inv::Bool=false)
    r == zeros(1,3) && return r
    ε̄  = obliquity(JD2JDC(JD))
    return inv ? rot1(r, ε̄) : rot1(r, -ε̄)
end


function mcl2mci(r::Array{Float64}; JD::Float64=2.45e6, inv::Bool=false)
    r == zeros(1,3) && return r
    ε̄  = meanObliq_Mars
    return inv ? rot1(r, ε̄) : rot1(r, -ε̄)
end


function ecl2EMsyn(r::Array{Float64}; JD::Float64=2.45e6, inv::Bool=false)
    r⃗☾ = αδr2ijk(moon(JD))#meanMoonCoords_ecl(JD)
    r☾ = norm(r⃗☾)
    θ  = atan(r⃗☾[2], r⃗☾[1])
    ϕ  = 0.5π - acos(r⃗☾[3]/r☾)

    if !inv
        r = (rot2mat(-ϕ) * rot3mat(θ) * r[:])'
        r[1] -= r☾ * μ☾ / (μ⨁ + μ☾)
        return r'
    else
        r[1] += r☾ * μ☾ / (μ⨁ + μ☾)
        return (rot3mat(-θ) * rot2mat(ϕ) * r[:])'
    end
end


function eci2ecef(r::Array{Float64}; JD::Float64=2.45e6, inv::Bool=false)
    r == zeros(1,3) && return r
    θG = LSTime(JD)[1] # [rad]
    return inv ? rot3(r, -θG) : rot3(r, θG)
end


hci2ecl(r::Array{Float64}; JD::Float64=2.45e6, inv::Bool=false) = (!inv) ? r - planCoords_hci(JD, "earth") : r + planCoords_hci(JD, "earth")

hci2mcl(r::Array{Float64}; JD::Float64=2.45e6, inv::Bool=false) = (!inv) ? r - planCoords_hci(JD,  "mars") : r + planCoords_hci(JD,  "mars")


function hci2ESsyn(r::Array{Float64}; JD::Float64=2.45e6, inv::Bool=false)
    r⃗⨁ = -αδr2ijk(sun(JD))#planCoords_hci(JD, "EARTH")
    r⨁ = norm(r⃗⨁)
    θ  = atan(r⃗⨁[2], r⃗⨁[1])
    ϕ  = 0.5π - acos(r⃗⨁[3]/r⨁)

    if !inv
        r = (rot2mat(-ϕ) * rot3mat(θ) * r[:])'
        r[1] -= r⨁ * μ⨁ / (μ⨁ + μ☉)
        return r'
    else
        r[1] += r⨁ * μ⨁ / (μ⨁ + μ☉)
        return (rot3mat(-θ) * rot2mat(ϕ) * r[:])'
    end
end

# Convert ecef geocentric equatorial position vector into [lon (rad), lat (rad), h (km)]. Source: Vallado ijk2ll. h = height above ellipsoid [km]
function ecef2gdc(r::Array{Float64}; JD::Float64=NaN, inv::Bool=false, max_i::Int64=20, tol::Float64=1e-6)
    if !inv
        i, Nᵩ = 1, 0.0
        dist = norm(r[1:2])

        λ = (abs(dist) < tol) ? 0.5π*sign(r[3]) : atan(r[2], r[1])
        λ = putInRange(λ, -π, π)

        ϕgd  = asin(r[3]/norm(r))
        ϕgd₀ = ϕgd + 10
        while abs(ϕgd₀ - ϕgd) ≥ tol && i < max_i
            ϕgd₀ = ϕgd
            Nᵩ   = R⨁/√(1 - e1⨁2*sin(ϕgd)^2)
            ϕgd  = atan((r[3] + Nᵩ*e1⨁2*sin(ϕgd))/dist)
            i += 1
        end

        h = (0.5π - abs(ϕgd) > deg2rad(1)) ? dist*sec(ϕgd) - Nᵩ : r[3]*csc(ϕgd) - Nᵩ*(1 - e1⨁2)
        return [λ ϕgd h] # [rad, rad, km]
    else # from Vallado Site algorithm
        sϕ  = sin(r[2])
        C⨁  = R⨁/√(1-e1⨁2*sϕ^2)
        S⨁  = C⨁*(1-e1⨁2)

        r_δ = (C⨁ + r[3])*cos(r[2])
        r_K = (S⨁ + r[3])*sϕ
        return [r_δ*cos(r[1]) r_δ*sin(r[1]) r_K] # [km, km, km]
    end
end


# Converts between geodetic and geocentric [lat, lon, h] for positions on the Earth's surface. Source: Vallado. Using degrees.
function gdc2gcc(r::Array{Float64}; JD::Float64=NaN, inv::Bool=false)
    r == zeros(1,3) && return r
    r[2] = inv ? atan(tan(r[2])/(1 - e1⨁2)) : atan((1 - e1⨁2)*tan(r[2]))
    return r
end

function ecef2aer(r::Array{Float64}; JD::Float64=NaN, inv::Bool=false) # From Vallado
    if !inv
        global azElGSECEF # [km km km]
        global azElGSGDC # [rad rad km]

        # tol = 1e-8
        # ------- find ecef range vector from site to satellite -------
        # @show r 
        # @show azElGSECEF
        ρEcef = r - azElGSECEF
        # drhoEcef = vEcef

        # ------------- convert to sez for calculations ---------------
        # tempvec = rot3(rhoEcef, lon)
        # rhoSez  = rot2(tempvec, pi/2 - ϕgdgd)

        # tempvec = rot3(drhoEcef, lon)
        # drhoSez = rot2(tempvec, pi/2 - ϕgdgd)

        ρSez = ecef2sez(ρEcef, azElGSGDC[1], azElGSGDC[2])

        # ------------- calculate azimuth and elevation ---------------
        # temp = norm(ρSez[1:2])
        # if temp < tol                   # directly over the north pole
            # el = sign(rhoSez[3])*pi/2   # +- 90 deg
            # az = atan(drhoSez[2], -drhosez[1])
        # else
            el = asin(ρSez[3]/norm(ρSez))
            # el = putInRange(el, 0, π/2)
            az = atan(ρSez[2], -ρSez[1])
        # end
        return [mod2pi(az) el norm(ρEcef)] # [rad, rad, km]
    else
        # get spacecraft ecef coords from az el r measurement
        # azl2radc.m
    end
end


"""
    ecef2aer(X::Array{T}, X_site::Array{T}, λ_site::T, ϕgd_site::T) where T <: Float64

Transform object's position vector `X` from body-centered, body fixed (ECEF) frame
to Azimuth/Elevation/Range frame as seen from position `X_site`.

### Arguments (Inputs):
| Parameter | Description                         | Units        |
| --------: | :------------------------------------------------------------- | :---- |
| `r⃗`       | Object's ECI/IJK state vector `[x y z]` | DU (e.g. km, km/s) |
| `r⃗_site`  | Site's ECI/IJK state vector `[x y z]` | DU |
| `λ_site` | Site's longitude | rad |
| `ϕgd_site` | Site's geodetic latitude | rad |


### Value (Output)
| Parameter | Description                         | Units        |
| --------: | :------------------------------ | :---- |
| `X`  | Az/El/R state vector `[β el ρ]` where | |
| `β`  | Azimuth | rad |
| `el` | Elevation | rad |
| `ρ`  | Range | DU |

Author: James Spicer (2017.12.08) <br/>
Source: Vallado 'Razel' algorithm without derivatives or θ_LST.  
"""
function ecef2aer(r⃗::Array{T}, r⃗_site::Array{T}, λ_site::T, ϕgd_site::T) where T <: Float64
    ρ⃗_ECEF = r⃗ - r⃗_site
    ρ⃗_SEZ  = ecef2sezmat(λ_site, ϕgd_site)*ρ⃗_ECEF[:]

    ρ      = norm(ρ⃗_SEZ)
    el     = asin(ρ⃗_SEZ[3]/ρ)
    t      = hypot(ρ⃗_SEZ[1], ρ⃗_SEZ[2])
    az     = atan(ρ⃗_SEZ[2]/t, -ρ⃗_SEZ[1]/t)
    return [mod2pi(az) el ρ] # [rad, rad, km]
end

# function aer2eci
#     else
#         # get spacecraft ecef coords from az el r measurement
#         # site-track (algorithm 48)
#         # azl2radc.m
#     end
# end


# Converts range, azimuth & elevation into topocentric horizon (sez) system.
# Source: Vallado raz2rvs.m & rvs2raz.m

function aer2sez(r::Array{Float64}; JD::Float64=NaN, inv::Bool=false) # Input aer in [rad]
    if !inv # az, el, range to SEZ:
        return r[3]*[-cos(r[2])*cos(r[1])
                      cos(r[2])*sin(r[1])
                      sin(r[2])          ]
    else # SEZ to az, el, range:
        # el = asin(r[3] / r[4])
        # az = atan(r[2], -r[1])
        # return
    end
end


"""
    aer2sez(X::Array{Float64}) 

Transform object's state vector `X` from azimuth, elevation, range (AER) frame to topocentric-horizon
system centered at a body-fixed site (SEZ). <br/> N.B. The topocentric-horizon frame is parallel to
the local body surface, NOT the same as the topocentric RaDec frame which is parallel to the equator.

### Arguments (Inputs):
| Parameter  | Description                              | Units       |
| ---------: | :--------------------------------------- | :---------- |
| `X`        | AER state vector `[β el ρ β̇ ėl ρ̇]` where |             |
| `β`, `β̇`   | Topocentric-horizon (rate of) azimuth    | rad, rad/TU |
| `el`, `ėl` | Topocentric-horizon (rate of) elevation  | rad, rad/TU |
| `ρ`, `ρ̇`   | Range and range rate                     | DU, DU/TU   |

### Value (Output)
| Parameter | Description                               | Units                     |
| --------: | :---------------------------------------- | :------------------------ |
| `Y`       | Object's SEZ state vector `[x y z ẋ ẏ ż]` | DU, DU/TU (e.g. km, km/s) |

Author: James Spicer (2018.11.07) <br/>
Source: Vallado p.250-251.  
"""
function aer2sez2(X::Array{Float64}) # Need to rename to not be confused with the other aer2sez
    cβ, sβ, cel, sel = cos(X[1]), sin(X[1]), cos(X[2]), sin(X[2])
    Y    = similar(X)
    Y[1] = -X[3]*cel*cβ
    Y[2] =  X[3]*cel*sβ
    Y[3] =  X[3]*sel
    length(X) == 3 && return Y
    Y[4] = -X[6]*cel*cβ + X[3]*sel*cβ*X[5] + X[3]*cel*sβ*X[4]
    Y[5] =  X[6]*cel*sβ - X[3]*sel*sβ*X[5] + X[3]*cel*cβ*X[4]
    Y[6] =  X[6]*sel    + X[3]*cel   *X[5]
    return Y
end


"""
    sez2aer(X::Array{T}; tol::T=1e-3) where T <: Float64

Transform object's state vector `X` from topocentric-horizon frame (SEZ) to azimuth,
elevation, range (AER). <br/> N.B. The topocentric-horizon frame is parallel to the
local body surface, NOT the same as the topocentric RaDec frame which is parallel to the equator.

### Arguments (Inputs):
| Parameter | Description                              | Units                     |
| --------: | :--------------------------------------- | :------------------------ |
| `X`       | Object's SEZ state vector `[x y z ẋ ẏ ż]` | DU, DU/TU (e.g. km, km/s) |

#### Optional Keyword Arguments:
| Parameter  | Description            | Units |
| ---------: | :--------------------- | :---- |
| `tol=1e-3` | Tolerance for el = 90° | rad   |

### Value (Output)
| Parameter  | Description                              | Units       |
| ---------: | :--------------------------------------- | :---------- |
| `Y`        | AER state vector `[β el ρ β̇ ėl ρ̇]` where |             |
| `β`, `β̇`   | Topocentric-horizon (rate of) azimuth    | rad, rad/TU |
| `el`, `ėl` | Topocentric-horizon (rate of) elevation  | rad, rad/TU |
| `ρ`, `ρ̇`   | Range and range rate                     | DU, DU/TU   |

Author: James Spicer (2018.11.07) <br/>
Source: Vallado RAZEL algorithm 27, p.252-255.  
"""
function sez2aer(X::Array{T}; tol::T=1e-3) where T <: Float64
    Y    = similar(X)
    Y[3] = norm(X[1:3])                # ρ
    Y[2] = asin(X[3]/Y[3])             # el
    if abs(Y[2] - 0.5π) > tol 
        Y[1] = atan(X[2], -X[1])       # β
    elseif length(X) == 6
        Y[1] = atan(X[5], -X[4])       # β
    else 
        Y[1] = NaN                     # β
    end
    Y[1] = mod2pi(Y[1])
    length(X) == 3 && return Y
    Y[6] = (X[1:3]⋅X[4:6])/Y[3]        # ρ̇
    t    =  hypot(X[1], X[2])
    Y[4] = (X[4]*X[2] - X[5]*X[1])/t^2 # β̇
    Y[5] = (X[6] - Y[6]*X[3]/Y[3])/t   # ėl
    return Y
end


"""
    eci2topo(X::Array{T}, X_site::Array{T}) where T <: Float64

Transform object's state vector `X` from body-centered intertial (ECI/IJK) frame to topocentric
RaDec frame centered at a body-fixed site. <br/> N.B. The topocentric frame is parallel to the 
equator, NOT the same as the SEZ/RAzEl/AER frame which is parallel to the local body surface. 

### Arguments (Inputs):
| Parameter | Description                         | Units        |
| --------: | :------------------------------------------------------------- | :---- |
| `X`       | Object's ECI/IJK state vector `[x y z ẋ ẏ ż]` | DU, DU/TU (e.g. km, km/s) |
| `X_site`  | Site's ECI/IJK state vector `[x y z ẋ ẏ ż]` | DU, DU/TU |


### Value (Output)
| Parameter | Description                         | Units        |
| --------: | :------------------------------------------------------------- | :---- |
| `Y`       | Topocentric state vector `[α_t δ_t ρ α̇_t δ̇_t ρ̇]` where | |
| `α_t`, `α̇_t` | Topocentric (rate of) right ascension | rad, rad/TU |
| `δ_t`, `δ̇_t` | Topocentric (rate of) declination | rad, rad/TU |
| `ρ`, `ρ̇` | Topocentric range and range rate | DU, DU/TU |

Author: James Spicer (2017.11.14) <br/>
Source: Vallado 'topocentric' algorithm.  
"""
eci2topo(X::Array{T}, X_site::Array{T}) where T <: Float64 = eci2geo(X - X_site)


"""
    eci2geo(X::Array{Float64}) 

Transform object's state vector `X` from body-centered intertial (ECI/IJK) frame to geocentric
RaDec frame.

### Arguments (Inputs):
| Parameter | Description                         | Units        |
| --------: | :------------------------------------------------------------- | :---- |
| `X`       | Object's ECI/IJK state vector `[x y z ẋ ẏ ż]` | DU, DU/TU (e.g. km, km/s) |

### Value (Output)
| Parameter | Description                         | Units        |
| --------: | :------------------------------------------------------------- | :---- |
| `Y`      | Geocentric state vector `[α δ r α̇ δ̇ ṙ]` where | |
| `α`, `α̇` | Geocentric (rate of) right ascension | rad, rad/TU |
| `δ`, `δ̇` | Geocentric (rate of) declination | rad, rad/TU |
| `r`, `ṙ` | Radius vector and derivative | DU, DU/TU |

Author: James Spicer (2017.11.15) <br/>
Source: Vallado 'Geocentric Radec' algorithm.  
"""
function eci2geo(X::Array{Float64})
    Y    = similar(X)
    Y[3] = norm(X[1:3])
    t    = hypot(X[1], X[2])
    Y[2] = atan(X[3]/Y[3], t/Y[3]) # only need sine expression (δ = asin(rK/r)) for declinations within ±90°.
    if t ≉ 0.0 
        Y[1] = atan(X[2]/t, X[1]/t)
    else
        u    = hypot(X[4], X[5])
        Y[1] = atan(X[5]/u, X[4]/u)
    end
    Y[6] = (X[1:3]⋅X[4:6])/Y[3]
    Y[4] = (X[1]*X[5]-X[2]*X[4])/t^2
    Y[5] = (X[6]-Y[6]*X[3]/Y[3])/t
    return Y
end

### Site-Track (gdc_site and aer => ijk)
# ϕgd, θLST, h_ellp, ρ, β, el, ρ̇, β̇, ėl => r⃗ijk, v⃗ijk, r⃗site_ijk
# 1. use ϕgd, h_ellp, and θLST to get r⃗site_ijk
# 2. convert ρ, β, el, ρ̇, β̇, ėl to ρ_SEZ and ρ̇_SEZ.
# 3. use ϕgd, θLST to get IJK/SEZ matrix
# 4. use IJK/SEZ matrix and r⃗site_ijk to get r⃗ijk, v⃗ijk.
# 
#
# Split into several methods: 
# gdc_site, time => X_siteIJK
# X_SEZ = aer2sez(X_AER)
# X_IJK = sez2ijk(X_SEZ, X_site, time)


function ecef2sez(r::Array{T}, λ::T, ϕgd::T; inv::Bool=false) where T <: Float64 # Source: Vallado.
    r == zeros(1,3) && return r
    return (!inv) ? (rot2mat(0.5π-ϕgd) * rot3mat(λ) * r[:])' : (rot3mat(-λ) * rot2mat(ϕgd-0.5π) * r[:])'
end

function eci2sez(r::Array{T}, θ_LST::T, ϕgd::T; inv::Bool=false) where T <: Float64 # Source: Vallado.
    r == zeros(1,3) && return r
    return (!inv) ? (rot2mat(0.5π-ϕgd) * rot3mat(θ_LST) * r[:])' : (rot3mat(-θ_LST) * rot2mat(ϕgd-0.5π) * r[:])'
end

# Source: Vallado. "Covariance Transformations for Satellite Flight Dynamics Operations", p.7
# Also Razel algorithm
function ecef2sezmat(λ::T, ϕgd::T) where T <: Float64
    M                      = zeros(Float64,3,3)
    sλ, cλ                 = sin(λ  ), cos(λ  ) 
    sϕ, cϕ                 = sin(ϕgd), cos(ϕgd)
    M[1,1], M[1,2], M[1,3] = sϕ*cλ, sϕ*sλ, -cϕ
    M[2,1], M[2,2]         =   -sλ,    cλ
    M[3,1], M[3,2], M[3,3] = cϕ*cλ, cϕ*sλ,  sϕ
    return M
end

# ------------------------------------------------------------------------------
# function coe2peri
#
# This function finds the position and velocity vectors in perifocal coordinate
# system (pqw) system given the classical orbit elements.
#
#  author        : James Spicer                   650-999-0331   22 mar 2016
#
#  revisions
#
#
#  inputs          description                    range / units
#    p           - semilatus rectum               [km]
#    e           - eccentricity
#    μ           - body's gravitational parameter 0.0  to  pi rad
#    ν           - true anomaly                   0.0  to 2pi rad
#
#  outputs         description                    range / units
#    r⃗           - pqw position vector            [km]
#    v⃗           - pqw velocity vector            [km / s]
#
#  references    :
#    vallado       2007, 126, alg 10, ex 2-5 (coe2rv.m)
# ------------------------------------------------------------------------------
function coe2peri(p::Real, e::T, ν::T; μ::T=μ⨁) where T <: Float64
    cν, sν = cos(ν), sin(ν)
    r⃗ = p / (1 + e*cν) * [cν;
                          sν;
                          0.0]
    v⃗ = √(μ/p) * [-sν    ;
                   cν + e;
                   0.0   ]
    return r⃗, v⃗
end

# ------------------------------------------------------------------------------
# function peri2bci
#
# This function finds the position and velocity vectors in body-centered inertial
# coordinate system (ijk) system from the perifocal (pqw) coordinate system.
#
#  author        : James Spicer                   650-999-0331   22 mar 2016
#
#  revisions
#
#
#  inputs          description                    range / units
#    r           - pqw (perifocal) vector
#    ι           - inclination                    0.0  to  pi rad
#    Ω           - longitude of ascending node    0.0  to 2pi rad
#    ω           - argument of perigee            0.0  to 2pi rad
#
#  outputs         description                    range / units
#    v2          - ijk (BCI) vector
#
#  references    :
#    vallado       2007, 126, alg 10, ex 2-5 (coe2rv.m)
# ------------------------------------------------------------------------------
function peri2bci(r::Array{T}, ω::T, ι::T, Ω::T) where T <: Float64
    temp = rot3(r   , -ω)
    temp = rot1(temp, -ι)
    return rot3(temp, -Ω)
end


# ------------------------------------------------------------------------------
# function peri2ecimat
#
# This function returns the matrix to transform to the body-centered inertial
# coordinate system (ijk) system from the perifocal (pqw) coordinate system.
#
#  author        : James Spicer                   650-999-0331   22 mar 2016
#
#  revisions
#
#
#  inputs          description                    range / units
#    ι           - inclination                    0.0  to  pi rad
#    Ω           - longitude of ascending node    0.0  to 2pi rad
#    ω           - argument of perigee            0.0  to 2pi rad
#
#  outputs         description                    range / units
#    v2          - ijk (BCI) vector
#
#  references    :
#    vallado       2007, 126, alg 10, ex 2-5 (coe2rv.m)
# ------------------------------------------------------------------------------
peri2bcimat(ω::T, ι::T, Ω::T) where T <: Float64 = rot3mat(-Ω) * rot1mat(-ι) * rot3mat(-ω)


"""
    rot1(v⃗::Array{T}, ϕ::T) where T <: Float64

Rotate vector about the 1st axis 

### Arguments (Inputs):
| Parameter | Description       | Units        |
| --------: | :---------------- | :----------- |
| `v⃗`       | Input vector      | DU (e.g. km) |
| `ϕ`       | Angle of rotation | rad          |

### Value (Output)
| Parameter | Description    | Units |
| --------: | :------------- | :---- |
| `w⃗`       | Rotated vector | DU    |

Author: James Spicer (2016.03.22) <br/>
Source: Vallado 
"""
function rot1(v⃗::Array{T}, ϕ::T) where T <: Float64
    sϕ, cϕ = sin(ϕ), cos(ϕ)
    w⃗      = similar(v⃗)
    w⃗[1]   =    v⃗[1]
    w⃗[2]   = cϕ*v⃗[2] + sϕ*v⃗[3]
    w⃗[3]   = cϕ*v⃗[3] - sϕ*v⃗[2]
    return w⃗
end


"""
    rot2(v⃗::Array{T}, ϕ::T) where T <: Float64

Rotate vector about the 2nd axis 

### Arguments (Inputs):
| Parameter | Description       | Units        |
| --------: | :---------------- | :----------- |
| `v⃗`       | Input vector      | DU (e.g. km) |
| `ϕ`       | Angle of rotation | rad          |

### Value (Output)
| Parameter | Description    | Units |
| --------: | :------------- | :---- |
| `w⃗`       | Rotated vector | DU    |

Author: James Spicer (2016.03.22) <br/>
Source: Vallado 
"""
function rot2(v⃗::Array{T}, ϕ::T) where T <: Float64
    sϕ, cϕ = sin(ϕ), cos(ϕ)
    w⃗      = similar(v⃗)
    w⃗[1]   = cϕ*v⃗[1] - sϕ*v⃗[3]
    w⃗[2]   =    v⃗[2]
    w⃗[3]   = cϕ*v⃗[3] + sϕ*v⃗[1]
    return w⃗
end


"""
    rot3(v⃗::Array{T}, ϕ::T) where T <: Float64

Rotate vector about the 3rd axis 

### Arguments (Inputs):
| Parameter | Description       | Units        |
| --------: | :---------------- | :----------- |
| `v⃗`       | Input vector      | DU (e.g. km) |
| `ϕ`       | Angle of rotation | rad          |

### Value (Output)
| Parameter | Description    | Units |
| --------: | :------------- | :---- |
| `w⃗`       | Rotated vector | DU    |

Author: James Spicer (2016.03.22) <br/>
Source: Vallado 
"""
function rot3(v⃗::Array{T}, ϕ::T) where T <: Float64
    sϕ, cϕ = sin(ϕ), cos(ϕ)
    w⃗      = similar(v⃗)
    w⃗[1]   = cϕ*v⃗[1] + sϕ*v⃗[2]
    w⃗[2]   = cϕ*v⃗[2] - sϕ*v⃗[1]
    w⃗[3]   =    v⃗[3]
    return w⃗
end


"""
    rot1mat(ϕ::Float64) 

Returns rotation matrix to rotate an angle about the 1st axis. 

### Arguments (Inputs):
| Parameter | Description       | Units |
| --------: | :---------------- | :---- |
| `ϕ`       | Angle of rotation | rad   |

### Value (Output)
| Parameter | Description     | Units |
| --------: | :-------------- | :---- |
| `M`       | Rotation matrix |       |

Author: James Spicer (2016.03.22) <br/>
Source: Vallado 
"""
function rot1mat(ϕ::Float64)
    M                      = zeros(Float64,3,3)
    sϕ, cϕ                 = sin(ϕ), cos(ϕ) 
    M[1,1], M[2,2], M[3,3] = 1.0, cϕ, cϕ
    M[2,3], M[3,2]         = sϕ, -sϕ
    return M
end


"""
    rot2mat(ϕ::Float64) 

Returns rotation matrix to rotate an angle about the 2nd axis. 

### Arguments (Inputs):
| Parameter | Description       | Units |
| --------: | :---------------- | :---- |
| `ϕ`       | Angle of rotation | rad   |

### Value (Output)
| Parameter | Description     | Units |
| --------: | :-------------- | :---- |
| `M`       | Rotation matrix |       |

Author: James Spicer (2016.03.22) <br/>
Source: Vallado 
"""
function rot2mat(ϕ::Float64)
    M                      = zeros(Float64,3,3)
    sϕ, cϕ                 = sin(ϕ), cos(ϕ) 
    M[1,1], M[2,2], M[3,3] =  cϕ, 1.0, cϕ
    M[1,3], M[3,1]         = -sϕ, sϕ
    return M
end


"""
    rot3mat(ϕ::Float64) 

Returns rotation matrix to rotate an angle about the 3rd axis. 

### Arguments (Inputs):
| Parameter | Description       | Units |
| --------: | :---------------- | :---- |
| `ϕ`       | Angle of rotation | rad   |

### Value (Output)
| Parameter | Description     | Units |
| --------: | :-------------- | :---- |
| `M`       | Rotation matrix |       |
Author: James Spicer (2016.03.22) <br/>
Source: Vallado 
"""
function rot3mat(ϕ::Float64)
    M                      = zeros(Float64,3,3)
    sϕ, cϕ                 = sin(ϕ), cos(ϕ) 
    M[1,1], M[2,2], M[3,3] = cϕ,  cϕ, 1.0
    M[1,2], M[2,1]         = sϕ, -sϕ
    return M
end


# Converts geocentric (not topocentric) α, δ, r coordinates to ijk.
# S = state vector [r₁ r₂ r₃ ṙ₁ ṙ₂ ṙ₃] e.g. [α δ r α̇ δ̇ ṙ]
function αδr2ijk(S::Array{Float64}; inv::Bool=false)
    # From Vallado 2001 p.247-248 'Geocentric Radec'
    if inv # find α, δ, r, α̇, δ̇, ṙ given r and v
        r, q = norm(S[1:3]), norm(S[1:2])
        δ    = atan(S[3]/r, q/r) # only need sine expression (δ = asin(S[3]/r)) for declinations within ±90°.
        α    = (q ≉ 0) ? atan(S[2], S[1]) : atan(S[5], S[4])
        (length(S) == 3) && return [mod2pi(α) δ r]

        ṙ    = (S[1:3]⋅S[4:6])/r
        α̇    = (S[5]*S[1]-S[4]*S[2])/(S[2]^2+S[1]^2)
        δ̇    = (S[6]-ṙ*S[3]/r)/q

        return [mod2pi(α) δ r α̇ δ̇ ṙ]

    else # find rijk, vijk given α, δ, r, α̇, δ̇, ṙ

        r⃗ = S[3]*[cos(S[2])*cos(S[1]);
                  cos(S[2])*sin(S[1]);
                  sin(S[2])]
        (length(S) == 3) && return r⃗'

        v⃗ = [S[6]*cos(S[2])*cos(S[1]) - S[3]*sin(S[2])*cos(S[1])*S[5] - S[3]*cos(S[2])*sin(S[1])*S[4];
             S[6]*cos(S[2])*sin(S[1]) - S[3]*sin(S[2])*sin(S[1])*S[5] + S[3]*cos(S[2])*cos(S[1])*S[4];
             S[6]*sin(S[2])                                           + S[3]*cos(S[2])          *S[5]]

        return [r⃗' v⃗']
    end
end

# Converts ecliptic longitude & latitude to equatorial lon & lat. r = [α δ r] or [ϕ λ r]
# Source: Vallado (2001), p.258-261
function λϕr2αδr(r::Array{Float64}; JD::Float64=2.45e6, inv::Bool=false)
    ε̄ = obliquity(JD2JDC(JD))
    if inv
        α = r[1]
        δ = r[2]
        ϕ = asin( -cos(δ)*sin(α)*sin(ε̄)+sin(δ)*cos(ε̄))
        λ = atan( cos(δ)*sin(α)*cos(ε̄)+sin(δ)*sin(ε̄), cos(δ)*cos(α))
        return [λ ϕ r[3]]
    else
        λ = r[1]
        ϕ = r[2]
        α = atan(-sin(ϕ)*sin(ε̄)+cos(ϕ)*cos(ε̄)*sin(λ), cos(ϕ)*cos(λ))
        δ = asin(  sin(ϕ)*cos(ε̄)+cos(ϕ)*sin(ε̄)*sin(λ))

        return [mod2pi(α) δ r[3]]
    end
end


"""
    RANDV(p::T, e::T, ι::T, Ω::T, ω::T, ν::T; μ::T=μ⨁) where T <: Float64

Converts classical orbital elements to body-centered inertial coordinates.

### Arguments (Inputs):
| Parameter | Description                                 | Units |
| --------: | :------------------------------------------ | :---- |
| `p`      | Semilatus rectum                             | DU    |
| `e`      | Eccentricity                                 | -     |
| `ι`      | Inclination                                  | rad   |
| `Ω`      | Right ascension of the ascending node (RAAN) | rad   |
| `ω`      | Argument of periapsis                        | rad   |
| `ν`      | True anomaly                                 | rad   |

#### Optional Keyword Arguments:
| Parameter | Description                          | Units   |
| --------: | :----------------------------------- | :------ |
| `μ=μ⨁`    | Central body gravitational parameter | DU³/TU² |

### Values (Outputs):
| Parameter | Description               | Units     |
| --------: | :------------------------ | :-------- |
| `X⃗`       | BCI state ([x y z ẋ ẏ ż]) | DU, DU/TU |

Author: James Spicer (2016.11.06) <br/>
Source: Vallado (2001), Algorithm 10, p.125<br/>
Source: GMAT Math Spec 2013a, Section 3.1.3, p.49-50.
"""
function RANDV(p::T, e::T, ι::T, Ω::T, ω::T, ν::T; μ::T=μ⨁) where T <: Float64 
    # if e < tol
    #     if ι < tol || abs(ι-π) < tol    # circular equatorial
    #         ω, Ω, ν = 0.0, 0.0, λtrue
    #     else                            # circular inclined
    #         ω, ν    = 0.0, u
    #     end
    # else                                # elliptical equatorial
    #     if ι < tol || abs(ι-π) < tol
    #         Ω, ω    = 0.0, ω̃true
    #     end
    # end

    # (abs(p) < 1e-4) && (p = 1e-4)

    # r, v = coe2peri(p, e, ν, μ=μ) # Perifocal coordinates
    # r    = peri2bci(r, ω, ι, Ω)   # Body-centered inertial coordinates
    # v    = peri2bci(v, ω, ι, Ω)   # Body-centered inertial coordinates

    # return r', v'
    
    e > 1 && abs(abs(ν) - π) ≤ asec(e) && error("ν, e input values not physically possible.")

    sν, cνe, sΩ, cΩ, sι, cι, sω, cω, sων, cων = sin(ν), cos(ν)+e, sin(Ω), cos(Ω), sin(ι), cos(ι), sin(ω), cos(ω), sin(ω+ν), cos(ω+ν)

    r = p/(1+e*cos(ν))
    v = √(μ/p)
    
    x = r*(cων*cΩ - cι*sων*sΩ)
    y = r*(cων*sΩ + cι*sων*cΩ)
    z = r*(            sων*sι)

    ẋ = v*(cνe*(-sω*cΩ-cι*sΩ*cω) - sν*(cω*cΩ-cι*sΩ*sω))
    ẏ = v*(cνe*(-sω*sΩ+cι*cΩ*cω) - sν*(cω*sΩ+cι*cΩ*sω))
    ż = v*(cνe*  sι*cω           - sν*sι*sω           )

    return [x y z ẋ ẏ ż]
end


"""
    RANDV(X::Array{T}; μ::T=μ⨁) where T <: Float64 

Converts classical orbital elements to body-centered inertial coordinates.

### Arguments (Inputs):
| Parameter | Description                                  | Units     |
| --------: | :------------------------------------------- | :-------- |
| `X`       | COE state where X = [p e ι Ω ω ν]            |           |
| `p`       | Semilatus rectum                             | DU        |
| `e`       | Eccentricity                                 | -         |
| `ι`       | Inclination                                  | rad       |
| `Ω`       | Right ascension of the ascending node (RAAN) | rad       |
| `ω`       | Argument of periapsis                        | rad       |
| `ν`       | True anomaly                                 | rad       |

#### Optional Keyword Arguments:
| Parameter | Description                                  | Units     |
| --------: | :------------------------------------------- | :-------- |
| `μ=μ⨁`    | Central body gravitational parameter         | DU³/TU²   |

### Values (Outputs):
| Parameter | Description                                  | Units     |
| --------: | :------------------------------------------- | :-------- |
| `Y`       | BCI state ([x y z ẋ ẏ ż])                    | DU, DU/TU |

Author: James Spicer (2016.11.06) <br/>
Source: GMAT Math Spec 2013a, Section 3.1.3, p.49-50.
"""
function RANDV(X::Array{T}; μ::T=μ⨁) where T <: Float64 
    X[2] > 1 && abs(abs(X[6]) - π) ≤ asec(X[2]) && error("ν, e input values not physically possible.")

    sν, cνe, sΩ, cΩ, sι, cι, sω, cω, sων, cων = sin(X[6]), cos(X[6])+X[2], sin(X[4]), cos(X[4]), sin(X[3]), cos(X[3]), sin(X[5]), cos(X[5]), sin(X[5]+X[6]), cos(X[5]+X[6])

    r = X[1]/(1+X[2]*cos(X[6]))
    v = √(μ/X[1])
    
    Y = similar(X)
    Y[1] = r*(cων*cΩ - cι*sων*sΩ)
    Y[2] = r*(cων*sΩ + cι*sων*cΩ)
    Y[3] = r*(            sων*sι)

    Y[4] = v*(cνe*(-sω*cΩ-cι*sΩ*cω) - sν*(cω*cΩ-cι*sΩ*sω))
    Y[5] = v*(cνe*(-sω*sΩ+cι*cΩ*cω) - sν*(cω*sΩ+cι*cΩ*sω))
    Y[6] = v*(cνe*  sι*cω           - sν*sι*sω           )

    return Y
end

# R = radial direction
# S = in direction of (but not parallel to) velocity and perpendicular to radial
# W = normal (to orbital plane) direction
function rsw2ijkmat(r⃗::AbstractVector, v⃗::AbstractVector; inv::Bool=false)
    R̂ = vec(normalize(r⃗))
    Ŵ = vec(normalize(r⃗ × v⃗))
    Ŝ = vec(Ŵ × R̂)
    return inv ? [R̂ Ŝ Ŵ]' : [R̂ Ŝ Ŵ]
end

# N = normal (to velocity) direction
# T = tangential (to orbit) direction
# W = normal (to orbital plane) direction
function ntw2ijkmat(r⃗::AbstractVector, v⃗::AbstractVector; inv::Bool=false)
    T̂ = vec(normalize(v⃗))
    Ŵ = vec(normalize(r⃗ × v⃗))
    N̂ = vec(T̂ × Ŵ)
    return inv ? [N̂ T̂ Ŵ]' : [N̂ T̂ Ŵ]
end


"""
    sph2cart(r⃗::Array{Float64}) 

Transform vector from spherical [az el r] into cartesian [x y z] coordinates.

### Arguments (Inputs):
| Parameter | Description                         | Units        |
| --------: | :------------------------------------------------------------- | :---- |
| `az`      | Azimuth, CCW angle in the x-y plane measured from the +x-axis. | rad |
| `el`      | Elevation angle from the x-y plane. | rad | 
| `r`       | Radius | DU (e.g. km) |    

Author: James Spicer (2017.10.30) <br/>  
"""
function sph2cart(r⃗::Array{Float64})
    cel  = cos(r⃗[2])
    v⃗    = similar(r⃗)
    v⃗[1] = cos(r⃗[1]) .* cel
    v⃗[2] = sin(r⃗[1]) .* cel
    v⃗[3] = sin(r⃗[2])
    return r⃗[3].*v⃗
end

function sph2cart(az::T, el::T, r::T) where T <: Float64
    v⃗ = sph2cart([az el r])
    return v⃗[1], v⃗[2], v⃗[3]
end


"""
    cart2sph(r⃗::Array{Float64}) 

Transform vector from cartesian [x y z] into spherical [az el r] coordinates.

## Values (Outputs):
| Parameter | Description                         | Units        |
| --------: | :------------------------------------------------------------- | :---- |
| `az`      | Azimuth, CCW angle in the x-y plane measured from the +x-axis. | rad |
| `el`      | Elevation angle from the x-y plane. | rad | 
| `r`       | Radius | DU (e.g. km) |    

Author: James Spicer (2017.10.30) <br/>  
"""
function cart2sph(r⃗::Array{Float64})
    v⃗    = similar(r⃗)
    v⃗[1] = atan(r⃗[2], r⃗[1])
    v⃗[2] = atan(r⃗[3], norm(r⃗[1:2]))
    v⃗[3] = norm(r⃗)
    return v⃗
end

function cart2sph(x::T, y::T, z::T) where T <: Float64
    v⃗ = cart2sph([x y z])
    return v⃗[1], v⃗[2], v⃗[3]
end