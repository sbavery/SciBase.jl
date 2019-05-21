export
    dms2rad,
    rad2dms,
    polarTranslate,
    polarDistance,
    learningCurve,
    antennaRDECost,
    antennaFUCost,
    antennaMass,
    num2si,
    rotMat,
    GSEqArc,
    horizonAngle,
    azElPol,
    putInRange,
    fix,
    boundInd,
    biLinearInterpSq,
    biLinearInterpTrap,
    linearInterp,
    biCubicInterp,
    sanitizeSysInputs!,
    randRange,
    quadratic,
    rgb,
    rangeInRange,
    nansum,
    nanmean,
    inRange,
    line2segments,
    fresnelS,
    fresnelC




fix(x) = trunc(x)


# Source: Vallado DMStoRad
function dms2rad(d, m, s)
    all([d m s] .> 0.0) && return deg2rad(d + m/60 + s/3600)
    return -deg2rad(abs(d) + abs(m)/60 + abs(s)/3600)
end

# Source: Vallado radToDMS
function rad2dms(α)
    t = rad2deg(α)
    d = trunc(t)
    m = trunc(60(t - d))
    s = 3600(t - d - m/60)
    return [d m s]
end


# Translate (r,t) by (r₀, t₀) in polar coordinates
function polarTranslate(r, t, r₀, t₀)
    r₂ = polarDistance(r, t, r₀, t₀)
    t₂ = atan(r*sin(t) - r₀*sin(t₀), r*cos(t) - r₀*cos(t₀))
    return r₂, t₂
end

polarDistance(r₁, t₁, r₂, t₂) = √(r₁^2 + r₂^2 - 2r₁*r₂*cos(t₁-t₂))


# Cost of R&DE antenna cost in FY92$K (SMAD '92 p. 726)
antennaRDECost(mass) = 1015(mass^0.59)

# Cost of antenna first unit cost in FY92$K (SMAD '92 p. 727)
function antennaFUCost(item, X)
    item == "antenna" && return 20 + 230(X^0.59)
    item == "bus" && return 16253 + 110(X^1)
end


# Parabolic antenna mass from diameter
antennaMass(D) = 7.08D - 6 # From Harris Corp.


# Learning curve: Cost reduction factor for producing N units (SMAD '92 Eq 20-4)
function learningCurve(N)
    if N < 10
        S = 0.95
    elseif N < 50
        S = 0.9
    else
        S = 0.85
    end

    B = 1 + log(S)/log(2)
    N^B
end


# Inputs (*==default):
#  val = NumericScalar, the value to be converted to string <str>.
#  sgf = NumericScalar, the significant figures in the coefficient, *3.
#  pfx = LogicalScalar, true/false* -> select SI prefix as name/symbol.
#  trz = LogicalScalar, true*/false -> select if decimal trailing zeros are required.
#
# Output:
#  str = Input <val> as a string: coefficient + space character + SI prefix.
#  str = num2si(val,*sgf,*pfx,*trz)
function num2si(val::Real; sgf::Int64=3, pfx::T=false, trz::T=true) where T <: Bool
    val = Float64(val)

    fmt = (trz && sgf > 1) ? "%#.$(sgf)g %s" : "%.$(sgf)g %s"

    !isfinite(val) && return @sprintf("%f ", val)

    # Calculate coefficient value:
    xpt = rem.(min.(9, max.(-9, [0; 1] + floor.(log10.(abs.(val))/3))), 9)
    cof = val.*1000 .^-xpt
    # Round coefficient value:
    ord = 1 + floor.(log10.(abs.(cof)))
    (val ≠ 0) && (cof = 10 .^(ord-sgf).*round.(cof.*10 .^(sgf-ord)))
    # Select prefix symbol/name:
    pfc = ["yocto" "zepto" "atto" "femto" "pico" "nano" "micro" "milli" "" "kilo" "mega" "giga" "tera" "peta" "exa" "zetta" "yotta";
           "y"     "z"     "a"    "f"     "p"    "n"    "μ"     "m"     "" "k"    "M"    "G"    "T"    "P"    "E"   "Z"     "Y"]
    idx = 1 + any(abs.(cof) == [1000; 1])
    pfs = pfc[2 - Int(pfx), 9 + Int(xpt[idx])]
    # Convert to string (without prefix || whole part<digits):
    if 4 < abs.(ord[idx]) || floor.(1 + log10.(abs.(cof[idx]))) < sgf
        @eval f = (x, y) -> @sprintf($fmt, x, y)
        return Base.invokelatest(f, cof[idx], pfs)
    end # (whole part>=digits)
    return @sprintf("%.0f %s", cof[idx], pfs)
end


# Returns marix for a rotation by an angle θ around axis in the direction of UNIT vector û.
function rotMat(θ::T, û::Array{T}) where T <: Float64
    @assert norm(û) ≈ 1.0 "`û` must be a unit direction vector."
    ux = [ 0.0  -û[3] û[2];
           û[3]  0.0 -û[1];
          -û[2]  û[1] 0.0] # 'Cross product matrix'
    cosθ = cos(θ)
    a  = [cosθ 0.0 0.0; 0.0 cosθ 0.0; 0.0 0.0 cosθ] # Could be replaced by cos(θ)*eye(3), but currently 10× slower. 

    return a + sin(θ)*ux + (1-cosθ)*(û ⊗ û')
end


function GSEqArc(R; lat=0, halfAngle=45, r=R⨁)
    (halfAngle ≤ lat) && error("Error: No equatorial coverage at latitudes greater than the view half-angle.")

    a   = (cosd(lat)*tand(halfAngle))^2 - sind(lat)^2
    b   = r*(tand(halfAngle)^2)*cosd(lat)

    x₀  = b/a

    ϕ   = atan(√a)
    arc = ϕ - asin(x₀*sin(ϕ)/R)
    return 2arc
end

function horizonAngle(ϕ, λ, R) # [deg], [deg], [km]. Returns angle [rad] from the horizon of a ray from radius R at (lat=0, lon=0) to reference relative lat, lon.
    # Angle sweeping out horizon to relay
    q = acos(cosd(ϕ)*cosd(λ))

    # Distance from horizon to relay
    d = √(R⨁^2 + R^2 - 2R⨁*R*cos(q)) # [km]

    # Angle sweeping out center of the Earth to the relay
    z = acos(0.5(R⨁^2 + d^2 - R^2)/(R⨁*d))

    # Angle from the horizon to the relay
    a = z - 0.5π # [rad]
    return d, a
end

# Finds az, el, pol of antenna A pointing at satellite S. Lats, lons, az & el all radians, r's must be in same unit. Assumes spherical Earth. Source: http://www.eutelsat.com/files/contributed/support/pdf/azimuth-elevation-polarization.pdf
function azElPol(ϕA, λA, λS; ϕS=0, rS=r_E_GEO, rA=R⨁, polS=0) 
    w   = λA - λS
    r   = √((rS*cos(ϕS)*cos(w) - rA*cos(ϕA))^2 + (rS*cos(ϕS)*sin(w))^2 + (rS*sin(ϕS) - rA*sin(ϕA))^2)
    el  = -asin(-(rS*cos(ϕA)*cos(ϕS)*cos(w) + rS*sin(ϕA)*sin(ϕS) - rA)/r)
    az  = atan(-rS*cos(ϕS)*sin(w)/(r*cos(el)), -rS*(sin(ϕA)*cos(ϕS)*cos(w) - cos(ϕA)*sin(ϕS))/(r*cos(el)))

    y   = cos(az)*sin(el)*(sin(ϕA)*sin(ϕS)*sin(polS)*cos(w) - sin(ϕA)*cos(polS)*sin(w) + cos(ϕA)*cos(ϕS)*sin(polS)) + sin(az)*sin(el)*(sin(ϕS)*sin(polS)*sin(w) + cos(polS)*cos(w)) + cos(el)*(cos(ϕA)*sin(ϕS)*sin(polS)*cos(w) - cos(ϕA)*cos(polS)*sin(w) - sin(ϕA)*cos(ϕS)*sin(polS))
    x   = -sin(az)*(sin(ϕA)*sin(ϕS)*sin(polS)*cos(w) - sin(ϕA)*cos(polS)*sin(w) + cos(ϕA)*cos(ϕS)*sin(polS)) + cos(az)*(sin(ϕS)*sin(polS)*sin(w) + cos(polS)*cos(w))
    pol = atan(y, x)

    az, el, pol
end


function putInRange(x::Real, min::Real, max::Real; delta=2π)
    while x > max; x -= delta; end
    while x < min; x += delta; end
    return x
end

function putInRange(A::Array{Float64,1}, min::Real, max::Real; delta=2π)
    for i in eachindex(A)
        while A[i] > max; A[i] -= delta; end
        while A[i] < min; A[i] += delta; end
    end
    return A
end

randRange(min::Real, max::Real) = rand()*(max-min) + min


# Calculates I(r,c), where r is a fractional row number between R and R + 1 and c is a
# fractional column number between C and C + 1, using bi-linear interpolation.
# function biLinearInterpSq(r, c, R, C, I) # From ITU-R P.1144-7, Annex 1, Sec. 1b
    # return I[R,  C  ]*((R+1-r)*(C+1-c)) +
    #        I[R+1,C  ]*((r-R  )*(C+1-c)) +
    #        I[R,  C+1]*((R+1-r)*(c-C  )) +
    #        I[R+1,C+1]*((r-R  )*(c-C  ))
# end


# Find the value of f(x, y) assuming value of f is known at four points (Q₁₁ = f(x₁, y₁), etc.)
function biLinearInterpSq(x, y, x₁, x₂, y₁, y₂, Q₁₁, Q₂₁, Q₁₂, Q₂₂)
    den = (x₂-x₁)*(y₂-y₁)
    num = Q₁₁*(x₂-x)*(y₂-y) +
          Q₂₁*(x-x₁)*(y₂-y) +
          Q₁₂*(x₂-x)*(y-y₁) +
          Q₂₂*(x-x₁)*(y-y₁)
    return num/den
end

# Find the value of f(x, y) assuming value of f is known at four points (2 lats but 4 lons) (QC1 = f(xC, y1), etc.)
# ITU-R P.1144-7, Annex 1, Sec. 1a
function biLinearInterpTrap(x, y, xA, xB, xC, xD, y0, y1, QA0, QB0, QC1, QD1)
    t = (y-y0)/(y1-y0)
    s = (x - xA + t*(xA-xC))/(xB - xA + t*(xA-xC+xD-xB))
    return  (1-s)*(1-t)*QA0 +
            (1-s)*t    *QC1 +
            s    *(1-t)*QB0 +
            s    *t    *QD1
end

# Calculates I(r,c), where r is a fractional row number between R+1 and R+2 and c is a
# fractional column number between C+1 and C+2, using bi-cubic interpolation.
# From ITU-R P.1144-7, Annex 1, Sec. 2
biCubicInterp(r, c, R, C, I) = sum([sum([I[i,j]*K(c-j) for j = C:C+3])*K(r-i) for i = R:R+3])

function K(δ) # Helper function for bicubic interpolation
    a = -0.5
    δ = abs(δ)
        (δ ≤ 1) && return (a + 2)*δ^3 - (a + 3)*δ^2 + 1
    (1 ≤ δ ≤ 2) && return a*δ^3 - 5a*δ^2 + 8a*δ - 4a
                   return 0
end


function boundInd(x, arr) # Finds index of array right before target value x.
    @assert minimum(arr) ≤ x ≤ maximum(arr) "Linear interpolation: target is out of bounds."
    x == maximum(arr) && return length(arr)
    x == minimum(arr) && return 1
    i = 1
    while sign(x - arr[i]) == sign(x - arr[i+1]); i += 1; end
    return i
end

linearInterp(x, C, arr) = (x - arr[C])/(arr[C+1]-arr[C]) + C # Returns the fractional index value c in (C, C+1) corresponding to x in array arr.


function sanitizeSysInputs!(sys)
    if length(sys["refFrame"]) < 5; sys["refFrame"] = lowercase(sys["refFrame"]); end
    if !in(sys["refFrame"], ["mci", "mcl", "aer", "hci", "eci", "ecef", "ecl", "EMsyn", "gcc", "gdc", "ESsyn"])
        error("ERROR: Please specify a valid reference frame")
    end

    if !(2 ≤ length(sys["view"]) ≤ 3); error("'view' must be a 2- or 3-letter combination of x, y, or z.");  end
    for c in sys["view"]
        if !(Int('x') ≤ UInt32(c) ≤ Int('z') || UInt32(c) == Int('t')); error("'view' must be some combination of t, x, y, or z."); end
    end
    sys["view"] = lowercase(sys["view"])

    if length(sys["plot"]) ≠ length(intersect(sys["plot"], ["position", "range", "atmoLoss", "pathLoss", "capacity", "totalGroundLinkCapacity", "clientRates"]))
        error("Invalid item in 'plot' array.")
    end

    if length(sys["report"]) ≠ length(intersect(sys["report"], ["latency", "coverage", "eclipse", "range", "atmoLoss", "pathLoss", "capacity", "totalGroundLinkCapacity", "clientRates"]))
        error("Invalid item in 'report' array.")
    end

    for p in sys["points"]
        if !in(p, ["EML1", "EML2", "EML3", "EML4", "EML5", "ESL1", "ESL2", "ESL3", "ESL4", "ESL5"])
            error("Invalid point specification.")
        end
    end
end


# Returns all positive roots of equation ax^2 + bx + c = 0
function quadratic(a::T, b::T, c::T) where T <: Float64
    d = b^2 - 4a*c
    (d  > 0.0) && return 0.5(-b .+ [1 -1]*√d)/a
    (d == 0.0) && return [-0.5b/a]
                  return []
end


rgb(r, g, b) = (r/255.0, g/255.0, b/255.0)


# Returns range of n evenly-spaced integers between a and b.
function rangeInRange(a, b; n::Int64=100)
    q = ceil((b-a+1)/n)
    c = collect(a:b)
    return filter(c -> c .% q == 0.0, c)
end

nansum(X) = sum(X[!isnan(X)])

nanmean(X) = mean(X[!isnan(X)])

inRange(i::Int64,N::Int64) = ((i+N-1)%N)+1


"""
line2segments(X::Array{Float64,2})

Returns array of n-1 line segment coordinates. Used when plotting different segments of a line
different colors.

### Arguments (Inputs):
| Parameter | Description                          | Units |
| --------: | :----------------------------------- | :---- |
| `X`  | Line defined by `[x1 y1; x2 y2; ...; xn yn]` | DU (e.g. m) |

### Value (Output)
| Parameter | Description       | Units |
| --------: | :---------------- | :---- |
| `Y`  | Line segments: `[[x1 y1; x2 y2], [x2 y2; x3 y3], ..., [xn-1 yn-1; xn yn]]` | DU    |

Author: James Spicer (2018.03.02)
"""
line2segments(X::Array{Float64,2}) = [X[i:i+1,:] for i = 1:size(X,1)-1]


"""Returns Fresnel sine integral of `z`. Tests [here](https://github.com/ebertolazzi/Clothoids/blob/master/src/Fresnel.cc)"""
function fresnelS(z::Float64)
    k = 0.5abs(z)*√π
    S = real(0.25(1+im)*(erf(k*(1+im)) - im*erf(k*(1-im))))
    return sign(z)*S
end


"""Returns Fresnel cosine integral of `z`. Tests [here](https://github.com/ebertolazzi/Clothoids/blob/master/src/Fresnel.cc)"""
function fresnelC(z::Float64)
    k = 0.5abs(z)*√π
    C = real(0.25(1-im)*(erf(k*(1+im)) + im*erf(k*(1-im))))
    return sign(z)*C
end