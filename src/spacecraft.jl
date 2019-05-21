export
    spacecraft,
    addSpacecraft,
    loadSpacecraft!,
    updateSCLib,
    listSC,
    TLE2OrbEl,
    parseTLE,
    orbEl2TLE,
    loadGEO!

# TLEurls = [ "http://www.celestrak.com/NORAD/elements/resource.txt",
#             "http://www.celestrak.com/NORAD/elements/stations.txt",
#             "http://www.celestrak.com/NORAD/elements/geo.txt",
#             "http://www.celestrak.com/NORAD/elements/intelsat.txt",
#             "http://www.celestrak.com/NORAD/elements/science.txt",
#             "http://www.celestrak.com/NORAD/elements/orbcomm.txt",
#             "http://www.celestrak.com/NORAD/elements/other-comm.txt",
#             "http://www.celestrak.com/NORAD/elements/visual.txt",
#             "http://ephemerides.planet-labs.com/planet_mc.tle"]



mutable struct Spacecraft
    name::String
    center::String
    t₀::Float64         # Epoch [JD]
    X⃗₀::Array{Float64}  # State at epoch [km, km/s] in bci.
    children            # Array to be filled with subcomponents, e.g. antennas, etc.
end

# '0' = 'at epoch'
# u     = Argument of latitude [rad]
# ω̃true = True Longitude of Perigee [rad]
# λtrue = True Longitude [rad]
function addSpacecraft(name::String; center::String="EARTH", e::Float64=0.0, ra::Float64=0.0, rp::Float64=0.0, ha::Float64=0.0, hp::Float64=0.0, a::Float64=0.0, p::Float64=0.0, i::Float64=0.0, epoch=DateTime(2000), u::Float64=0.0, λtrue::Float64=0.0, ω̃true::Float64=0.0, RAAN::Float64=0.0, AP::Float64=0.0, T::Float64=0.0, n::Float64=0.0, MA0::Float64=0.0, TA0::Float64=0.0, EA0::Float64=0.0, updatingLib::Bool=false, tol::Float64=1e-6)
    # global dName

    name   = uppercase(name)
    name   = replace(name, "/" => "_")
    center = uppercase(center)
    μ      = getMu(center)                          # [km3/s2]
    t₀     = datetime2julian(epoch)                 # day
    (ha  ≠ 0.0) && (ra = ha + R⨁)                   # Todo: make this generic for any planet
    (hp  ≠ 0.0) && (rp = hp + R⨁)
    (a  == 0.0) && (a  = findA(μ, n, T, ra, rp, e)) # [km]
    (n  == 0.0) && (n  = meanMotionA(a, μ))         # [rad/s]
    (e  == 0.0) && (e  = findE(ra, rp, a))
    # (RAAN == 0.0 && LAN ≠ 0.0) && (RAAN = LAN)
    # if (MA0  == 0.0) && [TA0, EA0, AL, TAN, TP] ≠ [0.0, 0.0, 0.0, 0.0, 0.0]
    #     MA0  = findMA0(e, TA0, EA0, AL, TAN, TP, AP, t₀, n)
    # end
    p  = semilatusRectum(a, e)


    if (TA0  == 0.0) && [MA0, EA0] ≠ [0.0, 0.0, 0.0, 0.0, 0.0]
        TA0  = findTA0(e, MA0, EA0)
    end

    # Special cases
    # if e < tol
    #     if i < tol || abs(i-π) < tol    # circular equatorial
    #         (λtrue == 0.0) && (λtrue = TA0)
    #     else                            # circular inclined
    #         (u     == 0.0) && (u     = TA0)
    #     end
    # else                                # elliptical equatorial
    #     if i < tol || abs(i-π) < tol
    #         (ω̃true == 0.0) && (ω̃true = AP)
    #     end
    # end 

    # println("hi")
    X⃗₀ = RANDV(p, e, i, RAAN, AP, TA0, μ=μ)#, u=u, λtrue=λtrue, ω̃true=ω̃true)
    # println("hu")
    # @show r⃗₀, v⃗₀
    SC = Spacecraft(name, center, t₀, X⃗₀, Any[])
    # println("ho")
    if !updatingLib
        jldopen(joinpath(datapath, "scLib.jld"), "r+") do scLib
            if !exists(scLib, name)
                write(scLib, name, SC)
            else
                o_delete(scLib, name)
                write(scLib, name, SC)
            end
        end
    #     # writecsv(dName * "/data/scLib/$(name).txt", [getfield(SC, i) for i in 1:6])#10])
    else
        return SC
    end
    # println("hum")
end

function loadSpacecraft!(sys)
    # global dName
    for key in [sys["clients"]; sys["constellation"]]
        SC = load(joinpath(datapath, "scLib.jld"), key)

        X⃗₀, t₀, μ = SC.X⃗₀, SC.t₀, getMu(SC.center) # [km, km/s, JD]
        X⃗ = zeros(sys["tSteps"], 6)
        for i = 1:sys["tSteps"]
            Δt = (sys["t"][i] - t₀) * solSecsEarth # [sec]
            X⃗[i,:] = kepler(X⃗₀, Δt, μ=μ)
        end
        sys[key]["coords"] = X⃗[:,1:3]
        sys[key]["vel"]    = X⃗[:,4:6]

        (SC.center == "EARTH") && (sys[key]["frame"] = "eci")
        (SC.center == "MARS" ) && (sys[key]["frame"] = "mci")
        (SC.center == "SUN"  ) && (sys[key]["frame"] = "hci")

        if !isempty(SC.children)
            sys[key]["antennas"] = SC.children
        end
    end
end


function findE(ra::T, rp::T, a::T) where T <: Float64
    (ra       == rp)       && return 0.0
    (ra ≠ 0.0 && rp ≠ 0.0) && return (ra - rp)/(ra + rp)
    (a  ≠ 0.0 && rp ≠ 0.0) && return 1 - rp/a
    (a  ≠ 0.0 && ra ≠ 0.0) && return ra/a - 1
                              return 0.0
end

function findA(μ, n, T, ra, rp, e)
    (n  ≠ 0.0)              && return semimajorAxisN(n, μ)
    (T  ≠ 0.0)              && return semimajorAxisT(T, μ)
    (rp ≠ 0.0 && ra == 0.0) && return rp/(1-e)
    (ra ≠ 0.0 && rp == 0.0) && return ra/(1+e)
    (rp ≠ 0.0 && ra  ≠ 0.0) && return 0.5(ra + rp)
                               return 0.0
end

# function findMA0(e, TA0, EA0, AL, TAN, TP, AP, epoch, n)
#     if TA0 ≠ 0.0
#         EA0 = eccentricAnomalyNu(TA0, e)
#         return meanAnomaly(EA0, e)
#     elseif EA0 ≠ 0.0
#         return meanAnomaly(EA0, e)
#         return 0.0
#     end
# end

# Find ν given MA or A(EA/PA/HA). 
function findTA0(e::T, MA::T, A::T) where T <: Float64
    (MA ≠ 0.0) && (A = anomaly(e, M=MA))
    ( A ≠ 0.0) && return A2ν(A, e)
    return 0.0
end


# Converts perifocal coordinates to body-centered inertial coordinates
# function peri2bci(t, μ, epoch, ω, ι, Ω, e, M₀, p)
#     n = meanMotionP(p, μ, e) # [rad/s]
#     peri_C_eci = rot3mat(-Ω) * rot1mat(-ι) * rot3mat(-ω)
#     r = zeros(length(t), 3)
#     # v = zeros(length(t), 3)
#     for i = 1:length(t)
#         ∆t = (t[i]-epoch)*solSecsEarth   # [sec]
#         M  = mod2pi(M₀ + n*∆t)           # [rad]
#         E  = eccentricAnomalyM(M, e)     # [rad]
#         ν  = trueAnomaly(E, e)           # [rad]
#         O_r_S_peri = p*[cos(ν); sin(ν); 0.0]/(1 + e*cos(ν))
#         r[i,:]     = peri_C_eci * O_r_S_peri

#         # r, v = coe2peri(p, e, μ, ν) # Classical orbital elements -> perifocal
#         # r[i,:], v[i,:] = peri2bci(r, ω, ι, Ω)', peri2bci(v, ω, ι, Ω)' # perifocal to body-centered inertial
#     end
#     return r#, v
# end


function readTLEurl(TLEurl::String)
    page = HTTP.request("GET", TLEurl)

    lines = split(String(page.body), "\r\n")[1:end-1]
    sc = String[]
    for i in 1:3:length(lines); push!(sc, join(lines[i:i+2], "\n")); end
    return sc
end


function TLE2OrbEl(s::String)
    lines = split(s, "\n")
    name = (split(lines[1])[1] == "0") ? join(split(lines[1])[2:end]) : join(split(lines[1]))
    t = split(join(lines[2:3], " "))

    return name,
        parse(Int, t[4][1:2]) + 2000,           # Epoch year, e.g. 2016
        parse(Float64, t[4][3:end]),            # day of year + fraction [day]
        deg2rad(parse(Float64, t[12])),         # ι [rad]
        deg2rad(parse(Float64, t[13])),         # Ω [rad]
        parse(Float64, "0." * t[14]),           # e [-]
        deg2rad(parse(Float64, t[15])),         # ω [rad]
        deg2rad(parse(Float64, t[16])),         # M [rad]
        parse(Float64, t[17])*2π/solSecsEarth   # n [rad/sec]
end


"""
    parseTLE(line1::String, line2::String) -> Int64, String, String, Float64, Float64, Float64, Float64, Int64,
    Float64, Float64, Float64, Float64, Float64, Float64, Int64

Parse two-line element sets (TLEs) into object parameters.

### Arguments (Inputs)
  * `line1::String` : First line of TLE.
  * `line2::String` : Second line of TLE.

### Values (Outputs)
  * `satnum1`  : Object number                                                      e.g.  25544
  * `class`    : Classification (U=Unclassified)                                    e.g.  "U"
  * `intldesg` : International Designator                                           e.g.  98067A
  * `epoch`    : Object epoch                                          [Julian day] e.g.  2.4576921763149304e6
  * `ṅ`        : Mean motion 1st derivative ÷ 2 (Ballistic Coeffient)  [rev/day²]   e.g. −2.182e-5
  * `n̈`        : Mean motion 2nd derivative ÷ 6                        [rev/day³]   e.g.  0.000
  * `b★`       : 'BSTAR' drag term                                     [-]          e.g. -1.1606e-4
  * `TLEnum`   : Element set number.                                                e.g.  292
  * `i`        : Inclination                                           [rad]        e.g.  0.90131
  * `Ω`        : Right ascension of the ascending node                 [rad]        e.g.  4.31904
  * `e`        : Eccentricity                                          [-]          e.g.  6.703e-4
  * `ω`        : Argument of perigee                                   [rad]        e.g.  2.278283
  * `M`        : Mean Anomaly                                          [rad]        e.g.  5.672823
  * `n`        : Mean Motion                                           [rev/day]    e.g. 15.72125391
  * `revnum`   : Revolution number at epoch                            [rev]        e.g. 56353

Author: James Spicer (2016.10.30)
"""
function parseTLE(line1::String, line2::String)
    (length(line1) ≠ 69 || length(line2) ≠ 69) && error("TLE input error: lines should both be 69 characters long.")

    for j = 11:16
        (line1[j] == ' ') && (line1 = line1[1:j-1] * "_" * line1[j+1:end])
    end

    (line1[45]  ≠ ' ') && (line1 = line1[1:43] * string(line1[45]) * line1[45:end])
    line1 = line1[1:44] * "." * line1[46:end]

    (line1[ 8] == ' ') && (line1 = line1[1:7] * "U" * line1[9:end])

    (line1[10] == ' ') && (line1 = line1[1:9] * "." * line1[11:end])

    for j = 46:50
        (line1[j] == ' ') && (line1 = line1[1:j-1] * "0" * line1[j+1:end])
    end

    (line1[52] == ' ') && (line1 = line1[1:51] * "0" * line1[53:end])

    (line1[54]  ≠ ' ') && (line1 = line1[1:52] * string(line1[54]) * line1[54:end])

    line1 = line1[1:53] * "." * line1[55:end]
    line2 = line2[1:25] * "." * line2[27:end]

    for j = 27:33
        (line2[j] == ' ') && (line2 = line2[1:j-1] * "0" * line2[j+1:end])
    end

    line1[63] == ' ' && (line1 = line1[1:62] * "0" * line1[64:end])

    if length(line1) < 68 || line1[68] == ' '
        line1 = line1[1:67] * "0" * line1[69:end]
    end

    # parse first line
    linenum1    = parse(Int64,   line1[ 1: 1])
    satnum1     = parse(Int64,   line1[ 3: 7])
    class       = line1[8]
    intldesg    = line1[10:17]
    epochyr     = parse(Int64,   line1[19:20])
    epochdays   = parse(Float64, line1[21:32])
    ṅ           = parse(Float64, line1[34:43])
    n̈           = parse(Float64, line1[45:50]) *
               10^parse(Float64, line1[51:52])
    b★          = parse(Float64, line1[53:59]) *
               10^parse(Float64, line1[60:61])
    zero        = parse(Int64,   line1[63:63])
    TLEnum      = parse(Int64,   line1[65:68])

    # parse second line
    linenum2    = parse(Int64,   line2[ 1: 1])
    satnum2     = parse(Int64,   line2[ 3: 7])
    i           = parse(Float64, line2[ 8:16])
    Ω           = parse(Float64, line2[17:25])
    e           = parse(Float64, line2[26:33])
    ω           = parse(Float64, line2[34:42])
    M           = parse(Float64, line2[43:51])
    n           = parse(Float64, line2[52:63])
    revnum      = parse(Int64,   line2[64:68])

    (linenum1 ≠ 1 || linenum2 ≠ 2) && error("Invalid TLE, line numbers not correct.")

    (satnum1 ≠ satnum2) && error("Invalid TLE, satellite numbers on lines 1 and 2 are not identical.")

    epochyr = (epochyr < 57) ? epochyr + 2000 : epochyr + 1900

    mon, day, hr, minute, sec = days2mdh(epochyr, epochdays, output="sec")

    epoch = jday(epochyr, mon, day, hr, minute, sec) # [julian days]

    return satnum1, class, intldesg, epoch, ṅ, n̈, b★, TLEnum, deg2rad(i), deg2rad(Ω), e, deg2rad(ω), deg2rad(M), n, revnum
end

function orbEl2TLE(satnum1::I, class::S, intldesg::S, epoch::F, ṅ::F, n̈::F, b★::F, TLEnum::I, i::F, Ω::F, e::F, ω::F, M::F, n::F, revnum::I) where {F <: Float64, S <: String, I <: Int64}
    # convert epoch (JDay) to year, day of year + fraction    
    dt = julian2datetime(epoch)
    yr = Dates.year(dt)
    yr -= (yr < 2000) ? 1900 : 2000
    days = mdh2days(dt)
    
    @assert length(class) == 1
    @assert length(intldesg) == 8
    @assert yr < 100 # must be 2-digit
        
    bstar = TLEexpStr(b★)
    ndot = getNdot(ṅ)
    nddot = TLEexpStr(n̈)
    
    @assert length(ndot)  == 10
    @assert length(nddot) ==  8
    @assert length(bstar) ==  8
    
    line1 = @sprintf "1 %05d%s %s %02d%3.8f %s %s %s 0 %04d"  satnum1 class intldesg yr days ndot nddot bstar TLEnum
    line1 *= string(TLEchecksum(line1))
    
    e_str = @sprintf "%0.7f" e
    line2 = @sprintf "2 %05d %08.4f %08.4f %s %08.4f %08.4f %11.8f%05d" satnum1 rad2deg(i) rad2deg(Ω) e_str[3:end] rad2deg(ω) rad2deg(M) n revnum
    line2 *= string(TLEchecksum(line2))
    
    @assert length(line1) == 69
    @assert length(line2) == 69
    return line1, line2
end

# letters, blanks, periods, plus signs = 0, minus signs = 1. Source: Celestrak.
function TLEchecksum(line)
    sum = 0
    for c in line
        if isalpha(c) || c == ' ' || c == '.' || c == '+' || c == '_'
            continue
        elseif c == '-'
            sum += 1
        else
            sum += parse(Int64, c)
        end
    end    
    return sum % 10
end

# Turn -0.11606x10^-4 to "-11606-4"
function TLEexpStr(n)
    n == 0 && return "-00000-0"
    exponent = fix(log10(abs(n)))
    temp = n*10^-exponent
    str = @sprintf "%0.5f" temp
    signstr = (n < 0) ? "-" : "+"
    return @sprintf "%s%s%+d" signstr str[4:end] exponent
end

function getNdot(ṅ)
    str = @sprintf "%+10.8f"  ṅ
    return string(str[1]) * str[3:end] # remove leading 0
end


function updateSCLib()
    # global dName

    scNames = Set()
    jldopen(joinpath(datapath, "scLib.jld"), "w") do scLib
        addrequire(scLib, AstroSciKit)
        for TLEurl in TLEurls
            pageSC = readTLEurl(TLEurl)
            for scTLE in pageSC
                lines = split(scTLE, "\n")
                name = (split(lines[1])[1] == "0") ? join(split(lines[1])[2:end]) : join(split(lines[1]))
                if !in(name, scNames)
                    name, yr, t₀, ι, Ω, e, ω, M₀, MM = TLE2OrbEl(scTLE)
                    SC = addSpacecraft(name, epoch=days2mdh(yr, t₀), i=ι, RAAN=Ω, e=e, AP=ω, MA0=M₀, n=MM, updatingLib=true)
                    write(scLib, name, SC)
                    push!(scNames, name)
                end
            end
        end
    end
    @info "Spacecraft library updated"
end


function listSC()
    # global dName
    scDict = load(joinpath(datapath, "scLib.jld"))
    for key in sort(collect(keys(scDict))); println(key); end
    # for file in readdir(dName * "/data/scLib/"); println(file[1:end-4]); end
end


function loadGEO!(sys)
    θ = linspace(0, 2π, sys["tSteps"])
    sys["GEO"]["coords"] = r_E_GEO*[cos.(θ) sin.(θ) zeros(length(θ))]
    sys["GEO"]["frame"] = "ecef"
end