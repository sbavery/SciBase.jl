export
	JD2JDC,
    # JD2GMST,
    # GMST2GAST,
    LSTime,
    days2mdh,
    mdh2days,
    hms2rad,
    rad2hms,
    hms2τ,
    τ2hms,
    dayOfWeek,
    jday,
    lightTime




#=
# Based on Wolfram algorithm, works for all CE Gregorian dates
function juliandate(y, mo, d, h, mi, s)

    UT = h + (mi/60) + (s/3600);

    JD0 = 367*y - floor(7*(y + floor((mo+9)/12))/4) - floor(3*(floor((y + (mo-9)/7)/100)+1)/4) + floor(275*mo/9) + d + 1721028.5;

    JD = JD0 + UT/24;

    JD, JD0

end
=#

"""
    JD2JDC(JD::Float64) → Float64

Compute Julian centuries from a Julian date.
"""
JD2JDC(JD::Float64) = (JD - 2451545)/36525


"""
    LSTime(JD_UT1::T; λ::T=0.0) where T <: Float64

Compute Greenwich Meridian sidereal time angle (θ_GMST) and 
local sidereal time angle at longitude λ (θ_LST) at a UT1 Julian date. 
"""
function LSTime(JD_UT1::T; λ::T=0.0) where T <: Float64
    T_UT1  = JD2JDC(JD_UT1) # Centuries
    θ_GMST =   -6.2e-6*T_UT1^3 + 
             9.3104e-2*T_UT1^2 + 
            (876600*3600 + 8640184.812866)*T_UT1 + 
            67310.54841 # [sec]

    θ_GMST = mod2pi(θ_GMST*2π/solSecsEarth) # [rad]

    return θ_GMST, θ_GMST + λ
end

# Find Greenwich Mean Sidereal Time based on USNO algorithm
# function JD2GMST(JD::Float64)
#     D  = JD - 2451545
#     D₀ = datetime2julian(trunc(julian2datetime(JD), Day)) - 2451545

#     T  = JD2JDC(JD)
#     H  = 24.0*(D-D₀)

#     GMST = 6.697374558 + 0.06570982441908*D₀ + 1.00273790935*H + 0.000026*(T^2) # [Hrs]

#     while GMST > 24; GMST -= 24; end
#     GMST, D
# end


# # Find Greenwich Apparent Sidereal Time based on USNO algorithm
# function GMST2GAST(GMST, D)
#     ω = 125.04   - 0.052954 *D
#     L = 280.4665 + 0.98565  *D
#     e =  23.4393 - 0.0000004*D

#     ω = deg2rad(ω)
#     L = deg2rad(L)
#     e = deg2rad(e)

#     deltaPsi = -0.000319*sin(ω) - 0.000024*sin(2*L)
#     eqeq = deltaPsi*cos(e)
#     GAST = GMST + eqeq

#     while GAST > 24.0; GAST -= 24.0; end
#     return GAST
# end

# Takes in a year and day of year + fraction, returns a DateTime. Source: Vallado
function days2mdh(year::Int64, days::Float64; output::String="DT")
    # --------------- set up array of days in month  --------------
    lmonths = [31, (rem(year-1900, 4) == 0) ? 29 : 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

    dayofyr = floor(days)

    # ----------------- find month and day of month ---------------
    inttemp, mon = 0, 1
    while dayofyr > inttemp + lmonths[mon] && mon < 12
        inttemp += lmonths[mon]
        mon += 1
    end

    day = dayofyr - inttemp

    # ----------------- find hours minutes and seconds ------------
    temp = 24(days - dayofyr)
    H   = trunc(temp)
    temp = 60(temp - H)
    M  = trunc(temp)
    temp = 60(temp - M)
    (output ≠ "DT") && (return Int(mon), Int(day), Int(H), Int(M), temp)
    S  = trunc(temp)
    temp = 1000(temp - S)
    ms   = trunc(temp)

    return DateTime(Int(year), Int(mon), Int(day), Int(H), Int(M), Int(S), Int(ms))
end

function mdh2days(D::DateTime)
    lmonths = [31, (rem(Dates.year(D) - 1900, 4) == 0) ? 29 : 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

    m  = Dates.month(D)
    d  = Dates.day(D)
    h  = Dates.hour(D)
    M  = Dates.minute(D)
    s  = Dates.second(D)
    ms = Dates.millisecond(D)

    return sum(lmonths[1:m-1]) + d + h/24 + M/(24*60) + s/(24*60*60) + ms/(24*60*60*1000)
end



# Source: Vallado HMStoRad
hms2rad(h::T, m::T, s::Float64) where T <: Int64 = deg2rad(15(h + m/60 + s/3600))


# Source: Vallado radToHMS
function rad2hms(τ::Float64)
    t = rad2deg(τ)/15
    h = trunc(t)
    m = trunc(60(t - h))
    s = 3600(t - h - m/60)
    return [h m s]
end


# Source: Vallado HMStoτ (time of day)
hms2τ(h::T, m::T, s::Float64) where T <: Int64 = 3600h + 60m + s


# Source: Vallado τ(time of day)ToHMS
function τ2hms(τ::Float64)
    t = τ/3600
    h = trunc(t)
    m = trunc(60(t - h))
    s = 3600(t - h - m/60)
    return [h m s]
end


function dayOfWeek(y, m, d) # Returns 0 for Sunday, 6 for Saturday, etc.
    t  = (0, 3, 2, 5, 0, 3, 5, 1, 4, 6, 2, 4)
    y -= (m < 3)
    return (y + floor(0.25y) - floor(0.01y) + floor(0.0025y) + t[m] + d) % 7
end


function jday(yr::T, mon::T, day::T, hr::T, min::T, sec::Float64) where T <: Int64
    jd = 367 * yr -
         floor( (7 * (yr + floor( (mon + 9) / 12) ) ) * 0.25 ) +
         floor( 275 * mon / 9 ) +
         day + 1721013.5 +
         ( (sec/60 + min ) / 60 + hr ) / 24
    #  - 0.5 * sign(100 * yr + mon - 190002.5) + 0.5;
    # return jd
end

"""Calculate time taken for light to travel distance `d` (km) through a vacuum."""
lightTime(d) = d/c₀