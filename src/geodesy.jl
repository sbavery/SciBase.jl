export
    googleDist,
    googleElev,
    finalCoord,
    bearing,
    ellipsoidDistance








# A, B = (lon, lat) [deg, deg]. Returns distance between A, B in m
function googleDist(A, B; key="AIzaSyAccIcYQuXjjJWTJR2VXVki9aJYY6bAK-A", units="m")
    page = get("https://maps.googleapis.com/maps/api/distancematrix/json?units=metric&origins=$(A[2]),$(A[1])&destinations=$(B[2])%2C$(B[1])%7C&key=$(key)")
    csv = readcsv(IOBuffer(page.data))
    return parse(Int, csv[10, 1][29:end])
end

# lon [deg], lat [deg]. Returns elevation in m. Does not take into account building height
function googleElev(lon, lat; key="AIzaSyAOZjBoSfW2z2i4bJNQ2ZDiCMud5EkQEzc")
    page = get("https://maps.googleapis.com/maps/api/elevation/json?locations=$(lat),$(lon)&key=$(key)")
    csv = readcsv(IOBuffer(page.data))
    return parse(Float64, split(csv[4])[end])
end




# Turns "N 45°0'0" E" to "45.0" etc. See http://www.georgialandsales.com/blogspot/?p=266
# function bearing2deg(NS, num, EW)
#     if length(num) == 1
#         num = (num, 0, 0)
#     end

#     bearing = degMinSec2Deg(num[1], num[2], num[3])

#     if EW == "W"; bearing = 360 - bearing; end
#     if NS == "S"; bearing = 180 - bearing; end

#     return putInRange(bearing, 0, 360, delta=360)
# end

# d in km, λ₁, ϕ₁, θ in rad
# see http://www.movable-type.co.uk/scripts/latlong.html
# Vincenty Source: https://www.ngs.noaa.gov/PUBS_LIB/inverse.pdf
function finalCoord(λ₁::T, ϕ₁::T, θ::T, d::T; R::T=R⨁, f::T=f⨁, shape::String="sphere", max_i::Int64=10, tol::T=1e-12, a::T=R⨁_a, b::T=R⨁_b) where T <: Float64
	if shape == "sphere"
	    δ  = d/R
	    ϕ₂ = asin(sin(ϕ₁)*cos(δ) + cos(ϕ₁)*sin(δ)*cos(θ))
	    λ₂ = λ₁ + atan(sin(θ)*sin(δ)*cos(ϕ₁), cos(δ) - sin(ϕ₁)*sin(ϕ₂))
    	return λ₂, ϕ₂

	elseif shape == "ellipsoid"
        s, α₁  = d, θ 

        tU₁ = (1-f)*tan(ϕ₁)
		U₁  = atan(tU₁)
		sU₁, cU₁ = sin(U₁), cos(U₁)
        sα₁, cα₁ = sin(α₁), cos(α₁) 
		sα  = cU₁*sα₁                                             # Eqn. 2
		c²α = 1 - sα^2
		u²  = c²α*(a^2/b^2 - 1)

		A   = u²*(4096 + u²*(u²*(320 - 175u²)-768))/16384 + 1     # Eqn. 3
		B   = u²*( 256 + u²*(u²*( 74 -  47u²)-128))/ 1024         # Eqn. 4

        σ₀, σ₁ = s/(b*A), atan(tU₁, √c²α)
	    σ, i, sσ, cσ, c2σm = σ₀, 1, 0.0, 0.0, 0.0
	    while true
            c2σm, sσ, cσ = cos(2σ₁ + σ), sin(σ), cos(σ)           # Eqn. 5
            Δσ = B*sσ*(c2σm + 0.25B*(cσ*(-1+2c2σm^2) - B*c2σm*(-3+4sσ^2)*(-3+4c2σm^2)/6)) # Eqn. 6
            (abs(Δσ) ≤ tol || i > max_i) && (σ = σ₀ + Δσ; break)
	        σ  = σ₀ + Δσ                                          # Eqn. 7
            i += 1
		end

		ϕ₂ = atan(sU₁*cσ + cU₁*sσ*cα₁, (1-f)*hypot(sα, sU₁*sσ - cU₁*cσ*cα₁)) # Eqn. 8

		λ  = atan(sσ*sα₁, cU₁*cσ - sU₁*sσ*cα₁)                   # Eqn. 9
		C  = f*c²α*(0.25 + f*(0.25 - 0.1875c²α))                  # Eqn. 10
		L  = λ - (1-C)*f*sα*(σ + C*sσ*(c2σm + C*cσ*(-1+2c2σm^2))) # Eqn. 11

		return λ₁ + L, ϕ₂
	end
end


# R, L in km
# Delta, lon1, lat1 in rad
# function platCurve(R, Delta, lon1, lat1, bearing1; n=10)
#     dist = 2*R*sin(Delta*0.5/n)

#     lon, lat = lon1, lat1
#     bearing = bearing1
#     for i = 1:n
#         lon, lat = finalCoordEllipse(lon, lat, bearing + i*Delta/n, dist)
#         print(rad2deg(lon), ",", rad2deg(lat), ",0 ")
#     end
#     return (lon, lat)
# end


# Finds bearing from (λ₁, ϕ₁) to (λ₂, ϕ₂) [all in rad]. Source: http://www.movable-type.co.uk/scripts/latlong.html
bearing(λ₁, ϕ₁, λ₂, ϕ₂) = atan(sin(λ₂ - λ₁)*cos(ϕ₂), cos(ϕ₁)*sin(ϕ₂) - sin(ϕ₁)*cos(ϕ₂)*cos(λ₂ - λ₁))


# Uses Vincenty's algorithm to calculate the distance between two points on an ellipsoid. Coordinate inputs in [rad]
function ellipsoidDistance(λ₁::T, ϕ₁::T, λ₂::T, ϕ₂::T; max_i::Int64=200, tol::T=1e-12, R::T=R⨁, f::T=f⨁, a::T=R⨁_a, b::T=R⨁_b) where T <: Float64
    λ₀       = λ₂ - λ₁
    U        = atan.((1-f)*tan.([ϕ₁ ϕ₂]))
    sU₁, sU₂ = sin.(U)
    cU₁, cU₂ = cos.(U)

    λ, i, c²α, c2σm, sσ, cσ, σ = λ₀, 1, 0.0, 0.0, 0.0, 0.0, 0.0
    while true
        cλ, sλ = cos(λ), sin(λ)
        sσ     = hypot(cU₂*sλ, cU₁*sU₂ - sU₁*cU₂*cλ)
        cσ     = sU₁*sU₂ + cU₁*cU₂*cλ
        σ      = atan(sσ, cσ)
        sα     = csc(σ)*(cU₁*cU₂*sλ)
        c²α    = 1 - sα^2
        c2σm   = cσ - 2sU₁*sU₂/c²α
        C      = 0.0625f*c²α*(4 + f*(4 - 3c²α))
        λ′, λ  = λ, λ₀ + (1-C)*f*sα*(σ + C*sσ*(c2σm + C*cσ*(2c2σm^2-1)))
        (abs(λ - λ′) ≤ tol || i > max_i) && break
        i += 1
    end

    u² = c²α*((a/b)^2 - 1)
    A  = 1 + u²*(4096 + u²*(u²*(320-175u²) - 768))/16384
    B  =     u²*( 256 + u²*(u²*( 74- 47u²) - 128))/ 1024
    Δσ = B*sσ*(c2σm + 0.25B*(cσ*(2c2σm^2-1) - B*c2σm*(4sσ^2-3)*(4c2σm^2-3)/6))

    return b*A*(σ-Δσ)
end