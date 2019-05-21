export
	A2ν,
	ν2A,
	kepEqtnE,
	kepEqtnP,
	kepEqtnH,
	anomaly,
	c₂c₃,
	kepler,
    pKepler,
    ELORB,
	SGPwrap,
	SGP,
    predict,
    visibilityAngles,
	vCirc,
    deltaVcirc,
    orbitalEnergy,
    ellipseRadius,
    ellipseSpeed,
    period,
    semimajorAxisT,
    semimajorAxisN,
    meanMotionA,
    meanMotionP,
    deltaV,
	ve,
	deltaVEllipse,
	v0,
	trueAnomaly,
    meanAnomaly,
    eccentricAnomalyM,
    eccentricAnomalyNu,
    semilatusRectum,
    e2c,
    c2e,
    e2rv,
    rv2e,
    me2c,
    c2me,
    me2rv,
    eqMo_ME



"""
	A2ν(A::T, e::T; tol::T=1e-6) where T <: Float64

Compute true anomaly from eccentric (E), parabolic (B), or (H) hyperbolic anomaly.

Author: James Spicer (2016.10.30)
Source: Vallado (2001), Algorithm 6, p. 85
"""
function A2ν(A::T, e::T; tol::T=1e-6) where T <: Float64
    # (e ≤ tol)        && return NaN # Commenting out as it messes with MA → ν for circular orbits.
    if e ≤ 1 - tol
        cA = cos(A)
        return mod2pi(atan(  sin(A)*√(1-e^2)/(1-e*cA), (cA-e)/(1-e*cA)))
    elseif abs(e - 1) ≤ tol
        return mod2pi(2atan(A))
    else
        cA = cosh(A)
        return        atan(-sinh(A)*√(e^2-1)/(1-e*cA), (cA-e)/(1-e*cA)) # No mod2pi!
    end
end


"""
	ν2A(ν::T, e::T; tol::T=1e-6) where T <: Float64

Compute eccentric (E), parabolic (B), or hyperbolic (H) anomaly from true anomaly.

Author: James Spicer (2016.10.30)
Source: Vallado (2001), Algorithm 5, p. 85
"""
function ν2A(ν::T, e::T; tol::T=1e-6) where T <: Float64
    e ≤ tol && return NaN
    if e ≤ 1 - tol
        cν = cos(ν)
        return mod2pi(atan(sin(ν)*√(1-e^2)/(1+e*cν), (e+cν)/(1+e*cν)))
    elseif abs(e - 1) ≤ tol
        abs(abs(ν)/π % 2 - 1) ≤ tol && @error "ν input value of ±π not possible for parabola."
        return tan(0.5ν) # Don't use mod2pi!
    else
        νM = ν
        νM >  π && (νM -= 2π)
        νM < -π && (νM += 2π)
        abs(νM) ≥ π - asec(e) && @error "ν, e input values not physically possible."
        return asinh(sin(ν)*√(e^2-1)/(1+e*cos(ν))) # Don't use mod2pi!
    end
end


"""
	kepEqtnE(M::T, e::T; tol::T=1e-8, max_i::Int64=20) where T <: Float64

Solve Kepler's equation for elliptical orbits.

Author: James Spicer (2016.10.30)
Source: Vallado (2001), Algorithm 2, p.73
"""
function kepEqtnE(M::T, e::T; tol::T=1e-8, max_i::Int64=20) where T <: Float64
    E, i = (-π < M < 0.0 || M > π) ? M - e : M + e, 1
    while true
        ΔE = (M - E + e*sin(E)) / (1 - e*cos(E))
        (abs(ΔE) ≤ tol || i > max_i) && (E += ΔE; break)
        E += ΔE
        i += 1
    end
    return mod2pi(E)
end


"""
	kepEqtnP(∆t, p; μ=μ⨁)

Solve Kepler's equation for parabolic orbits

Author: James Spicer (2016.10.30)
Source: Vallado (2001), Algorithm 3, p.77
"""
function kepEqtnP(∆t, p; μ=μ⨁)
	n = 2*√(μ/p^3)
	s = 0.5acot(1.5n*∆t)
	w = atan(∛tan(s))
	B = 2cot(2w)
end


"""
	kepEqtnH(M::T, e::T; tol::T=1e-8, max_i::Int64=20) where T <: Float64

Solve Kepler's equation for hyperbolic orbits

Author: James Spicer (2016.10.30)
Source: Vallado (2001), Algorithm 4, p.79
"""
function kepEqtnH(M::T, e::T; tol::T=1e-8, max_i::Int64=20) where T <: Float64
	if e < 1.6
		H = (-π < M < 0.0 || M > π) ? M - e  		: M + e
	else
		H = (e < 3.6 && abs(M) > π) ? M - sign(M)*e : M/(e-1)
	end

	i = 1
    while true
		ΔH = (M - e*sinh(H) + H) / (e*cosh(H) - 1)
        (abs(ΔH) ≤ tol || i > max_i) && (H += ΔH; break)
        H += ΔH
		i += 1
	end
	return mod2pi(H)
end

"Compute eccentric, parabolic, or hyperbolic anomaly for any eccentricity."
function anomaly(e::T; M::T=0.0, Δt::T=0.0, p::T=0.0, μ::T=μ⨁, tol::T=1e-8, max_i::Int64=20) where T <: Float64
	(e == 0.0) && return M
	(e  < 1.0) && return kepEqtnE(M, e, tol=tol, max_i=max_i)
	(e == 1.0) && return kepEqtnP(Δt, p, μ=μ)
	(e  > 1.0) && return kepEqtnH(M, e, tol=tol, max_i=max_i)
    @error "Invalid eccentricity value"
end


"""
    c₂c₃(ψ::T; tol::T=1e-6) where T <: Float64

Compute c₂ and c₃ values for Kepler algorithm.

Author: James Spicer (2016.10.30) <br/> 
Source: Vallado (2001), Algorithm 1, p.71
"""
function c₂c₃(ψ::T; tol::T=1e-6) where T <: Float64
    if ψ > tol
		sqψ = √ψ
		return (1-cos(sqψ))/ψ, (sqψ-sin(sqψ))/sqψ^3
	elseif ψ < -tol
		sqψ = √-ψ
		return (1-cosh(sqψ))/ψ, (sinh(sqψ)-sqψ)/sqψ^3
	end
    return 0.5, 1/6
end


"""
    kepler(X⃗₀::Array{T}, ∆t::T; μ::T=μ⨁, tol::T=1e-6, max_i::Int64=20) where T <: Float64

Solve Kepler's problem for all orbit types using universal variables.

Author: James Spicer (2016.11.06) <br/>
Source: Vallado (2001), Algorithm 8, p.101-2
"""
function kepler(X⃗₀::Array{T}, ∆t::T; μ::T=μ⨁, tol::T=1e-6, max_i::Int64=20) where T <: Float64
    r₀, v₀  = norm(X⃗₀[1:3]), norm(X⃗₀[4:6])
    # ξ       =  v₀^2*0.5 - μ/r₀
    α       = -v₀^2/μ   + 2/r₀

    sqμ = √μ
    if α > tol                          # Ellipse
        χ = sqμ*∆t*α
        (abs(α-1) < tol) && (χ *= 0.97) # Circle
    end

    if abs(α) < tol                     # Parabola
        h⃗ = X⃗₀[1:3] × X⃗₀[4:6]
        p = norm(h⃗)^2/μ
        s = 0.5acot(3sqμ*∆t/√(p^3))
        w = atan(∛tan(s))
        χ = 2*√p*cot(2w)
    end

    rdv = X⃗₀[1:3]⋅X⃗₀[4:6]
    if α < -tol                         # Hyperbola
        a = √-inv(α)
        χ = sign(∆t)*a*log(-2μ*α*∆t / (rdv + sign(∆t)*sqμ*a*(1-r₀*α)))
    end

    c₂, c₃, χ², χ³, k, r, i = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1
    while true
        ψ         = α*χ^2
        c₂, c₃    = c₂c₃(ψ)
        χ², χ³, k = χ^2, χ^3, χ*(1-ψ*c₃)
        r         = c₂*χ² + rdv*k/sqμ + r₀*(1-ψ*c₂)
        Δχ        = (sqμ*∆t - c₃*χ³ - rdv*c₂*χ²/sqμ - r₀*k)/r
        (abs(Δχ) ≤ tol || i > max_i) && (χ += Δχ; break)
        χ        += Δχ
        i        += 1
    end

    (f, ġ) = 1 .- c₂*χ²./[r₀, r]
    g      = ∆t - c₃*χ³/sqμ
    ḟ      = -sqμ*k/(r*r₀)

    X⃗      = similar(X⃗₀)
    X⃗[1:3] = f*X⃗₀[1:3] + g*X⃗₀[4:6]
    X⃗[4:6] = ḟ*X⃗₀[1:3] + ġ*X⃗₀[4:6]

    return X⃗
end


"""
    pKepler(X⃗₀::Array{T}, Δt::T, ṅ₀::T, n̈₀::T; R::T=R⨁, J₂::T=J₂⨁, μ::T=μ⨁, tol::T=1e-6) where T <: Float64

Propogate an orbit under the influence of secular gravitational
 and drag peturbations.

### Arguments (Inputs):
| Parameter | Description                        | Units     |
| --------: | :--------------------------------- | :-------  |
| `X⃗₀`      | ECI state ([x y z ẋ ẏ ż]) at epoch | DU, DU/TU |
| `Δt`      | Time since epoch                   | TU        |
| `ṅ₀`      | Mean-motion rate                   | rad/TU²   |
| `n̈₀`      | Mean-motion acceleration           | rad/TU³   |

#### Optional Keyword Arguments:
| Parameter | Description                          | Units   |
| --------: | :----------------------------------- | :------ |
| `R=R⨁`    | Central body mean radius             | DU      |
| `J₂=J₂⨁`  | Central body 2nd dynamic form factor | -       |
| `μ=μ⨁`    | Central body gravitational parameter | DU³/TU² |
| `tol=1e-6`| Ecc. condition for circular orbits   | -       |

### Values (Outputs):
| Parameter | Description                                | Units     |
| --------: | :----------------------------------------- | :-------- |
| `X⃗`       | ECI state ([x y z ẋ ẏ ż]) `Δt` after epoch | DU, DU/TU |

Author: James Spicer (2016.11.22) <br/>  
Source: Vallado (2001), Algorithm 60, p.646-7
"""
function pKepler(X⃗₀::Array{T}, Δt::T, ṅ₀::T, n̈₀::T; R::T=R⨁, J₂::T=J₂⨁, μ::T=μ⨁, tol::T=1e-6) where T <: Float64
    p₀, a₀, e₀, ι₀, Ω₀, ω₀, ν₀ = ELORB(X⃗₀, μ=μ)

    E₀   = abs(e₀) ≥ tol ? ν2A(ν₀, e₀) : ν₀ # ELORB returns u = ν₀ if necessary.
    M₀   = E₀ - e₀*sin(E₀)
    n₀   = √(μ/a₀^3)

    n₀Δt, ṅ₀Δt = n₀*Δt, ṅ₀*Δt
    α, β = 0.75J₂*n₀Δt*R^2/p₀^2, 2ṅ₀Δt/(3n₀)

    a    = a₀*(1 - β)
    e    = e₀ - β*(1-e₀)
    Ω    = Ω₀ - 2α*    cos(ι₀)
    ω    = ω₀ +  α*(4-5sin(ι₀)^2)
    M    = M₀ + n₀Δt + 0.5ṅ₀Δt*Δt + n̈₀*Δt^3/6
    p    = a*(1-e^2)
      
    E    = kepEqtnE(M, e)
  
    ν    = abs(e) ≥ tol ? A2ν(E, e) : E

    return RANDV(p, e, ι₀, Ω, ω, ν, μ=μ)
end


"""
    ELORB(X⃗::Array{T}; μ::T=μ⨁, tol::T=1e-6) where T <: Float64

Compute Keplerian elements given Cartesian (IJK) state vector ([r⃗ v⃗]).

Author: James Spicer (2016.11.06) <br/>
Source: Vallado (2001), Algorithm 9, p.120-121 <br/>
Special case treatment from GMAT Math Spec, p.48-49, eqns. 3.11-23.
"""
function ELORB(X⃗::Array{T}; μ::T=μ⨁, tol::T=1e-6) where T <: Float64
    r = norm(X⃗[1:3])
    v = norm(X⃗[4:6])

    h⃗ = X⃗[1:3] × X⃗[4:6]
    h = norm(h⃗)

    n⃗ = [0, 0, 1] × h⃗
    n = norm(n⃗)

    e⃗ = ((v^2-μ/r)*X⃗[1:3] - (X⃗[1:3]⋅X⃗[4:6])*X⃗[4:6])/μ
    e = norm(e⃗)
    
    ι = acos(max(min(h⃗[3]/h, 1), -1))

    circular   = e ≤ tol
    parabolic  = abs(e - 1) ≤ tol
    equatorial = ι ≤ tol

    if parabolic
        a = Inf
        p = h^2/μ
    else 
        ξ = 0.5v^2 - μ/r
        a = -0.5μ/ξ
        p = a*(1-e^2)
    end

    if circular && equatorial
        Ω = 0.0
        ω = 0.0
        ν = acos(max(min(X⃗[1]/r     , 1), -1)) # Vallado calls this 'λtrue'

        (X⃗[2] < 0) && (ν = 2π - ν)
    elseif equatorial # & elliptical
        Ω = 0.0
        ω = acos(max(min(e⃗[1]/e          , 1), -1)) # Vallado calls this 'ωtrue'
        ν = acos(max(min((e⃗⋅X⃗[1:3])/(e*r), 1), -1))

        (e⃗[2]           < 0) && (ω = 2π - ω)
        (X⃗[1:3]⋅X⃗[4:6]  < 0) && (ν = 2π - ν)
    elseif circular # & inclined
        Ω = acos(max(min(n⃗[1]/n          , 1), -1))
        ω = 0.0
        ν = acos(max(min((n⃗⋅X⃗[1:3])/(n*r), 1), -1)) # Vallado calls this 'u'

        (n⃗[2] < 0) && (Ω = 2π - Ω)
        (X⃗[3] < 0) && (ν = 2π - ν)
    else 
        Ω = acos(max(min(n⃗[1]/n          , 1), -1))
        ω = acos(max(min((e⃗⋅n⃗     )/(e*n), 1), -1))
        ν = acos(max(min((e⃗⋅X⃗[1:3])/(e*r), 1), -1))

        (n⃗[2]          < 0) && (Ω = 2π - Ω)
        (e⃗[3]          < 0) && (ω = 2π - ω)
        (X⃗[1:3]⋅X⃗[4:6] < 0) && (ν = 2π - ν)
    end

    return p, a, e, ι, Ω, ω, ν
end


# lines 1 & 2 are TLE strings, τ is propogation time [sec]
function SGPwrap(line1, line2, τ; R=R⨁, μ=μ⨁, J₂=J₂⨁, J₃=J₃⨁)
    satnum, class, intldesg, epoch, ṅ₀, n̈₀, b☆, TLEnum, i₀, Ω₀, e₀, ω₀, M₀, n₀, revnum = parseTLE(line1, line2)

    n₀  *= 2π / 1440    # [rev/min]
    ṅ₀  *= 2π / 1440^2
    n̈₀  *= 2π / 1440^3

    τ   /= 60 # [min]
    kₑ   = 60/√(R^3/μ)
    r⃗, v⃗ = SGP(τ, n₀, ṅ₀, n̈₀, e₀, i₀, M₀, ω₀, Ω₀, kₑ, J₂, J₃, 1.0) # [DU], [DU/min]

    return r⃗*R, v⃗*R/60 # [km], [km/sec]
end

# τ = t - t₀ in [min], n, ṅ, n̈ in [rev/min^x], angles in [rad].
# Source: Spacetrak Report #3 p3-6 (https://celestrak.com/NORAD/documentation/spacetrk.pdf)
function SGP(τ, n₀, ṅ₀, n̈₀, e₀, i₀, M₀, ω₀, Ω₀, kₑ, J₂, J₃, aE; max_i=20, tol=1e-8)
    ci₀, si₀, β, e₀² = cos(i₀), sin(i₀), 0.25J₂*aE^2, e₀^2
    ci₀²   = ci₀^2

    a₁     = ∛(kₑ/n₀)^2 # [DU]
    δ₁     = 3β*(3ci₀^2 - 1) / (a₁^2*√(1 - e₀²)^3)
    a₀     = a₁*(1 - δ₁/3 - δ₁^2 - 134δ₁^3/81)
    p₀     = a₀*(1 - e₀²)
    q₀     = a₀*(1 - e₀ )
    L₀     = M₀ + ω₀ + Ω₀
 
    p₀⁻²   = 1/p₀^2
    ζ      = 3β*n₀*p₀⁻²
    Ω̇      = -2ζ*  ci₀       
    ω̇      =   ζ*(5ci₀² - 1)
 
    τ²     = τ^2
    a      = a₀*∛(n₀/(n₀ + 2ṅ₀*τ + 3n̈₀*τ²))^2
    e      = (a > q₀) ? 1 - q₀/a : 1e-6
    p      = a*(1 - e^2)

    Ω̇τ, ω̇τ = Ω̇*τ, ω̇*τ
    Ωs₀    = Ω₀ + Ω̇τ
    ωs₀    = ω₀ + ω̇τ
    Ls     = L₀ + ω̇τ + Ω̇τ + n₀*τ + ṅ₀*τ² + n̈₀*τ^3
 
    θ      = 0.25J₃*aE*si₀/(J₂*p)
    ayNSL  = e*sin(ωs₀) - 2θ
    axNSL  = e*cos(ωs₀)
    L      = Ls - θ*axNSL*(3 + 5ci₀)/(1 + ci₀)

    # solve Kepler's equation for (E⁺ω)
    U      = L - Ωs₀
    i, E⁺ω, cE⁺ω, sE⁺ω  = 1, U, 0.0, 0.0
    while true
        cE⁺ω, sE⁺ω = cos(E⁺ω), sin(E⁺ω)
        ΔE⁺ω = (U - ayNSL*cE⁺ω + axNSL*sE⁺ω - E⁺ω) / (1 - ayNSL*sE⁺ω - axNSL*cE⁺ω)
        (abs(ΔE⁺ω) ≤ tol || i > max_i) && (E⁺ω += ΔE⁺ω; break)
        E⁺ω += ΔE⁺ω
        i   += 1
    end

    # Intermediate (partially-osculating) quantities
    ecosE = axNSL*cE⁺ω + ayNSL*sE⁺ω
    esinE = axNSL*sE⁺ω - ayNSL*cE⁺ω

    eL²   = axNSL^2 + ayNSL^2
    pL    = a*(1 - eL²  )
    r     = a*(1 - ecosE)
    ṙ     = kₑ*√a*esinE/r
    rv̇    = kₑ*√pL     /r

    γ     = esinE/(1 + √(1 - eL²))
    su    = sE⁺ω - ayNSL - axNSL*γ#*(a/r)
    cu    = cE⁺ω - axNSL + ayNSL*γ#*(a/r)
    u     = atan(su, cu)

    # Short-period peturbations
    βc2u, βs2u, pL⁻² = β*cos(2u), β*sin(2u), 1/pL^2
    rₖ    = r   +    βc2u/pL  *  si₀^2    
    uₖ    = u   - 0.5βs2u*pL⁻²*(7ci₀² - 1)
    Ωₖ    = Ωs₀ +   3βs2u*pL⁻²*  ci₀      
    iₖ    = i₀  +   3βc2u*pL⁻²*  ci₀*si₀

    cΩₖ, sΩₖ, cuₖ, suₖ, ciₖ = cos(Ωₖ), sin(Ωₖ), cos(uₖ), sin(uₖ), cos(iₖ)

    M     = [-sΩₖ*ciₖ;
              cΩₖ*ciₖ;
              sin(iₖ)]

    N     = [ cΩₖ;
              sΩₖ;
              0.0]

    U     = M*suₖ + N*cuₖ
    V     = M*cuₖ - N*suₖ

    r⃗     = rₖ*U        # [DU]
    v⃗     = ṙ*U + rv̇*V  # [DU/min]
    return r⃗, v⃗
end


"""
    predict(X⃗₀::Array{T}, JD₀::T, Δt::T, JD_start::T, JD_end::T, λ::T, ϕ_gd::T, h_ellp::T; ṅ::T=0, n̈::T=0, R::T=R⨁, e2::T=e2⨁, J₂::T=J₂⨁, μ::T=μ⨁, tol::T=1e-6) where T <: Float64

Predict a satellite's look angles and visibility. 

### Arguments (Inputs):
| Parameter  | Description                        | Units     |
| ---------: | :--------------------------------- | :-------- |
| `X⃗₀`       | ECI state ([x y z ẋ ẏ ż]) at epoch | DU, DU/TU |
| `JD₀`      | Epoch                              | days      |
| `Δt`       | Simulation timestep                | TU        |
| `JD_start` | Simulation start time              | days      |
| `JD_end`   | Simulation end time                | days      |
| `λ`        | Look site's longitude              | rad       |
| `ϕ_gd`     | Look site's geodetic latitude      | rad       |
| `h_ellp`   | Look site's height above ellipsoid | DU        |

#### Optional Keyword Arguments:
| Parameter | Description                          | Units   |
| --------: | :----------------------------------- | :------ |
| `ṅ=0`     | 1st derivative of mean motion        | rad/s²  |
| `n̈=0`     | 2nd derivative of mean motion        | rad/s³  |
| `R=R⨁`    | Central body mean radius             | DU      |
| `e2=e2⨁`  | Body’s eccentricity squared          | -       |
| `J₂=J₂⨁`  | Central body 2nd dynamic form factor | -       |
| `μ=μ⨁`    | Central body gravitational parameter | DU³/TU² |
| `tol=1e-6`| Ecc. condition for circular orbits   | -       |

### Values (Outputs):
| Parameter | Description                      | Units |
| --------: | :------------------------------- | :---- |
| `JDs`     | Array of simulation Julian dates | days  |
| `Vis`     | Array of satellite's visibility from the site at each time. <br/>
                `Visible`: Satellite is sunlight, site is dark.  <br/>
                `Radar Night`: Satellite and site both dark.  <br/>
                `Radar Sun`: Satellite and site both sunlit.  <br/>
                `Not visible`: Satellite below site horizon.  <br/>
| `Vecs`    | Array of look angles in the form `[az el ρ]` at each time. | rad, rad, km |


Author: James Spicer (2016.11.24) <br/>
Source: Vallado (2001), Algorithm 68, p.828
"""
function predict(X⃗₀::Array{T}, JD₀::T, Δt::T, JD_start::T, JD_end::T, λ::T, ϕ_gd::T, h_ellp::T; ṅ::T=0, n̈::T=0, R::T=R⨁, e2::T=e2⨁, J₂::T=J₂⨁, μ::T=μ⨁, tol::T=1e-6) where T <: Float64
    (Δt == 0) && @error "Δt cannot be 0."
    (Δt  > 0 && JD_end < JD_start) && @error "Impossible propogation conditions."
    (Δt  < 0 && JD_end > JD_start) && @error "Impossible propogation conditions."
    
    JDs  = Float64[]
    Vis  = String[]
    Vecs = []

    X⃗ = pKepler(X⃗₀, 86400(JD_start-JD₀), ṅ, n̈, R=R, J₂=J₂, μ=μ, tol=tol)
    JD   = JD_start
    while true
        vis, look = visibilityAngles(X⃗[1:3], JD, λ, ϕ_gd, h_ellp, R=R, e2=e2)

        push!(JDs, JD)
        push!(Vis, vis)
        Vecs = isempty(Vecs) ? look : [Vecs; look]
        
        JD  += Δt/86400
        abs(JD_start-JD) > abs(JD_start-JD_end) && break
        X⃗ = pKepler(X⃗, Δt, ṅ, n̈, R=R, J₂=J₂, μ=μ, tol=tol)
    end
    return JDs, Vis, Vecs
end


"""
    visibilityAngles(r⃗::Array{T}, JD::T, λ::T, ϕ_gd::T, h_ellp::T; R::T=R⨁, e2::T=e2⨁) where T <: Float64

Predict a satellite's look angles and visibility. 

### Arguments (Inputs):
| Parameter  | Description                         | Units |
| ---------: | :---------------------------------- | :---- |
| `r⃗`        | Satellite's current position vector | DU    |
| `JD`       | Current Julian date                 | days  |
| `λ`        | Site's longitude                    | rad   |
| `ϕ_gd`     | Site's geodetic latitude            | rad   |
| `h_ellp`   | Site's heigh above ellipsoid        | DU    |

#### Optional Keyword Arguments:
| Parameter | Description                   | Units  |
| --------: | :---------------------------- | :----- |
| `R=R⨁`    | Radius of site’s central body | DU     |
| `e2=e2⨁`  | Body’s eccentricity squared   | -      |

### Values (Outputs):
| Parameter | Description                      | Units |
| --------: | :------------------------------- | :---- |
| `vis`     | Satellite's visibility from the site: <br/>
                `Visible`: Satellite is sunlight, site is dark.  <br/>
                `Radar Night`: Satellite and site both dark.  <br/>
                `Radar Sun`: Satellite and site both sunlit.  <br/>
                `Not visible`: Satellite below site horizon.  <br/>
| `vec`     | Look angles in the form `[az el ρ]`. | rad, rad, DU |


Author: James Spicer (2016.11.24) <br/>
Source: Vallado (2001), Algorithm 68, p.828
"""
function visibilityAngles(r⃗::Array{T}, JD::T, λ::T, ϕ_gd::T, h_ellp::T; R::T=R⨁, e2::T=e2⨁) where T <: Float64
    θ_LST  = LSTime(JD, λ=λ)[2]                        # [rad]
    r⃗_site = site_eci(θ_LST, ϕ_gd, h_ellp, R=R, e2=e2) # [DU]
    ρ⃗_ijk  = r⃗' - r⃗_site'                              # [DU]
    ρ⃗_sez  = rot2mat(0.5π - ϕ_gd)*rot3mat(θ_LST)*ρ⃗_ijk # [DU]
    
    if ρ⃗_sez[3] > 0
        r⃗☉ = sun(JD)*km_in_AU # [km]

        if r⃗☉⋅r⃗_site > 0
            vis  = "Radar Sun"
        else
            ζ    = asin(norm(r⃗☉ × vec(r⃗))/(norm(r⃗☉)*norm(r⃗)))
            Dist = norm(r⃗)*cos(ζ - 0.5π)
            vis  = (Dist > R⨁) ? "Visible" : "Radar Night"
        end
        β  = mod2pi(atan(ρ⃗_sez[2], -ρ⃗_sez[1]))
        ρ  = norm(ρ⃗_sez)
        el = asin(ρ⃗_sez[3]/ρ)
        return vis, [β el ρ]
    else
        return "Not visible", [NaN NaN NaN]
    end
end


vCirc(μ, r) = √(μ/r)

deltaVcirc(μ, r₁, r₂) = √(μ*norm(inv(r₁) - inv(r₂)))

orbitalEnergy(μ, r, v) = 0.5v^2 - μ/r

ellipseRadius(e, a, ν) = a*(1-e^2) / (1+e*cos(ν))

ellipseSpeed(μ, a, r) = √(μ*(2/r - inv(a)))

period(μ, a) = 2π*√((a^3)/μ) # [sec]

semimajorAxisT(T, μ) = ∛(μ*(0.5T/π)^2) # T in [s], μ in [km3/s2]

semimajorAxisN(n, μ) = ∛(μ/(n^2)) # n in [rad/s], μ in [km3/s2]

meanMotionA(a, μ) = √(μ/(a^3)) # [rad/s]

meanMotionP(p, μ, e) = √(μ*((1-e^2)/p)^3) # [rad/s]

deltaV(ve, mi, mf) = ve*log(mi/mf) # Tsiolkovsky rocket equation: ve = effective exhaust velocity, mi/f = initial, final mass

ve(Isp) = Isp*g₀ # Effective exhaust velocity = Specific impulse * 9.81 m/s

deltaVEllipse(r₀, rp, ra, lat, i) = 2μ⨁*(r₀-rp+ra) - v0(lat, i) # Delta-V to achieve elliptical orbit. openAerospace.org

v0(ϕ, ι) = 2π*ω⨁*R⨁*cos(ϕ)*(cos(ι-ϕ)-sin(ι-ϕ)) # Velocity due to the rotation of the earth at latitude ϕ and inclination ι.

# trueAnomaly(E, e) = 2*atan(sqrt(1+e)*sin(E/2), sqrt(1-e)*cos(E/2)) # nu, E in radians

trueAnomaly(E::Float64, e::Float64) = mod2pi(atan(√(1-e^2)*sin(E), cos(E)-e)) # ν, E, in [rad]. Source: Vallado

meanAnomaly(E::Float64, e::Float64) = mod2pi(E - e*sin(E)) # M, E in [rad]

function eccentricAnomalyM(M::T, e::T; tol::T=1e-10, max_i::Int64=20) where T <: Float64 # M, E in [rad]
    E, i = M + e, 1
    while true
        ΔE = (M + e*sin(E) - E)/(1 - e*cos(E))
        (abs(ΔE) ≤ tol || i > max_i) && (E += ΔE; break)
        E += ΔE
        i += 1
    end
    return mod2pi(E)
end

eccentricAnomalyNu(ν::Real, e::Real) = mod2pi(atan(√(1-e^2)*sin(ν), e+cos(ν))) # E, ν in [rad]. Source: Vallado

semilatusRectum(a::Real, e::Real) = a*(1 - e^2)


"""
    c2e(X::Array{Float64})

Transforms classical into equinoctial orbital elements.
Inputs must be in order a, e, i, Ω, ω, M, row vector if array.
Outputs are in order a, h, k, p, q, λ.
Author: James Spicer (2017.03.15) <br/>
Source: Danielson (1995), Eq 2.1.2 (1-2), p.6
"""
function c2e(X::Array{Float64})
    j = 1#(0 ≤ i ≤ 0.5π) ? 1 : -1 # Retrograde factor
    return [X[1]
            X[2]*sin(X[5] + j*X[4])
            X[2]*cos(X[5] + j*X[4])
            tan(0.5X[3])^j*sin(X[4])
            tan(0.5X[3])^j*cos(X[4])
            mod2pi(X[6] + X[5] + j*X[4])]'
end

function c2e(a::T, e::T, i::T, Ω::T, ω::T, M::T) where T <: Float64
    X = c2e([a e i Ω ω M])
    return X[1], X[2], X[3], X[4], X[5], X[6]
end


"""
    e2c(a::T, h::T, k::T, p::T, q::T, λ::T) where T <: Float64
    e2c(X::Array{Float64})

Transforms equinoctial into classical orbital elements.
Inputs must be in order a, h, k, p, q, λ, row vector if array.
Outputs are in order a, e, i, Ω, ω, M.
Author: James Spicer (2017.03.15) <br/>
Source: Danielson (1995), Eq 2.1.3 (1-2), p.7
"""
function e2c(X::Array{Float64})
    j = 1#(0 ≤ i ≤ 0.5π) ? 1 : -1 # Retrograde factor, how to calculate?
    e = hypot(X[2]  , X[3]  )
    ζ = atan(X[2]/e, X[3]/e)
    ϕ = hypot(X[4]  , X[5]  )
    Ω = atan(X[4]/ϕ, X[5]/ϕ)

    return [X[1]
            e 
            0.5π*(1-j) + 2j*atan(ϕ)
            mod2pi(Ω)
            mod2pi(ζ - j*Ω)
            mod2pi(X[6] - ζ)       ]'
end

function e2c(a::T, h::T, k::T, p::T, q::T, λ::T) where T <: Float64
    X = e2c([a h k p q λ])
    return X[1], X[2], X[3], X[4], X[5], X[6]
end


"""
    e2rv(a::T, h::T, k::T, p::T, q::T, λ::T; μ::T=μ⨁, tol::T=1e-10, max_i::Int64=20) where T <: Float64
    e2rv(X::Array{T}; μ::T=μ⨁, tol::T=1e-10, max_i::Int64=20) where T <: Float64

Transforms equinoctial orbital elements into Cartesian state.
Inputs must be in order a, h, k, p, q, λ (row vector if array).
Author: James Spicer (2017.03.15) <br/>
Source: GMAT Math. Spec. (2015), Section 3.1.4, p.51-52
"""
function e2rv(X::Array{T}; μ::T=μ⨁, tol::T=1e-10, max_i::Int64=20) where T <: Float64
    F, i, cF, sF = X[6], 1, 0.0, 0.0
    while true
        cF, sF = cos(F), sin(F)
        ΔF = (F + X[2]*cF - X[3]*sF - X[6])/(1 - X[2]*sF - X[3]*cF)
        (abs(ΔF) ≤ tol || i > max_i) && (F -= ΔF; break)
        F -= ΔF
        i += 1
    end

    X₂², X₃², X₄², X₅², X₄X₅ = X[2]^2, X[3]^2, X[4]^2, X[5]^2, X[4]*X[5]

    β          = inv(1 + √(1-X₂²-X₃²))
    a, b, c, d = √(μ/X[1])/(1 - X[3]*cF - X[2]*sF), X[2]*X[3]*β, 1-X₂²*β, 1-X₃²*β
    bcF, bsF   = b*cF, b*sF
    
    X₁         = X[1]*( c*cF + bsF - X[3])
    Y₁         = X[1]*( d*sF + bcF - X[2])
        
    Ẋ₁         =    a*(-c*sF + bcF       )
    Ẏ₁         =    a*( d*cF - bsF       )
        
    j          = 1 # Retrograde factor, how to calculate?
    ϕ          = inv(1+X₄²+X₅²)
    Q          = ϕ*[1-X₄²+X₅²  2j*X₄X₅               2X[4];
                    2X₄X₅      j*(1+X₄²-X₅²)        -2X[5];
                      -2j*X[4] 2X[5]         j*(1-X₄²-X₅²)]
        
    f̂          = Q[:,1]'
    ĝ          = Q[:,2]'

    return [X₁*f̂+Y₁*ĝ Ẋ₁*f̂+Ẏ₁*ĝ]
end

function e2rv(a::T, h::T, k::T, p::T, q::T, λ::T; μ::T=μ⨁, tol::T=1e-10, max_i::Int64=20) where T <: Float64
    X = e2rv([a h k p q λ], μ=μ, tol=tol, max_i=max_i)
    return X[1:3], X[4:6]
end


"""
    rv2e(r⃗::Array{T}, v⃗::Array{T}; μ::T=μ⨁) where T <: Float64
    rv2e(X::Array{T}; μ::T=μ⨁) where T <: Float64

Transforms Cartesian state into equinoctial orbital elements.
Inputs must be row vector(s).
Output order is a, h, k, p, q, λ.
Author: James Spicer (2017.03.15) <br/>
Source: GMAT Math. Spec. (2015), Section 3.1.5, p.52-53
"""
function rv2e(X::Array{T}; μ::T=μ⨁) where T <: Float64
    r⃗, v⃗   = X[1:3], X[4:6]
    r, v   = norm(r⃗), norm(v⃗)
    γ, v²  = μ/r, v^2
  
    e⃗      = ((v² - γ)*r⃗ - (r⃗⋅v⃗)*v⃗)/μ
    ξ      = 0.5v² - γ 
    a      = -0.5μ/ξ
    ĥ      = normalize(r⃗ × v⃗)
  
    j      = 1 # Retrograde factor
    m      = inv(1+ĥ[3]*j)
    p, q   = m*ĥ[1], -m*ĥ[2]
     
    fx     = 1-p*ĥ[1]
    fy     =   q*ĥ[1]
    fz     =  -j*ĥ[1]
    f̂      = [fx, fy, fz]
    
    ĝ      = ĥ × f̂
    h, k   = e⃗⋅ĝ, e⃗⋅f̂
    h², k² = h^2, k^2
  
    X₁, Y₁ = r⃗⋅f̂, r⃗⋅ĝ

    b      = √(1-h²-k²)
    β      = 1/(1 + b)
    c, d   = h*k*β, inv(a*b)
  
    cF     = k + d*((1-k²*β)*X₁ - c*Y₁)
    sF     = h + d*((1-h²*β)*Y₁ - c*X₁)
    F      = atan(sF, cF)
    λ      = F + h*cF - k*sF
    
    return [a h k p q mod2pi(λ)]
end

function rv2e(r⃗::Array{T}, v⃗::Array{T}; μ::T=μ⨁) where T <: Float64
    X = rv2e([r⃗ v⃗])
    return X[1], X[2], X[3], X[4], X[5], X[6]
end


"""
    c2me(a::T, e::T, i::T, Ω::T, ω::T, M::T) where T <: Float64
    c2me(X::Array{Float64})

Transforms classical into modified equinoctial orbital elements.
Inputs must be in order a, e, i, Ω, ω, ν, row vector if array.
Outputs are in order p, f, g, h, k, L (ESA notation).
Author: James Spicer (2017.03.16) <br/>
Source: Walker (1985), Eqns. (2), p.2.
"""
c2me(X::Array{Float64}) = [
    X[1]*(1-X[2]^2)     
    X[2]*cos(X[4]+X[5])
    X[2]*sin(X[4]+X[5])
    tan(0.5X[3])*cos(X[4])
    tan(0.5X[3])*sin(X[4])
    mod2pi(X[4]+X[5]+X[6])
]'

function c2me(a::T, e::T, i::T, Ω::T, ω::T, ν::T) where T <: Float64
    X = c2me([a e i Ω ω ν])
    return X[1], X[2], X[3], X[4], X[5], X[6]
end


"""
    me2c(p::T, f::T, g::T, h::T, k::T, L::T) where T <: Float64
    me2c(X::Array{Float64})

Transforms modified equinoctial into classical orbital elements.
Inputs must be in order p, f, g, h, k, L, row vector if array.
Outputs are in order a, e, i, Ω, ω, ν.
Author: James Spicer (2017.03.16) <br/>
Source: Walker (1985)
"""
function me2c(X::Array{Float64})
    e = hypot(X[2]  , X[3]  )
    ζ = atan(X[3]/e, X[2]/e)
    ϕ = hypot(X[4]  , X[5]  )
    Ω = atan(X[5]/ϕ, X[4]/ϕ)

    return [X[1]*inv(1-X[2]^2-X[3]^2)
            e 
            2atan(ϕ)
            mod2pi(Ω)
            mod2pi(ζ - Ω)
            mod2pi(X[6] - ζ)        ]'
end

function me2c(p::T, f::T, g::T, h::T, k::T, L::T) where T <: Float64
    X = me2c([p f g h k L])
    return X[1], X[2], X[3], X[4], X[5], X[6]
end


"""
    me2rv(p::T, f::T, g::T, h::T, k::T, L::T; μ::T=μ⨁) where T <: Float64
    me2rv(X::Array{T}; μ::T=μ⨁) where T <: Float64

Transforms modified equinoctial orbital elements into Cartesian state.
Inputs must be in order p, f, g, h, k, L (row vector if array).
Author: James Spicer (2017.03.16) <br/>
Source: Low Thrust Transfer Between Non-coplanar Circular Orbits, p.6-7.
"""
function me2rv(X::Array{T}; μ::T=μ⨁) where T <: Float64
    sL, cL = sin(X[6]), cos(X[6])
    h², k² = X[4]^2, X[5]^2
    α²     = h² - k²
    p      = 2X[4]*X[5]

    A      = inv(1 + h² + k²)
    B      = A*X[1]/(1 + X[2]*cL + X[3]*sL)
    C      = A*√(μ/X[1])

    return [ B*(       cL*(1+α²) + p*sL       )
             B*(       sL*(1-α²) + p*cL       )
            2B*(       sL*X[4]   - X[5]*cL    )
            -C*((sL+X[3])*(α²+1) - p*(cL+X[2]))
            -C*((cL+X[2])*(α²-1) + p*(sL-X[3]))
            2C*( cL*X[4]      + X[5]*(sL+X[3]))]'
end

function me2rv(p::T, f::T, g::T, h::T, k::T, L::T; μ::T=μ⨁) where T <: Float64
    X = me2rv([p f g h k L], μ=μ)
    return X[1:3], X[4:6]
end


"""
    eqMo_ME(Y::Array{T}; μ::T=μ⨁, P::Array{T}=zeros(1,3)) where T <: Float64
    eqMo_ME(p::T, f::T, g::T, h::T, k::T, L::T; μ::T=μ⨁, P::Array{T}=zeros(1,3)) where T <: Float64

Returns equinotial velocity Ẏ given an equinoctial state Y.
Inputs must be in order p, f, g, h, k, L (row vector if array).
Column array P is peturbing vector in RSW (radial, tangential, normal) frame.
Author: James Spicer (2017.03.16) <br/>
Source: Low Thrust Transfer Between Non-coplanar Circular Orbits, p.7-8.
"""
function eqMo_ME(Y::Array{T}; μ::T=μ⨁, P::Array{T}=zeros(3,1)) where T <: Float64
    sL, cL = sin(Y[6]), cos(Y[6])
    w      = 1 + Y[2]*cL + Y[3]*sL
    b      = [0 0 0 0 0 √μ*w^2/√Y[1]^3]
    P == zeros(3,1) && return b

    p, q = Y[4]*sL - Y[5]*cL, √(Y[1]/μ)
    s²   = 0.5(1 + Y[4]^2 + Y[5]^2)
    A    = q/w*[ 0     2Y[1]            0
                 sL*w  ((w+1)*cL+Y[2]) -Y[3]*p;
                -cL*w  ((w+1)*sL+Y[3])  Y[2]*p; 
                 0     0                 s²*cL;
                 0     0                 s²*sL;
                 0     0                     p]

    return (A*P + b')'
end

function eqMo_ME(p::T, f::T, g::T, h::T, k::T, L::T; μ::T=μ⨁, P::Array{T}=zeros(1,3)) where T <: Float64
    X = eqMo_ME([p f g h k L], μ=μ, P=P)
    return X[1], X[2], X[3], X[4], X[5], X[6]
end