export
	hohmann,
    bielliptic,
    incOnlyCirc,
    lanOnly,
    incLan


# For coplanar, circular orbits.
function hohmann(r_i::T, r_f::T; μ::T=μ⨁) where T <: Float64
    a_t   = 0.5(r_i + r_f)

    v_i   = √(μ/r_i)
    v_t_a = √(2μ/r_i - μ/a_t)
    v_f   = √(μ/r_f)
    v_t_b = √(2μ/r_f - μ/a_t)

    Δv_a  = v_t_a - v_i
    Δv_b  = v_f - v_t_b

    Δv    = abs(Δv_a) + abs(Δv_b)
    τ     = π*√(a_t^3/μ)

    return Δv, τ
end


# For coplanar, circular orbits
function bielliptic(r_i::T, r_f::T, r_b::T; μ::T=μ⨁) where T <: Float64
    a_t_1  = 0.5(r_i + r_b)
    a_t_2  = 0.5(r_b + r_f)

    v_i    = √(μ/r_i)
    v_t_1a = √(2μ/r_i - μ/a_t_1)
    v_t_1b = √(2μ/r_b - μ/a_t_1)
    v_t_2b = √(2μ/r_b - μ/a_t_2)
    v_t_2c = √(2μ/r_f - μ/a_t_2)

    Δv_a   = v_t_1a - v_i
    Δv_b   = v_t_2b - v_t_1b
    Δv_c   = v_f - v_t_2c

    Δv     = abs(Δv_a) + abs(Δv_b) + abs(Δv_c)
    τ      = π*√(a_t_1^3/μ) + π*√(a_t_2^3/μ)

    return Δv, τ
end


# One-Tangent Burn


# For circular orbits
function incOnlyCirc(Δι::T, r::T; μ::T=μ⨁) where T <: Float64
    v  = √(μ/r)
    Δv = 2v*sin(0.5Δι)
    return Δv
end


# For elliptic orbits


# For circular orbits
function lanOnly(ΔΩ::T, ι_i::T, r::T; μ::T=μ⨁) where T <: Float64
    v  = √(μ/r)
    θ  = acos(1 + sin(ι_i)^2*(cos(ΔΩ)-1))
    Δv = 2v*sin(0.5θ)
end


# For circular orbits
function incLan(ι_i::T, ι_f::T, ΔΩ::T, r::T; μ::T=μ⨁) where T <: Float64
    v  = √(μ/r)
    θ  = acos(cos(ι_i)cos(ι_f) + sin(ι_i)sin(ι_f)cos(ΔΩ))
    Δv = 2v*sin(0.5θ)
end