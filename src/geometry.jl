export
	linePlaneIntersectC,
    circleLineIntersectC,
    circleCircleIntersectC,
    lineSphereIntersectV,
    ellipsoidLineIntersectV,
    lineTriangleIntersectV,
    lineSegmentIntersectC,
    PIPoly,
    perimPoly,
    areaPoly,
    PICone,
    areaUnderCurve,
    lineLineIntersectC,
    pointOnPlaneV,
    circleOverlapArea,
    shadow,
    sight,
    slantRange,
    angularDiameter,
    angularRadius,
    hav,
    ahav,
    haversine,
    hexagonal,
    centeredHexagonal,
    incenter,
    polygonPatternAngle,
    polygonPatternEdge,
    minRadius,
    VsphericalSector,
    circumRadius,
    inradius,
    sideLength,
    vectorAngle,
    vectorDistance,
    vEllipsoid,
    vSphere,
    footprint,
    limb,
    normalize,
    atan3

import LinearAlgebra.normalize, LinearAlgebra.cross

normalize(A::Array{Float64,2}) = A/norm(A)
normalize(A::Array{Int64  ,2}) = A/norm(A)


cross(u::Array{T}, v::Array{T}) where T <: Float64 = [u[2]v[3] - u[3]v[2] u[3]v[1] - u[1]v[3] u[1]v[2] - u[2]v[1]]
cross(u::Array{T}, v::Array{T}) where T <: Int64   = [u[2]v[3] - u[3]v[2] u[3]v[1] - u[1]v[3] u[1]v[2] - u[2]v[1]]
cross(u::Array{Float64}, v::Array{Int64})          = [u[2]v[3] - u[3]v[2] u[3]v[1] - u[1]v[3] u[1]v[2] - u[2]v[1]]
cross(u::Array{Int64}, v::Array{Float64})          = [u[2]v[3] - u[3]v[2] u[3]v[1] - u[1]v[3] u[1]v[2] - u[2]v[1]]
cross(u::Array{T}, v::Array{T}) where T <: Real    = [u[2]v[3] - u[3]v[2] u[3]v[1] - u[1]v[3] u[1]v[2] - u[2]v[1]]


### CARTESIAN GEOMETRY

# Intersection of a plane defined by p1 = (x1, y1, z1), p2, p3 with line passing through points p4 & p5. Source: Wikipedia
function linePlaneIntersectC(p1, p2, p3, p4, p5)
    p1, p2, p3, p4, p5 = p1[:], p2[:], p3[:], p4[:], p5[:] # Converts all to column vectors
    return p4 + (p5-p4)*(inv([p4-p5 p2-p1 p3-p1])*(p4-p1))[1]
end

# Returns point(s) of intersection (if any) of an infinite line (determined by p1 and p2) with circle of radius r and center (0,0).
# Source: http://mathworld.wolfram.com/Circle-LineIntersection.html
function circleLineIntersectC(p₁, p₂, R::Real)
    (dx, dy) = p₂ - p₁
    dr = hypot(dx, dy)
    D  = p₁[1]*p₂[2] - p₁[2]*p₂[1]

    ∆  = R^2*dr^2 - D^2
    (∆ <  0.0) && return false
    (∆ == 0.0) && return D*[dy -dx]/dr^2
    x1 =  D*dy + sign(dy)*dx*√∆
    y1 = -D*dx +  abs(dy)   *√∆
    x2 =  D*dy - sign(dy)*dx*√∆
    y2 = -D*dx -  abs(dy)   *√∆
    return [x1 y1]/dr^2, [x2 y2]/dr^2
end

# Returns points of intersection of circles with centers pᵢ and radii Rᵢ.
function circleCircleIntersectC(p₁, R₁::Real, p₂, R₂::Real)
    d = norm(p₁-p₂)
    (d >     R₁+R₂)         && return false         # Circles don't overlap
    (d < abs(R₁-R₂))        && return false         # Once circle contained within the other
    (d == 0.0 && R₁ == R₂)  && return false         # Circles are coincident
    (d ==    R₁+R₂)         && return 0.5(p₁+p₂)    # Single intersection
    (δx, δy) = p₂ - p₁
    l  = 0.5(R₁^2 - R₂^2 + d^2)*inv(d)
    hd = √(R₁^2 - l^2)/d
    ld = l/d
    return (ld*[δx δy] + hd*[δy -δx]) + p₁, (ld*[δx δy] + hd*[-δy δx]) + p₁
end


"""
    PIPoly(X::Array{T,2}, P::Array{T}) where T <: Float64

Returns `true` if test point is inside polygon, `false` otherwise. Inaccurate for
test points on the polygon boundary (edge).

### Arguments (Inputs):
| Parameter | Description                                                | Units |
| --------: | :--------------------------------------------------------- | :---- |
| `X`       | N × 2 array of polygon coordinates `[x₁ y₁; ...; xₙ yₙ]` | DU    |
| `P`       | 1 × 2 test coordinate `[xp yp]`                            | DU    |

Author: James Spicer (2018.02.28) <br/>
Source: https://wrf.ecse.rpi.edu//Research/Short_Notes/pnpoly.html  
"""
function PIPoly(X::Array{T,2}, P::Array{T}) where T <: Float64
    c, n = false, size(X,1)
    j = n
    for i = 1:n
        if ((X[i,2] > P[2]) ≠ (X[j,2] > P[2])) && (P[1] < (P[2]-X[i,2])*(X[j,1]-X[i,1])/(X[j,2]-X[i,2]) + X[i,1])
            c = !c
        end
        j = i
    end
    return c
end


"""
    perimPoly(X::Array{Float64,2}) 

Returns perimeter of closed polygon defined by coordinate array `X`.

### Arguments (Inputs):
| Parameter | Description                                                       | Units |
| --------: | :---------------------------------------------------------------- | :---- |
| `X`       | N × 2 array of closed polygon coordinates `[x₁ y₁; ...; xₙ yₙ]` | DU    |

### Value (Output)
| Parameter | Description       | Units |
| --------: | :---------------- | :---- |
| `P`       | Polygon perimeter | DU    |

Author: James Spicer (2018.02.28) <br/>
"""
function perimPoly(X::Array{Float64,2})
    P, n = 0.0, size(X,1)
    j = n
    for i = 1:n
        P += norm(X[i,:] - X[j,:])
        j = i
    end
    return P 
end


"""
    areaPoly(X::Array{Float64,2}) 

Returns perimeter of polygon defined by coordinate array `X`.

### Arguments (Inputs):
| Parameter | Description                                                | Units |
| --------: | :--------------------------------------------------------- | :---- |
| `X`       | N × 2 array of polygon coordinates `[x₁ y₁; ...; xₙ yₙ]` | DU    |

### Value (Output)
| Parameter | Description       | Units |
| --------: | :---------------- | :---- |
| `A`       | Polygon area      | DU²   |

Author: James Spicer (2018.02.28) <br/>
"""
function areaPoly(X::Array{Float64,2})
    A, n = 0.0, size(X,1)
    j = n
    for i = 1:n
        A += (X[j,1] + X[i,1])*(X[j,2] - X[i,2])
        j = i
    end
    return 0.5abs(A)
end


"""
    PIPoly{T <: Float64}(X::Array{T,2}, P::Array{T}) 

Returns `true` if test point is inside polygon, `false` otherwise. Inaccurate for
test points on the polygon boundary (edge).

### Arguments (Inputs):
| Parameter | Description                                           | Units |
| --------: | :---------------------------------------------------- | :---- |
| `P⃗₀`      | Coordinates of cone tip `[x₀ y₀ z₀]`                  | DU    |
| `n̂`       | Cone's axis direction unit vector (tip towards base). | -     |
| `secθ`    | sec() of cone half-angle `θ`.                         | rad   |
| `p⃗`       | Point to test `[xₚ yₚ zₚ]`                         | DU    |

#### Optional Keyword Arguments:
| Parameter | Description                          | Units |
| --------: | :----------------------------------- | :---- |
| `h=Inf`   | Cone height (max. tip-base distance) | DU    |

Author: James Spicer (2019.03.18)
"""
function PICone(P⃗₀::Array{T}, n̂::Array{T}, secθ::T, p⃗::Array{T}; h::T=Inf) where T <: Float64
    @assert norm(n̂) ≈ 1.0 "`n̂` must be a unit direction vector."
    d = n̂⋅(p⃗ - P⃗₀) # distance from tip to point on axis closest to test point.
    !(0 ≤ d ≤ h) && return false # if point is behind the tip or beyond the cone end.
    return norm(p⃗ - P⃗₀) < d*secθ # compare shortest distance from p to axis with cone radius at distance 'd' from the tip. Simplify.
end


# Uses trapezium method to find area under curve defined by arrays X & Y.
function areaUnderCurve(X, Y)
    w = (maximum(X)-minimum(X))/(length(X)-1)
    return 0.5w*(Y[1] + 2sum(Y[2:end-1])+ Y[end])
end


# finds point of intersection of lines y = ax + c and y = bx + d
function lineLineIntersectC(a, c, b, d)
    if a == b
        c ≠ d && @error "Lines do not intersect."
        @error "Lines are identical."
    end
    return (d-c)/(a-b), (a*d-b*c)/(a-b)
end


"""
    lineSegmentIntersectC(x₁::T,y₁::T,x₂::T,y₂::T,x₃::T,y₃::T,x₄::T,y₄::T) where T <: Float64

Returns `true` if line segment A (defined by (x₁,y₁) and (x₂,y₂)) and line segment B 
(defined by (x₃,y₃) and (x₄,y₄)) intersect.

### Arguments (Inputs):
| Parameter | Description                                | Units |
| --------: | :----------------------------------------- | :---- |
| `x₁`      | Line segment A starting point x-coordinate | DU    |
| `y₁`      | Line segment A starting point y-coordinate | DU    |
| etc.      | |

* Source: http://bryceboe.com/2006/10/23/line-segment-intersection-algorithm/
* Source: http://jeffe.cs.illinois.edu/teaching/373/notes/x06-sweepline.pdf
* Source: http://jeffe.cs.illinois.edu/teaching/373/notes/x05-convexhull.pdf
Author: James Spicer
"""
ccw(x₁::T,y₁::T,x₂::T,y₂::T,x₃::T,y₃::T) where T <: Float64 = (y₃-y₁)*(x₂-x₁) > (y₂-y₁)*(x₃-x₁)
lineSegmentIntersectC(x₁::T,y₁::T,x₂::T,y₂::T,x₃::T,y₃::T,x₄::T,y₄::T) where T <: Float64 = ccw(x₁,y₁,x₃,y₃,x₄,y₄) ≠ ccw(x₂,y₂,x₃,y₃,x₄,y₄) && ccw(x₁,y₁,x₂,y₂,x₃,y₃) ≠ ccw(x₁,y₁,x₂,y₂,x₄,y₄)



## VECTOR GEOMETRY
# Planes defined by (P - P₀)·n⃗ = 0, where P₀ is a point on the plane and n⃗ is a normal vector.
# Lines defined by P = L₀ + dl⃗, where L₀ is a point on the line and l⃗ is a parallel vector.
# Spheres defined by radius r centered at C.
# Ellipsoids defined by major axes half-lengths a, b, c, centered at C.
# Triangles defined by vertices V₀, V₁, V₂
# Points are column vectors, e.g. [1; 2; 3]


# Returns a random point on the plane.
function pointOnPlaneV(P₀, n⃗)
    v⃗ = rand(3)
    return [v⃗[1:2]; P₀[3] - inv(n⃗[3])*(n⃗[1:2]⋅(v⃗[1:2] - P₀[1:2]))]
end


# Returns values of d.
# Source: https://www.cs.oberlin.edu/~bob/cs357.08/VectorGeometry/VectorGeometry.pdf
function lineSphereIntersectV(L₀::Array{T}, l⃗::Array{T}; r::T=R⨁, C::Array{T}=zeros(3)) where T <: Float64
    A = C[:] - L₀[:]
    return quadratic(norm(l⃗)^2, -2l⃗⋅A, norm(A)^2 - r^2)
end


# Returns values of d.
# Source: https://www.cs.oberlin.edu/~bob/cs357.08/VectorGeometry/VectorGeometry.pdf
function ellipsoidLineIntersectV(L₀, l⃗; a::Float64=R⨁_a, b::Float64=R⨁_a, c::Float64=R⨁_b, C=zeros(3))
    L₀, l⃗, C = reshape(L₀, (1, 3)), reshape(l⃗, (1, 3)), reshape(C, (1, 3)) # Convert to row matrices
    M        = diagm(1 ./[a, b, c])
    l⃗₁       = l⃗*M
    L₁       = L₀*M - C*M
    return quadratic(norm(l⃗₁)^2, 2vecdot(L₁, l⃗₁), norm(L₁)^2 - 1)
end


function linePlaneIntersectV(P₀, n⃗, L₀, l⃗)
    P₀, n⃗, L₀, l⃗ = P₀[:], n⃗[:], L₀[:], l⃗[:] # Convert all to column vectors
    if l⃗⋅n⃗ == 0.0
        ((P₀-L₀)⋅n⃗ == 0.0) && error("No p.o.i as line is contained in the plane.")
        error("No p.o.i as line and plane are parallel.")
    end
    return L₀ + l⃗*(P₀-L₀)⋅n⃗ / (l⃗⋅n⃗)
end


"""
    lineTriangleIntersectV(V₀::Array{F}, V₁::Array{F}, V₂::Array{F}, L₀::Array{F}, l⃗::Array{F}; ε::F=1e-6, culling::Bool=false) where F <: Float64

Returns value of d for intersection of ray L₀ + dl⃗ with triangle described by vertices V₀, V₁, V₂. NaN if no intersection.
If culling is true, rays intersecting triangles from behind are discarded.
Source: Möller, Trumbore (1997) https://www.cs.virginia.edu/~gfx/Courses/2003/ImageSynthesis/papers/Acceleration/Fast%20MinimumStorage%20RayTriangle%20Intersection.pdf. 
"""
function lineTriangleIntersectV(V₀::Array{F}, V₁::Array{F}, V₂::Array{F}, L₀::Array{F}, l⃗::Array{F}; ε::F=1e-6, culling::Bool=false) where F <: Float64
    E₁  = V₁ - V₀
    E₂  = V₂ - V₀

    P   = l⃗ × E₂
    D   = P⋅E₁
    
    if culling
        (D < ε)        && return NaN # Ray and triangle are parallel
        
        # Calculate u, v, t parameters and test bounds
        T   = L₀ - V₀
        u   = P⋅T
        !(0 ≤ u ≤ D  ) && return NaN

        Q   = T × E₁
        v   = Q⋅l⃗
        !(0 ≤ v ≤ D-u) && return NaN
        
        return (Q⋅E₂)/D
        
    else # no culling
        
        (abs(D) < ε)   && return NaN # Ray and triangle are parallel
        D⁻¹ = 1/D

        # Calculate u, v, t parameters and test bounds
        T   = L₀ - V₀
        u   = D⁻¹*(P⋅T)
        !(0 ≤ u ≤ 1  ) && return NaN

        Q   = T × E₁
        v   = D⁻¹*(Q⋅l⃗)
        !(0 ≤ v ≤ 1-u) && return NaN
        
        return D⁻¹*(Q⋅E₂)       
    end  
end


A(R′::Real, d′::Real) = R′^2*acos(d′/R′) - d′*√(R′^2-d′^2)
"Compute overlap area of two circles of radii `R₁` and `R₂`, with centers separated by distance `d`. Source: Wolfram"
function circleOverlapArea(R₁::Real, R₂::Real, d::Real)
    (d < abs(R₁-R₂)) && return π*min(R₁, R₂)^2  # One circle inside other
    (d ≥     R₁+R₂ ) && return 0.0              # Circles don't overlap
    (d₁, d₂) = 0.5inv(d)*(d^2 + [-1.0, 1.0]*(R₂^2 - R₁^2))
    return A(R₁, d₁) + A(R₂, d₂)
end


# r⃗☉, r⃗ in eci or ecl [km] Source: Vallado (n.b. errors in early editions)
# (αᵤ, αₚ) = atan((R☉ + R⨁*[-1.0, 1.0])/km_in_AU)  # [rad] Umbra & penumbra angles
# function shadow(r⃗☉::Array{Float64}, r⃗::Array{Float64}) # maybe change output to 0.0, 0.5, 1.0
#     (r⃗☉⋅r⃗ ≥ 0.0)          && return "none"
#     ζ = vectorAngle(-r⃗☉, r⃗) # [rad]
#     (hₛ, vₛ) = norm(r⃗)*[cos(ζ), sin(ζ)]     # Satellite horizon & vertex
#     vₛ ≤ (R⨁ - tan(αᵤ)*hₛ)  && return "umbra"    # Umbra vertex
#     vₛ ≤ (R⨁ + tan(αₚ)*hₛ)  && return "penumbra" # Penumbra vertex
#     return "none"
# end

"""
    shadow(r⃗☆::Array{T}, r⃗::Array{T}; r⃗●::Array{T}=zeros(3), re●::T=R⨁, re☆::T=R☉) where T <: Float64 -> Float64

Compute the fraction of a star visible from an observer with a blocking body.

### Arguments (Inputs)
Required Arguments
  * `r⃗☆::Array{Float64}` : Location of star.
  * `r⃗::Array{Float64}`  : Location of observer (usually orbiting satellite).

Optional Keyword Arguments
  * `r⃗●::Array{Float64}` : Location of blocking body (defaults to origin).
  * `re●::Real`          : Radius of blocking body (defaults to R⨁ [km]).
  * `re☆::Real`          : Radius of star (defaults to R☉ [km]).

### Value (Output)
Fraction of star visible `[0.0, 1.0]`. <br/>
  * `1.0` implies full illumination (no shadow).<br/>
  * `0.0` implies umbra (complete shadow).<br/>
  * Anything else implies penumbra or antumbra (incomplete shadow).

!!! note
    `r⃗☆`, `r⃗`, `r⃗●`, `re●`, and `re☆` must be in the same distance units (ER, km, etc.)<br/>
    `r⃗☆`, `r⃗`, and `r⃗●` must be in the same reference frame.

Author: James Spicer (2016.10.28)
TODO: Add generalized penumbra & umbra conditions so that circle overlap bit only happens if r⃗ is definitely in penumbra.
"""
function shadow(r⃗☆::Array{T}, r⃗::Array{T}; r⃗●::Array{T}=zeros(3), re●::T=R⨁, re☆::T=R☉) where T <: Float64
#     if the satellite is closer to the star than the body is, you're good.
    norm(r⃗-r⃗☆) < norm(r⃗●-r⃗☆) && return 1.0

#     get angular diameter (radius) and spherical coordinates as viewed from the satellite of both body and star.
    re●′ = angularRadius(re●, norm(r⃗-r⃗●))
    r⃗●′  = αδr2ijk(r⃗●-r⃗, inv=true)

    re☆′ = angularRadius(re☆, norm(r⃗-r⃗☆))
    r⃗☆′  = αδr2ijk(r⃗☆-r⃗, inv=true)

    d    = norm(r⃗☆′[1:2]-r⃗●′[1:2])

    A☆′  = π*re☆′^2
    A◐′  = circleOverlapArea(re●′, re☆′, d)

    return 1 - A◐′/A☆′
end


## Could do this using shortest distance between line and point?
"""
    sight(r⃗₁::Array{T}, r⃗₂::Array{T}; R::T=R⨁, bodyShape::String="e", f::T=f⨁) where T <: Float64

Check two points have line of sight around blocking body at origin.

### Arguments (Inputs):
| Parameter | Description | Units |
| --------: | :---------- | :---- |
| `r⃗₁`      | [x₁ y₁ z₁]  | DU    |
| `r⃗₂`      | [x₂ y₂ z₂]  | DU    |

#### Optional Keyword Arguments:
| Parameter       | Description                                           | Units |
| --------------: | :---------------------------------------------------- | :---- |
| `R=R⨁`          | Blocking body mean radius                             | DU    |
| `bodyShape="e"` | If "e" then z₁ and z₂ are scaled by flattening factor | -     |
| `f=f⨁`          | Body flattening factor                                | -     |

### Values (Outputs):
| Parameter      | Description                      | Units |
| -------------: | :------------------------------- | :---- |
| `true`/`false` | If two points have line of sight | -     |

Author: James Spicer (2016.11.22) <br/>  
Source: Vallado (2001), Algorithm 35, p.294 (simplifed by JPWS)
"""
function sight(r⃗₁::Array{T}, r⃗₂::Array{T}; R::T=R⨁, bodyShape::String="e", f::T=f⨁) where T <: Float64 # From Vallado 'sight.m'.
    a, b = copy(r⃗₁), copy(r⃗₂)

    # --------------------- scale z component ---------------------
    (bodyShape == "e") && ((a[3], b[3]) = inv(1-f)*[a[3], b[3]])

    ȧb = a⋅b
    sign(norm(a)^2-ȧb) ≠ sign(norm(b)^2-ȧb) && return true
    norm(a × b)/norm(a - b) ≥ R && return true
    return false
end

## LOS between points A and B, with blocking sphere radius R at point C.
function sight(A::Array{T}, B::Array{T}, C::Array{T}, R::T) where T <: Float64
    A², ȦB, ȦC, ḂC = norm(A)^2, A⋅B, A⋅C, B⋅C
    τ = (A² + ḂC - ȦB - ȦC)/(A² + norm(B)^2 - 2ȦB)
    !(0 ≤ τ ≤ 1) && return true
    norm(C)^2 + (1-τ)*(A² - 2ȦC) + τ*(ȦB - ȦC - ḂC) ≥ R^2 && return true    
    return false
end

"Compute slant range [km] from ground station given orbit height [km] and elevation angle [rad]."
slantRange(h::T, ε::T; R::T=R⨁) where T <: Float64 = -R*sin(ε) + √((R+h)^2 - (R*cos(ε))^2)


"Compute angular diameter [rad] of sphere of radius r observed distance R from its center."
angularDiameter(r::T, R::T) where T <: Float64 = 2asin(r./R)


"Compute angular radius [rad] of sphere of radius r observed distance R from its center."
angularRadius(r::T, R::T) where T <: Float64   =  asin(r./R)


hav(θ::Float64)  = sin(0.5θ).^2
ahav(θ::Float64) = 2asin(√θ)


"Compute central angle between spherical lat/lon vectors [ϕ1, λ1] & [ϕ2, λ2] with same center"
haversine(ϕ₁::T, λ₁::T, ϕ₂::T, λ₂::T) where T <: Float64 = ahav(hav(ϕ₂-ϕ₁) + cos(ϕ₁).*cos(ϕ₂).*hav(λ₂-λ₁))


"Compute `n`th centered hexagonal number."
centeredHexagonal(n::Int64) = 3n*(n-1) + 1


"Compute `n`th hexagonal number."
hexagonal(n::Int64) = n*(2n-1)


"Compute distance from center to vertex of `n`-sided polygon with inradius `r₀`."
minRadius(r₀::Real, n::Int64) = r₀*sec(π/n)


"Compute side length of an `n`-sided polygon with circumradius `R`."
sideLength(R::Real, n::Int64) = 2R*sinpi(1/n)


"Compute inradius of regular `n`-sided polygon with side length `s`."
inradius(s::Real, n::Int64) = 0.5s*cot(π/n)


"Compute circumradius of an `n`-sided polygon with side length `s`."
circumRadius(s::Real, n::Int64) = 0.5s*csc(π/n)

"Compute center coordinates of 2D triangle's incircle."
function incenter(xa,ya,xb,yb,xc,yc)
    a = hypot(xb-xc, yb-yc)
    b = hypot(xa-xc, ya-yc)
    c = hypot(xa-xb, ya-yb)
    x = (a*xa + b*xb + c*xc)/(a+b+c)
    y = (a*ya + b*yb + c*xc)/(a+b+c)
    return x, y
end


VsphericalSector(ϕ::Real, r::Real) = 2π*(r^3)*(1-cos(ϕ))/3 # ϕ is half the cone angle


"Return `k`x2 array of `(x,y)` coordinates of `k` evenly-spaced (by central angle) points along the edge of a regular `n`-sided polygon."
function polygonPatternAngle(n::Int64, k::Int64; r=1.0)
    θ  = (k == 1) ? π : range(0, stop=2π, length=k+1)
    r *= cos(-π/n) * sec(θ % (2π/n) - π/n)
    return r.*[cos(θ) sin(θ)]
end


"Return `k`×2 array of `(x,y)` coordinates of `k` evenly-spaced (along perimeter) points along the edge of a regular `n`-sided polygon."
function polygonPatternEdge(n::Int64, k::Int64; r=1)
    (k % n ≠ 0) && error("Polygon must have an integer number of points per side, i.e. k % n = 0.")

    j = Int(k/n) # central angles per edge
    α = 2π/n # Central angle of each isosceles triangle
    s = 2r*sin(π/n) # Side length opposite α

    β = 0.5(π-α)
    R = r
    subAngles = zeros(j)
    for i = 1:j
        sβ = s*sin(β)
        γ  = atan(sβ, R*j - s*cos(β)) # j'th central angle
        R  = sβ*csc(γ)/j  # Update R
        β += γ                  # Update β
        subAngles[i] = γ
    end

    θ = zeros(k)
    for u = 0:n-1, v = 1:j
        θ[u*j + v] = u*α + sum(subAngles[1:v])
    end

    r *= cos(-π/n) * sec(θ % (2π/n) - π/n)
    return r.*[cos(θ) sin(θ)]
end


"Compute angle between vectors v₁ and v₂."
vectorAngle(v₁, v₂) = acos((v₁⋅v₂)/(norm(v₁)*norm(v₂)))

vectorDistance(v₁, v₂) = √sum([(v₁[i]-v₂[i])^2 for i = 1:length(v₁)])

"Compute volume of ellipsoid with major axis half-lengths a, b, and c."
vEllipsoid(a::Real, b::Real, c::Real) = 4π*(a.*b.*c)/3

"Compute volume of a sphere of radius `r`."
vSphere(r::Real) = 4π*(r.^3)/3


# ------------------------------------------------------------------------------
#                           function footprint
#
#  Calculates the Earth coverage "footprint" of a satellite-borne antenna.
#  Antenna beam is conical with half-angle α pointed at target point P at
#  (λT, ϕT).
#
#  author        : James Spicer                   650-999-0331   22 mar 2016
#
#  revisions
#                -
#
#  inputs          description                    range / units
#    v           - input vector
#    ϕ           - angle of rotation              [rad]
#    λS          - satellite's longitude/azimuth  [rad]
#    ϕS          - satellite's latitude/azimuth   [rad]
#    RS          - satellite's orbit radius       [DU] (e.g. [km])
#    λT          - target's longitude             [rad]
#    ϕT          - target's latitude              [rad]
#    α           - beam's half-angle              [rad]
#    n           - number of points in footprint  [-] (Default n=500)
#    RP          - radius of target's planet      [DU] (e.g. [km], default RP=R⨁)
#
#  outputs       :
#    λF          - Array of footprint boundary
#                  latitudes                      [rad]
#    ϕF          - Array of footprint boundary
#                  latitudes                      [rad]
#
# λF, ϕF = footprint(λS, ϕS, RS, λT, ϕT, α; n=500, RP=R⨁)
#
# references     :
#   Zolnay, 1971
# -----------------------------------------------------------------------------
# function footprint(λS, ϕS, RS, λT, ϕT, α; n=500, RP=R⨁)
#     λd  = λT - λS

#     # Range vector of satellite
#     # Rd  = RP*[ cos(ϕT)*sin(λd),
#     #           -cos(ϕT)*cos(λd),
#     #            sin(ϕT)        ]
#     # Rd  = RP*[0, -cos(ϕT), sin(ϕT)]
#     Rd  = -[0, RP, 0]
#     Rd  = rot1(Rd,  ϕT)
#     Rd  = rot3(Rd, -λd)
#     Rd += RS*[0, cos(ϕS), sin(ϕS)]
#     R   = rot1(Rd,  ϕS)

#     pd  = asin(Rd[1]/norm(Rd[1:2]))     # Eqn. (6)
#     ρd  = atan(Rd[3], norm(Rd[1:2]))   # Eqn. (7)

#     modu = norm(R)*tan(α)
#     bmax = asin(RP/RS)

#     β  = linspace(0, 2π, n)
#     λF = fill(NaN, length(β))
#     ϕF = fill(NaN, length(β))
#     for n = 1:length(β)
#         # undd = modu*[cos(β[n]), 0, sin(β[n])]
#         undd = [modu, 0, 0]
#         undd = rot2(undd, β[n])
#         und  = rot1(undd, -ρd)
#         und  = rot3(und,   pd)

#         Mnd  = Rd + und
#         Mn   = rot1(Mnd, ϕS)

#         bn   = acos(Mn[2]/norm(Mn))
#         if bn ≤ bmax
#             cn  = asin(RS*sin(bn)/RP)
#             dn  = cn - bn
#             Mn  = √(RP^2 + RS^2 - 2*RS*RP*cos(dn))
#             pn  = asin( Mnd[1]/ norm(Mnd[1:2]))
#             ρn  = atan(Mnd[3], norm(Mnd[1:2]))
#             # rn  = Mn*[cos(ρn)*sin(pn),
#             #           cos(ρn)*cos(pn),
#             #           sin(ρn)        ]
#             rn  = Mn*[0, cos(ρn), sin(ρn)]
#             rn  = rot3(rn, pn)
#             rn -= RS*[0, cos(ϕS), sin(ϕS)]
#         else
#             τn  = acos(Mn[3]/√(Mn[1]^2 + Mn[3]^2))
#             # rn  = RP*[cos(bmax)*sin(τn),
#             #           sin(bmax)        ,
#             #           cos(bmax)*cos(τn)]
#             rn  = RP*[0, sin(bmax), cos(bmax)]
#             rn  = rot2(rn, -τn)
#             rn  = rot1(rn,  ϕS)
#         end

#         λdn   = asin( rn[1]/ norm(rn[1:2]))
#         ϕF[n] = atan(rn[3], norm(rn[1:2]))
#         λF[n] = λS + λdn
#     end
#     return λF, ϕF
# end


# ------------------------------------------------------------------------------
#                           function footprint
#
#  Calculates the Earth coverage "footprint" of a satellite-borne antenna.
#  Antenna beam is conical with half-angle α pointed at a target point.
#
#  author        : James Spicer                   650-999-0331   15 jun 2016
#
#  revisions
#                -
#
#  inputs          description                    range / units
#    O_r_S_ECEF  - satellite's position vector
#                  in ECEF coordinates            [km]
#    O_r_T_ECEF  - target's position vector
#                  in ECEF coordinates            [km]
#    α           - beam's half-angle              [rad]
#    n           - number of points in footprint  [-] (Default n=500)
#    a, b, c     - axes of footprint's ellipse    [km]
#
#  outputs       :
#    coords      - array of ECEF coordinates
#                  describing footprint (NaN if
#                  no Earth intersection)         [km]
#
# coords = footprint(O_r_S_ECEF, O_r_T_ECEF, α; n=500)
# TODO: Add Earth limb (last point of contact with ellipsoid)
# function footprint(O_r⃗_S_ECEF::Array{Float64}, O_r⃗_T_ECEF::Array{Float64}, α::Float64; n::Int64=500, a::Float64=R⨁_a, b::Float64=R⨁_a, c::Float64=R⨁_b)
#     coords = fill(NaN, n, 3)
#     Φ      = linspace(0, 2π, n)

#     if sight(O_r⃗_S_ECEF, O_r⃗_T_ECEF)
#         # Rotate + translate satellite position into T's SEZ frame (T becomes origin).
#         O_r⃗_T_GDC = ecef2gdc(O_r⃗_T_ECEF)
#         O_r⃗_S_SEZ = ecef2sez(O_r⃗_S_ECEF - O_r⃗_T_ECEF, O_r⃗_T_GDC[1], O_r⃗_T_GDC[2])

#         for i = 1:n # Step round T's local horizon (x-y plane in SEZ frame) 0 to 2π
#             # Locate point P where SrP is angle α from SrT line
#             v          = [cos(Φ[i]) sin(Φ[i]) 0]
#             β          = vectorAngle(O_r⃗_S_SEZ, v) # angle point vector makes with satellite vector
#             l          = norm(O_r⃗_S_SEZ)*sin(α)/sin(α+β) # distance to α-point on line.
#             O_r⃗_P_SEZ  = v*l

#             # rotate P back into ECEF frame
#             O_r⃗_P_ECEF = ecef2sez(O_r⃗_P_SEZ, O_r⃗_T_GDC[1], O_r⃗_T_GDC[2], inv=true) + O_r⃗_T_ECEF
#             S_r⃗_P_ECEF = O_r⃗_P_ECEF - O_r⃗_S_ECEF

#             # Find intersection values of t of line O_r⃗_S_ECEF + t*S_r⃗_P_ECEF with ellipsoid Earth
#             T = ellipsoidLineIntersectV(O_r⃗_S_ECEF, S_r⃗_P_ECEF, a=a, b=b, c=c)

#             !isempty(T) && (coords[i,:] = O_r⃗_S_ECEF + minimum(T)*S_r⃗_P_ECEF)
#         end
#     end
#     return coords
# end


# function footprint{T <: Float64}(O_r⃗_S::Array{T}, O_r⃗_T::Array{T}, α::T; bodyShape::String="e", n::Int64=500, a::T=R⨁_a, b::T=R⨁_a, c::T=R⨁_b, R::T=R⨁)
#     coords = fill(NaN, n, 3)

#     !sight(O_r⃗_S, O_r⃗_T) && return coords

#     S_r⃗_T = O_r⃗_T - O_r⃗_S            # Satellite → target vector
#     S_r̂_T = S_r⃗_T[:]/norm(S_r⃗_T)     # Satellite → target direction

#     n̂₀    = (S_r⃗_T[1] == S_r⃗_T[2] == 0) ? normalize([1 0 0] × S_r⃗_T)' : normalize([0 0 1] × S_r⃗_T)' # X/Z × S_r⃗_T

#     # O_r_S = norm(O_r⃗_S)
#     # βmax     = asin(R/O_r_S) # ∠ between the satellite position vector and the satellite's horizon.
#     # S_r_Pmax = √(O_r_S^2 - R^2)

#     θ     = linspace(0, 2π, n)
#     for i = 1:n
#         n̂     = rotMat(θ[i], S_r̂_T)*n̂₀ # Rotate normal vector θ[i] around satellite → target axis.
#         S_r̂_P = rotMat(α, n̂)*S_r̂_T     # Satellite → point direction

#         # β = acos(-(S_r̂_P⋅O_r⃗_S)/O_r_S) # ∠ between satellite position vector and satellite → point vector.

#         if β < βmax # This ray intersects the Earth, find p.o.i.
        
#             # Find intersection values (t) of line O_r⃗_S + t*S_r̂_P with Earth
#             t = (bodyShape == "e") ? ellipsoidLineIntersectV(O_r⃗_S, S_r̂_P, a=a, b=b, c=c) : lineSphereIntersectV(O_r⃗_S, S_r̂_P, r=R)
#             # if !isempty(t)
#             coords[i,:] = O_r⃗_S + minimum(t)*S_r̂_P'
#         else  # This ray does not intersect Earth, find point of last intersection.
#             # S_r̂_P       = rotMat(βmax, n̂)*-O_r̂_S
#             # coords[i,:] = O_r⃗_S + S_r_Pmax*S_r̂_P
#         end
#     end
#     return coords
# end

# DONE: Add Earth limb. 
# TODO: Make work if body is not at origin.
# TODO: Use normalize()
# TODO: Make work for an ellipsoid Earth.
# TODO: Gallery of footprint plots alongside satellite pov plots of test scenarios including edge cases.
# DONE: Check if beam intersects body. 
"""
    footprint(O_r⃗_S::Array{T}, O_r⃗_T::Array{T}, α::T; n::Int64=500, R::T=R⨁) where T <: Float64

Compute coordinates of conical beam intersection with spherical body.

### Arguments (Inputs):
| Parameter | Description                           | Units        |
| --------: | :------------------------------------ | :----------- |
| `O_r⃗_S`   | Body-centered satellite coordinates   | DU (e.g. km) |
| `O_r⃗_T`   | Body-centered beam target coordinates | DU           |
| `α`       | Beam cone half-angle                  | rad          |

#### Optional Keyword Arguments:
| Parameter       | Description                                     | Units  |
| --------------: | :---------------------------------------------- | :----- |
| `n=500`         | Number of points in footprint (accuracy)        | -      |
| `R=R⨁`          | Radius of spherical body                        | DU     |     

### Values (Outputs):
| Parameter | Description                                      | Units |
| --------: | :----------------------------------------------- | :---- |
| `coords`  | n×3 array of body-centered footprint coordinates | DU    |

Author: James Spicer (2017.01.30) <br/>  
"""
function footprint(O_r⃗_S::Array{T}, O_r⃗_T::Array{T}, α::T; n::Int64=500, R::T=R⨁) where T <: Float64
    coords = fill(NaN, n, 3)

    !sight(O_r⃗_S, O_r⃗_T, R=R, bodyShape="s") && return coords         # Body blocks Target, no footprint.

    O_r⃗_T, O_r⃗_S = O_r⃗_T[:], O_r⃗_S[:]
    S_r⃗_T = O_r⃗_T - O_r⃗_S                                             # Satellite → Target vector
    S_r̂_T = S_r⃗_T/norm(S_r⃗_T)                                         # Satellite → Target direction
                    
    O_r_S = norm(O_r⃗_S)                                               # Body → Satellite distance
    O_r̂_S = O_r⃗_S/O_r_S                                               # Body → Satellite direction
                        
    βmax  = (O_r_S < R) ? Inf : asin(R/O_r_S)                         # Apparent ∠ between Body center and its Horizon (`Inf` if Satellite inside Body).
    ∠OST  = acos(-max(min(O_r̂_S⋅S_r̂_T, 1), -1))                       # Apparent ∠ between Body center and Target.
                    
    if (βmax > α && ∠OST < βmax - α)                                  # Beam entirely intersects Body.
        partial = false                 
    elseif βmax < α && ∠OST < α - βmax                                # Beam entirely envelops Body.
        return limb(O_r⃗_S, n=n, R=R)                   
    elseif ∠OST < βmax + α                                            # Beam partially intersects Body.
        partial   = true                    
        S_r_Pmax  = √(O_r_S^2 - R^2)                                  # Satellite-Horizon distance.
        S_r⃗_T_sph = cart2sph( S_r⃗_T)                                  # Angular position of Target from Satellite.
        S_r⃗_O_sph = cart2sph(-O_r⃗_S - O_r⃗_S)                          # Angular position of Body from Satellite.
                    
        if α < βmax < ∠OST                                            # If not all of beam's "points of last intersection" are in the footprint.
            issues = true 
            OdT    = norm(S_r⃗_T_sph[1:2]-S_r⃗_O_sph[1:2])              # Angular Body ↔ Target distance.
            y      = √(βmax^2 - (0.5(OdT^2 - α^2 + βmax^2)/OdT)^2)    # Angular distance from line OdT to poi between beam and body. From intersection of two circles (Wolfram).
            ϕ      = asin(y/βmax)                                     # Apparent angle between line OdT and poi.
            argT   = abs(atan((S_r⃗_T_sph[2]-S_r⃗_O_sph[2])/(S_r⃗_T_sph[1]-S_r⃗_O_sph[1]))) # Argument of the Target's angular position.
        else
            issues = false
        end
    else                                            # Beam does not intersect Body at all.
        ∠OST > α + βmax && return coords
    end

    # Create vector normal to Satellite ↔ Target direction to be rotated in the loop.
    n̂₀    = (S_r⃗_T[1] == S_r⃗_T[2] == 0) ? normalize([1 0 0] × S_r⃗_T)' : normalize([0 0 1] × S_r⃗_T)' # X/Z × S_r⃗_T
       
    θ  = range(0, stop=2π, length=n)
    for i = 1:n        
        n̂     = rotMat(θ[i], S_r̂_T)*n̂₀ # Rotate normal vector θ[i] around satellite → target axis.
        S_r̂_P = rotMat(α, n̂)*S_r̂_T     # Satellite → Point direction

        partial && (∠OSP = acos(-O_r̂_S⋅S_r̂_P)) # ∠ between Point and Body.
        
        # Find direction S_r̂_P and distance S_r_P of each Point on the footprint.
        if partial && ∠OSP > βmax # This ray does not intersect Body, find point of last intersection.
            # Line defined by P and O, circle defined by βmax and O. New P is poi of line and circle.
            S_r̂_P_sph = cart2sph(S_r̂_P) # Angular position of Point.
            
            if issues # Check the point of last intersection would be within footprint.
                argP = abs(atan((S_r̂_P_sph[2]-S_r⃗_O_sph[2])/(S_r̂_P_sph[1]-S_r⃗_O_sph[1]))) # Argument of Point.
                abs(argT-argP) > ϕ && continue # If Point argument is outside that of lens-shaped beam-body intersection, ignore.
            end
            
            # Find intersection point of Target-Body line with horizon and update P.
            p1         = S_r̂_P_sph - S_r⃗_O_sph  # Put body at origin.
            p2         = S_r⃗_O_sph - S_r⃗_O_sph  # Put body at origin.
            poi1, poi2 = circleLineIntersectC(p1[1:2], p2[1:2], βmax)
            poi        = sign(p1[1]) == sign(poi1[1]) ? poi1 : poi2
            S_r̂_P_sph  = [poi + S_r⃗_O_sph[1:2]' 1.0] # New angular position of Point.
            
            S_r_P      = S_r_Pmax
            S_r̂_P      = sph2cart(S_r̂_P_sph)[:] # Cartesian direction of Point.
        else # This ray intersects the Body, find intersection values (S_r_P) of line O_r⃗_S + S_r_P*S_r̂_P with Body.
            arr        = lineSphereIntersectV(O_r⃗_S, S_r̂_P, r=R) # Get distances of points of intersection.
            S_r_P      = minimum(arr[arr .> 0]) # Pick smallest (nearest) positive distance.
        end
        coords[i,:] = O_r⃗_S[:]' + S_r_P*S_r̂_P' 
    end
    return coords
end


# TODO: Make work for an ellipsoid Earth.
# TODO: Make work if body is not at origin.
"""
    limb(O_r⃗_S::Array{T}; n::Int64=500, R::T=R⨁, h::T=0.0) where T <: Float64    

Compute ring of coordates at which a viewer at altitude `h` above a blocking
spherical body of radius `R` cannot see a satellite at position `O_r⃗_S`.

### Arguments (Inputs):
| Parameter | Description                         | Units        |
| --------: | :---------------------------------- | :----------- |
| `O_r⃗_S`   | Body-centered satellite coordinates | DU (e.g. km) |

#### Optional Keyword Arguments:
| Parameter | Description                            | Units  |
| --------: | :------------------------------------- | :----- |
| `n=500`   | Number of points in horizon (accuracy) | -      |
| `R=R⨁`    | Body's radius                          | DU     |
| `h=0.0`   | Altitude of viewer above body          | DU     |     

### Values (Outputs):
| Parameter | Description                                 | Units |
| --------: | :------------------------------------------ | :---- |
| `coords`  | n×3 array of body-centered limb coordinates | DU    |

Author: James Spicer (2017.01.30) <br/>  
"""
function limb(O_r⃗_S::Array{T}; n::Int64=500, R::T=R⨁, h::T=0.0) where T <: Float64
    n̂₀     = (O_r⃗_S[1] == O_r⃗_S[2] == 0) ? normalize([1 0 0] × O_r⃗_S)' : normalize([0 0 1] × O_r⃗_S)' # X/Z × O_r⃗_S consider nullspace?

    O_r_S  = norm(O_r⃗_S)
    O_r̂_S  = O_r⃗_S[:]/O_r_S
    β      = asin(R/O_r_S) # Angle between the satellite position vector and the satellite's horizon
    S_r_P  = √(O_r_S^2 - R^2) + √(2R*h + h^2) # Simplified cosine rule.

    θ      = range(0, stop=2π, length=n)
    coords = zeros(n, 3)
    for i = 1:n
        n̂           = rotMat(θ[i], O_r̂_S)*n̂₀
        S_r̂_P       = rotMat(β, n̂)*-O_r̂_S
        coords[i,:] = O_r⃗_S[:] + S_r_P*S_r̂_P
    end
    return coords
end


""" 
    atan3{T <: Float64}(sinϕ::T, cosϕ::T; ϵ::T=1e-10)

Four-quadrant inverse tangent. Returns radian angle 0 ≤ ϕ ≤ 2π.
"""
function atan3(sinϕ::T, cosϕ::T; ϵ::T=1e-10) where T <: Float64
    (abs(sinϕ) < ϵ) && return 0.5π*(1 - sign(cosϕ))
    c = π*(1 - 0.5sign(sinϕ))
    (abs(cosϕ) < ϵ) && return c
    return c + sign(sinϕ)*sign(cosϕ)*(abs(atan(sinϕ/cosϕ)) - 0.5π)
end