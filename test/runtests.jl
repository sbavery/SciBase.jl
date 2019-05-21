using AstroSciKit, Base.Test, Base.Dates

println("Running tests:")

include("test_geometry.jl")


# Testing dms2rad & rad2dms:
d1, m1, s1 = 35, 15, 53.63 # [deg min sec]
d = 0.6154886 # [rad]

@test_approx_eq_eps dms2rad( d1,  m1,  s1)  d 1e-6
@test_approx_eq_eps dms2rad(-d1,  m1,  s1) -d 1e-6
@test_approx_eq_eps dms2rad(-d1, -m1, -s1) -d 1e-6
@test_approx_eq_eps dms2rad( d1, -m1,  s1) -d 1e-6

@test all(abs(rad2dms(-d) - [-d1 -m1 -s1]) .< 1e-2)
@test all(abs(rad2dms( d) - [ d1  m1  s1]) .< 1e-2)

d2, m2, s2 = rad2dms(dms2rad(d1, m1, s1))
@test all(abs([d1 m1 s1] - [d2 m2 s2]) .< 1e-10)

d2, m2, s2 = rad2dms(dms2rad(-d1, -m1, -s1))
@test all(abs([-d1 -m1 -s1] - [d2 m2 s2]) .< 1e-10)



# Testing hms2rad & rad2hms:
h1, m1, s1 = 15, 15, 53.63 # [hr min sec]
τ = 3.996341 # [rad]

@test_approx_eq_eps hms2rad(h1, m1, s1) τ 1e-6
@test all(abs(rad2hms(τ) - [h1 m1 s1]) .< 1e-2)

h2, m2, s2 = rad2hms(hms2rad(h1, m1, s1))
@test all(abs([h1 m1 s1] - [h2 m2 s2]) .< 1e-10)



# Testing hms2τ & τ2hms:
h1, m1, s1 = 13, 22, 45.98 # [hr min sec]
τ = 48165.98 # [sec] (time of day)

@test_approx_eq_eps hms2τ(h1, m1, s1) τ 1e-6
@test all(abs(τ2hms(τ) - [h1 m1 s1]) .< 1e-2)

h2, m2, s2 = τ2hms(hms2τ(h1, m1, s1))
@test all(abs([h1 m1 s1] - [h2 m2 s2]) .< 1e-10)


# Test shadow:
println("Testing Shadow")
# solns = ["umbra", "penumbra", "penumbra", "none", "none", "penumbra"]
rSun  = [-km_in_AU, 0, 0]

# rSats = ([R⨁*1.5, 0, 0],
# 		 [R⨁*1.5, 0, R⨁],
# 		 [R⨁*1.5, R⨁, 0],
# 		 [0, R⨁*1.5, 0],
# 		 [-R⨁*1.5, 0, 0],
# 		 [1.001*km_in_AU*R⨁/(R☉ - R⨁), 0, 0])

# for i = 1:length(solns)
# 	@test shadow(rSun, rSats[i]) == solns[i]
# end

@test shadow(rSun, [R⨁*1.5, 0, 0])  == 0.0
@test 0.0 < shadow(rSun, [R⨁*1.5, 0, R⨁]) < 1.0
@test 0.0 < shadow(rSun, [R⨁*1.5, R⨁, 0]) < 1.0
@test shadow(rSun, [0, R⨁*1.5, 0])  == 1.0
@test shadow(rSun, [-R⨁*1.5, 0, 0]) == 1.0
@test 0.0 < shadow(rSun, [1.001*km_in_AU*R⨁/(R☉ - R⨁), 0, 0]) < 1.0


# Test sight:
println("Testing Sight")
r1 = [0, -4464.696, -5102.509]
r2 = [0, 5740.323, 3189.068]

@test sight(r1, r2) == false



# Test Lagrange point location
# Source: http://descanso.jpl.nasa.gov/monograph/series12/LunarTraj--Overall.pdf, p.44
println("Testing Lagrange Points")
μ = 4.90280058e3 / (4.90280058e3+3.98600432897e5)
@test isapprox(lagrangePoint(μ, 1)[1],  0.8369151324)
@test isapprox(lagrangePoint(μ, 2)[1],  1.1556821603)
@test isapprox(lagrangePoint(μ, 3)[1], -1.0050626453)
@test isapprox(lagrangePoint(μ, 4)[1],  0.4878494157)
@test isapprox(lagrangePoint(μ, 4)[2],  0.8660254038)
@test isapprox(lagrangePoint(μ, 5)[1],  0.4878494157)
@test isapprox(lagrangePoint(μ, 5)[2], -0.8660254038)

μ = (3.98600432897e5+4.90280058e3) / (3.98600432897e5+1.32712440e11+4.90280058e3)
@test isapprox(lagrangePoint(μ, 1)[1],  0.9899859823)
@test isapprox(lagrangePoint(μ, 2)[1],  1.0100752000)
@test isapprox(lagrangePoint(μ, 3)[1], -1.0000012670)
@test isapprox(lagrangePoint(μ, 4)[1],  0.4999969596)
@test isapprox(lagrangePoint(μ, 4)[2],  0.8660254038)
@test isapprox(lagrangePoint(μ, 5)[1],  0.4999969596)
@test isapprox(lagrangePoint(μ, 5)[2], -0.8660254038)

rₗ₁, rₗ₂, rₗ₃, rₗ₄, rₗ₅, r₁, r₂ = LCoords(4.90280058e3, 3.98600432897e5, 384400.0)
@test isapprox(rₗ₁[1],  321710.177)
@test isapprox(rₗ₂[1],  444244.222)
@test isapprox(rₗ₃[1], -386346.081)
@test isapprox(rₗ₄[1],  187529.315)
@test isapprox(rₗ₄[2],  332900.165)
@test isapprox(rₗ₅[1],  187529.315)
@test isapprox(rₗ₅[2], -332900.165)

@test isapprox(r₁[1], -  4670, rtol=1e-3)
@test isapprox(r₂[1],  379729, rtol=1e-3)

rₗ₁, rₗ₂, rₗ₃, rₗ₄, rₗ₅, r₁, r₂ = LCoords(3.98600432897e5+4.90280058e3, 1.32712440e11, km_in_AU)
@test isapprox(rₗ₁[1],  148099795.0)
@test isapprox(rₗ₂[1],  151105099.2)
@test isapprox(rₗ₃[1], -149598060.2)
@test isapprox(rₗ₄[1],   74798480.5)
@test isapprox(rₗ₄[2],  129555556.4)
@test isapprox(rₗ₅[1],   74798480.5)
@test isapprox(rₗ₅[2], -129555556.4)


# Test Lagrange equations. Source: wikipedia.org/wiki/Lagrangian_point
for (μ₁, μ₂, R) in [(μ⨁, μ☾, r_E_M), (μ☉, μ⨁, km_in_AU)]
	rₗ₁, rₗ₂, rₗ₃, rₗ₄, rₗ₅, r₁, r₂ = LCoords(μ₁, μ₂, R)

	r   = abs(r₂[1] - rₗ₁[1])
	LHS = μ₁ /(r-R)^2
	RHS = μ₂/r^2 + μ₁/R^2 - r*(μ₁+μ₂)/R^3
	@test isapprox(LHS, RHS)

	r   = abs(r₂[1] - rₗ₂[1])
	LHS = μ₁ /(r+R)^2 + μ₂/r^2
	RHS = μ₁/R^2 + r*(μ₁+μ₂)/R^3
	@test isapprox(LHS, RHS)

	r   = abs(abs(rₗ₃[1] - r₁[1]) - abs(r₂[1] - r₁[1]))
	LHS = μ₁/(R-r)^2 + μ₂/(2*R-r)^2
	RHS = (μ₂*R/(μ₁+μ₂) + R - r)*(μ₁+μ₂)/R^3
	@test isapprox(LHS, RHS)
end


# Test Sun
# Vallado (2001) p.266-7
println("Testing Sun")
t = DateTime(1994, 4, 2)
JD = datetime2julian(t)

rsun = [sun(JD) 0.0 0.0 0.0] # [λecl ϕecl r] [rad rad km]

@test isapprox(rsun[1], deg2rad(12.022725), rtol=1e-6)
@test isapprox(rsun[3], 0.9994850*km_in_AU, rtol=1e-6)

rsun = αδr2ijk(rsun)[1:3] # XYZ
rsun = [ecl2eci(rsun') 0.0 0.0 0.0] # IJK

@test isapprox(rsun[1], 146241097, rtol=1e-4)
@test isapprox(rsun[2],  28574940, rtol=1e-4)
@test isapprox(rsun[3],  12389196, rtol=1e-4)

rsun = αδr2ijk(rsun, inv=true)

@test isapprox(rsun[1], deg2rad(11.056071 ), rtol=1e-4)
@test isapprox(rsun[2], deg2rad( 4.7529393), rtol=1e-4)


# Test Moon
# Vallado (2001) p.274-5
println("Testing Moon")
t = DateTime(1994, 4, 28)
JD = datetime2julian(t)

rmoon = [moon(JD) 0.0 0.0 0.0] # [λecl ϕecl r] [rad rad km]

@test isapprox(rad2deg(rmoon[1]), 248.2370233, rtol=1e-6)
@test isapprox(rad2deg(rmoon[2]),   1.2185048, rtol=1e-6)
@test isapprox(rmoon[3]/reEarth,   56.7790575, rtol=1e-6)

rmoon = αδr2ijk(rmoon)[1:3]
rmoon = [ecl2eci(rmoon') 0.0 0.0 0.0]

@test isapprox(rmoon[1], -134241.192, rtol=1e-4)
@test isapprox(rmoon[2], -311571.349, rtol=1e-4)
@test isapprox(rmoon[3], -126693.681, rtol=1e-4)

rmoon = αδr2ijk(rmoon, inv=true)

@test isapprox(rad2deg(rmoon[1]), 246.691102, rtol=1e-4)
@test isapprox(rad2deg(rmoon[2]), -20.477711, rtol=1e-4)



println("Testing Kepler's Equation Solvers")
# Vallado (2001) p. 74-75
M = deg2rad(235.4)
e = 0.4

@test kepEqtnE(M, e, max_i=0) ≈ 3.84697110600265
@test kepEqtnE(M, e, max_i=1) ≈ 3.84866146081460
@test kepEqtnE(M, e, max_i=2) ≈ 3.84866174509716
@test kepEqtnE(M, e, max_i=3) ≈ 3.84866174509717
@test kepEqtnE(M, e)          ≈ deg2rad(220.512074767522)

# Vallado (2001), Example 2-1, p. 284
M = deg2rad(-150.443142)
e = 0.048486

E = kepEqtnE(M, e)
@test isapprox(A2ν(E, e), deg2rad(206.95453), rtol=1e-6)

# Vallado (2001), Example 2-2, p.77-8
∆t = 4 # [TU]
p  = 4 # [ER]
@test isapprox(kepEqtnP(∆t, p, μ=1.0), 0.8177316, rtol=1e-6)

∆t = 53.7874*60 # [sec]
p  = 25512 # [km]
@test isapprox(kepEqtnP(∆t, p, μ=398600.4415), 0.8177316, rtol=1e-4)

# Vallado (2001) p. 74-75
M = deg2rad(235.4)
e = 2.4

@test kepEqtnH(M, e, max_i=0) ≈ 1.6074355626
@test kepEqtnH(M, e, max_i=1) ≈ 1.6013962815
@test kepEqtnH(M, e, max_i=2) ≈ 1.60137614515
@test kepEqtnH(M, e, max_i=3) ≈ 1.6013761449
@test kepEqtnH(M, e)          ≈ 1.601376144

# Vallado (2001) p. 74-75
M = deg2rad(235.4)
e = 0.4
@test anomaly(e, M=M) ≈ deg2rad(220.512074767522)

# # Vallado (2001), Example 2-2, p.77-8
∆t = 4.0 # [TU]
p  = 4.0 # [ER]
e  = 1.0
@test isapprox(anomaly(e, Δt=∆t, p=p, μ=1.0), 0.8177316, rtol=1e-6)

∆t = 53.7874*60 # [sec]
p  = 25512.0 # [km]
@test isapprox(anomaly(e, Δt=∆t, p=p, μ=μ⨁), 0.8177316, rtol=1e-4)

# # Vallado (2001) p. 74-75
M = deg2rad(235.4)
e = 2.4
@test anomaly(e, M=M) ≈ 1.601376144





println("Testing Kepler's Problem")
# Vallado (2001), Example 2-6, p.125-6 (See errata)
p = 11067.79 # [km]
e = 0.83285
ι = deg2rad( 87.87)
Ω = deg2rad(227.89)
ω = deg2rad( 53.38)
ν = deg2rad( 92.335)

r, v = RANDV(p, e, ι, Ω, ω, ν)

@test isapprox(r[1],  6525.344, rtol=1e-4)
@test isapprox(r[2],  6861.535, rtol=1e-6)
@test isapprox(r[3],  6449.125, rtol=1e-6)
@test isapprox(v[1],  4.902276, rtol=1e-6)
@test isapprox(v[2],  5.533124, rtol=1e-4)
@test isapprox(v[3], -1.975709, rtol=1e-6)





println("Testing Planets")
# Vallado (2001), Example 5-5, p. 283-5
t = DateTime(1994, 5, 20, 20)
JD = datetime2julian(t)
TTDB = JD2JDC(JD)
a, e, ι, Ω, ω̃, λM = planetEphemerides("jupiter", TTDB)

@test isapprox(a, 5.202603, rtol=1e-6)
@test isapprox(e, 0.048486, rtol=1e-4)
@test isapprox(rad2deg(ι),     1.303382, rtol=1e-6)
@test isapprox(rad2deg(Ω),   100.454519, rtol=1e-6)
@test isapprox(rad2deg(ω̃),    14.319203, rtol=1e-6)
@test isapprox(rad2deg(λM), 360-136.12394,  rtol=1e-6)

@test isapprox(rad2deg(obliquity(TTDB)), 23.440021, rtol=1e-6)

r, v = planetRV("jupiter", JD) # (xyz)
r    = ecl2eci(r, JD=JD) # [km]     (xyz(fk5))
v    = ecl2eci(v, JD=JD) # [km/sec] (xyz(fk5))

@test isapprox(r[1], -609750815, rtol=1e-6)
@test isapprox(r[2], -497437968, rtol=1e-6)
@test isapprox(r[3], -198394488, rtol=1e-6)
@test isapprox(v[1],   0.145626, rtol=1e-3)
@test isapprox(v[2],  -0.144349, rtol=1e-3)
@test isapprox(v[3],  -0.065424, rtol=1e-3)





println("Testing Propogators")

# Vallado (2001), p.102-103, Example 2-4
r₀ = [1131.34, -2282.343, 6672.423]
v₀ = [-5.64305, 4.30333, 2.42879]
Δt = 40.0*60.0

r, v = kepler(r₀, v₀, Δt, μ=398600.4415)

@test isapprox(r[1], -4219.853, rtol=1e-4)
@test isapprox(r[2],  4363.116, rtol=1e-4)
@test isapprox(r[3], -3958.789, rtol=1e-4)
@test isapprox(v[1],  3.689736, rtol=1e-4)
@test isapprox(v[2], -1.916620, rtol=1e-4)
@test isapprox(v[3], -6.112528, rtol=1e-4)


r₀ /= 6378.1363
s   = √(6378.1363^3/398600.4415)
v₀ *= s/6378.1363
Δt /= s

r, v = kepler(r₀, v₀, Δt, μ=1.0)

@test isapprox(r[1], -0.6616125, rtol=1e-4)
@test isapprox(r[2],  0.6840739, rtol=1e-4)
@test isapprox(r[3], -0.6206809, rtol=1e-4)
@test isapprox(v[1],  0.466738,  rtol=1e-4)
@test isapprox(v[2], -0.2424455, rtol=1e-4)
@test isapprox(v[3], -0.7732126, rtol=1e-4)





println("Testing Coordinate Frame Conversions")

# Test ecl2EMsyn
JD = 2.45e6 + rand(-1000:1000)

rE = zeros(3)
rM = αδr2ijk(moon(JD))
rS = r_E_M*0.5*normalize(vec(rM))

rE = ecl2EMsyn(rE, JD=JD)
rS = ecl2EMsyn(rS, JD=JD)
rM = ecl2EMsyn(rM, JD=JD)

@test rE[1] < 0.0
@test -1.0 ≤ -rE[1]/R⨁ ≤ 1.0
@test isapprox(rE[2], 0.0, atol=1e-9)
@test isapprox(rE[3], 0.0, atol=1e-9)

@test rS[1] > 0.0
@test 0 ≤ rS[1] ≤ r_E_M
@test isapprox(rS[2], 0.0, atol=1e-9)
@test isapprox(rS[3], 0.0, atol=1e-9)

@test rM[1] > 0.0
@test 0.8 ≤ rM[1]/r_E_M ≤ 1.2
@test isapprox(rM[2], 0.0, atol=1e-9)
@test isapprox(rM[3], 0.0, atol=1e-9)

# Testing hci2ESsyn
JD = 2.45e6 + rand(-1000:1000)

rS = zeros(3)
rE = -αδr2ijk(sun(JD))
rA = km_in_AU*0.5*normalize(vec(rE))

rS = hci2ESsyn(rS, JD=JD)
rE = hci2ESsyn(rE, JD=JD)
rA = hci2ESsyn(rA, JD=JD)

@test rS[1] < 0.0
@test -1.0 ≤ -rS[1]/R☉ ≤ 1.0
@test isapprox(rS[2], 0.0, atol=1e-9)
@test isapprox(rS[3], 0.0, atol=1e-9)

@test rA[1] > 0.0
@test 0 ≤ rA[1] ≤ km_in_AU
@test isapprox(rA[2], 0.0, atol=1e-7)
@test isapprox(rA[3], 0.0, atol=1e-7)

@test rE[1] > 0.0
@test 0.9 ≤ rE[1]/km_in_AU ≤ 1.1
@test isapprox(rE[2], 0.0, atol=1e-7)
@test isapprox(rE[3], 0.0, atol=1e-7)




println("Testing SGP")
line1 = "1 88888U          80275.98708465  .00073094  13844-3  66816-4 0     8"
line2 = "2 88888  72.8435 115.9689 0086731  52.6988 110.5714 16.05824518  105 "

t = linspace(0, 1440*60, 5)
rTest, vTest = zeros(5, 3), zeros(5, 3)
for i = eachindex(t)
    rTest[i,:], vTest[i,:] = SGPwrap(line1, line2, t[i], R=6378.135, μ=398600.79964, J₂=1.082616e-3, J₃=-2.53881e-6)
end

rSoln = [2328.96594238 -5995.21600342 1719.97894287;
    	 2456.00610352 -6071.94232177 1222.95977784;
    	 2567.39477539 -6112.49725342  713.97710419;
    	 2663.03179932 -6115.37414551  195.73919105;
    	 2742.85470581 -6079.13580322 -328.86091614]

vSoln = [   2.91110113    -0.98164053   -7.09049922;
    	    2.67852119    -0.44705850   -7.22800565;
    	    2.43952477     0.09884824   -7.31889641;
    	    2.19531813     0.65333930   -7.36169147;
    	    1.94707947     1.21346101   -7.35499924]

for i = 1:5, j=1:3
    @test isapprox(rTest[i,j], rSoln[i,j], rtol=1e-4)
end


# Vallado (2001), Example 3-5, page 192.
JD = jday(1992, 8, 20, 12, 14, 0.0)
@test JD ≈ 2448855.009722

θ_GMST, θ_LST = LSTime(JD, λ=-deg2rad(104))
@test θ_GMST ≈ deg2rad(152.578787886)
@test θ_LST  ≈ deg2rad( 48.578787886)



# Vallado (2001), Example 3-14, page 223-4.
JD_UTC = jday(1991, 4, 6, 7, 51, 28.386009)

rfk5   =  [5102.5096, 6123.01152, 6378.1363]
vfk5   =  [-4.7432196, 0.7905366, 5.53375619]

# θGMST = LSTime(JD_UTC)[1]

recef  = eci2ecef(rfk5, JD=JD_UTC)#rot3(rfk5, θGMST)
vecef  = eci2ecef(vfk5, JD=JD_UTC) - ω⃗⨁ × recef #rot3(vfk5, θGMST)

@test isapprox(recef[1], -1120.598506, rtol=1e-1)
@test isapprox(recef[2],  7894.483204, rtol=1e-1)
@test isapprox(recef[3],  6374.079611, rtol=1e-1)

@test isapprox(vecef[1], -3.18701800 , rtol=1e-1)
@test isapprox(vecef[2], -2.90527125 , rtol=1e-1)
@test isapprox(vecef[3],  5.53765280 , rtol=1e-1)






println("All tests passed!")