export
    arrayInterp,
    rainRateCalc,
    rainAttenuation,
    specificRainAttenuation,
    cloudAttenuation,
    totalColumnWater,
    specificCloudAttenuation,
    cloudKl,
    scintAttenuation,
    specificAttenuationWet,
    specificAttenuationDry,
    waterVaporContent,
    waterVaporAttenuation,
    gasAttenuation,
    totalAttenuation,
    galacticNoise,
    brightTempSat,
    brightTempEarth,
    meanSurfaceTemp,
    antennaTempSat,
    antennaTempGround,
    rainHeightCalc,
    altCalc



#### ITU-R Rain Attenuation Model #### (ITU-R P.618-12)

const h_H = [ -5.33980 -0.10008 1.13098;
              -0.35351  1.26970 0.45400;
              -0.23789  0.86036 0.15354;
              -0.94158  0.64552 0.16817]
const m_k_H = -0.18961
const c_k_H =  0.71147
"""Get k_H coefficient for ITU-R P.838-3. f′ = log10(signal frequency [GHz])."""
getk_H(f′::Float64) = 10^(sum(h_H[:,1].*exp.(-((f′-h_H[:,2])./h_H[:,3]).^2)) + m_k_H*f′ + c_k_H)

const h_V = [ -3.80595  0.56934 0.81061;
              -3.44965 -0.22911 0.51059;
              -0.39902  0.73042 0.11899;
               0.50167  1.07319 0.27195]
const m_k_V = -0.16398
const c_k_V =  0.63297
"""Get k_V coefficient for ITU-R P.838-3. f′ = log10(signal frequency [GHz])."""
getk_V(f′::Float64) = 10^(sum(h_V[:,1].*exp.(-((f′-h_V[:,2])./h_V[:,3]).^2)) + m_k_V*f′ + c_k_V)

const α_H =[  -0.14318  1.82442 -0.55187;
               0.29591  0.77564  0.19822;
               0.32177  0.63773  0.13164;
              -5.37610 -0.96230  1.47828;
              16.1721  -3.29980  3.43990]
const m_α_H =  0.67849
const c_α_H = -1.95537
"""Get α_H coefficient for ITU-R P.838-3. f′ = log10(signal frequency [GHz])."""
getα_H(f′::Float64) = sum(α_H[:,1].*exp.(-((f′-α_H[:,2])./α_H[:,3]).^2)) + m_α_H*f′ + c_α_H

const α_V = [ -0.07771 2.33840 -0.76284 ;
               0.56727 0.95545  0.54039 ;
              -0.20238 1.14520  0.26809 ;
             -48.2991  0.791669 0.116226;
              48.5833  0.791459 0.116479]
const m_α_V = -0.053739
const c_α_V =  0.83433
"""Get α_V coefficient for ITU-R P.838-3. f′ = log10(signal frequency [GHz])."""
getα_V(f′::Float64) = sum(α_V[:,1].*exp.(-((f′-α_V[:,2])./α_V[:,3]).^2)) + m_α_V*f′ + c_α_V

# Loads all ESA data files specifying rainfall rate and rain height for any given latitude and longitude.
global rain_data = Dict{String,Any}()
function loadRainData(variable::String)
    global rain_data#, dName
    path = joinpath(datapath, "rain_data")#dName * "/data/rain_data"

    if !haskey(rain_data, variable)
        file = open(joinpath(path, "ESA$(uppercase(variable)).txt"))#$(datapath)/ESA$(uppercase(variable)).txt")
        rain_data[lowercase(variable)] = readdlm(file)
        close(file)
    end
    return rain_data[variable]
end


"""
    rainRateCalc(ϕ::T, λ::T; p::T=0.01) where T <: Float64

Return rainfall rate exceeded for p% of the year at earth station location.

### Arguments (Inputs):
| Parameter | Description             | Units |
| --------: | :---------------------- | :---- |
| ϕ         | Earth station latitude  | °     |
| λ         | Earth station longitude | °     |

### Values (Outputs):
| Parameter | Description                            | Units |
| --------: | :------------------------------------- | :---- |
| R_p       | Rain height at earth station location  | mm/hr |

Author: Sam Avery (2016), updated James Spicer (2018.03.07) <br/>
Source: ITU-R P.837-6. <br/>
Validation Examples: https://www.itu.int/en/ITU-R/study-groups/rsg3/Pages/iono-tropo-spheric.aspx  <br/>
MATLAB routines (in zip): https://www.itu.int/rec/R-REC-P.837-6-201202-S/en
"""
function rainRateCalc(ϕ::T, λ::T; p::T=0.01) where T <: Float64
    β    = arrayInterp(λ, ϕ, loadRainData("rain_lon")[1,:], loadRainData("rain_lat")[:,1], loadRainData("rain_beta"))
    M_T  = arrayInterp(λ, ϕ, loadRainData("rain_lon")[1,:], loadRainData("rain_lat")[:,1], loadRainData("rain_mt"))
    P_r6 = arrayInterp(λ, ϕ, loadRainData("rain_lon")[1,:], loadRainData("rain_lat")[:,1], loadRainData("rain_pr6"))

    M_c  = β*M_T
    M_s  = M_T - M_c

    P₀   = P_r6*(1-exp(-0.0079(M_s/P_r6)))

    a    = 1.09
    b    = M_T/(21797P₀)
    c    = 26.02b
    A    = a*b
    C    = log(p/P₀)
    B    = a + c*C

    return 0.5(-B + √(B^2 - 4A*C))/A
end


# Finds value of f(x,y) given x in array X, y in array Y, and f(X,Y) in array arr.
function arrayInterp(x::T, y::T, X::Array{T}, Y::Array{T}, arr::Array{T,2}) where T <: Float64
    C = boundInd(x, X)
    c = linearInterp(x, C, X)
    R = boundInd(y, Y)
    r = linearInterp(y, R, Y)
    return biLinearInterpSq(c, r, C, C+1, R, R+1, arr[R,C], arr[R,C+1], arr[R+1,C], arr[R+1,C+1])
end


"""
    rainHeightCalc(ϕ::T, λ::T) where T <: Float64

Return rain height at earth station location.

### Arguments (Inputs):
| Parameter | Description             | Units |
| --------: | :---------------------- | :---- |
| ϕ         | Earth station latitude  | °     |
| λ         | Earth station longitude | °     |

### Values (Outputs):
| Parameter | Description                            | Units |
| --------: | :------------------------------------- | :---- |
| H         | Rain height at earth station location  | km    |

Author: Sam Avery (2016), updated James Spicer (2018.03.06) <br/>
Source: ITU-R P.839-4.  
"""
rainHeightCalc(ϕ::T, λ::T) where T <: Float64 = 0.36 + arrayInterp(λ, ϕ, loadRainData("lon")[1,:], loadRainData("lat")[:,1], loadRainData("0height"))


"""
    specificRainAttenuation(f::T, R::T; el::T=90.0, pol::T=45.0) where T <: Float64

Return specific attenuation due to rainfall of a satellite ↔ earth station signal.

### Arguments (Inputs):
| Parameter | Description       | Units |
| --------: | :---------------- | :---- |
| f         | Frequency         | GHz   |
| R         | Rain rate         | mm/hr |

#### Optional Keyword Arguments:
| Parameter | Description                                                     | Default value | Units  |
| --------: | :-------------------------------------------------------------- | :------------ | :----- |
| el        | Signal elevation above earth station's local horizon            | 90.0          | °      |
| pol       | Polarization tilt angle relative to horizontal (45° = circular) | 45.0          | °      |

### Values (Outputs):
| Parameter | Description                                     | Units |
| --------: | :---------------------------------------------- | :---- |
| Y         | Specific attenuation of signal due to rainfall  | dB/km |

Author: Sam Avery (2016), updated James Spicer (2018.03.06) <br/>
Source: ITU-R P.838-3.  
"""
function specificRainAttenuation(f::T, R::T; el::T=90.0, pol::T=45.0) where T <: Float64
    # !(1.0 ≤ f ≤ 1e3) && warn("ITU-R P.838-3 only accurate for frequencies between 1 and 1,000 GHz.")
    f′  = log10(f)
    k_H = getk_H(f′)
    k_V = getk_V(f′)
    α_H = getα_H(f′)
    α_V = getα_V(f′)

    ce² = cosd(el)^2
    c2p = cosd(2pol)

    k = 0.5(k_H     + k_V     + (k_H     - k_V    )*ce²*c2p)
    α = 0.5(k_H*α_H + k_V*α_V + (k_H*α_H - k_V*α_V)*ce²*c2p)/k

    return k*R^α
end

# Check results here: https://logiciels.cnes.fr/sites/default/files/validation_examples.pdf
"""
    rainAttenuation(ϕ::T, λ::T, f::T; p::T=0.01, el::T=90.0, h_s::T=0.0, pol::T=45.0) where T <: Float64

Return attenuation due to rainfall of a satellite ↔ earth station signal.

### Arguments (Inputs):
| Parameter | Description             | Units |
| --------: | :---------------------- | :---- |
| ϕ         | Earth station latitude  | °     |
| λ         | Earth station longitude | °     |
| f         | Frequency               | GHz   |

#### Optional Keyword Arguments:
| Parameter | Description                                                     | Default value | Units  |
| --------: | :-------------------------------------------------------------- | :------------ | :----- |
| p         | Percentage of time result will be exceeded                      | 0.01          | %      |
| el        | Signal elevation above earth station's local horizon            | 90.0          | °      |
| h_s       | Earth station height above sea level                            | 0.0           | km     |
| pol       | Polarization tilt angle relative to horizontal (45° = circular) | 45.0          | °      |

### Values (Outputs):
| Parameter | Description                            | Units |
| --------: | :------------------------------------- | :---- |
| A_p       | Attenuation of signal due to rainfall  | dB    |

Author: Sam Avery (2016), updated James Spicer (2018.03.06) <br/>
Source: ITU-R P.618-12.  
"""
function rainAttenuation(ϕ::T, λ::T, f::T; p::T=0.01, el::T=90.0, h_s::T=0.0, pol::T=45.0) where T <: Float64 # p = percent
    # f > 55.0 && warn("ITU-R P.618-7 2.2.1.1 only accurate for frequencies below 55 GHz.")
    # if p < 0.001 || p > 5; warn("ITU-R P.618-7 2.2.1.1 only accurate for time percents between 0.001% and 5%."); end
    λ = putInRange(λ,   0, 360, delta=360)
    ϕ = putInRange(ϕ, -90,  90, delta=360)

    sel, cel = sind(el), cosd(el)

    # Step 1: Calculate the effective rain height h_R [km]
    h_R = rainHeightCalc(ϕ, λ)

    (h_R - h_s ≤ 0) && return 0.0

    # Step 2: Calculate the slant path length [km] below the freezing rain height
    Re = 8500 # effective radius of the Earth [km]
    L_s = (el ≥ 5) ? (h_R-h_s)/sel : 2(h_R-h_s) / (√(sel^2 + 2(h_R-h_s)/Re) + sel) # [km]

    # Step 3: Horizontal projection of slant path length
    L_G = L_s*cel

    # Step 4: Obtain the rainfall rate R_0.01 only exceeded for 0.01% of the average year (ITU)
    R0_01 = max(rainRateCalc(ϕ, λ), 0) # rainfall rate percent value [mm/hr]
    R0_01 == 0 && return 0.0

    # Step 5: Specific attenuation [dB/km] based on ITU-R P.838-3
    Y = specificRainAttenuation(f, R0_01, pol=pol, el=el)

    # Step 6: Calculate the horizontal adjustment factor
    r0_01 = inv(1 + 0.78*√(L_G*Y/f) - 0.38(1-exp(-2L_G)))

    # Step 7: Calculate the vertical adjustment factor
    ζ = atand((h_R - h_s)/(L_G*r0_01)) # degrees

    L_R = (ζ > el) ? L_G*r0_01/cel : (h_R-h_s)/sel

    χ = (abs(ϕ) < 36) ? 36-abs(ϕ) : 0.0

    v0_01 = inv(1 + √(sel)*(31*(1-exp(-el/(1+χ)))*√(L_R*Y)/(f^2)-0.45))

    # Step 8: Calculate the effective path length
    L_E = L_R*v0_01 # [km]

    # Step 9: Predicted attenuation exceeded for 0.01% of an average year
    A0_01 = Y*L_E # [dB]

    # Step 10: Estimated attenuation to be exceeded for other percentages of the year (0.001% to 5%)
    if p ≥ 1 || abs(ϕ) ≥ 36
        β = 0
    elseif p < 1 && abs(ϕ) < 36 && el ≥ 25
        β = -0.005(abs(ϕ)-36)
    else
        β = -0.005(abs(ϕ)-36) + 1.8 - 4.25sel
    end

    return A0_01*(100p)^(-(0.655 + 0.033log(p) - 0.045log(A0_01) - β*(1-p)*sel)) #[dB]
end





#### ITU-R Cloud & Fog Attenuation Model #### (ITU-R P.840-6)

# Load Cloud Data
# ---------------
# Inputs:
# ** variable: string name of the text file to be opened
#
# Outputs:
# ** Array from the input file
global cloud_data = Dict{String,Any}()
function loadCloudData(variable)
    global cloud_data#, dName
    path = joinpath(datapath, "cloud_data")#dName * "/data/cloud_data"

    if !haskey(cloud_data, variable)
        file = open(joinpath(path, "$(variable).txt"))#$(datapath)/$(variable).txt")
        cloud_data[variable] = readdlm(file)
        close(file)
    end
    return cloud_data[variable]
end

# Specific Attenuation Coefficient (ITU-R P.840-6 2)
# --------------------------------
# Inputs:
# ** frequency [GHz] (<1000 GHz)
# ** Cloud temperature [K]
# Outputs:
# ** K_l: Specific attenuation coefficient
function cloudKl(f, T)
    # f > 1e3 && warn("Rayleigh scattering model only valid for frequencies below 1000 GHz.")
    th  = 300/T
    e₀  = 77.66 + 103.3(th-1)
    e₁  = 0.0671e₀
    e₂  = 3.52

    f_p = 20.2 - 146(th-1) + 316(th-1)^2 # GHz
    f_s = 39.8f_p # GHz

    # Complex dielectric permittivity of water
    d_ε =    (e₀ - e₁)/     (1 + (f/f_p)^2)  +   (e₁ - e₂)/     (1 + (f/f_s)^2) + e₂
    d2_ε = f*(e₀ - e₁)/(f_p*(1 + (f/f_p)^2)) + f*(e₁ - e₂)/(f_s*(1 + (f/f_s)^2))

    η = (2 + d_ε)/d2_ε

    K_l = 0.819f/(d2_ε*(1+η^2)) # [dB/km * g/m3]
end

# Cloud Specific Attenuation (ITU-R P.840-6 1)
# --------------------------
# Inputs:
# ** frequency [GHz] (<200 GHz)
# ** Cloud temperature [K]
# ** M = Cloud liquid water density [g/m^3] (0.05 - 0.5 g/m^3 medium to thick fog)
# Outputs:
# ** gamma_c: Cloud specific attenuation [dB/km]
function specificCloudAttenuation(f, T, M)
    # f > 200 && warn("Rayleigh approximation is only valid for frequencies below 200 GHz.")
    cloudKl(f, T)*M # [dB/km]
end

# Cloud Attenuation
# -----------------
# Calculates the total attenuation through the cloud layer with probability p.
# Inputs:
# ** frequency [GHz] (<200 GHz)
# ** Cloud temperature [K]
# ** el: Elevation angle [deg] (90 >= el >= 5 deg)
# ** L: Total columnar content of liquid [kg/m^2].
# Outputs:
# ** A: Cloud attenuation [dB]
function cloudAttenuation(ϕ, λ, f; p=1, el=90, T=273.15)
    # f > 200        && warn("ITU-R P.840-6 only valid for frequencies below 200 GHz.")
    # !(5 ≤ el ≤ 90) && warn("ITU-R P.840-6 only valid for elevation angles between 90° and 5°.")
    λ = putInRange(λ,   0, 360, delta=360)
    ϕ = putInRange(ϕ, -90,  90, delta=360)

    L = totalColumnWater(ϕ, λ, p=p)
    return L*cloudKl(f, T)/sind(el)
end

# Total Columnar Content of Liquid Water
# --------------------------------------
# Inputs:
# ** p: Desired probability percent [%]
# Outputs:
# ** L: Total columnar content of reduced cloud liquid water [kg/m^2]
function totalColumnWater(ϕ, λ; p=1)
    # !(0.1 ≤ p ≤ 99) && warn("ITU-R P.840-6 only valid for probabilities between 0.1% and 99%.")
    λ = putInRange(λ,   0, 360, delta=360)
    ϕ = putInRange(ϕ, -90,  90, delta=360)

    cloudLat = loadCloudData("ESALAT1dot125")
    cloudLon = loadCloudData("ESALON1dot125")

    pSet = [0.1, 0.2, 0.3, 0.5, 1, 2, 3, 5, 10, 20, 30, 50, 60, 70, 80, 90, 95, 99] # Annual
    ind = boundInd(p, pSet)
    pBelow = pSet[ind]
    pAbove = pSet[ind+1]

    L = Float64[]
    for q in (pAbove, pBelow)
        str = (q < 1) ? replace(string(q), ".", "") : string(Int(q))
        L_p = loadCloudData("Lred_" * str * "_v4")
        l = arrayInterp(λ, ϕ, cloudLon[1,:], cloudLat[:,1], L_p)
        push!(L, l)
    end

    return L[1] + (L[2]-L[1])*log(p/pBelow)/log(pAbove/pBelow) # [kg/m2]
end

## Scintillation Loss
# Saturation Vapor Pressure (ITU-R P.453)
# -------------------------
# Inputs:
# **temp: Temperature (deg C) measured at ground level for 1 month or longer
#     Water: -40 < temp < 50 °C     Ice:  -80 < temp < 0 °C
# **press: Pressure [hPa]
# **isWater: boolean of water (true) or ice (false)
function satVapPress(T, P; isWater=true)
    if isWater # Liquid water - valid from -40 < T < 50 °C
        EF =   1.0 + 1e-4(7.2 + P*(3.20e-3 + 5.9e-7T^2))
        a  =   6.1121
        b  =  18.678
        c  = 257.14
        d  = 234.5
    else         # Ice - valid from -80 < T < 0 °C
        EF =   1.0 + 1e-4(2.2 + P*(3.82e-3 + 6.4e-7T^2))
        a  =   6.1115
        b  =  23.036
        c  = 279.82
        d  = 333.7
    end
    e_s = EF*a*exp((b-T/d)*T/(T+c))
end

# Radio Refractivity Wet Term (ITU-R P.453)
# ---------------------------
# Inputs:
# **temp: Temperature (deg C) measured at ground level for 1 month or longer
#     Water: -40 < temp < 50 °C     Ice:  -80 < temp < 0 °C
# **press: Pressure [hPa]
# **relHum: Average surface relative humidity (%) at the site for a period of 1 month or longer
# **isWater: boolean of water (true) or ice (false)
# Outputs:
# N_wet: radio refractivity wet term [ppm]
function refractWet(T, P, relHum; isWater=true)
    e_s   = satVapPress(T, P, isWater=isWater)
    T    += 273.15 # temp [K]
    e     = 0.01relHum*e_s
    N_wet = 72e/T + 3.75e5e/T^2
end

# Load Refractivity Data
# ----------------------
# Inputs:
# ** variable: string name of the text file to be opened
# Outputs:
# ** Array from the input file
global refract_data = Dict{String,Any}()
function loadRefractData(variable)
    global refract_data#, dName
    path = joinpath(datapath, "refract_data")#dName * "/data/refract_data"

    if !haskey(refract_data, variable)
        file = open(joinpath(path, "ESA$(variable).txt"))#$(datapath)/ESA$(variable).txt")
        refract_data[variable] = readdlm(file)
        close(file)
    end
    return refract_data[variable]
end

# Radio Refractivity Wet Term Data
# Inputs:
# **p: Desired probability percent [%]
# **ϕ: latitude [deg] (-90° < lat < +90°)
# **λ: longitude [deg] (0° < lon < 360°)
# Outputs:
# **N_wet: radio refractivity wet term [ppm]
refractWetData(ϕ, λ) = arrayInterp(λ, ϕ, loadRefractData("LON")[1,:], loadRefractData("LAT")[:,1], loadRefractData("NWET"))

# Scintillation Loss (ITU-R P.618-12 2.4.1)
# Inputs:
# **frequency [GHz] (4≤f≤20 GHz)
# **refractWet: radio refractivity wet term [ppm]
# **el: Elevation angle [deg] Free-space elevation angle (θ >= 5°)
# **p: Time percent in a year [%] 0.01% < p <= 50%
# **D: Antenna Diameter [m]
# **n: Antenna Efficiency
# Output:
# **A_p: Scintillation Loss [dB]
function scintAttenuation(ϕ, λ, f; p=1, el=90, D=1, n=0.5)
    # if !(4 ≤ f ≤ 20); warn("ITU-R P.618-12 2.4.1 only valid for frequencies between 4 and 20 GHz."); end
    # el < 5 && warn("ITU-R P.618-12 2.4.1 only valid for elevation angles greater than 5°.")
    λ = putInRange(λ,   0, 360, delta=360)
    ϕ = putInRange(ϕ, -90,  90, delta=360)

    # Step 1: Calculate the saturation water vapor pressure
    ## Calculated in N_wet

    # Step 2: Calculate radio refractivity N_wet
#     N_wet = refractWet
    N_wet = refractWetData(ϕ, λ)

    # Step 3: Calculate standard deviation of reference signal amplitude
    σ_ref = 3.6e-3 + 1e-4N_wet # [dB]

    # Step 4: Calculate the effective path length, L
    h_L = 1000.0 # Height [m] of the turbulence layer (by ITU)
    L = 2h_L/(sind(el) + √(sind(el)^2 + 2.35e-4)) # Length of the path through the turbulence [m]

    # Step 5: Estimate the effective antenna diameter
    D_eff = D*√n # Effective antenna diameter [m]

    # Step 6: Calculate the antenna averaging factor
    x = 1.22D_eff^2*f/L

    x ≥ 7 && return 0.0

    g = √(3.86(x^2 + 1)^(11/12)*sin(11acot(x)/6) - 7.08x^(5/6))

    # Step 7: Calculate the standard deviation of the signal for the applicable period and propagation path
    σ = σ_ref*f^(7/12)*g/(sind(el)^1.2)

    # Step 8: Calculate the time percentage factor for 0.01<p<=50 %
    a_p = -0.061log10(p)^3 + 0.072log10(p)^2 - 1.71log10(p) + 3

    # Step 9: Calculate the fade depth exceeded for p% of the time
    A_p = max(a_p*σ, 0.0)
end


phiCalc(r_p, r_t, a, b, c, d) = (r_p^a)*(r_t^b)*exp(c*(1-r_p) + d*(1-r_t)) # Helper function for specificAttenuationDry

# Dry Air Specific Attenuation (ITU-R P.676-10, Annex 2, Sec. 1)
# ----------------------------
# Inputs:
# ** f [GHz] (1 - 350 GHz)
# ** T: Temperature [°C] Mean temperature values can be obtained from maps
# given in Recommendation ITU-R P.1510, when no adequate temperature data are available
# ** p_tot: Total Air Pressure [hPa]
#
# Outputs:
# ** g_0: dry air specific attenuation [dB/km]
#
# Note: Only applies up to an altitude of 10 km (assuming this is ground station altitude)
# ACCURATE TO ITU VALIDATION
function specificAttenuationDry(f, T, p_tot)
    # !(1 ≤ f ≤ 350) && warn("ITU-R P.676-10, Annex 2, Sec. 1 only valid for frequencies between 1 and 350 GHz.")
    r_p = p_tot/1013.0
    r_t = 288.0/(273.0+T) # 273.15?

    e_1 =         phiCalc(r_p, r_t,  0.0717, -1.8132,  0.0156, -1.6515)
    e_2 =         phiCalc(r_p, r_t,  0.5146, -4.6368, -0.1921, -5.7416)
    e_3 =         phiCalc(r_p, r_t,  0.3414, -6.5851,  0.2130, -8.5854)
    e_4 =         phiCalc(r_p, r_t, -0.0112,  0.0092, -0.1033, -0.0009)
    e_5 =         phiCalc(r_p, r_t,  0.2705, -2.7192, -0.3016, -4.1033)
    e_6 =         phiCalc(r_p, r_t,  0.2445, -5.9191,  0.0422, -8.0719)
    e_7 =         phiCalc(r_p, r_t, -0.1833,  6.5589, -0.2402,  6.131 )

    g_54 =  2.192*phiCalc(r_p, r_t,  1.8286, -1.9487,  0.4051, -2.8509)
    g_58 = 12.59 *phiCalc(r_p, r_t,  1.0045,  3.5610,  0.1588,  1.2834)
    g_60 = 15.00 *phiCalc(r_p, r_t,  0.9003,  4.1335,  0.0427,  1.6088)
    g_62 = 14.28 *phiCalc(r_p, r_t,  0.9886,  3.4176,  0.1827,  1.3429)
    g_64 =  6.819*phiCalc(r_p, r_t,  1.4320,  0.6258,  0.3177, -0.5914)
    g_66 =  1.908*phiCalc(r_p, r_t,  2.0717, -4.1404,  0.4910, -4.8718)

    δ  = -3.06e-3*phiCalc(r_p, r_t,  3.211, -14.94,    1.583, -16.37  )

    f ≤ 54  && return (7.2(r_t^2.8)/(f^2+0.34(r_p^2)*(r_t^1.6)) + 0.62*e_3/((54-f)^(1.16*e_1) + 0.83*e_2))*(f^2)*(r_p^2)*1e-3
    f ≤ 60  && return exp(log(g_54)/24.0*(f-58.0)*(f-60.0) - log(g_58)/8.0*(f-54.0)*(f-60.0) + log(g_60)/12.0*(f-54.0)*(f-58.0))
    f ≤ 62  && return g_60 + (g_62-g_60)*(f-60.0)*0.5
    f ≤ 66  && return exp(log(g_62)/8.0*(f-64.0)*(f-66.0) - log(g_64)/4.0*(f-62.0)*(f-66.0) + log(g_66)/8.0*(f-62.0)*(f-64.0))
    f ≤ 120 && return (3.02e-4*r_t^3.5 + 0.283*r_t^3.8/((f-118.75)^2 + 2.91*r_p^2*r_t^1.6) + 0.502*e_6*(1.0-0.0163*e_7*(f-66))/((f-66)^(1.4346*e_4) + 1.15*e_5))*f^2*r_p^2*1e-3
               return (3.02e-4/(1.0+1.9e-5*f^1.5) + 0.283*r_t^0.3/((f-118.75)^2 + 2.91*r_p^2*r_t^1.6))*f^2*r_p^2*r_t^3.5*1e-3 + δ
end

gCalc(f, f_i) = 1 + ((f-f_i)/(f+f_i))^2 # Helper function for specificAttenuationWet

# Water Vapor Specific Attenuation (ITU-R P.676-10, Annex 2, Sec. 1)
# --------------------------------
# Inputs:
# ** frequency [GHz] (1 - 350 GHz)
# ** T: Temperature [°C]
# ** p_tot: Total Air Pressure [hPa]
# ** ρ: Water-vapor density [g/m^3]
#
# Outputs:
# ** g_w: water vapor specific attenuation [dB/km]
#
# Note: Only applies up to an altitude of 10 km (assuming this is ground station altitude)
# ACCURATE TO ITU VALIDATION
function specificAttenuationWet(f, T, p_tot, ρ)
    # !(1 ≤ f ≤ 350) && warn("ITU-R P.676-10, Annex 2, Sec. 1 only valid for frequencies between 1 and 350 GHz.")
    r_p = p_tot/1013
    r_t = 288/(273+T) # 273.15?

    n₁ = 0.955r_p*r_t^0.68 + 0.006ρ
    n₂ = 0.735r_p*r_t^0.5  + 0.0353ρ*r_t^4

    g_w = ( 3.98    *n₁*exp(2.23(1-r_t))/((f- 22.235)^2 +  9.42n₁^2)*gCalc(f, 22)
        +  11.96    *n₁*exp(0.7*(1-r_t))/((f-183.31 )^2 + 11.14n₁^2)
        +   0.081   *n₁*exp(6.44(1-r_t))/((f-321.226)^2 +  6.29n₁^2)
        +   3.66    *n₁*exp(1.6*(1-r_t))/((f-325.153)^2 +  9.22n₁^2)
        +  25.37    *n₁*exp(1.09(1-r_t))/((f-380    )^2            )
        +  17.4     *n₁*exp(1.46(1-r_t))/((f-448    )^2            )
        + 844.6     *n₁*exp(0.17(1-r_t))/((f-557    )^2            )*gCalc(f, 557)
        + 290       *n₁*exp(0.41(1-r_t))/((f-752    )^2            )*gCalc(f, 752)
        +   8.3328e4*n₂*exp(0.99(1-r_t))/((f-1780   )^2            )*gCalc(f, 1780))*f^2*r_t^2.5*ρ*1e-4
end

# Approximate Gaseous Attenuation (ITU-R P.676-10, Annex 2, Sec. 2.2)
# -------------------------------
# Inputs:
# ** f: [GHz] (1 - 350 GHz)
# ** T: Temperature [°C]
# ** p_tot: Total Air Pressure [hPa]
# ** rho: Water-vapor density [g/m^3]
# ** ϕ: latitude [deg] (+90° to -90°)
# ** λ: longitude [deg] (0° to 360°)
# ** alt: Altitude of the ground station (?) [km]
# ** el: Elevation Angle [degrees] (min of 5°)
# ** p: Time Percent [%]
#
# Outputs:
# ** A_g: total approximate atmospheric gaseous attenuation [dB]
#
# ACCURATE TO ITU MODEL
function gasAttenuation(ϕ, λ, f; T=0, p_tot=1013, rho=7.5, alt=0, el=90, p=1)
    # (alt > 10)     && warn("ITU-R P.676-10 Annex 2 Sec. 1 only valid for ground altitudes below 10 km.")
    # !(5 ≤ el ≤ 90) && warn("ITU-R P.676-10 Annex 2 Sec. 2.2.1.1 only valid for elevation angles between 90° and 5°.")
    λ = putInRange(λ,   0, 360, delta=360)
    ϕ = putInRange(ϕ, -90,  90, delta=360)

    r_p = p_tot/1013.0 # p_tot is the total air pressure
    # r_t = 288.0/(273.0+T) # <- not sure this is needed anywhere

    t₁ = 4.64/(1+0.066(r_p^-2.3))*exp(-((f-59.7)/(2.87 + 12.4*exp(-7.9r_p)))^2)
    t₂ = 0.14*exp(2.12r_p)/((f-118.75)^2 + 0.031*exp(2.2r_p))
    t₃ = 0.0114/(1+0.14r_p^(-2.6))*f*(-0.0247 + 0.0001f + 1.61e-6f^2)/(1-0.0169f + 4.1e-5f^2 + 3.2e-7f^3)

    h₀ = 6.1/(1+0.17r_p^-1.1)*(1+t₁+t₂+t₃)

    (f < 70) && (h₀ = min(h₀, 10.7r_p^0.3))

    s_w = 1.013/(1 + exp(-8.6(r_p-0.57)))

    h_w = 1.66(1 + 1.39s_w/((f-22.235)^2 + 2.56s_w) +
                   3.37s_w/((f-183.31)^2 + 4.69s_w) +
                   1.58s_w/((f-325.1 )^2 + 2.89s_w))

    g_0 = specificAttenuationDry(f, T, p_tot)
#     g_w = specificAttenuationWet(f, T, p_tot, rho)

#     A_g = r0*(g_0 + g_w) # [dB] Path attenuation for horizontal, terrestrial paths. r0 = path length [km]
#     A_g = (g_0*h₀ + g_w*h_w)/sind(el) # [dB] for path attenuation based on surface meteorological data

    A_g = (h₀*g_0 + waterVaporAttenuation(ϕ, λ, f, alt=alt, p=p))/sind(el) # [dB] for path attenuation based on integrated water vapour content
end

# Zenith Path Water Vapor Attenuation (ITU-R P.676-10 Annex 2, Sec. 2.3)
# -----------------------------------
# Inputs:
# ** f: frequency [GHz] (1 - 350 GHz)
# ** ϕ: latitude [deg] (+90° to -90°)
# ** λ: longitude [deg] (0° to 360°)
# ** alt: Altitude of the ground station (?) [km]
# ** p: Time Percent [%]
#
# Outputs:
# ** A_w: Approximate water vapor attenuation [dB]
#
# ACCURATE TO ITU MODEL
function waterVaporAttenuation(ϕ, λ, f; alt=0, p=1) # [GHz]
    # alt > 10 && warn("ITU-R P.676-10 Annex 2 Sec. 2.2 only valid for ground altitudes below 10 km.")
    λ = putInRange(λ,   0, 360, delta=360)
    ϕ = putInRange(ϕ, -90,  90, delta=360)

    f_ref = 20.6 # reference frequency [GHz]

    V_t = waterVaporContent(ϕ, λ, alt=alt, p=p) # [kg/m2 or mm]
    t_ref = 14log(0.055V_t) + 3 #temperature (°C)
    p_ref = 780 # p_tot is the total air pressure [hPa]
    ρ_ref = 0.25V_t # water-vapour density [g/m^3]

    A_w = 0.0173V_t*specificAttenuationWet(f, t_ref, p_ref, ρ_ref)/specificAttenuationWet(f_ref, t_ref, p_ref, ρ_ref)
end

# Load Water Vapor Data
# ---------------------
# Inputs:
# ** variable: string name of the text file to be opened
# Outputs:
# ** Array from the input file
global vapor_data = Dict{String,Any}()
function loadVaporData(variable)
    global vapor_data#, dName
    path = joinpath(datapath, "vapor_data")#dName * "/data/vapor_data"

    if !haskey(vapor_data, variable)
        file = open(joinpath(path, "$(variable).txt"))#$(datapath)/$(variable).txt")
        vapor_data[variable] = readdlm(file)
        close(file)
    end
    return vapor_data[variable]
end


# Water Vapor Content (ITU-R P.836-5, Annex 2)
# -------------------
# Inputs:
# ** ϕ: latitude [deg] (+90° to -90°)
# ** λ: longitude [deg] (0° to 360°)
# ** alt: Altitude of the ground station (?) [km]
# ** p: Time Percent [%]
#
# Outputs:
# ** V: Water vapor content [g/m^3 ????] (listed as kg/m^2 and g/m^3)
#
# Note: Only applies up to an altitude of 10 km (assuming this is ground station altitude)
# ACCURATE TO ITU MODEL
function waterVaporContent(ϕ, λ; alt=0, p=1) # from ITU-R P.836-5, Annex 2
    λ = putInRange(λ,   0, 360, delta=360)
    ϕ = putInRange(ϕ, -90,  90, delta=360)

    vLat = loadVaporData("LAT1dot125")
    vLon = loadVaporData("LON1dot125")

    C = boundInd(λ, vLon[1,:])
    R = boundInd(ϕ, vLat[:,1])
    c = linearInterp(λ, C, vLon[1,:])
    r = linearInterp(ϕ, R, vLat[:,1])

    pSet   = [0.1, 0.2, 0.3, 0.5, 1, 2, 3, 5, 10, 20, 30, 50, 60, 70, 80, 90, 95, 99] # Annual
    ind    = boundInd(p, pSet)
    pBelow = pSet[ind]
    pAbove = pSet[ind+1]

    V = Float64[]
    for q in (pAbove, pBelow)
        str = (q < 1.0) ? replace(string(q), ".", "") : string(Int(q))
        V_p = loadVaporData("V_" * str * "_v4")
        vsch_p = loadVaporData("VSCH_" * str * "_v4")

        vis = Float64[]
        for (P,Q) in [(R,C), (R+1,C), (R,C+1), (R+1,C+1)] # 4 closest points. Maybe switch 2 & 3
            vi    = V_p[P,Q]
            vschi = vsch_p[P,Q]
            alti  = altCalc(vLat[P,Q], vLon[P,Q])
            vi   *= exp(-(alt-alti)/vschi)
            push!(vis, vi)
        end

        v = biLinearInterpSq(r, c, R, R+1, C, C+1, vis[1], vis[2], vis[3], vis[4])
        push!(V, v)
    end
    V[1] + (V[2]-V[1])*log(p/pBelow)/log(pAbove/pBelow) # [kg/m2] or [g/m2]
end

# Load Altitude Data
# ------------------
# Inputs:
# ** variable: string name of the text file to be opened
# Outputs:
# ** Array from the input file
global altitude_data = Dict{String,Any}()
function loadAltitudeData(variable)
    global altitude_data#, dName
    path = joinpath(datapath, "altitude_data")#dName * "/data/altitude_data"

    if !haskey(altitude_data, variable)
        file = open(joinpath(path, "$(variable).txt"))#$(datapath)/$(variable).txt")
        altitude_data[variable] = readdlm(file)
        close(file)
    end
    return altitude_data[variable]
end

# Altitude Calculation (ITU-R P.1511-1)
# --------------------
# Inputs:
# ** ϕ: latitude [deg] (+90° to -90°)
# ** λ: longitude [deg] (0° to 360°)
#
# Outputs:
# ** altitude: topographical altitude to 0.5 degree [degrees]
#
# ACCURATE TO ITU MODEL
function altCalc(ϕ, λ)
    λ = putInRange(λ,   0, 360, delta=360)
    ϕ = putInRange(ϕ, -90,  90, delta=360)

    altLat = loadAltitudeData("TOPOLAT")[:,1]
    altLon = loadAltitudeData("TOPOLON")[2,:]
    alt = loadAltitudeData("TOPO_0DOT5")

    C = boundInd(λ, altLon)
    c = linearInterp(λ, C, altLon)
    R = boundInd(ϕ, altLat)
    r = linearInterp(ϕ, R, altLat)
    biCubicInterp(r, c, R-1, C-1, alt) # altitude at (r,c) [km]
end


# Estimation of total attenuation (ITU-R P.618-12, Sec. 2.5) [dB]
# Inputs: ϕ [deg], λ [deg], f [GHz], p [%], el [deg] (elevation angle),
# D [m] (Rx antenna diameter), n [-] (Rx antenna efficiency), h_s [km]
function totalAttenuation(ϕ, λ, f; p=1.0, el=90.0, D=1.0, n=0.5, h_s=0.0, T=273.15, pol=45.0)
    # if !(0.001 ≤ p ≤ 50); warn("ITU-R P.618-12 Sec. 2.5 only accurate for time percents between 0.001% and 50%"); end
    p = min(p, 98.9)
    p = max(p, 0.1)
    λ = putInRange(λ,   0, 360, delta=360)
    ϕ = putInRange(ϕ, -90,  90, delta=360)

    A_G = gasAttenuation(  ϕ, λ, f, p=p,         el=el, T=T, alt=h_s)
    A_R = rainAttenuation( ϕ, λ, f, p=p,         el=el, h_s=h_s, pol=pol)
    A_C = cloudAttenuation(ϕ, λ, f, p=max(p, 1), el=el, T=T)
    A_S = scintAttenuation(ϕ, λ, f, p=p,         el=el, D=D, n=n)
    return A_G + hypot(A_R + A_C, A_S) # [dB]
end



# Brightness Temperature for Ground to Satellite
# ----------------------------------------------
# ITU-R P.372-12 Sec. 4.1
#
# Inputs:
# ** ϕ: latitude [deg] (+90° to -90°)
# ** λ: longitude [deg] (0° to 360°)
# ** f: frequency [GHz] (0.1 - 100 GHz)
#
# Outputs:
# ** T_b: Earth's brightness temperature [K]
function brightTempEarth(ϕ, λ, f) # Estimates the brightness temperature in the direction of propagation.
    λ    = putInRange(λ,   0, 360, delta=360)
    ϕ    = putInRange(ϕ, -90,  90, delta=360)

    A    = totalAttenuation(ϕ, λ, f)
    T_s  = meanSurfaceTemp(ϕ, λ)                 # [K] Surface temperature at the ground site
    T_mr = 37.34 + 0.81T_s                       # [K] atmospheric mean radiating temperature
    T_b  = T_mr*(1 - db2pow(-A)) + 2.7db2pow(-A)
end

# Load Earth Surface Temperature Data
# -----------------------------------
# Inputs:
# ** variable: string name of the text file to be opened
# Outputs:
# ** Array from the input file
#
global temp_data = Dict{String,Any}()
function loadTempData(variable)
    global temp_data, dName
    path = joinpath(datapath, "temp_data")#dName * "/data/temp_data"

    if !haskey(temp_data, variable)
        file = open(joinpath(path, "ESA$(uppercase(variable)).txt"))#$(datapath)/ESA$(uppercase(variable)).txt")
        temp_data[variable] = readdlm(file)
        close(file)
    end
    return temp_data[variable]
end

# Earth Mean Surface Temperature
# ------------------------------
# Inputs:
# ** lat: latitude [deg] (+90° to -90°)
# ** λ: longitude [deg] (0° to 360°)
#
# Outputs:
# ** T: Earth's mean surface temperature [K]
meanSurfaceTemp(ϕ, λ) = arrayInterp(λ, ϕ, loadTempData("lon")[1,:], loadTempData("lat")[:,1], loadTempData("temp"))

# Brightness Temperature for Satellite to Ground
# ----------------------------------------------
# ITU-R P.372-12 Sec. 4.2
#
# Inputs:
# ** ϕ: latitude [deg] (+90° to -90°)
# ** λ: longitude [deg] (0° to 360°)
# ** f: frequency [GHz] ((0.1 - 100 GHz)
# ** em: emissivity (~0.3, typically smaller than reflectivity)
# ** p: reflectivity (~0.7, larger than emissivity)
#
# Outputs:
# ** T: brightness temperature from the Earth to the Satellite [K]
function brightTempSat(ϕ, λ, f; em=0.3, p=0.7)
    # This calculation involves integration of downwelling radiation
    # over all angles and includes atmospheric attenuation.
    λ = putInRange(λ,   0, 360, delta=360)
    ϕ = putInRange(ϕ, -90,  90, delta=360)
    T_surf = meanSurfaceTemp(ϕ, λ)          # [K] physical temperature of the Earth's surface
    while isnan(T_surf)
        ϕ += 0.01
        λ += 0.01
        T_surf = meanSurfaceTemp(ϕ, λ)
    end
    T_atm = brightTempEarth(ϕ, λ, f)        # [K] weighted average of the sky brightness temperature
    T = em*T_surf + p*T_atm                 # [K] simplified brightness tem
end

# Brightness Temperature for Galactic Sources
# -------------------------------------------
# ITU-R P.372-12 Sec. 6
#
# Inputs:
# ** f: frequency [GHz] (0.1 - 100 GHz)
# ** sun: bool indicating if the antenna is pointing at the Sun
# ** moon_new: bool indicating if the antenna is pointing at a new moon
# ** moon_full: bool indicating if the antenna is pointing at a full moon
#
# Outputs:
# ** T: brightness temperature of the galactic background [K]
function galacticNoise(f; sun=false, moon_new=false, moon_full=false)
    T_b_0 = 200         # [K] reference temperature
    f₀    = 0.408       # [GHz] reference frequency

    sun       && return 1e4      # [K] under quiet sun conditions at 10 GHz
    moon_new  && return 140      # [K] brightness temperature of new Moon
    moon_full && return 280      # [K] brightness temperature of full Moon
                 return T_b_0*(f/f₀)^(-2.75) + 2.7  # Galactic background temp ~ 2.7 K
end

# Latitude and Longitude Calculator
# ---------------------------------
# Calculates the latitude and longitude on the ground given a satellite at
# a starting latitude and longitude with a given beam direction.
#
# Inputs:
# ** latInit: initial latitude of the satellite [deg]
# ** lonInit: initial longitude of the satellite [deg]
# ** θ: rotation angle of a beam [rad]
# ** ϕ: offset angle from perpendicular face of satellite to the ground [rad]
# ** reSat: radius of the satellite orbit [km]
#
# Outputs:
# ** lat: adjusted latitude [deg]
# ** lon: adjusted longitude [deg]
function latLonCalc(latInit, lonInit, θ, ϕ; reSat=26562)
    η   = (180 - asind(reSat/R⨁*sin(ϕ)))
    lat = (180 - ϕ - η)*cos(θ) + latInit
    lon = (180 - ϕ - η)*sin(θ) + lonInit

    lat = max(lat, -88.5)
    lat = min(lat,  88.5)
    lon = putInRange(lon, 0, 360, delta=360)
    lon = max(lon, 1.5)
    lon = min(lon, 358)

    return lat, lon
end

# Brightness Temperature for Satellite Calculator
# -----------------------------------------------
# Helper function to Satellite Antenna Temperature that ensures the temperature output
# is a real number value.
#
# Inputs:
# ** lat: latitude [deg] (+90° to -90°)
# ** lon: longitude [deg] (0° to 360°)
# ** f: frequency [GHz]
#
# Outputs:
# ** TB: satellite estimated brightness temperature [K]
function brightTempSatCalc(lat, lon, f)
    TB = brightTempSat(lat, lon, f)
    while isnan(TB)
        lat += 0.1
        lon += 0.1
        TB = brightTempSat(lat, lon, f)
    end
    return TB
end

# Brightness Temperature for Ground Calculator
# --------------------------------------------
# Helper function to Ground Antenna Temperature that ensures the temperature output
# is a real number value.
#
# Inputs:
# ** lat: latitude [deg] (+90° to -90°)
# ** lon: longitude [deg] (0° to 360°)
# ** f: frequency [GHz] (0.1 - 100 GHz)
#
# Outputs:
# ** TB: ground estimated brightness temperature [K]
function brightTempGroundCalc(lat, lon, f)
    TB = brightTempEarth(lat, lon, f)
    while isnan(TB)
        lat += 0.1
        lon += 0.1
        TB = brightTempEarth(lat, lon, f)
    end
    return TB
end

# Satellite Antenna Temperature
# -----------------------------
# Antenna temperature integration at the satellite.
# Antenna Theory: Analysis and Design, Constantine A. Balanis Eqn. 2-138
#
# Inputs:
# ** offsetAngle: Angle offset from the center of the beam [rad]
# ** f: frequency [GHz] (0.1 - 100 GHz)
# ** Gmax: maximum gain of the antenna [dB]
# ** lat: latitude of the satellite [deg]
# ** lon: longitude of the satellite [deg]
# ** n: number of integration steps
# ** reSat: radius of the satellite [km]
# ** TB_space: brightness temperature of the galactic background [K]
#
# Output:
# ** TA_num/TA_den: averaged temperature at the satellite antenna [K]
function antennaTempSat(latSat, lonSat, f, offsetAngle, Gmax; n=10, reSat=26562, TB_space=65)
    L  = -30 # ITU antenna model
    BW = BWfromG(Gmax)
    Θ  = linspace(0, π, n)
    Φ  = linspace(0, π, n)

    minEarthAngle = 0.5 - asin(R⨁/reSat)
    maxEarthAngle = 0.5 + asin(R⨁/reSat)

    TA_num, TA_den = 0, 0
    for θ in Θ, ϕ in Φ
        if minEarthAngle < ϕ < maxEarthAngle
            lat, lon = latLonCalc(latSat, lonSat, θ, 0.5π - ϕ)
            TB = brightTempSat(lat, lon, f)
            G = ITUrefRadPat672(Gmax, L, BW, abs(ϕ-offsetAngle-0.5π))
        else
            TB = TB_space
            G = max(ITUrefRadPat672(Gmax, L, BW, abs(ϕ-offsetAngle-0.5π)), 0)
        end
        TA_num += TB*G*sin(θ)
        TA_den += G*sin(θ)
    end
    return TA_num / TA_den
end

# Ground Antenna Temperature
# --------------------------
# Antenna temperature integration at the ground station.
# Antenna Theory: Analysis and Design, Constantine A. Balanis Eqn. 2-138
#
# Inputs:
# ** offsetAngle: Angle offset from the center of the beam [rad]
# ** f: frequency [GHz] (0.1 - 100 GHz)
# ** Gmax: maximum gain of the antenna [dB]
# ** lat: latitude of the ground dish [deg]
# ** lon: longitude of the ground dish [deg]
# ** n: number of integration steps
#
# Output:
# ** TA_num/TA_den: averaged temperature at the satellite antenna [K]
function antennaTempGround(latGround, lonGround, f, offsetAngle, Gmax; n=10)
    L = -30 # ITU antenna model
    BW = BWfromG(Gmax)
    Θ = linspace(0, π, n)
    Φ = linspace(0, π, n)

    TA_num, TA_den = 0, 0
    for θ in Θ, ϕ in Φ
        TB = brightTempGroundCalc(latGround, lonGround, f)
        G = ITUrefRadPat672(Gmax, L, BW, abs(ϕ-offsetAngle-0.5π))
        TA_num += TB*G*sin(θ)
        TA_den += G*sin(θ)
    end
    return TA_num / TA_den
end