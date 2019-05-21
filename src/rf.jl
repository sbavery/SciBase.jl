export
	wavelength,
    frequency,
    db2pow,
	pow2db,
    dbm2pow,
    pow2dbm,
	GfromDL,
    GfromGRMS,
    DfromGL,
	hpbw,
	pathLoss,
	BWfromRDN,
	BWfromD,
	DfromBW,
	demod,
	GfromBW,
	BWfromG,
    farField,
	BER2EbN0,
	EbN02BER,
	EbN02BER_PSK,
	EbN02BER_QAM,
	pointLoss,
	pointError,
	PFDcalc,
	EPFDcalc,
    ITUrefRadPat1528,
	ITUrefRadPat1428,
	ITUrefRadPat672,
	ITUrefRadPat465,
	ITUrefRadPat580,
    ITUrefRadPat1336,
    antenna,
    addAntenna,
    ITUPFDLims,
    FCCPFDLims,
    apparentElevation



wavelength(f) = 1e3c₀/f # f in Hz, returns λ in m.

frequency(λ) = 1e3c₀/λ # λ in m, returns f in Hz.

db2pow(x) = 10 .^(0.1x)

function pow2db(x)
    in(true, x .≤ 0.0) && @error "Negative power values cannot be converted to decibels."
    return 10log10.(x)
end

dbm2pow(x) = db2pow(x-30)

pow2dbm(x) = pow2db(x) + 30

GfromDL(D, λ; η=0.6) = pow2db(η*((π*D/λ)^2))

DfromGL(G, λ; η=0.6) = √(db2pow(G)/η)*(λ/π) # G in dBi

BWfromD(D, f) = deg2rad(21e9/(f*D)) # [Rad]

DfromBW(BW, f) = 21e9/(f*rad2deg(BW)) # BW in Rad

pathLoss(S, λ) = 2pow2db(4π*S./λ) # [dB]. S, λ must have same distance units

BWfromRDN(r, d, n) = 2atan(r*sin(π/n)/(d-r*cos(π/n))) # rad, r = radius of coverage, d = relay-spot distance, n = number of nodes

GfromBW(BW; η=0.6) = pow2db(η.*((21e9π./(1e3c₀*rad2deg(BW))).^2)) # [dB, BW in rad], from SMAD '92 p. 524 (c0 in km/s)

BWfromG(G; η=0.6) = deg2rad(√(η/db2pow(G))*21e9π/(1e3c₀)) # [rad, G in dB], from SMAD '92 p. 524 (c0 in km/s)

function GfromGRMS(ε, G₀, λ; phased_array=false) # Ruze's equation. ε is reflector surface's RMS errors, G₀ is ideal gain [dB].
    phased_array && return G₀ - pow2db(exp((2π)^2))*(ε/λ)^2
    return G₀ - pow2db(exp((4π)^2))*(ε/λ)^2 
end

pointLoss(BW, error) = 12((error/BW)^2) # error and BW in same units, result in dB

pointError(BW, loss) = BW*√(loss/12) # loss in dB, error in same units as BW

farField(D, λ) = 2D^2/λ # Minimum far field range (Fraunhofer distance), in same units as D and λ.

PFDcalc(EIRP, BW, BWref, S) = EIRP + pow2db(BWref/(BW*4π*S^2)) # [dB(W/m2)] BW, BWref in Hz, EIRP in dbW, S in m.

EPFDcalc(EIRP, BW, BWref, Gref, GmaxRef, S) = PFDcalc(EIRP, BW, BWref, S) + Gref - GmaxRef # [dB(W/m2)] BW & BWref in Hz, Gref & GmaxRef in dB, EIRP in dbW, S in m.


function demod(modulation)
	modulation = replace(modulation, "-", "")

	contains(modulation, "PSK") && (keying = "PSK")
	contains(modulation, "QAM") && (keying = "QAM")

	m = 0
	if modulation[1] == 'B'
		m = 1
	elseif modulation[1] == 'Q'
		m = 2
	else
		a = matchall(r"\d", modulation)
		m = log2(parseint(join(a)))
	end

	m, keying
end


function BER2EbN0(m, keying, BER)
    M = 2^m
    keying == "PSK" && return pow2db(((erfcinv(m*BER)/sinpi(1/M))^2)/m)
    keying == "QAM" && return pow2db((2(M-1)/(3m)) * erfcinv(m*BER*√M/(2(√M-1)))^2)
    # return 9.0
end


function EbN02BER(m, keying, EbN₀)
    keying == "PSK" && return EbN02BER_PSK(EbN₀, m)
    keying == "QAM" && return EbN02BER_QAM(EbN₀, m)
end


# TODO: Check units of EbN0
# in AWGN Channel. Eq 2.2 from "Performance of Digital Communication Over Fading Channels"
EbN02BER_PSK(EbN₀::Float64, m::Int64) = sum([erfc(√(2EbN₀*m)*sinpi((2k-1)/2^(0.5+m))) for k = 1:max(2^(m-2), 1)])/max(m, 2)


function EbN02BER_QAM(EbN₀::Float64, m::Int64) # EbN₀ in dB. From "On the General BER Expression of One- and Two-Dimensional Amplitude Modulations"
    if m % 2 == 0
        q, H = √(1.5m*db2pow(EbN₀)/(2^m-1)), 2^(0.5m)
        return 2EbN02BER_QAM_helper(H, q)/(H*m)
    end

    I, J = 2^ceil(0.5m), 2^floor(0.5m)
    q    = √(3m*db2pow(EbN₀)/(I^2+J^2-2))
    return (EbN02BER_QAM_helper(I, q)/I + EbN02BER_QAM_helper(J, q)/J)/m # Eq. 20, 21, 22
end

EbN02BER_QAM_helper(I::T, q::T) where T <: Float64 = sum([sum([(-1)^(floor(i*2^(j-1)/I)) * erfc((2i+1)*q) *  # Eq. 20/1 & 22
                                                          (2^(j-1) - floor(i*2^(j-1)/I+0.5))
                                                for i = 0:(1-2^(-j))*I-1]) for j = 1:log2(I)])


#Reference radiation pattern of antenna diameter D [m], freq f [Hz], major/minor axis ratio z, HALF-3-dB BW ψb ([rad], if known), and max gain Gm, returns gain ([dBi]) at off-axis angle ψ [rad].
function ITUrefRadPat1528(D, f, Gm; z=1, Ln=-15, ψ=0.0, ψb=NaN)  # ITU-R S.1528-0 "Satellite antenna radiation patterns for non-geostationary orbit satellite antennas operating in the fixed-satellite service below 30 GHz".
#     !(f ≤ 30e9) && @warn "ITU-R S.1528 only accurate for frequencies less than 30 GHz."
    
    (Ln == -15) && (a = 1.4)
    (Ln == -20) && (a = 1.0)
    (Ln == -25) && (a = 0.6)
    (Ln == -30) && (a = 0.4)

#     !isdefined(:a) && @error "ITU-R S.1528 only accepts values for L_N of -15, -20, -25, or -30 dB."
    a    = 2.58*√(1 - a*log(z))

    b, α = 6.32, 1.5

    λ    = wavelength(f) # [m]
    ψb   = isnan(ψb) ? 20λ*√3/D : rad2deg(ψb) # [deg]
    Lf   = 0 # Far-out side-lobe level [dBi]
    X    = Gm + Ln + 25log(b*ψb)
    Y    = b*ψb*10^(0.04(Gm + Ln - Lf))
    
    ψ    = rad2deg(ψ)

          (0 ≤ ψ ≤ a*ψb)    && return Gm - 3(ψ/ψb)^α
       (a*ψb < ψ ≤ 0.5b*ψb) && return Gm + Ln + 20log(z)
    (0.5b*ψb < ψ ≤ b*ψb)    && return Gm + Ln
       (b*ψb < ψ ≤ Y)       && return X - 25log(ψ)
          (Y < ψ ≤ 90)      && return Lf
         (90 < ψ ≤ 180)     && return max(15 + Ln + 0.25Gm + 5log(z), 0) # Back-lobe level Lb [dBi]
    @error "ψ must be between 0 and π."
end


#Reference radiation pattern of antenna diameter D [m] and freq f [Hz], returns gMax and gain (both [dBi]) at off-axis angle ϕ [rad].
function ITUrefRadPat1428(D, f; ϕ=0.0) # ITU-R S.1428-1 "Reference FSS earth-station radiation patterns for use with non-GSO satellites".
	!(0 ≤ ϕ ≤ π) && @error "ϕ must be a radian value between 0 and π."
	ϕ = rad2deg(ϕ)

	!(10.7 ≤ f/1e9 ≤ 30) && @warn "ITU-R S.1428-1 only accurate for frequencies between 10.7 and 30 GHz."
	λ = wavelength(f)                  # Wavelength [m]

    Dλ, λD = D/λ, λ/D

	Gmax = 20log10(Dλ) + 7.7       # [dBi]
	G₁ = 29 - 25log10(95λD)   # [dBi]
	ϕₘ = 20λD*√(Gmax - G₁)       # [deg]

	if 20 ≤ Dλ ≤ 25
		if ϕ < ϕₘ
			G = Gmax - 2.5e-3(Dλ*ϕ)^2
		elseif ϕ < 95λD
			G = G₁
		elseif ϕ < 33.1
			G = 29 - 25log10(ϕ)
		elseif ϕ ≤ 80
			G = -9
		else
			G = -5
		end
	elseif 25 < Dλ ≤ 100
		if ϕ < ϕₘ
			G = Gmax - 2.5e-3(Dλ*ϕ)^2
		elseif ϕ < 95λD
			G = G₁
		elseif ϕ < 33.1
			G = 29 - 25log10(ϕ)
		elseif ϕ ≤ 80
			G = -9
		elseif ϕ ≤ 120
			G = -4
		else
			G = -9
		end
	elseif Dλ > 100
		Gmax = 20log10(Dλ) + 8.4     # [dBi]
		G₁ = -1 + 15log10(Dλ)
		ϕₘ = 20λD*√(Gmax - G₁)        # [deg]
		ϕᵣ = 15.85Dλ^-0.6             # [deg]
		if ϕ < ϕₘ
			G = Gmax - 2.5e-3(Dλ*ϕ)^2
		elseif ϕ < ϕᵣ
			G = G₁
		elseif ϕ < 10
			G = 29 - 25log10(ϕ)
		elseif ϕ < 34.1
			G = 34 - 30log10(ϕ)
		elseif ϕ ≤ 80
			G = -12
		elseif ϕ ≤ 120
			G = -7
		else
			G = -12
		end
	else
		@error "ERROR: The D, λ inputs were not recognized by the ITU standard."
	end
	G, Gmax
end

function ITUrefRadPat672(Gₘ, Lₛ, ψ₀, ψ) # ITU-R S.672-4 Reference radiation pattern for GEO space station antenna with max gain Gm, half 3-dB beamwidth ψ0 [rad] and near-in-side-lobe level Ls [dB], returns G [dBi] at off-axis angle ψ [rad].
	if Lₛ == -10 # From ITU-R S.1433 p.3
		a, b = 1.83, 6.32
	elseif Lₛ == -20
		a, b = 2.58, 6.32
	elseif Lₛ == -25
		a, b = 2.88, 6.32
	elseif Lₛ == -30
		a, b = 3.16, 6.32
	else
		@error "The value of Lₛ is not recognized by ITU-R S.672-4."
	end

	ψ, ψ₀ = rad2deg(ψ), rad2deg(ψ₀)
	ψ₁    = ψ₀*10^(0.04(20+Gₘ+Lₛ))

	(  0  ≤ ψ <   ψ₀) && return Gₘ - 3(ψ/ψ₀)^2 # This condition not in ITU standard.
	(  ψ₀ ≤ ψ ≤ a*ψ₀) && return Gₘ - 3(ψ/ψ₀)^2
	(a*ψ₀ < ψ ≤ b*ψ₀) && return Gₘ + Lₛ
	(b*ψ₀ < ψ ≤   ψ₁) && return Gₘ + Lₛ + 20 - 25log10(ψ/ψ₀)
	(  ψ₁ < ψ       ) && return 0

	@error "Could not be processed."
end


function ITUrefRadPat465(D, f, Gm, ϕ) # ITU-R S.465-5 "Reference earth-station radiation pattern". Parabolic antenna with max gain Gm [dB] and diameter D [m], at frequency f [Hz]. Returns G [dBi] at off-axis angle ϕ [rad].
    !(2 ≤ f/1e9 ≤ 30) && @warn "ITU-R S.465-5 only accurate for frequencies between 2 and 30 GHz."
    λ = wavelength(f)                  # Wavelength [m]
    ϕ = rad2deg(ϕ)
    if     D/λ < 33.3
        ϕmin = 2.5
    elseif D/λ < 50
        ϕmin = max(2, 114(D/λ)^-1.09)
    else
        ϕmin = max(1, 100λ/D)
    end

    ( 0   ≤ ϕ <  ϕmin) && return (ϕ^2)*(32 - 25log10(ϕmin) - Gm)/(ϕmin^2) + Gm # This condition not in ITU standard
    (ϕmin ≤ ϕ <  48  ) && return        32 - 25log10(ϕ)
    (48   ≤ ϕ ≤ 180  ) && return       -10
    @error "ϕ must be a radian angle between 0 and π."
end


function ITUrefRadPat580(D, f, Gm, ϕ) # ITU-R S.580-6 "Radiation diagrams for antennas of earth stations operating with geostationary satellites". Parabolic antenna with max gain Gm [dB] and diameter D [m] at frequency f [Hz]. Returns G [dBi] at off-axis angle phi [rad].
    λ = wavelength(f) # Wavelength [m]
    ϕ = rad2deg(ϕ)
    ϕmin = max(1, 100λ/D)
    ( 0.0 ≤ ϕ <  ϕmin) && return (ϕ^2)*(29 - 25log10(ϕmin) - Gm)/(ϕmin^2) + Gm # This condition not in ITU standard
    (ϕmin ≤ ϕ <  48.0) && return        29 - 25log10(ϕ)
    (48.0 ≤ ϕ ≤ 180.0) && return       -10
    @error "ϕ must be a radian angle between 0 and π."
end


# ITU-R F.1336-1, 2.3: Low-cost, low-gain antenna in the 1-3 GHz range with circular symmetry about the 3 dB beamwidth and with a gain less than about 20 dBi.
# Returns gain (dBi) at off-angle θ (rad).
function ITUrefRadPat1336(G₀; θ=0)
    ϵ = 1e-6
    !(0 ≤ abs(θ - ϵ) ≤ π) && @error "θ must be a radian value between 0 and π."
    θ  = rad2deg(θ)             # [°]
    ϕ₃ = √(2.7*10^(4 - 0.1G₀))  # [°]
    0      ≤ θ < 1.08ϕ₃ && return G₀ - 12(θ/ϕ₃)^2
    ϕ₁ = 1.9ϕ₃                  # [°]
    1.08ϕ₃ ≤ θ < ϕ₁     && return G₀ - 14
    ϕ₂ = ϕ₁*10^((G₀-6)/32)      # [°]
        ϕ₁ ≤ θ < ϕ₂     && return G₀ - 14 - 32log(θ/ϕ₁)
                           return -8    
end


mutable struct Antenna{S <: String, F <: Float64}
    parent::S # station or spacecraft name, e.g. "SFO" or "ISS(ZARYA)"
    target::S # points at a "station", "constellation", or "client"
    G::F # Gain [dB]
    GT::F # G/T (Gain/Temperature) [dB/K]
    P::F # Power [dBW]
    f::F # Frequency [Hz]
    # T::F # Antenna noise temperature [K]
    L::F # Array of antenna losses [dB] (Don't include path losses)
    EbN0_req::F # Required EbN0 to close a link at target BER [dB] (Rx ants only)
    FEC::F # Coding rate [-]
    margin::F # Link margin [dB]
end

function addAntenna(parent, target; G=0.0, P_dBW=0.0, P_dBm=0.0, P_W=0.0, f=0.0, T=0.0, L=[0.0], GT=0.0, D=0.0, η=0.6, EbN0_req=0.0, FEC=1.0, margin=0.0)
    global dName
    parent = uppercase(parent)
    target = lowercase(target)
    !in(target, ["station", "client", "constellation"]) && @error "Invalid antenna target."

    (P_dBW == 0.0 && P_W   ≠ 0.0) && (P_dBW = pow2db(P_W))
    (P_dBW == 0.0 && P_dBm ≠ 0.0) && (P_dBW = P_dBm - 30.0)

    if G == 0.0 && GT == 0.0 && GT ≠ 0.0
        λ = wavelength(f)
        G = GfromDL(D, λ, η=η)
    end

    if GT == 0.0 && G ≠ 0.0 && T ≠ 0.0
        Rx_NF = 3.0 # dB
        Temp_cont = 290.0  # 290K Connection temperature
        Temp_eff = Temp_cont*0.1db2pow(Rx_NF) # K effective input noise temperature of receiver
        noiseTemp = T/db2pow(Rx_ant.L) + Temp_cont*(1-db2pow(-Rx_ant.L)) + Temp_eff # dBK System noise
        GT = G - pow2db(noiseTemp)
    end
    ant = Antenna(parent, target, G, GT, P_dBW, f, abs(sum(L)), EbN0_req, FEC, margin)

    try
        jldopen(scLib_path, "r+") do scLib
            SC = read(scLib, parent)
            o_delete(scLib, SC.name)
            push!(SC.children, ant)
            write(scLib, SC.name, SC)
        end
    catch
        if isdefined(parse(parent))
            push!(eval(parse(parent * ".children")), ant)
        else
            @error "$parent not found."
        end
    end
end




function ITUPFDLims(f, δ; GSO=true, i=0.0, N=1, a=8000) # ITU RR21-6 Limits of PFD from space stations
    !(0.0 ≤ δ ≤ 90.0) && @error "Angle of arrival δ must be between 0° and 90°."

    if 10.7e9 ≤ f < 11.7e9
        if GSO
                if δ ≤ 5.0
                    return -150.0, 4e3
                elseif δ ≤ 25.0
                    return -150.0 + 0.5*(δ-5.0), 4e3
                else
                    return -140.0, 4e3
                end
        else # NGSO
            if 35.0 < i < 145.0 && a > 18000.0+R⨁
                if δ ≤ 5.0
                    return -129.0, 1e6
                elseif δ ≤ 25.0
                    return -129.0 + 0.75*(δ-5.0), 1e6
                else
                    return -114.0, 1e6
                end
            else
                if δ ≤ 5.0
                    return -129.0, 1e6
                elseif δ ≤ 25.0
                    return -129.0 + 0.75*(δ-5.0), 1e6
                else
                    return -114.0, 1e6
                end
            end
        end
    elseif 15.43e9 ≤ f < 15.63e9
                if δ ≤ 20.0
                    return -127.0, 1e6
                elseif δ ≤ 25.0
                    return -127.0 + 0.56*(δ-20.0), 1e6
                elseif δ ≤ 29.0
                    return -113.0, 1e6
                elseif δ ≤ 31.0
                    return -136.9 + 25.0*log10(δ-20.0), 1e6
                else
                    return -111, 1e6
                end
    elseif 17.7e9 ≤ f < 19.3e9
                if δ ≤ 5.0
                    return -115.0, 1e6
                elseif δ ≤ 25.0
                    return -115.0 + 0.5*(δ-5.0), 1e6
                else
                    return -105.0, 1e6
                end
    elseif 19.30e9 ≤ f < 19.70e9 ||
           22.55e9 ≤ f < 23.55e9 ||
           24.45e9 ≤ f < 24.75e9 ||
           25.25e9 ≤ f < 27.50e9 ||
           27.50e9 ≤ f < 27.501e9
                if δ ≤ 5.0
                    return -115.0, 1e6
                elseif δ ≤ 25.0
                    return -115.0 + 0.5*(δ-5.0), 1e6
                else
                    return -105.0, 1e6
                end
    elseif 31.00e9 ≤ f < 31.30e9 ||
           34.70e9 ≤ f < 35.20e9
                if δ ≤ 5.0
                    return -115.0, 1e6
                elseif δ ≤ 25.0
                    return -115.0 + 0.5*(δ-5.0), 1e6
                else
                    return -105.0, 1e6
                end
    elseif 31.8e9 ≤ f < 32.3e9
                if δ ≤ 5.0
                    return -120.0, 1e6
                elseif δ ≤ 25.0
                    return -120.0 + 0.75*(δ-5.0), 1e6
                else
                    return -105.0, 1e6
                end
    elseif 32.3e9 ≤ f < 33e9
                if δ ≤ 5.0
                    return -135.0, 1e6
                elseif δ ≤ 25.0
                    return -135.0 + (δ-5.0), 1e6
                else
                    return -115.0, 1e6
                end
    elseif 37.5e9 ≤ f < 40e9
        if GSO
                if δ ≤ 5.0
                    return -127.0, 1e6
                elseif δ ≤ 20.0
                    return -127.0 + (4.0/3.0)*(δ-5.0), 1e6
                elseif δ ≤ 25.0
                    return -107.0 + 0.4*(δ-20.0), 1e6
                else
                    return -105.0, 1e6
                end
        else # NGSO
                if δ ≤ 5.0
                    return -120.0, 1e6
                elseif δ ≤ 25.0
                    return -120.0 + 0.75*(δ-5.0), 1e6
                else
                    return -105.0, 1e6
                end
        end
    elseif 40e9 ≤ f < 40.5e9
                if δ ≤ 5.0
                    return -115.0, 1e6
                elseif δ ≤ 25.0
                    return -115.0 + 0.5*(δ-5.0), 1e6
                else
                    return -105.0, 1e6
                end
    elseif 40.5e9 ≤ f < 42e9
        if GSO
                if δ ≤ 5.0
                    return -120.0, 1e6
                elseif δ ≤ 15.0
                    return -120.0 + (δ-5.0), 1e6
                elseif δ ≤ 25.0
                    return -110.0 + 0.5*(δ-15.0), 1e6
                else
                    return -105.0, 1e6
                end
        else # NGSO
                if δ ≤ 5.0
                    return -115.0, 1e6
                elseif δ ≤ 25.0
                    return -115.0 + 0.5*(δ-5.0), 1e6
                else
                    return -105.0, 1e6
                end
        end
    elseif 42e9 ≤ f < 42.5e9
        if GSO
                if δ ≤ 5.0
                    return -127.0, 1e6
                elseif δ ≤ 20.0
                    return -127.0 + (4.0/3.0)*(δ-5.0), 1e6
                elseif δ ≤ 25.0
                    return -107.0 + 0.4*(δ-20.0), 1e6
                else
                    return -105.0, 1e6
                end
        else # NGSO
                if δ ≤ 5.0
                    return -120.0, 1e6
                elseif δ ≤ 25.0
                    return -120.0 + 0.75*(δ-5.0), 1e6
                else
                    return -105.0, 1e6
                end
        end
    else
        return NaN, NaN
    end
end


function getXITU(N)
    N ≤  50 && return 0.0
    N ≤ 288 && return (5/119)*(N- 50)
               return (1/ 69)*(N+402)
end


function FCCPFDLims(f; δ=90.0, GSO=true)#, i=0.0, N=1, a=8000)
    !(0.0 ≤ δ ≤ 90.0) && @error "Angle of arrival δ must be between 0° and 90°."

    if 17.70e9 ≤ f < 17.80e9 ||
       18.30e9 ≤ f < 18.80e9 ||
       19.30e9 ≤ f < 19.70e9 ||
       22.55e9 ≤ f < 23.55e9 ||
       24.45e9 ≤ f < 24.75e9
        if δ ≤ 5.0
            return -115.0, 1e6
        elseif δ ≤ 25.0
            return -115.0 + 0.5*(δ-5.0), 1e6
        else
            return -105.0, 1e6
        end
    elseif 37.5e9 ≤ f < 40e9
        if GSO
#                 if δ ≤ 5.0
#                     return -127.0, 1e6
#                 elseif δ ≤ 20.0
#                     return -127.0 + (4.0/3.0)*(δ-5.0), 1e6
#                 elseif δ ≤ 25.0
#                     return -107.0 + 0.4*(δ-20.0), 1e6
#                 else
#                     return -105.0, 1e6
#                 end
        else # NGSO
            if δ ≤ 5.0
                return -132.0, 1e6
            elseif δ ≤ 25.0
                return -132.0 + 0.75*(δ-5.0), 1e6
            else
                return -117.0, 1e6
            end
        end
    elseif 40.0e9 ≤ f < 40.5e9
        if δ ≤ 5.0
            return -115.0, 1e6
        elseif δ ≤ 25.0
            return -115.0 + 0.5*(δ-5.0), 1e6
        else
            return -105.0, 1e6
        end
    elseif 40.5e9 ≤ f < 42.0e9
        if GSO
#                 if δ ≤ 5.0
#                     return -127.0, 1e6
#                 elseif δ ≤ 20.0
#                     return -127.0 + (4.0/3.0)*(δ-5.0), 1e6
#                 elseif δ ≤ 25.0
#                     return -107.0 + 0.4*(δ-20.0), 1e6
#                 else
#                     return -105.0, 1e6
#                 end
        else # NGSO
            if δ ≤ 5.0
                return -115.0, 1e6
            elseif δ ≤ 25.0
                return -115.0 + 0.5*(δ-5.0), 1e6
            else
                return -105.0, 1e6
            end
        end
    elseif 54.25e9 ≤ f < 56.9e9 ||
            57.0e9 ≤ f < 58.2e9 ||
            59.0e9 ≤ f < 59.3e9
#         return -147, 1e8
        return -167.0, 1e6
    else
        return NaN, NaN
    end
end



# Effects of Tropospheric Refraction on Radiowave Propogation (ITU-R P.834-7)
# --------------------------------
# Inputs:
# ** θ₀: Elevation angle of a space station in free-space propagation conditions [rad]
# ** h: Altitude of Earth station [km] (default=0)
# ** r: Earth's radius
# Outputs:
# ** θ: Apparent elevation angle [rad]
function apparentElevation(θ₀; h=0, r=R⨁)
    θₘ = -acos(r/(r+h) * n(0)/n(h))
    (θₘ - τ(h, θₘ) ≤ θ₀) && return θ₀ + τs(h, θ₀) # Eqn. (11)
    return 0.0
end

function n(x)
    a = 3.15e-4
    b = 0.1361
    return 1 + a*exp(-b*x)
end

τ(h, θ)   = deg2rad(inv(1.314 + 0.6437θ  + 0.02869θ^2  + h*(0.2305 + 0.09428θ  + 0.01096θ ^2) + 0.008583h^2))
τs(h, θ₀) = deg2rad(inv(1.728 + 0.5411θ₀ + 0.03723θ₀^2 + h*(0.1815 + 0.06272θ₀ + 0.01380θ₀^2) + (h^2)*(0.01727 + 0.008288θ₀)))