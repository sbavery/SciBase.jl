export
	datapath,
	⊗,
	
	c0,
	c₀,
	h0,
	h₀,
	μ₀,
	ɛ₀,
	T0,
	S_E_M,
	r_E_M,
	S_E_LEO,
	S_E_GEO,
	r_E_GEO,

	reEarth,
	R⨁,
	reEarth_a,
	R⨁_a,
	reEarth_b,
	R⨁_b,
	muEarth,
	μ⨁,
	rotRateEarth,
	ω⨁,
	ω⃗⨁,
	solSecsEarth,
	sideSecsEarth,
	e1Earth,
	e1⨁,
	e1Earth2,
	e1⨁2,
	e2Earth,
	e2⨁,
	e2Earth2,
	e2⨁2,
	f⨁,
	J₂⨁,
	J₃⨁,

	 m_in_ft,
	ft_in_m,
	km_in_nm,
	in_in_m,
	 m_in_in,
	km_in_DU,
	km_in_AU,
 	 s_in_TU,

	reSun,
	R☉,
	muSun,
	μ☉,
	day_in_TU☉,
	s_in_TU☉,

	reMoon,
	R☾,
	muMoon,
	μ☾,
	muMercury,
	reMercury,
	muVenus,
	reVenus,
	muMars,
	reMars,
	meanObliq_Mars,
	muJupiter,
	reJupiter,
	muSaturn,
	reSaturn,
	muUranus,
	reUranus,
	muNeptune,
	reNeptune,

	DVBS2,
	DVBS2X,

	lightBlue,
	lightGreen,
	audacyBlue,
	audacyBlack,
	audacyCyan,
	audacyLightGray,
	audacyGreen,
	audacyOrange,

	groundstation,

	ALM,
	AYA,
	BKR,
	CAP,
	CNJ,
	CPJ,
	DLS,
	DON,
	DXB,
	EGE,
	GNP,
	GSP,
	HBK,
	HUX,
	HWD,
	ICN,
	KDK,
	KLN,
	KRO,
	KSC,
	LAX,
	LON,
	LUX,
	LYR,
	MAD,
	MEX,
	MHA,
	NPA,
	PER,
	PDX,
	PSK,
	SAN,
	SDN,
	SEA,
	SFO,
	SHE,
	SIN,
	SJC,
	STX,
	SYD,
	TMA,
	UIO,
	UNA,
	USH,
	VBG,
	VSB,
	VTY,
	WCN,
	WPS,
	WTX,
	XCN,

	RSA,
	HCH,
	YPG,
	BAB,
	CPK,
	NID,
	EDW,
	IRW,
	NXP,
	BKF,
	SRV,
	GDN,
	GUM,
	OHU,
	FDK,
	LSV,
	NTS,
	XSD,
	CVS,
	WSM,
	DYS,
	BIF,
	FSH,
	GOF,
	SKF,
	UTR,
	DAA,
	CSP,

	GBK,
	VLA,
	UAZ,
	CTC,
	AST,
	WFD,
	JCM,
	CMA,
	BSR,
	DVS,
	HCK,
	KPK,
	AMS,
	KEA,
	LBT,
	OWS,
	PTN,
	CRX,

	TLEurls

const datapath      = joinpath(@__DIR__, "..", "data")

const ⊗  			= kron

const c0 		    = 299792.458          	# [km/s]
const c₀ 			= c0
const h0 		    = 1.3806488e-23       	# [m²kg/s²K] Boltzmann's Constant
const h₀ 			= h0
const μ₀          	= 4e-7π	           		# [N/A²] Vacuum Permeability / permeability of free space / magnetic constant
const ɛ₀         	= inv(μ₀*(1000c₀)^2)  	# [F/m]  Vacuum Permittivity / Permittivity of free space / electric constant

const T0 		    = 290.0                	# [K]
const T₀ 			= T0

const reEarth	    = 6378.1366				# [km] Mean radius
const R⨁	    	= reEarth
const reEarth_a    	= 6378.1370    			# [km] 	Earth's semi-major axis (WGS84)
const R⨁_a 			= reEarth_a
const f⨁			= inv(298.257223563)	# Earth's flattening factor (WGS84)
const reEarth_b    	= reEarth_a*(1-f⨁)		# [km] Semi-minor axis (WGS84)
const R⨁_b			= reEarth_b
const e1Earth 		= √(2f⨁-f⨁^2)			# Earth's first eccentricity
const e1⨁     		= √(2f⨁-f⨁^2)			# Earth's first eccentricity
const e1Earth2      = 2f⨁-f⨁^2				# Earth's first eccentricity squared
const e1⨁2      	= 2f⨁-f⨁^2				# Earth's first eccentricity squared
const e2Earth 		= √(f⨁*(2-f⨁))/(1-f⨁)	# Earth's second eccentricity
const e2⨁           = √(f⨁*(2-f⨁))/(1-f⨁)	# Earth's second eccentricity
const e2Earth2 		= f⨁*(2-f⨁)/(1-f⨁)^2 	# Earth's second eccentricity squared
const e2⨁2 			= f⨁*(2-f⨁)/(1-f⨁)^2 	# Earth's second eccentricity squared
const muEarth      	= 3.986004418e5			# [km3/s²]
const μ⨁      		= muEarth
const rotRateEarth 	= 0.0000729211585530 	# [rad/s]
const ω⨁ 			= rotRateEarth
const ω⃗⨁ 			= [0.0, 0.0, ω⨁]
const g0 			= 9.8196 				# [m/s²]
const solSecsEarth  = 86400.0 				# [s] Secs in Earth's solar day
const sideSecsEarth	= 86164.0916			# [s] Secs in Earth's sidereal day
const J₂⨁			=  1.08262668e-3		# [-] J₂ - Earth's second dynamic form factor
const J₃⨁			= -2.53215306e-6		# [-] J₃ - Earth's third dynamic form factor


const S_E_M 	    = 377622.0           	# [km]
const r_E_M 		= 384400.0 				# [km]
const S_E_LEO      	= 400.0              	# [km] Between 160 and 2000km
const r_E_GEO 		= ∛(μ⨁*(0.5sideSecsEarth/π)^2) 	# [km]
const S_E_GEO      	= r_E_GEO - R⨁     		# [km]

const  m_in_in 		= 0.0254 				# meters in an inch
const in_in_m 		= inv(m_in_in) 			# inches in a meter
const  m_in_ft      = 0.3048				# meters in a foot
const ft_in_m 		= inv(m_in_ft) 			# feet in a meter
const km_in_nm      = 1.852 				# km in nautical mile
const km_in_DU 		= R⨁ 					# km in DU
const km_in_AU     	= 149597870.7          	# km in AU
const  s_in_TU 		= √(R⨁^3/μ⨁) 			# s in ⨁ TU [Vallado (2001), p. 230]

const muSun        	= 1.32712440018e11    	# [km³/s²]
const reSun 		= 696300.0 				# [km]
const μ☉        	= muSun 		    	# [km³/s²]
const R☉ 			= reSun 				# [km]
const day_in_TU☉    = 58.123440906 			# [solar days] Solar days in one Solar TU (Vallado (2001), p.905)
const   s_in_TU☉    = day_in_TU☉*solSecsEarth # [secs] Secs in a Solar TU (Vallado (2001))


const muMoon       	= 4902.8				# [km³/s²]
const reMoon 		= 1737.1				# [km] # Mean radius
const μ☾       		= muMoon				# [km³/s²]
const R☾	 		= reMoon				# [km] # Mean radius

const muMercury 	= 22032.0				# [km3/s²]
const reMercury 	= 2440.0				# [km] # Mean radius

const muVenus 		= 324859.0				# [km3/s²]
const reVenus 		= 6052.0				# [km] # Mean radius

const muMars 		= 42828.0				# [km³/s²]
const meanObliq_Mars = deg2rad(25.19)		# [rad]
const reMars		= 3396.2				# [km]

const muJupiter 	= 126686534.0			# [km³/s²]
const reJupiter 	= 71492.0				# [km] # Mean radius

const muSaturn 		= 37931187.0			# [km³/s²]
const reSaturn 		= 60268.0				# [km] # Mean radius

const muUranus 		= 579393913.0			# [km³/s²]
const reUranus 		= 25559.0				# [km] # Mean radius

const muNeptune 	= 6836529.0				# [km³/s²]
const reNeptune 	= 24766.0				# [km] # Mean radius



const colors       	= ["g", "r", "c", "m", "y", "k", "b"]
const lightBlue 	= rgb(135.0, 206.0, 250.0)
const lightGreen 	= rgb(206.0, 253.0, 206.0)
const audacyBlack 	= rgb(002.0, 006.0, 009.0)
const audacyBlue 	= rgb(023.0, 033.0, 104.0)
const audacyCyan 	= rgb(007.0, 209.0, 234.0)
const audacyLightGray = rgb(241.0, 241.0, 245.0)
const audacyOrange  = rgb(255.0, 109.0, 000.0)
const audacyGreen   = rgb(018.0, 199.0, 000.0)


const DVBS2 = Dict{String,Dict{String,Float64}}(
	  "QPSK 1/4"  => Dict("EbN0_req" => 0.75, "ρ" => 0.25, "Γ" => 0.49, "m" => 2.0),
	  "QPSK 1/3"  => Dict("EbN0_req" => 0.59, "ρ" => 0.33, "Γ" => 0.66, "m" => 2.0),
	  "QPSK 2/5"  => Dict("EbN0_req" => 0.73, "ρ" => 0.40, "Γ" => 0.79, "m" => 2.0),
	  "QPSK 1/2"  => Dict("EbN0_req" => 1.05, "ρ" => 0.50, "Γ" => 0.99, "m" => 2.0),
	  "QPSK 3/5"  => Dict("EbN0_req" => 1.48, "ρ" => 0.60, "Γ" => 1.19, "m" => 2.0),
	  "QPSK 2/3"  => Dict("EbN0_req" => 1.89, "ρ" => 0.67, "Γ" => 1.32, "m" => 2.0),
	  "QPSK 3/4"  => Dict("EbN0_req" => 2.31, "ρ" => 0.75, "Γ" => 1.49, "m" => 2.0),
	  "QPSK 4/5"  => Dict("EbN0_req" => 2.67, "ρ" => 0.80, "Γ" => 1.59, "m" => 2.0),
	  "QPSK 5/6"  => Dict("EbN0_req" => 2.99, "ρ" => 0.83, "Γ" => 1.65, "m" => 2.0),
	  "QPSK 8/9"  => Dict("EbN0_req" => 3.73, "ρ" => 0.89, "Γ" => 1.77, "m" => 2.0),
	  "QPSK 9/10" => Dict("EbN0_req" => 3.89, "ρ" => 0.90, "Γ" => 1.79, "m" => 2.0),
	  "8PSK 3/5"  => Dict("EbN0_req" => 3.00, "ρ" => 0.60, "Γ" => 1.78, "m" => 3.0),
	  "8PSK 2/3"  => Dict("EbN0_req" => 3.65, "ρ" => 0.67, "Γ" => 1.98, "m" => 3.0),
	  "8PSK 3/4"  => Dict("EbN0_req" => 4.43, "ρ" => 0.75, "Γ" => 2.23, "m" => 3.0),
	  "8PSK 5/6"  => Dict("EbN0_req" => 5.41, "ρ" => 0.83, "Γ" => 2.48, "m" => 3.0),
	  "8PSK 8/9"  => Dict("EbN0_req" => 6.46, "ρ" => 0.89, "Γ" => 2.65, "m" => 3.0),
	  "8PSK 9/10" => Dict("EbN0_req" => 6.70, "ρ" => 0.90, "Γ" => 2.68, "m" => 3.0),
	"16APSK 2/3"  => Dict("EbN0_req" => 4.76, "ρ" => 0.67, "Γ" => 2.64, "m" => 4.0),
	"16APSK 3/4"  => Dict("EbN0_req" => 5.49, "ρ" => 0.75, "Γ" => 2.97, "m" => 4.0),
	"16APSK 4/5"  => Dict("EbN0_req" => 6.03, "ρ" => 0.80, "Γ" => 3.17, "m" => 4.0),
	"16APSK 5/6"  => Dict("EbN0_req" => 6.42, "ρ" => 0.83, "Γ" => 3.30, "m" => 4.0),
	"16APSK 8/9"  => Dict("EbN0_req" => 7.42, "ρ" => 0.89, "Γ" => 3.52, "m" => 4.0),
	"16APSK 9/10" => Dict("EbN0_req" => 7.61, "ρ" => 0.90, "Γ" => 3.57, "m" => 4.0),
	"32APSK 3/4"  => Dict("EbN0_req" => 7.05, "ρ" => 0.75, "Γ" => 3.70, "m" => 5.0),
	"32APSK 4/5"  => Dict("EbN0_req" => 7.67, "ρ" => 0.80, "Γ" => 3.95, "m" => 5.0),
	"32APSK 5/6"  => Dict("EbN0_req" => 8.13, "ρ" => 0.83, "Γ" => 4.12, "m" => 5.0),
	"32APSK 8/9"  => Dict("EbN0_req" => 9.26, "ρ" => 0.89, "Γ" => 4.40, "m" => 5.0),
	"32APSK 9/10" => Dict("EbN0_req" => 9.57, "ρ" => 0.90, "Γ" => 4.45, "m" => 5.0)
)

const DVBS2X = Dict{String,Dict{String,Float64}}(
	   "QPSK 2/9"   => Dict("EbN0_req" =>  0.77, "ρ" =>  2/9 , "Γ" => 0.43, "m" => 2.0),
	   "QPSK 13/45" => Dict("EbN0_req" =>  0.43, "ρ" => 13/45, "Γ" => 0.57, "m" => 2.0),
	   "QPSK 9/20"  => Dict("EbN0_req" =>  0.73, "ρ" =>  9/20, "Γ" => 0.89, "m" => 2.0),
	   "QPSK 11/20" => Dict("EbN0_req" =>  1.08, "ρ" => 11/20, "Γ" => 1.09, "m" => 2.0),
	  "8APSK 5/9"   => Dict("EbN0_req" =>  2.56, "ρ" =>  5/9 , "Γ" => 1.65, "m" => 3.0),
	  "8APSK 26/45" => Dict("EbN0_req" =>  2.79, "ρ" => 26/45, "Γ" => 1.71, "m" => 3.0),
	   "8PSK 23/36" => Dict("EbN0_req" =>  3.34, "ρ" => 23/36, "Γ" => 1.90, "m" => 3.0),
	   "8PSK 25/36" => Dict("EbN0_req" =>  3.88, "ρ" => 25/36, "Γ" => 2.06, "m" => 3.0),
	   "8PSK 13/18" => Dict("EbN0_req" =>  4.18, "ρ" => 13/18, "Γ" => 2.15, "m" => 3.0),
	 "16APSK 1/2"   => Dict("EbN0_req" =>  3.02, "ρ" =>  1/2 , "Γ" => 1.97, "m" => 4.0),
	 "16APSK 8/15"  => Dict("EbN0_req" =>  3.32, "ρ" =>  8/15, "Γ" => 2.10, "m" => 4.0),
	 "16APSK 5/9"   => Dict("EbN0_req" =>  3.43, "ρ" =>  5/9 , "Γ" => 2.19, "m" => 4.0),
	 "16APSK 26/45" => Dict("EbN0_req" =>  3.93, "ρ" => 26/45, "Γ" => 2.28, "m" => 4.0),
	 "16APSK 3/5"   => Dict("EbN0_req" =>  4.05, "ρ" =>  3/5 , "Γ" => 2.37, "m" => 4.0),
	 "16APSK 3/5"   => Dict("EbN0_req" =>  3.66, "ρ" =>  3/5 , "Γ" => 2.37, "m" => 4.0),
	 "16APSK 28/45" => Dict("EbN0_req" =>  4.19, "ρ" => 28/45, "Γ" => 2.46, "m" => 4.0),
	 "16APSK 23/36" => Dict("EbN0_req" =>  4.36, "ρ" => 23/36, "Γ" => 2.52, "m" => 4.0),
	 "16APSK 2/3"   => Dict("EbN0_req" =>  4.22, "ρ" =>  2/3 , "Γ" => 2.64, "m" => 4.0),
	 "16APSK 25/36" => Dict("EbN0_req" =>  4.88, "ρ" => 25/36, "Γ" => 2.75, "m" => 4.0),
	 "16APSK 13/18" => Dict("EbN0_req" =>  5.15, "ρ" => 13/18, "Γ" => 2.86, "m" => 4.0),
	 "16APSK 7/9"   => Dict("EbN0_req" =>  5.77, "ρ" =>  7/9 , "Γ" => 3.08, "m" => 4.0),
	 "16APSK 77/90" => Dict("EbN0_req" =>  6.69, "ρ" => 77/90, "Γ" => 3.39, "m" => 4.0),
	 "32APSK 2/3"   => Dict("EbN0_req" =>  5.93, "ρ" =>  2/3 , "Γ" => 3.29, "m" => 5.0),
	 "32APSK 32/45" => Dict("EbN0_req" =>  6.30, "ρ" => 32/45, "Γ" => 3.51, "m" => 5.0),
	 "32APSK 11/15" => Dict("EbN0_req" =>  6.58, "ρ" => 11/15, "Γ" => 3.62, "m" => 5.0),
	 "32APSK 7/9"   => Dict("EbN0_req" =>  7.21, "ρ" =>  7/9 , "Γ" => 3.84, "m" => 5.0),
	 "64APSK 32/45" => Dict("EbN0_req" =>  7.74, "ρ" => 32/45, "Γ" => 4.21, "m" => 6.0),
	 "64APSK 11/15" => Dict("EbN0_req" =>  8.44, "ρ" => 11/15, "Γ" => 4.34, "m" => 6.0),
	 "64APSK 7/9"   => Dict("EbN0_req" =>  8.84, "ρ" =>  7/9 , "Γ" => 4.60, "m" => 6.0),
	 "64APSK 4/5"   => Dict("EbN0_req" =>  9.12, "ρ" =>  4/5 , "Γ" => 4.74, "m" => 6.0),
	 "64APSK 5/6"   => Dict("EbN0_req" =>  9.62, "ρ" =>  5/6 , "Γ" => 4.93, "m" => 6.0),
	"128APSK 3/4"   => Dict("EbN0_req" => 10.60, "ρ" =>  3/4 , "Γ" => 5.16, "m" => 7.0),
	"128APSK 7/9"   => Dict("EbN0_req" => 11.24, "ρ" =>  7/9 , "Γ" => 5.36, "m" => 7.0),
	"256APSK 29/45" => Dict("EbN0_req" =>  9.93, "ρ" => 29/45, "Γ" => 5.07, "m" => 8.0),
	"256APSK 2/3"   => Dict("EbN0_req" => 10.05, "ρ" =>  2/3 , "Γ" => 5.24, "m" => 8.0),
	"256APSK 31/45" => Dict("EbN0_req" => 10.76, "ρ" => 31/45, "Γ" => 5.42, "m" => 8.0),
	"256APSK 32/45" => Dict("EbN0_req" => 11.11, "ρ" => 32/45, "Γ" => 5.59, "m" => 8.0),
	"256APSK 11/15" => Dict("EbN0_req" => 11.23, "ρ" => 11/15, "Γ" => 5.77, "m" => 8.0),
	"256APSK 3/4"   => Dict("EbN0_req" => 11.86, "ρ" =>  3/4 , "Γ" => 5.90, "m" => 8.0)
)



mutable struct Groundstation
    name::String
    body::String # e.g. "EARTH", "MARS"
    lat::Float64 	# [deg] # Geodetic coordinates (gdc)
    lon::Float64 	# [deg] # Geodetic coordinates (gdc)
    alt::Float64 	# [m] above sea level
    children
end

const ALM = Groundstation("ALM", "Earth", +032.5014, -106.6136, 1447, Any[]) 	# White Sands, NM, USA
const AYA = Groundstation("AYA", "Earth", +069.2944, +016.0198, 0,    Any[])	# Andøya Space Center/Rocket Range, Norway.
const BKR = Groundstation("BKR", "Earth", +045.965 , +063.305 , 0,    Any[]) 	# Baikonur, Kazakhstan
const CAP = Groundstation("CAP", "Earth", -033.9253, +018.4239, 0,    Any[]) 	# Cape Town, South Africa
const CNJ = Groundstation("CNJ", "Earth", -020.7044, +140.5056, 0,    Any[])    # Cloncurry, Queensland, Australia
const CPJ = Groundstation("CPJ", "Earth", +025.30,   +091.70,   0,    Any[]) 	# Cherrapunji, India (Wettest place on Earth)
const DLS = Groundstation("DLS", "Earth", +054.1452, -004.4817, 0,    Any[]) 	# Douglas, Isle of Man, UK
const DON = Groundstation("DON", "Earth", -029.0458, +115.3487, 0,    Any[]) 	# Dongara, Australia (SSC/NEN)
const DXB = Groundstation("DXB", "Earth", +025.2048, +055.2708, 0,    Any[]) 	# Dubai
const EGE = Groundstation("EGE", "Earth", +067.889663108, +021.10416625, 0, Any[]) # Esrange, Kiruna, Sweden.
const EQU = Groundstation("EQU", "Earth", +000.0   , +000.0,    0,    Any[])    # Equator
const GNP = Groundstation("GNP", "Earth", +090.0   , +000.0,    0,    Any[]) 	# Geographic north pole.
const GSP = Groundstation("GSP", "Earth", -090.0   , +000.0,    0,    Any[]) 	# Geographic south pole.
const HBK = Groundstation("HBK", "Earth", -025.8900, +027.6853, 0,    Any[])    # Hartebeesthoek Radio Astronomy Observatory, South Africa
const HUX = Groundstation("HUX", "Earth", +015.8340, -096.3199, 0,    Any[])    # Huatulco, Oaxaca, Mexico
const HWD = Groundstation("HWD", "Earth", +037.6589, -122.1219, 0,    Any[])    # Hayward Executive Airport, CA, USA
const ICN = Groundstation("ICN", "Earth", +037.5485, +126.988,  0,    Any[]) 	# Seoul, South Korea
const KLN = Groundstation("KLN", "Earth", +009.1898, +167.4243, 0,    Any[]) 	# Kwajalein Atoll
const KDK = Groundstation("KDK", "Earth", +057.435222, -152.339500, 51,    Any[]) 	# Kodiak Launch Complex, AK, USA (Pacific Spaceport Complex - Alaska) Source: http://akaerospace.com/sites/default/files/download/PSCA%20Coordinate%20location.pdf
const KRO = Groundstation("KRO", "Earth", +005.2078, -052.7724, 0,    Any[]) 	# European Spacecraft, Korou, French Guiana
const KSC = Groundstation("KSC", "Earth", +028.5241, -080.6508, 0,    Any[]) 	# Kennedy Space Center, FL, USA
const LAX = Groundstation("LAX", "Earth", +034.0500, -118.2500, 0,    Any[]) 	# Los Angeles, CA, USA
const LON = Groundstation("LON", "Earth", +051.5072, -000.1275, 0,    Any[]) 	# London, UK
const LUX = Groundstation("LUX", "Earth", +049.6116, +006.1319, 0, 	  Any[])
const LYR = Groundstation("LYR", "Earth", +078.2461, +015.4656, 0,    Any[]) 	# Svalbard, Norway
const MAD = Groundstation("MAD", "Earth", +040.4000, -003.7167, 0,    Any[]) 	# Madrid, Spain
const MEX = Groundstation("MEX", "Earth", +019.4326, -099.1332, 0,    Any[]) 	# Mexico City
const MHA = Groundstation("MHA", "Earth", -039.1605, +177.9297, 0,    Any[]) 	# Mahia Peninsula, New Zealand
const NPA = Groundstation("NPA", "Earth", +038.2450, -122.2814, 10+2, Any[]) 	# Intelsat Napa, CA
const PER = Groundstation("PER", "Earth", -031.9522, +115.8589, 0,    Any[]) 	# Perth, Australia
const PDX = Groundstation("PDX", "Earth", +045.5200, -122.6819, 0,    Any[]) 	# Portland, OR, USA
const PSK = Groundstation("PSK", "Earth", +062.9279, +040.5748, 0,    Any[]) 	# Plesetsk Cosmodrome, Russia
const SAN = Groundstation("SAN", "Earth", +032.7150, -117.1625, 0,    Any[]) 	# San Diego, CA, USA
const SDN = Groundstation("SDN", "Earth", +013.7374, +080.2351, 0,    Any[]) 	# Satish Dhawan Space Center, India
const SEA = Groundstation("SEA", "Earth", +047.6097, -122.3331, 0,    Any[]) 	# Seattle, WA, USA
const SFO = Groundstation("SFO", "Earth", +037.7833, -122.4167, 0,    Any[]) 	# San Francisco, CA, USA
const SHE = Groundstation("SHE", "Earth", +039.6655, +124.7019, 0,    Any[]) 	# Sohae, China
const SIN = Groundstation("SIN", "Earth", +001.3000, +103.8000, 0,    Any[]) 	# Singapore
const SJC = Groundstation("SJC", "Earth", +037.3382, -121.8863, 0,    Any[]) 	# San Jose, CA, USA
const STX = Groundstation("STX", "Earth", +026.0703, -097.2288, 0,    Any[]) 	# SpaceX Texas, USA
const SYD = Groundstation("SYD", "Earth", -033.8650, +151.2094, 0,    Any[]) 	# Sydney, Australia
const TMA = Groundstation("TMA", "Earth", +030.4116, +130.9357, 0,    Any[]) 	# Tanegashima, Japan
const UIO = Groundstation("UIO", "Earth", -000.2333, -078.5167, 0,    Any[]) 	# Quito, Ecuador
const UNA = Groundstation("UNA", "Earth", +031.2503, +131.0726, 0,    Any[]) 	# Uchinoura, Japan
const USH = Groundstation("USH", "Earth", -054.8   , -068.3   , 0,    Any[]) 	# Ushuaia, Argentina
const VBG = Groundstation("VBG", "Earth", +034.7417, -120.5722, 0,    Any[]) 	# Vandenburg Air Force Base, CA, USA
const VSB = Groundstation("VSB", "Earth", +008.5314, +076.8690, 0,    Any[]) 	# Vikram Sarabhai Space Center, India
const VTY = Groundstation("VTY", "Earth", +051.7673, +128.1200, 0,    Any[]) 	# Vostochny, Russia
const WCN = Groundstation("WCN", "Earth", +019.5434, +110.7970, 0,    Any[]) 	# Wenchang, China
const WPS = Groundstation("WPS", "Earth", +037.9402, -075.4664, 0,    Any[]) 	# Wallops Flight Facility, VA, USA
const WTX = Groundstation("WTX", "Earth", +031.4504, -104.7622, 0,    Any[]) 	# Blue Origin Landing Facility, West Texas
const XCN = Groundstation("XCN", "Earth", +027.8931, +102.2448, 0,    Any[]) 	# Xichang, China


# Military Installations for US389
const RSA = Groundstation("RSA", "Earth", +034.6841, -086.6487, 0,    Any[]) 	# Redstone Arsenal							AL		Huntsville
const HCH = Groundstation("HCH", "Earth", +031.5552, -110.3499, 0,    Any[]) 	# Fort Huachuca								AZ		Sierra Vista
const YPG = Groundstation("YPG", "Earth", +033.0178, -114.2525, 0,    Any[]) 	# Yuma Proving Ground						AZ		Yuma
const BAB = Groundstation("BAB", "Earth", +039.1115, -121.3599, 0,    Any[]) 	# Beale AFB									CA		Marysville
const CPK = Groundstation("CPK", "Earth", +037.7261, -121.8977, 0,    Any[]) 	# Camp Parks Reserve Forces Training Area	CA		Dublin
const NID = Groundstation("NID", "Earth", +035.6890, -117.6833, 0,    Any[]) 	# China Lake Naval Air Weapons Station		CA		Ridgecrest
const EDW = Groundstation("EDW", "Earth", +034.9240, -117.8912, 0,    Any[]) 	# Edwards AFB								CA		Rosamond
const IRW = Groundstation("IRW", "Earth", +035.2628, -116.6846, 0,    Any[]) 	# Fort Irwin								CA		Barstow
const NXP = Groundstation("NXP", "Earth", +034.2317, -116.0617, 0,    Any[]) 	# Marine Corps Air Ground Combat Center		CA		Twentynine Palms
const BKF = Groundstation("BKF", "Earth", +039.7035, -104.7618, 0,    Any[]) 	# Buckley AFB								CO		Aurora (Denver)
const SRV = Groundstation("SRV", "Earth", +038.8027, -104.5218, 0,    Any[]) 	# Schriever AFB								CO		Colorado Springs
const GDN = Groundstation("GDN", "Earth", +033.4189, -082.1403, 0,    Any[]) 	# Fort Gordon								GA		Augusta
const GUM = Groundstation("GUM", "Earth", +013.5833, +144.8494, 0,    Any[]) 	# Naval Satellite Operations Center			GU		Finegayan (Guam)
const OHU = Groundstation("OHU", "Earth", +021.5200, -157.9947, 0,    Any[]) 	# Naval Computer and Telecommunications Area Master Station, Pacific	HI	Wahiawa (Oahu Is.)
const FDK = Groundstation("FDK", "Earth", +039.4383, -077.4231, 0,    Any[]) 	# Fort Detrick								MD		Frederick
const LSV = Groundstation("LSV", "Earth", +036.2414, -115.0508, 0,    Any[]) 	# Nellis AFB								NV		Las Vegas
const NTS = Groundstation("NTS", "Earth", +037.1164, -116.1889, 0,    Any[]) 	# Nevada Test Site							NV		Amargosa Valley
const XSD = Groundstation("XSD", "Earth", +037.7947, -116.7786, 0,    Any[]) 	# Tonapah Test Range Airfield				NV		Tonapah
const CVS = Groundstation("CVS", "Earth", +034.3898, -103.3183, 0,    Any[]) 	# Cannon AFB								NM		Clovis
const WSM = Groundstation("WSM", "Earth", +033.2385, -106.3464, 0,    Any[]) 	# White Sands Missile Range					NM		White Sands
const DYS = Groundstation("DYS", "Earth", +032.4224, -099.8419, 0,    Any[]) 	# Dyess AFB									TX		Abilene
const BIF = Groundstation("BIF", "Earth", +031.8124, -106.4213, 0,    Any[]) 	# Fort Bliss								TX		El Paso
const FSH = Groundstation("FSH", "Earth", +029.4527, -098.4500, 0,    Any[]) 	# Fort Sam Houston							TX		San Antonio
const GOF = Groundstation("GOF", "Earth", +031.4328, -100.3993, 0,    Any[]) 	# Goodfellow AFB							TX		San Angelo
const SKF = Groundstation("SKF", "Earth", +029.3820, -098.5785, 0,    Any[]) 	# Kelly AFB									TX		San Antonio
const UTR = Groundstation("UTR", "Earth", +040.6571, -113.4383, 0,    Any[]) 	# Utah Test and Training Range				UT
const DAA = Groundstation("DAA", "Earth", +038.6884, -077.1464, 0,    Any[]) 	# Fort Belvoir								VA		Alexandria
const CSP = Groundstation("CSP", "Earth", +036.5600, -076.2678, 0,    Any[]) 	# Naval Satellite Operations Center			VA		Chesapeake


# Radio Astronomy Locations for US161
const GBK = Groundstation("GBK", "Earth", +038.4331, -079.8397, 0, 	  Any[]) 	# National Radio Astronomy Observatory (NRAO), Robert C. Byrd Telescope, Green Bank, WV
const VLA = Groundstation("VLA", "Earth", +034.0789, -107.6183, 0, 	  Any[]) 	# NRAO, Very Large Array, Socorro, 			NM
const UAZ = Groundstation("UAZ", "Earth", +031.9533, -111.6147, 0, 	  Any[]) 	# University of Arizona 12-m Telescope, Kitt Peak, AZ
const CTC = Groundstation("CTC", "Earth", +037.2317, -118.2933, 0, 	  Any[]) 	# Caltech Telescope, Owens Valley, 			CA
const AST = Groundstation("AST", "Earth", +042.3917, -072.3450, 0, 	  Any[]) 	# Five College Observatory, Amherst, 		MA
const WFD = Groundstation("WFD", "Earth", +042.6233, -071.4883, 0, 	  Any[]) 	# Haystack Observatory, Westford, 			MA
const JCM = Groundstation("JCM", "Earth", +019.8258, -155.4797, 0, 	  Any[]) 	# James Clerk Maxwell Telescope, Mauna Kea, HI
const CMA = Groundstation("CMA", "Earth", +037.2786, -118.1422, 0, 	  Any[]) 	# Combined Array for Research in Millimeter-wave Astronomy (CARMA), CA
const BSR = Groundstation("BSR", "Earth", +048.1311, -119.6833, 0, 	  Any[]) 	# Brewster,									WA
const DVS = Groundstation("DVS", "Earth", +030.6350, -103.9447, 0, 	  Any[]) 	# Fort Davis,								TX
const HCK = Groundstation("HCK", "Earth", +042.9336, -071.9867, 0, 	  Any[]) 	# Hancock, 									NH
const KPK = Groundstation("KPK", "Earth", +031.9564, -111.6125, 0, 	  Any[]) 	# Kitt Peak, 								AZ
const AMS = Groundstation("AMS", "Earth", +035.7750, -106.2456, 0, 	  Any[]) 	# Los Alamos, 								NM
const KEA = Groundstation("KEA", "Earth", +019.8014, -155.4556, 0, 	  Any[]) 	# Mauna Kea, 								HI
const LBT = Groundstation("LBT", "Earth", +041.7714, -091.5742, 0, 	  Any[]) 	# North Liberty, 							IA
const OWS = Groundstation("OWS", "Earth", +037.2317, -118.2769, 0, 	  Any[]) 	# Owens Valley, 							CA
const PTN = Groundstation("PTN", "Earth", +034.3011, -108.1192, 0, 	  Any[]) 	# Pie Town, 								NM
const CRX = Groundstation("CRX", "Earth", +017.7567, -064.5836, 0, 	  Any[]) 	# Saint Croix, 								VI



TLEurls = [ "visual", 
            # "amateur", 
            # "argos", 
            # "beidou", 
            # "2012-044", 
            # "cosmos-2251-debris", 
            # "cubesat", 
            # "dmc", 
            "resource", 
            # "education", 
            # "engineering", 
            # "x-comm", 
            # "1999-025", 
            # "galileo", 
            # "geodetic", 
            "geo", 
            # "globalstar", 
            # "glo-ops", 
            # "goes", 
            # "gorizont", 
            "gps-ops", 
            # "intelsat", 
            # "iridium", 
            # "iridium-33-debris", 
            # "tle-new", 
            # "military", 
            # "molniya", 
            # "nnss", 
            # "noaa", 
            # "orbcomm", 
            # "other", 
            # "other-comm", 
            # "radar", 
            # "raduga", 
            # "musson", 
            # "sbas", 
            # "sarsat", 
            # "science", 
            # "stations", 
            "tdrss", 
            "weather"]

TLEurls = ["http://www.celestrak.com/NORAD/elements/" * str * ".txt" for str in TLEurls]

push!(TLEurls, "http://ephemerides.planet-labs.com/planet_mc.tle")