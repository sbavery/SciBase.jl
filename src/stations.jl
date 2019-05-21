export
    loadStations!,
    site_eci


function loadStations!(sys)
    for key in sys["stations"]
        gs = eval(Meta.parse(key))
        sys[key]["coords"] = repeat([deg2rad(gs.lon) deg2rad(gs.lat) 1e-3gs.alt], sys["tSteps"]) # [lon, latgd, h], [rad, rad, km]
        # sys[key]["vel"] = zeros(3)
        sys[key]["frame"] = "gdc"
        if !isempty(gs.children) # Add antennas
        	sys[key]["antennas"] = gs.children
        end
    end
end


"""
   site_eci(θ_LST::T, ϕ_gd::T, h_ellp::T; R::T=R⨁, e²::T=e2⨁) where T <: Float64

Compute a site's IJK (ECI) position. N.B. Not as accurate as using full ITRF (ECEF) transformation.

### Arguments (Inputs):
| Parameter | Description                      | Units  |
| --------: | :------------------------------- | :----- |
| `θ_LST`   | Local sidereal time at site      | rad    |
| `ϕ_gd`    | Geodetic latitude of site        | rad    |
| `h_ellp`  | Altitude of site above ellipsoid | DU     |

#### Optional Keyword Arguments:
| Parameter | Description                   | Units  |
| --------: | :---------------------------- | :----- |
| `R=R⨁`    | Radius of site’s central body | DU     |
| `e²=e2⨁`  | Body’s eccentricity squared   | -      |

### Values (Outputs):
| Parameter     | Description                | Units |
| ------------: | :------------------------- | :---- |
| `r⃗_site_ijk`  | Site’s eci position vector | DU    |

Author: James Spicer (2016.11.13) <br/>
Source: Vallado (2001), Algorithm 48, p.408-9 (See errata)
"""
function site_eci(θ_LST::T, ϕ_gd::T, h_ellp::T; R::T=R⨁, e²::T=e2⨁) where T <: Float64
    sϕ  = sin(ϕ_gd)
    C   = R/√(1-e²*sϕ^2)
    S   = C*(1-e²)

    r_δ = (C + h_ellp)*cos(ϕ_gd)
    r_K = (S + h_ellp)*sϕ

   return [r_δ*cos(θ_LST) r_δ*sin(θ_LST) r_K]
end