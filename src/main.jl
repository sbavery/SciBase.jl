export
    system,
    setAzElStation

function system(; ti::DateTime=now(), tf::DateTime=now()+Dates.Day(1), tSteps::I=500, timeStep::F=300.0, bodies::Array{S}=["EARTH"], points::Array{S}=String[], constellation::Array{S}=String[], crosslinks::B=false, clients::Array{S}=String[], stations::Array{S}=String[], R::Array{F}=[1e3], report::Array{S}=String[], plot::Array{S}=["position"], showLegend::B=true, refFrame::S="eci", view::S="xy", animate::B=false, pyguiOn::B=false, azElStation::Array{S}=["SFO"], showAccess::B=true, AERpolar::B=false, showGEO::B=false) where {B <: Bool, I <: Int64, S <: String, F <: Float64}
    ti ≥ tf && error("ERROR: The start time must be before the end time.")
    t = (timeStep ≠ 300.0) ? collect(datetime2julian(ti) : timeStep/solSecsEarth : datetime2julian(tf)) : range(datetime2julian(ti), stop=datetime2julian(tf), length=tSteps)

    # global dName
    # include(dName * "/constants.jl") # Redefine constants & delete already-added antennas
    sys = Dict{String,Any}("constellation" => String[], "bodies" => String[], "clients" => String[], "stations" => String[], "points" => String[], "pairs" => String[], "links" => String[])

    for c in constellation
        sys[uppercase(c)] = Dict{String,Any}("type" => "constellation", "plot" => true)
        push!(sys["constellation"], uppercase(c))
    end
    for b in bodies
        sys[uppercase(b)] = Dict{String,Any}("type" => "body", "plot" => true)
        push!(sys["bodies"], uppercase(b))
    end
    for c in clients
        sys[uppercase(c)] = Dict{String,Any}("type" => "client", "plot" => true)
        push!(sys["clients"], uppercase(c))
    end
    for s in stations
        sys[uppercase(s)] = Dict{String,Any}("type" => "station", "plot" => true)
        push!(sys["stations"], uppercase(s))
    end
    for p in points
        sys[uppercase(p)] = Dict{String,Any}("type" => "point", "plot" => true)
        push!(sys["points"], uppercase(p))
    end

    sys["crosslinks"] = crosslinks
    sys["plot"]       = plot
    sys["report"]     = report

    sys["t"]          = t
    sys["tUnits"]     = "JD"
    sys["tSteps"]     = length(t)

    if showGEO
        sys["GEO"] = Dict{String,Any}("plot" => true)
        loadGEO!(sys)
    end

    sys["refFrame"]     = refFrame

    center = getCentralBody(sys)
    if !in(center, sys["bodies"])
        sys[center] = Dict{String,Any}("type" => "body", "plot" => false)
        push!(sys["bodies"], center)
    end

    sys["showLegend"] = showLegend
    sys["view"]       = in(sys["refFrame"], ["gdc", "gcc", "aer"]) ? "xy" : view
    sys["animate"]    = animate
    sys["pyguiOn"]    = pyguiOn
    sys["showAccess"] = showAccess
    sys["AERpolar"]   = AERpolar
    sys["AzElGS"]     = azElStation

    sanitizeSysInputs!(sys)

    setAzElStation(azElStation[1])
    loadBodies!(sys)
    loadPoints!(sys)
    loadStations!(sys)
    loadSpacecraft!(sys)
    calcStats!(sys)
    printReport(sys)
    plotCoords(sys)

    return sys
end

function setAzElStation(station)
    global azElGS = eval(Meta.parse(station))
    global azElGSGDC = [deg2rad(azElGS.lon) deg2rad(azElGS.lat) azElGS.alt/1e3] # [lon latgd h], [rad rad km]
    global azElGSECEF = ecef2gdc(azElGSGDC, inv=true) # Used in frames.jl
end