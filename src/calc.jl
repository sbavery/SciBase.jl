export
    calcBW!,
    calcCoverage!,
    calcEclipse!,
    calcLatency!,
    calcRanges!,
    calcAtmoLoss!,
    calcPathLoss!,
    calcCapacity!,
    calcTotalGroundLinkCapacity!,
    calcClientRates!,
    printReport,
    calcStats!,
    overall,
    printOverlap,
    printIndividual,
    printOverall,
    findNaNmean,
    printOverallCoverage,
    isBlackout,
    noCrosslinks,
    humanReadable!,
    findBestChain


function calcStats!(sys)
    in("range",     [sys["report"]; sys["plot"]]) && calcRanges!(sys)
    # in("beamwidth", [sys["report"]; sys["plot"]]) && calcBW!(sys)
    in("latency",   [sys["report"]; sys["plot"]]) && calcLatency!(sys)
    in("coverage",  [sys["report"]; sys["plot"]]) && calcCoverage!(sys)
    in("eclipse",   [sys["report"]; sys["plot"]]) && calcEclipse!(sys)
    in("atmoLoss",  [sys["report"]; sys["plot"]]) && calcAtmoLoss!(sys)
    if (sys["showAccess"] && sys["refFrame"] ≠ "aer" && in("position", sys["plot"])) || in("pathLoss",  [sys["report"]; sys["plot"]]); calcPathLoss!(sys);  end
    in("capacity", [sys["report"]; sys["plot"]]) && calcCapacity!(sys)
    in("totalGroundLinkCapacity", [sys["report"]; sys["plot"]]) && calcTotalGroundLinkCapacity!(sys)
    in("clientRates", [sys["report"]; sys["plot"]]) && calcClientRates!(sys)
end


function calcEclipse!(sys)
    if !haskey(sys, "SUN")
        sys["SUN"]   = Dict{Any,Any}("type" => "body")
        push!(sys["bodies"], "SUN")
    end
    if !haskey(sys, "EARTH")
        sys["EARTH"]   = Dict{Any,Any}("type" => "body")
        push!(sys["bodies"], "EARTH")
    end
    loadBodies!(sys)
    transformFrames!(sys, "eci")
    for key in [sys["clients"]; sys["constellation"]]
        e = falses(sys["tSteps"])
        for i = 1:sys["tSteps"]
            e[i] = !sight(sys[key]["coords"][i,:], sys["SUN"]["coords"][i,:])
        end
        sys[key]["eclipse"] = e
    end
end


function calcRanges!(sys)
    for A in [sys["clients"]; sys["constellation"]; sys["stations"]], B in [sys["clients"]; sys["constellation"]; sys["stations"]]
        if haskey(sys, "$(sort([A, B])[1]) <-> $(sort([A, B])[2])") && haskey(sys["$(sort([A, B])[1]) <-> $(sort([A, B])[2])"], "range"); continue; end
        if A == B || (A in sys["stations"] && B in sys["stations"]) || (A in sys["clients"] && B in sys["clients"]); continue; end

        transformFrames!(sys, "eci") # to use sight function

        d = zeros(sys["tSteps"])
        for i = 1:sys["tSteps"]
            rA = sys[A]["coords"][i,:]
            rB = sys[B]["coords"][i,:]
            d[i] = sight(rA, rB) ? norm(rA - rB) : NaN
        end
        sys["$(sort([A, B])[1]) <-> $(sort([A, B])[2])"] = Dict{String,Any}("type" => "pair", "range" => d)
        push!(sys["pairs"], "$(sort([A, B])[1]) <-> $(sort([A, B])[2])")
    end
end


function calcPathLoss!(sys)
    calcRanges!(sys)
    for A in [sys["clients"]; sys["constellation"]; sys["stations"]], B in [sys["clients"]; sys["constellation"]; sys["stations"]] # A = transmitter, B = receiver.
        if haskey(sys, "$A -> $B") && haskey(sys["$A -> $B"], "pathLoss"); continue; end
        if A == B || (A in sys["stations"] && B in sys["stations"]) || (A in sys["clients"] && B in sys["clients"]); continue; end
        if !haskey(sys[A], "antennas") || isempty(sys[A]["antennas"]); continue; end
        if !haskey(sys[B], "antennas") || isempty(sys[B]["antennas"]); continue; end
        for ant in sys[A]["antennas"]
            if ant.target == sys[B]["type"] && ant.P ≠ 0.0 # Find A's transmit antenna
                p = -pathLoss(sys["$(sort([A, B])[1]) <-> $(sort([A, B])[2])"]["range"], c0/ant.f) # Making resulting loss +ve
                if !haskey(sys, "$A -> $B")
                    sys["$A -> $B"] = Dict{String,Any}("type" => "link", "Tx" => A, "Rx" => B, "pathLoss" => p)
                    push!(sys["links"], "$A -> $B")
                else
                    sys["$A -> $B"]["pathLoss"] = p
                end
                break
            end
        end
    end
end


function calcLatency!(sys)
    calcRanges!(sys)

    for pair in sys["pairs"]
        sys[pair]["latency"] = sys[pair]["range"]./c0
    end
end


function calcCoverage!(sys)
    calcRanges!(sys)

    for pair in sys["pairs"]
        val = sys[pair]["range"]
        sys[pair]["coverage"] = (1.0 - (sum(isnan(val))/length(val)))*100.0
    end
end


function calcAtmoLoss!(sys)
    for A in [sys["clients"]; sys["stations"]; sys["constellation"]], B in [sys["clients"]; sys["stations"]; sys["constellation"]]
        if haskey(sys, "$A -> $B") && haskey(sys["$A -> $B"], "atmoLoss"); continue; end
        if A == B || (sys[A]["type"] == sys[B]["type"] && sys[A]["type"] ≠ "constellation"); continue; end
        if sort([sys[A]["type"], sys[B]["type"]]) == ["client", "station"] && !isempty(sys["constellation"]); continue; end
        if !haskey(sys[A], "antennas") || isempty(sys[A]["antennas"]); continue; end
        if !haskey(sys[B], "antennas") || isempty(sys[B]["antennas"]); continue; end
        for ant in sys[A]["antennas"]
            if ant.target == sys[B]["type"] && ant.P ≠ 0.0 # Find A's transmit antenna
                ls = zeros(sys["tSteps"])
                if A in sys["stations"] || B in sys["stations"]
                    gs = A in sys["stations"] ? A : B
                    sc = A in sys["stations"] ? B : A
                    setAzElStation(gs)
                    lon, lat, alt = transformItem(sys["t"][1], sys[gs]["coords"][1,:], sys[gs]["frame"],    "gdc") # [rad, rad, km]
                    p = transformItem(sys["t"],    sys[sc]["coords"],      sys[sc]["frame"],    "aer")
                    azSC, elSC = p[:,1], p[:,2]
                    if haskey(sys, "GEO")
                        q = transformItem(sys["t"], sys["GEO"]["coords"], sys["GEO"]["frame"], "aer")
                        azGEO, elGEO = q[:,1], q[:,2]
                    end
                    for i = 1:sys["tSteps"]
                        if 5.0 ≤ rad2deg(elSC[i]) ≤ 90.0
                            if haskey(sys, "GEO")
                                GEOangle = minimum(haversine(elSC[i], azSC[i], elGEO, azGEO)) # Smallest angle to GEO belt
                                if abs(rad2deg(GEOangle)) ≤ 5.0
                                    ls[i] = NaN
                                else
                                    ls[i] = totalAttenuation(rad2deg(lat), rad2deg(lon), ant.f/1e9; p=rand()*100, el=rad2deg(elSC[i]), h_s=alt)
                                end
                            else
                                ls[i] = totalAttenuation(rad2deg(lat), rad2deg(lon), ant.f/1e9; p=rand()*100, el=rad2deg(elSC[i]), h_s=alt)
                            end
                        else
                            ls[i] = NaN
                        end
                    end
                end
                if !haskey(sys, "$A -> $B")
                    sys["$A -> $B"] = Dict{String,Any}("type" => "link", "Tx" => A, "Rx" => B, "atmoLoss" => ls)
                    push!(sys["links"], "$A -> $B")
                    break
                else
                    sys["$A -> $B"]["atmoLoss"] = ls
                    break
                end
            end
        end
    end
end


function calcCapacity!(sys)
    calcRanges!(sys)
    calcPathLoss!(sys)
    calcAtmoLoss!(sys)

    for link in sys["links"]
        if haskey(sys[link], "capacity"); continue; end
        L = sys[link]["atmoLoss"] + sys[link]["pathLoss"] + 1.5 # Path losses. '1.5' is pol. loss + margin

        Tx = sys[link]["Tx"]
        Rx = sys[link]["Rx"]

        Rx_ant, Tx_ant = 0, 0
        for ant in sys[Tx]["antennas"] ## Find Tx transmit antenna
            if ant.target == sys[Rx]["type"] && ant.P ≠ 0
                Tx_ant = ant
                break
            end
        end
        for ant in sys[Rx]["antennas"] ## Find Rx receive antenna
            if ant.target == sys[Tx]["type"] && ant.P == 0 && ant.f == Tx_ant.f
                Rx_ant = ant
                break
            end
        end
        if Rx_ant == 0 || Tx_ant == 0; continue; end

        R_dB = zeros(sys["tSteps"])
        for i = 1:sys["tSteps"]
            R_dB[i] = isnan(L[i]) ? 0.0 : Tx_ant.G + Tx_ant.P + Rx_ant.GT - (L[i] + Tx_ant.L + Rx_ant.L) - Rx_ant.EbN0_req - pow2db(h0*Tx_ant.FEC) - Tx_ant.margin
        end
        sys[link]["capacity"] = db2pow(R_dB)
    end
end


function calcTotalGroundLinkCapacity!(sys)
    calcCapacity!(sys)
    downTot = ones(sys["tSteps"])
    upTot = ones(sys["tSteps"])

    for link in sys["links"]
        Tx = sys[link]["Tx"]
        Rx = sys[link]["Rx"]
        if sys[Tx]["type"] == "station" && (sys[Rx]["type"] == "client" || sys[Rx]["type"] == "constellation")
            upTot += sys[link]["capacity"]
        elseif sys[Rx]["type"] == "station" && (sys[Tx]["type"] == "client" || sys[Tx]["type"] == "constellation")
            downTot += sys[link]["capacity"]
        end
    end
    sys["totalUplinkCapacity"] = upTot
    sys["totalDownlinkCapacity"] = downTot
end


function calcClientRates!(sys)
    calcCapacity!(sys)

    for key in [sys["clients"]; sys["links"]]
        sys[key]["dataRate"] = ones(sys["tSteps"])
    end

    for i = 1:sys["tSteps"]
        for cli in sys["clients"] # Find client's best link at this timestep. Push data through it
            bestChain = findBestChain(sys, cli, i)
            for link in bestChain
                sys[link]["dataRate"][i] += sys[bestChain[1]]["capacity"][i]
            end
            if !isempty(bestChain); sys[cli]["dataRate"][i] = sys[bestChain[1]]["capacity"][i]; end
        end
    end
end

# TODO: Currently chains are max 3 links (4 nodes) long. Make recursive for longer chains?
function findBestChain(sys, cli, i) # Find a chain between the cli and the ground at timestep i
    chainNames = []
    chainCaps = Float64[]

    for link1 in sys["links"] # First check for direct client -> gs links
        if sys[link1]["Tx"] ≠ cli || sys[link1]["capacity"][i] == 1; continue; end

        if sys[link1]["Rx"] in sys["stations"]
            push!(chainNames, String[link1])
            push!(chainCaps, sys[link1]["capacity"][i])
        end
    end

    if !isempty(chainNames); return chainNames[indmax(chainCaps)]; end

    for link1 in sys["links"] # Now check for chains 2 links long
        if sys[link1]["Tx"] ≠ cli || sys[link1]["capacity"][i] == 1; continue; end

        for link2 in sys["links"]
            if sys[link2]["Tx"] ≠ sys[link1]["Rx"]; continue; end
            if sys[link2]["dataRate"][i] + sys[link1]["capacity"][i] > sys[link2]["capacity"][i]; continue; end # Check this link has room for new data rate

            if sys[link2]["Rx"] in sys["stations"]
                push!(chainNames, String[link1, link2])
                push!(chainCaps, min(sys[link1]["capacity"][i], sys[link2]["capacity"][i]))
            end
        end
    end

    if !isempty(chainNames); return chainNames[indmax(chainCaps)]; end

    for link1 in sys["links"] # Now check for chains 3 links long
        if sys[link1]["Tx"] ≠ cli || sys[link1]["capacity"][i] == 1; continue; end

        for link2 in sys["links"]
            if sys[link2]["Tx"] ≠ sys[link1]["Rx"]; continue; end
            if sys[link2]["dataRate"][i] + sys[link1]["capacity"][i] > sys[link2]["capacity"][i]; continue; end # Check this link has room for new data rate

            for link3 in sys["links"]
                if sys[link3]["Tx"] ≠ sys[link2]["Rx"]; continue; end
                if sys[link3]["dataRate"][i] + sys[link1]["capacity"][i] > sys[link3]["capacity"][i]; continue; end # Check this link has room for new data rate

                if sys[link3]["Rx"] in sys["stations"]
                    push!(chainNames, String[link1, link2, link3])
                    push!(chainCaps, min(sys[link1]["capacity"][i], sys[link2]["capacity"][i], sys[link3]["capacity"][i]))
                end
            end
        end
    end

    if !isempty(chainNames)
        return chainNames[indmax(chainCaps)] # Returns best chain tuple, e.g. ("SDO -> AUD1", "AUD1 -> SFO")
    else
        # println("No Chain")
        return ()
    end
end





function humanReadable!(sys)
    for key in keys(sys)
        if isa(sys[key], Dict) && haskey(sys[key], "frame") && in(sys[key]["frame"], ["gdc", "gcc", "aer"])
            if sys[key]["frame"] == "aer" && sys["AERpolar"]
                sys[key]["coords"][:,2] = 90.0 - rad2deg.(sys[key]["coords"][:,2])
                # sys[key]["coords"][:,1] = rad2deg(sys[key]["coords"][:,1])
            else
                sys[key]["coords"][:,1:2] = rad2deg.(sys[key]["coords"][:,1:2])
            end
        end
    end
    if sys["tUnits"] == "JD"
        sys["t"] = julian2datetime.(sys["t"])
        sys["tUnits"] = "dateTime"
    end
end




### Stuff below here doesn't work

function calcBW!(sys)
    sys["bw"] = Dict()
    for key in keys(sys["dists"])
        if in(key[1], sys["constellation"])
            if in(key[2], sys["clients"])
                sys["bw"][key] = rad2deg.(2*atan((R⨁+1000)./sys["dists"][key]))
            else
                sys["bw"][key] = rad2deg.(2*atan(250 ./sys["dists"][key]))
            end
        end
        if in(key[2], sys["constellation"])
            if in(key[1], sys["clients"])
                sys["bw"][key] = rad2deg.(2*atan((R⨁+1000)./sys["dists"][key]))
            else
                sys["bw"][key] = rad2deg.(2*atan(250 ./sys["dists"][key]))
            end
        end
    end
end


function inSameCategory(sys, A, B)
    for i in ["stations", "constellation", "clients"]
        if in(A, sys[i]) && in(B, sys[i])
            return true
        end
    end
    return false
end


# Returns true if there is a constellation specified (so don't report on client-station links).
function notInterested(sys, A, B)
    if !isempty(sys["constellation"])
        ((in(B, sys["clients"]) && in(A, sys["stations"])) || (in(A, sys["clients"]) && in(B, sys["stations"]))) && return true
    end
    return false
end


function printReport(sys)
    if !isempty(sys["report"])
        printIndividual(sys)
        println()
        overall(sys)
    end
end


function overall(sys)
    println("Overall:")
    for r in sys["report"]
        if r == "coverage"
            printOverallCoverage(sys)
            printOverlap(sys)
        elseif r == "range"
            # printOverall(sys["pairs"])
        elseif r == "latency"
            printOverall(sys["lats"])
        elseif r == "beamwidth"
            printOverall(sys["bw"])
        elseif r == "eclipse"
            continue
        else
            error("ERROR: '", r, "' not defined in 'report'.")
        end
    end
end


function printOverall(results)
    maxs, mins, avgs = Float64[], Float64[], Float64[]
    for key in keys(results)
        push!(maxs, maximum(results[key]))
        push!(mins, minimum(results[key]))
        push!(avgs, NaNmean(results[key]))
    end
    println(round(minimum(mins),1), " ", round(maximum(maxs),1), " ", round(NaNmean(avgs),1))
end


function printOverlap(sys)
    numReal = 0
    for key in keys(sys["dists"])
        for i = 1:sys["tSteps"]
            if !isnan(sys["dists"][key][i])
                numReal += 1
            end
        end
    end
    println("LoS fraction between all nodes: ", round(100*numReal/sys["tSteps"], 1), "%")
end


#=
function printGSEIRP(systemDists, constellation, stations, R)
    for a in constellation
        for s in stations
           key = sort([a, s])
            val = systemDists[key]
            max = maximum(val)
            min = minimum(val)
            avg = mean(val)

            theta = rad2deg(2*atan((reEarth+1000)/S_E_M))
            G = pow2db((70*pi/theta)^2)
            GB       = G; # [dB]

            EbN0_req = 8; # [dB]

            f        = 3.3e10          # [Hz]
            T_sys    = 60              # [K]

            tR = R./3                    # 1:10 JPEG compression applied
            tR = R./0.5                  # Encoding (Turbo / LDPC Codes)

            # EIRP = GA + pow2db(PA);
            println()
            println(key[1], " -> ", key[2], ":")
            for j in [max, avg, min]
                L_S = pathLoss(j, c0/f);    # Path loss [dB]
                println()
                for i = 1:length(R)
                    print("   ", round(R[i]), ": ")
                    EIRP = EbN0_req - L_S - GB + pow2db(h * T_sys * tR[i]); # [dB]
                    print(round(EIRP))
                end
            end
        end
    end
end
=#


function printIndividual(sys)
    if in("eclipse", sys["report"])
        for c in sys["constellation"]
            println(c, ":")
            println("   Eclipse: ", countnz(sys["eclipse"][c])*100/length(sys["eclipse"][c]), "%")
            println()
        end
        if sys["report"] == ["eclipse"]
            return
        end
    end

    for key in sys["pairs"]
        # println(key[1], " - ", key[2], ":")
        println(key, ":")
        for r in sys["report"]
            if r == "coverage"
                println("   Coverage: ", round(sys["covs"][key], 1), "%")
            elseif r == "range"
                println(round(minimum(sys[key]["range"]),1), " ", round(maximum(sys[key]["range"]),1), " ", round(NaNmean(sys[key]["range"]),1))
            elseif r == "latency"
                println(round(minimum(sys["lats"][key]),1), " ", round(maximum(sys["lats"][key]),1), " ", round(NaNmean(sys["lats"][key]),1))
            elseif r == "beamwidth"
                println(round(minimum(sys["bw"][key]),1), " ", round(maximum(sys["bw"][key]),1), " ", round(NaNmean(sys["bw"][key]),1))
            elseif r == "eclipse"
                continue
            else
                error("ERROR: '", r, "' not defined in 'report'.")
            end
        end
        println()
    end
end


function NaNmean(val)
    sum = 0
    nonNaNcount = 0;
    for v in val
        if !isnan(v)
            sum += v
            nonNaNcount += 1
        end
    end
    sum/nonNaNcount
end


function printOverallCoverage(sys)
    blackout = 0
    for i = 1:sys["tSteps"]
        if !isempty(sys["constellation"])
            if isempty(sys["stations"])
                if isBlackout(i, sys["dists"], sys["clients"], sys["constellation"])
                    blackout +=1
                end
            elseif isempty(sys["clients"])
                if isBlackout(i, sys["dists"], sys["constellation"], sys["stations"])
                    blackout +=1
                end
            else # Calculate coverage from clients, constellation, and stations
                if sys["crosslinks"]
                    if isBlackout(i, sys["dists"], sys["clients"], sys["constellation"]) || isBlackout(i, sys["dists"], sys["stations"], sys["constellation"])
                       blackout += 1
                    end
                else
                    if noCrosslinks(sys, i)
                        blackout += 1
                    end
                end
            end
        else
            if isBlackout(i, sys["dists"], sys["clients"], sys["stations"])
                blackout +=1
            end
        end
    end
    println("Coverage: ", round((1 - (blackout/sys["tSteps"]))*100, 1), "%")
end


function isBlackout(i, dists, arrA, arrB)
    (isempty(arrA) || isempty(arrB)) && return false
    for A in arrA
        for B in arrB
            !isnan(dists[sort([A, B])][i]) && return false
        end
    end
    return true
end


function noCrosslinks(sys, i)
    for c in sys["clients"]
        for a in sys["constellation"]
            for s in sys["stations"]
                if !isnan(sys["dists"][sort([c, a])][i]) && !isnan(sys["dists"][sort([a, s])][i])
                    return false
                end
            end
        end
    end
    return true
end


# function lineOfSight(x1, x2) # Formula: intersect of line x1+t*xd with ellipse. By JPWS
#     xd = x2 - x1
#     M = [reEarth_B^2, reEarth_B^2, reEarth_A^2]
#     a = dot(M, xd.^2)
#     b = 2 * dot(M, x1.*xd)
#     c = dot(M, x1.^2) - (reEarth_A*reEarth_B)^2
#     disc = b^2 - 4*a*c
#     if disc > 0
#         t1 = (-b + sqrt(disc))/(2a)
#         t2 = (-b - sqrt(disc))/(2a)
#         if t1 > 0 && t2 > 0
#             p = x1 + min(t1,t2)*(x2-x1)
#             if norm(x1-x2) > norm(x1-p)
#                 return false
#             end
#         end
#     end
#     return true
# end