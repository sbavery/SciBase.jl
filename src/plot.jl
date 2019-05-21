export
    plotCoords,
    plotPosn,
    plotRange,
    # plotBW,
    plotAtmoLoss,
    plotPathLoss,
    plotCapacity,
    plotTotalGroundLinkCapacity,
    plotClientRates,
    drawCircle,
    drawArc,
    planeTo3D,
    drawCircle3D,
    drawSquare,
    drawRectangle,
    roundedRectangle,
    drawSegment,
    drawSphere,
    drawPolygon,
    draw3DPolygon,
    drawCuboid,
    Arrow3D,
    drawAxes3D,
    verticalLine,
    horizontalLine,
    verticalBand,
    horizontalBand,
    drawCone,
    meshgrid,
    sanitizeLinestyle,
    maybeAxisLims,
    axisEqual3D,
    axisOff2,
    drawGSEqProj,
    drawBody,
    # drawGEO!,
    getCentralBody,
    splitPlot,
    plot_footprint2D,
    flattenTransparentImage,
    mimicAlpha


patch   = pyimport("matplotlib.patches")
proj3d  = pyimport("mpl_toolkits.mplot3d.proj3d")
art3d   = pyimport("mpl_toolkits.mplot3d.art3d")
ndimage = pyimport("scipy.ndimage")
mplot3d = pyimport("mpl_toolkits.mplot3d")


function plotCoords(sys::Dict{String,Any})
    isempty(sys["plot"]) && return
    in("position", sys["plot"]) && transformFrames!(sys, sys["refFrame"])
    humanReadable!(sys) # Convert radians to degrees, time to DateTime, etc.
    N = sys["animate"] ? rangeInRange(1, sys["tSteps"]) : [1] # Only export 100 frames no matter how high tSteps
    for n in N
        clf()
        pygui(sys["pyguiOn"])
        fig = figure(figsize=(8, 6length(sys["plot"])))
        for i = 1:length(sys["plot"])
            if sys["plot"][i] == "position"
                # spl = parse(Int, string(length(sys["plot"]))*"1"*string(i))
                spl = 100length(sys["plot"]) + 10 + i
                if !in(sys["refFrame"], ["gdc", "gcc", "aer"])
                    using3D()
                    PyPlot.PyObject(PyPlot.axes3D)
                    PyPlot.pyimport("mpl_toolkits.mplot3d.axes3d")
                    ax = subplot(spl, projection="3d")
                elseif sys["refFrame"] == "aer" && sys["AERpolar"]
                    ax = subplot(spl, polar="true")
                else
                    ax = subplot(spl)
                end
                plotPosn(sys, n, ax)
            else
                ax = subplot(parse(Int, string(length(sys["plot"]))*"1"*string(i)))
                sys["plot"][i] == "range"                   && plotRange(sys, n, ax)
                # sys["plot"][i] == "beamwidth" && plotBW(sys)
                sys["plot"][i] == "atmoLoss"                && plotAtmoLoss(sys, n, ax)
                sys["plot"][i] == "pathLoss"                && plotPathLoss(sys, n, ax)
                sys["plot"][i] == "capacity"                && plotCapacity(sys, n, ax)
                sys["plot"][i] == "totalGroundLinkCapacity" && plotTotalGroundLinkCapacity(sys, n, ax)
                sys["plot"][i] == "clientRates"             && plotClientRates(sys, n, ax)
            end
        end
        str  = "Plots/Orbit_Plots/" * replace(string(now()), ":" => "")
        str2 = length(N)==1 ? str * ".png" : str * "-" * string(n) * ".png"
        savefig(str2, dpi=300, bbox_inches="tight", pad_inches=0.0, transparent=false)
    end
end

# Used this version of plotCoords to get the 1x3 AzEl plot view with AzElStation=["SFO", "MAD", "SIN"]
# function plotCoords(sys)
#     if sys["animate"]
#         if sys["tSteps"] > 100 # Only export 100 frames no matter how high tSteps
#             q = ceil(sys["tSteps"]/100)
#             a = collect(1:sys["tSteps"])
#             N = filter(a -> a .% q == 0, a)
#         else
#             N = collect(1:sys["tSteps"])
#         end
#     else
#         N = [1]
#     end
#     for n in N
#         clf()
#         fig = figure(figsize=(8*length(sys["AzElGS"]), 6))
#         for i = 1:length(sys["AzElGS"])
#             loadSpacecraft!(sys)
#             setAzElStation(sys["AzElGS"][i])
#             transformFrames!(sys, sys["refFrame"])
#             # humanReadable!(sys) # Convert radians to degrees, etc.

#             ax = subplot(parse(Int, "13"*string(i)), polar="true")
#             plotPosn(sys, n, ax)
#         end
#         str = "Plots/Orbit_Plots/" * replace(string(now()), ":", "")
#         str2 = length(N)==1 ? str * ".png" : str * "-" * string(n) * ".png"
#         savefig(str2, dpi=150, bbox_inches="tight", pad_inches=0.0)
#     end
# end


function getPlotInds(view) # "txyz" => [-3, 1, 2, 3]
    inds = Int64[]
    for c in view; push!(inds, c); end
    sort!(inds) .- Int('w')
end


function drawBody(sys::Dict{String,Any}, body, center; zorder=10, ax=ax)
    body == 0 && return
    # println(body)

    if body == "EARTH"
        # if sys["view"] == "xy" && !sys["pyguiOn"]
        #     println("hi")
        #     global dName
        #     im = imread(dName * "/data/earthInfiniteTop.png")
        #     println(R⨁*[-1, 1, -1, 1])
        #     imshow(im, zorder=zorder, extent=R⨁*[-1, 1, -1, 1])
        #     return
        # elseif (sys["view"] == "yz" || sys["view"] == "xz") && !sys["pyguiOn"]
        #     global dName
        #     im = imread(dName * "/data/earthInfiniteFront.png")
        #     imshow(im, zorder=zorder, extent=R⨁*[-1, 1, -1, 1])
        #     return
        # else
            fc = "lightBlue"
            r = R⨁
            ec = "blue"
            if length(center) == 3
                drawSphere(r, color=fc, edgecolor=ec, ax=ax, zorder=zorder, linewidth=0.5)
            else
                drawCircle(r, color=fc, edgecolor=ec, ax=ax, zorder=zorder, linewidth=0.5)
            end
        # end
    else
        if      body == "MOON";     fc = "gray";      r = reMoon;     ec = "gray"
        elseif  body == "SUN";      fc = "yellow";    r = R☉;         ec = "black"
        elseif  body == "MERCURY";  fc = "gray";      r = reMercury;  ec = "black"
        elseif  body == "VENUS";    fc = "white";     r = reVenus;    ec = "orange"
        elseif  body == "MARS";     fc = "orange";    r = reMars;     ec = "k"
        elseif  body == "JUPITER";  fc = "yellow";    r = reJupiter;  ec = "red"
        elseif  body == "SATURN";   fc = "gold";      r = reSaturn;   ec = "yellow"
        elseif  body == "URANUS";   fc = "cyan";      r = reUranus;   ec = "blue"
        elseif  body == "NEPTUNE";  fc = "blue";      r = reNeptune;  ec = "lightBlue"
        end
        if length(center) == 3
            drawSphere(r, color=fc, edgecolor=ec, ax=ax, zorder=zorder, linewidth=0.5)
        else
            drawCircle(r, color=fc, edgecolor=ec, ax=ax, zorder=zorder, linewidth=0.5)
        end
    end
end


function getCentralBody(sys::Dict{String,Any})
    in(sys["refFrame"], ["eci", "ecl", "ecef", "EMsyn", "aer", "gdc", "gcc"]) && return "EARTH"
    in(sys["refFrame"], ["hci", "ESsyn"])                                     && return "SUN"
    in(sys["refFrame"], ["mci", "mcl"])                                       && return "MARS"
    @error "Could not resolve central body for this system."
end


function plotPosn(sys::Dict{String,Any}, n, ax)
    if in(sys["refFrame"], ["gdc", "gcc", "aer"])
        plotPosn2D(sys, n, ax)
    else
        plotPosn3D(sys, n, ax)
    end
end


function splitPlot(x, y, ∆; color="x", kwargs...)# color="x", hatch="", linewidth=1, linestyle="-", label="", alpha=1, zorder=0)
    d = [0.0; findall(abs.(diff(x)) .≥ ∆); length(x)] # Split at x discontinuities
    d = convert(Array{Int64,1}, d)
    if color == "x"
        p = plot(x[1], y[1], zorder=-1) # Get color so all plots in for loop are same color
        color = p[1].get_color()
    end
    obj = ""
    for k = 1:length(d)-1
        xp = x[collect(d[k]+1:d[k+1])]
        yp = y[collect(d[k]+1:d[k+1])]
        obj = plot(xp, yp; color=color, kwargs...)#, color=color, lw=linewidth, linestyle=linestyle, label=label, alpha=alpha, zorder=zorder)
    end
    return obj
end


function plotPosn2D(sys::Dict{String,Any}, n, ax)
    global azElGS

    # ax = sys["refFrame"] == "aer" && sys["AERpolar"] ? axes(polar = "true") : gca()
    # ax[:hold](true)
    inds = getPlotInds(sys["view"])
    legTexts = AbstractString[]
    legObjs = []

    for key in keys(sys)
        !isa(sys[key], Dict) && continue
        !haskey(sys[key], "plot") && continue
        !sys[key]["plot"] && continue
        sys["refFrame"] == "aer" && haskey(sys[key], "type") && sys[key]["type"] == "station" && continue
        key == getCentralBody(sys) && continue
        if key == "GEO" || sys[key]["type"] ≠ "station"
            x = in(-3, inds) ? sys["t"] : sys[key]["coords"][:,inds[1]]
            y = sys[key]["coords"][:,inds[2]]

            if key == "GEO"
                obj  = splitPlot(x, y, 240.0, color="r", linewidth=3, linestyle="--", zorder=3)
            else
                asdf = splitPlot(x, y, 180.0, zorder=3)#, color="k")
                obj  = scatter(x[n], y[n], color=ax.lines[end].get_color(), zorder=10, s=30)
            end
            push!(legTexts, key)
            legObjs = [legObjs; obj]
        else # If it's a station in a ground track plot
            va = (sys[key]["coords"][1,inds[2]] > 0.0) ? "bottom" : "top"
            scatter(sys[key]["coords"][1,inds[1]], sys[key]["coords"][1,inds[2]], color="k", zorder=10, s=30)
            text(sys[key]["coords"][1,inds[1]], sys[key]["coords"][1,inds[2]] + sign(sys[key]["coords"][1,inds[2]])*10, key, size=9, zorder=0, color="k", ha="center", va=va, bbox=Dict("facecolor"=>"wheat", "alpha"=>0.7, "boxstyle"=>"round"))
        end
    end

    if sys["refFrame"] == "aer" && sys["AERpolar"]
        ax.set_rgrids(collect(30:30:90))
        ax.set_yticklabels(["60°", "30°", "0°"])
        ylim([0, 90])
        ax.set_xticklabels(["N", "NE", "E", "SE", "S", "SW", "W", "NW"])
        ax.set_theta_zero_location("N")
        # ax[:spines]["polar"][:set_color]("none")
        ax.set_axisbelow(true)
        title("Azimuth-Elevation from $(azElGS.name)")

    elseif in(sys["refFrame"], ["gdc", "gcc"])
        # global dName
        im = imread(joinpath(datapath, "earthBackground.png"))#dName * "/data/earthBackground.png")
        imshow(im, zorder=-1, extent=[-180, 180, -90, 90], alpha=0.7)
        ylabel("Latitude [°]")
        xlabel("Longitude [°]")
        xticks(collect(-180:60:180))
        yticks(collect(- 90:30:90))
        xlim([-180, 180])
        ylim([-90, 90])
        plotAccess(sys, inds, n)

    elseif sys["refFrame"] == "aer"
        ylabel("Elevation [°]")
        # xlabel("Azimuth [°]")
        xlabel("Azimuth")
        # ax[:set_xticks](collect(0:60:360))
        xticks(collect(0:90:360))
        ax.set_xticklabels(["North", "East", "South", "West", "North"])
        yticks(collect(0:30:90))
        xlim([0, 360])
        ylim([0, 90])
        lightBlue = rgb(135.0, 206.0, 250.0)
        # ax[:set_axis_bgcolor](lightBlue)
        title("Azimuth-Elevation from $(azElGS.name)")
    end
    grid(true, ls="--")
    ax.set_aspect("equal")
    # ax[:hold](false)
    showLeg(sys["showLegend"], legTexts, legObjs=legObjs)
end


function plotAccess(sys, inds, n)
    if !in(4, inds) && sys["showAccess"]
        for key in keys(sys)
            if isa(sys[key], Dict) && haskey(sys[key], "type") && sys[key]["type"] == "link"
                Tx = sys[key]["Tx"]
                Rx = sys[key]["Rx"]

                col = "r"
                if (Tx in sys["clients"] && Rx in sys["constellation"]) || (Rx in sys["clients"] && Tx in sys["constellation"])
                    col = "SeaGreen"
                elseif Tx in sys["constellation"] && Rx in sys["constellation"]
                    col = "DarkViolet"
                elseif (Tx in sys["stations"] && Rx in sys["constellation"]) || (Tx in sys["constellation"] && Rx in sys["stations"])
                    col = "Teal"
                end

                if !isnan(sys["$(Tx) -> $(Rx)"]["pathLoss"][n]) # If they have line of sight
                    if length(inds) == 2
                        if sys[Tx]["coords"][n,inds[1]] < sys[Rx]["coords"][n,inds[1]]
                            x1 = sys[Tx]["coords"][n,inds[1]]
                            y1 = sys[Tx]["coords"][n,inds[2]]
                            x2 = sys[Rx]["coords"][n,inds[1]]
                            y2 = sys[Rx]["coords"][n,inds[2]]
                        else
                            x1 = sys[Rx]["coords"][n,inds[1]]
                            y1 = sys[Rx]["coords"][n,inds[2]]
                            x2 = sys[Tx]["coords"][n,inds[1]]
                            y2 = sys[Tx]["coords"][n,inds[2]]
                        end

                        if abs(x1 - x2) > 180
                            y3 = (180y1+180y2+x1*y2-x2*y1)/(x1-x2+360) # Find point at which access line intersects vertical axis
                            plot([-180.0, x1], [y3, y1], color=col, zorder=8)
                            plot([x2,  180.0], [y2, y3], color=col, zorder=8)
                        else
                            plot([x1, x2], [y1, y2], color=col, zorder=8)
                        end
                    else
                        plot3D([sys[Tx]["coords"][n,1],sys[Rx]["coords"][n,1]], [sys[Tx]["coords"][n,2],sys[Rx]["coords"][n,2]], [sys[Tx]["coords"][n,3],sys[Rx]["coords"][n,3]], color=col, zorder=8)
                    end
                end
            end
        end
    end
end


function plotPosn3D(sys, n, ax)
    using3D()

    centralBody = getCentralBody(sys)
    legTexts = String[]
    legObjs = []

    for key in keys(sys)
        # if !in(key, [sys["stations"]; sys["clients"]; sys["constellation"]; sys["bodies"]; sys["points"]]); continue; end
        !isa(sys[key], Dict) && continue
        !haskey(sys[key], "plot") && continue
        !sys[key]["plot"] && continue

        if key == centralBody
            drawBody(sys, key, [0, 0, 0], zorder=1, ax=ax)
            continue
        end

        x = sys[key]["coords"][:,1]
        y = sys[key]["coords"][:,2]
        z = sys[key]["coords"][:,3]
        if sys["refFrame"] == "ecef" && key in sys["stations"]
            p = plot3D([x[1]], [y[1]], [z[1]], zorder=3)
            obj = scatter3D(x[1], y[1], z[1], color=p[1].get_color(), s=30, zorder=10)
        else
            if key =="GEO"
                obj = plot3D(x, y, z, "r", lw=3, ls="--", zorder=3)
            else
                p = plot3D(x, y, z, zorder=3)#, "k")
                # if sys[key]["type"] ≠ "body"
                    # obj = scatter3D(x[n], y[n], z[n], color=ax.lines[end].get_color(), s=30, zorder=10)
                obj = scatter3D(x[n], y[n], z[n], color=p[1].get_color(), s=30, zorder=10)
                # else
                #     drawBody(sys, key, [x[n], y[n], z[n]], zorder=10, ax=ax)
                # end
            end
            scatter3D([minimum(x), maximum(x)], [minimum(y), maximum(y)], [minimum(z), maximum(z)], color="none", zorder=0, s=35) # Blank dots to make sure limits are set for worst cases.
        end
        push!(legTexts, key)
        legObjs = [legObjs; obj]
    end
    plotAccess(sys, [1, 2, 3], n)

    # legObjs, legTexts = drawGEO!(sys, [1, 2, 3], legObjs, legTexts)
    # ax[:hold](false)

    sys["view"] == "xy" && ax.view_init(elev=90, azim=-90)
    sys["view"] == "yz" && ax.view_init(elev=0 , azim=-90)
    sys["view"] == "xz" && ax.view_init(elev=0 , azim=0  )

    axis("off")
    axis("tight")
    axisEqual3D(ax)
    showLeg(sys["showLegend"], legTexts, legObjs=legObjs)
end


function plotRange(sys, n, ax)
    leg = AbstractString[]
    # ax=gca()
    # ax[:hold](true)
    for pair in sys["pairs"]
        push!(leg, pair)
        plot(sys["t"], sys[pair]["range"])
        scatter(sys["t"][n], sys[pair]["range"][n], color=ax[:lines][end][:get_color](), s=30)
        scatter([minimum(sys["t"]), maximum(sys["t"])], [minimum(sys[pair]["range"]), maximum(sys[pair]["range"])], color="none", zorder=0, s=35)
    end
    # ax[:hold](false)
    xlim([sys["t"][1], sys["t"][end]])
    ylim([max(0.0, ylim()[1]), ylim()[2]])
    verticalLine(sys["t"][n], ax=ax)
    showLeg(sys["showLegend"], leg, loc=4)
    grid(true, ls="--")
    ylabel("Range [km]")
end


function plotAtmoLoss(sys, n, ax)
    leg = AbstractString[]
    # ax=gca()
    # ax[:hold](true)
    for link in sys["links"]
        push!(leg, link)
        plot(sys["t"], sys[link]["atmoLoss"])
        scatter(sys["t"][n], sys[link]["atmoLoss"][n], color=ax.lines[end].get_color(), s=30)
        scatter([minimum(sys["t"]), maximum(sys["t"])], [minimum(sys[link]["atmoLoss"]), maximum(sys[link]["atmoLoss"])], color="none", zorder=0, s=35)
    end
    # ax[:hold](false)
    xlim([sys["t"][1], sys["t"][end]])
    ylim([max(0.0, ylim()[1]), min(40.0, ylim()[2])])
    verticalLine(sys["t"][n])
    showLeg(sys["showLegend"], leg)
    grid(true, ls="--")
    ylabel("Atmospheric Losses [dB]")
end


function plotCapacity(sys, n, ax)
    leg = AbstractString[]
    for link in sys["links"]
        push!(leg, link)
        semilogy(sys["t"], sys[link]["capacity"])
        scatter(sys["t"][n], sys[link]["capacity"][n], color=ax[:lines][end][:get_color](), s=30)
        scatter([minimum(sys["t"]), maximum(sys["t"])], [minimum(sys[link]["capacity"]), maximum(sys[link]["capacity"])], color="none", zorder=0, s=35)
    end
    xlim([sys["t"][1], sys["t"][end]])
    ylim([max(1e0, ylim()[1]), ylim()[2]])
    verticalLine(sys["t"][n])
    showLeg(sys["showLegend"], leg, loc=4)
    grid(true, ls="--")
    ylabel("Link Capacity [bps]")
end


function plotPathLoss(sys, n, ax)
    leg = AbstractString[]
    # ax=gca()
    # ax[:hold](true)
    for link in sys["links"]
        push!(leg, link)
        plot(sys["t"], sys[link]["pathLoss"])
        scatter(sys["t"][n], sys[link]["pathLoss"][n], color=ax[:lines][end][:get_color](), s=30)
        scatter([minimum(sys["t"]), maximum(sys["t"])], [minimum(sys[link]["pathLoss"]), maximum(sys[link]["pathLoss"])], color="none", zorder=0, s=35)
    end
    # ax[:hold](false)
    xlim([sys["t"][1], sys["t"][end]])
    ylim([max(0, ylim()[1]), ylim()[2]])
    verticalLine(sys["t"][n])
    showLeg(sys["showLegend"], leg)
    grid(true, ls="--")
    ylabel("Path Loss [dB]")
end


function plotTotalGroundLinkCapacity(sys, n, ax)
    leg = AbstractString[]
    # ax=gca()
    # ax[:hold](true)
    semilogy(sys["t"], sys["totalUplinkCapacity"], label="Total Uplink")
    scatter(sys["t"][n], sys["totalUplinkCapacity"][n], color=ax[:lines][end][:get_color](), s=30)
    semilogy(sys["t"], sys["totalDownlinkCapacity"], label="Total Downlink")
    scatter(sys["t"][n], sys["totalDownlinkCapacity"][n], color=ax[:lines][end][:get_color](), s=30)
    scatter([minimum(sys["t"]), maximum(sys["t"])], [minimum(sys["totalDownlinkCapacity"]), maximum(sys["totalDownlinkCapacity"])], color="none", zorder=0, s=35)
    scatter([minimum(sys["t"]), maximum(sys["t"])], [minimum(sys["totalUplinkCapacity"]), maximum(sys["totalUplinkCapacity"])], color="none", zorder=0, s=35)
    xlim([sys["t"][1], sys["t"][end]])
    ylim([max(1e0, ylim()[1]), ylim()[2]])
    verticalLine(sys["t"][n])
    # ax[:hold](false)
    grid(true, ls="--")
    if sys["showLegend"]; legend(fancybox="true", ncol=2, fontsize=8, loc=4); end
    ylabel("Total Ground Link Capacity [bps]")
end


function plotClientRates(sys, n, ax)
    leg = AbstractString[]
    # ax=gca()
    # ax[:hold](true)
    for cli in sys["clients"]
        push!(leg, cli)
        semilogy(sys["t"], sys[cli]["dataRate"])
        scatter(sys["t"][n], sys[cli]["dataRate"][n], color=ax[:lines][end][:get_color](), s=30)
        scatter([minimum(sys["t"]), maximum(sys["t"])], [minimum(sys[cli]["dataRate"]), maximum(sys[cli]["dataRate"])], color="none", zorder=0, s=35)
    end
    verticalLine(sys["t"][n])
    # ax[:hold](false)
    grid(true)
    ylim([max(1e0, ylim()[1]), ylim()[2]])
    xlim([sys["t"][1], sys["t"][end]])
    verticalLine(sys["t"][n])
    showLeg(sys["showLegend"], leg)
    ylabel("Client Data Rate [bps]")
end


# function plotBW(sys)
#     leg = AbstractString[]
#     ax=gca()
#     ax[:hold](true)
#     for A in keys(sys["bw"])
#         push!(leg, A[1] * " - " * A[2])
#         plot(julian2datetime(sys["t"]), sys["bw"][A])
#     end
#     ax[:hold](false)
#     showLeg(sys["showLegend"], leg)
# end


function showLeg(showLegend, legTexts; legObjs=0, loc=1)
    if showLegend && length(legTexts) ≠ 0
        if legObjs ≠ 0
            legend(legObjs, legTexts, scatterpoints=1, fancybox="true", ncol=2, fontsize=8, loc=loc)
        else
            legend(legTexts, fancybox="true", scatterpoints=1, ncol=2, fontsize=8, loc=loc)
        end
    end
end


### Stuff below here not part of orbit plotting tool.


function axisEqual3D(ax=gca())
    (xMn, xMx), (yMn, yMx), (zMn, zMx) = xlim(), ylim(), zlim()
    max_range = 0.5max(xMx-xMn, yMx-yMn, zMx-zMn)
    xlim([mean([xMx, xMn]) - max_range, mean([xMx, xMn]) + max_range])
    ylim([mean([yMx, yMn]) - max_range, mean([yMx, yMn]) + max_range])
    zlim([mean([zMx, zMn]) - max_range, mean([zMx, zMn]) + max_range])
    ax.set_aspect("equal")
end


function axisOff2(ax=gca())
    ax.set_xticklabels([""])
    xticks([])
    ax.set_yticklabels([""])
    yticks([])
    [s.set_visible(false) for s in values(ax.spines)]
end


function drawSphere(r::Real; center=[0, 0, 0], ax=axes(projection = "3d"), n=50, alter_axes=true, kwargs...)
    u = range(0, stop=2π, length=n)
    v = range(0, stop= π, length=n)

    x = r*cos.(u)*sin.(v)' .+ center[1]
    y = r*sin.(u)*sin.(v)' .+ center[2]
    z = r*ones(n)*cos.(v)' .+ center[3]

    plot_surface(x, y, z; kwargs...)

    if alter_axes
        maybeAxisLims("x", center[1]-r, center[1]+r)
        maybeAxisLims("y", center[2]-r, center[2]+r)
        maybeAxisLims("z", center[3]-r, center[3]+r)
        # in(color, ["none", "w", "white"]) && (color = edgecolor)
        # plot3D([center[1]], [center[2]], [center[3]])
        axisEqual3D(ax)
    end
end


function drawCircle(r::Real; center=[0, 0], ax=gca(), kwargs...)
    patch = pyimport("matplotlib.patches")
    c = patch.Circle(center, r; kwargs...)
    ax.add_patch(c)

    maybeAxisLims("x", center[1]-r, center[1]+r)
    maybeAxisLims("y", center[2]-r, center[2]+r)
    ax.set_aspect("equal")
end


# center: The center of the ellipse.
# w: The length of the horizontal axis.
# h: The length of the vertical axis.
# angle: Rotation of the ellipse in degrees (anti-clockwise).
# θ1, θ2 : Starting and ending angles of the arc in degrees. These values are relative to angle, .e.g. if angle = 45 and theta1 = 90 the absolute starting angle is 135. Default theta1 = 0, theta2 = 360, i.e. a complete ellipse.
function drawArc(w::Real, h::Real; ax=gca(), center=[0, 0], angle=0, θ1=0, θ2=360, edgecolor="none", linewidth=1, linestyle="-", label="", alpha=1, zorder=0)
    linestyle = sanitizeLinestyle(linestyle)
    patch = pyimport("matplotlib.patches")
    a = patch.Arc(center, w, h, angle=angle, theta1=θ1, theta2=θ2, edgecolor=edgecolor, linewidth=linewidth, label=label, linestyle=linestyle, alpha=alpha, zorder=zorder)
    ax.add_patch(a)

    wh = h*abs(cos(angle)) + w*abs(sin(angle))
    ww = h*abs(sin(angle)) + w*abs(cos(angle))
    l  = 0.5hypot(w, h)
    ϕ  = atan(h, w)
    xc = center[1] + l*cos(angle+ϕ)
    yc = center[2] + l*sin(angle+ϕ)

   maybeAxisLims("x", xc - 0.5ww, xc + 0.5ww)
   maybeAxisLims("y", yc - 0.5wh, yc + 0.5wh)
#     ax[:set_aspect]("equal")
end

# Take 2D array of points and projects them onto plane. P is n x 2 array of 2d points.
function planeTo3D(P, w; center=[0 0 0]) # w = nullspace(n̂) # where n̂ is normal vector to plane
    X = center[1] .+ w[1,1].*P[:,1] .+ w[1,2].*P[:,2] # Compute the corresponding cartesian coordinates
    Y = center[2] .+ w[2,1].*P[:,1] .+ w[2,2].*P[:,2] #   using the two vectors in w
    Z = center[3] .+ w[3,1].*P[:,1] .+ w[3,2].*P[:,2]
    return [X Y Z]
end


function drawCircle3D(r; n̂=[0 0 1], center=[0 0 0], ax=mplot3d.Axes3D(figure()), n=500, alter_axes=true, kwargs...)
    Θ = range(0, stop=2π, length=n)
    points2D = r.*[cos.(Θ) sin.(Θ)]
    w = nullspace(n̂)
    points3D = planeTo3D(points2D, w, center=center)
    draw3DPolygon(points3D, ax=ax; kwargs...)

    if alter_axes
        maybeAxisLims("x", center[1]-r, center[1]+r)
        maybeAxisLims("y", center[2]-r, center[2]+r)
        maybeAxisLims("z", center[3]-r, center[3]+r)
        ax.set_aspect("equal")
    end
end

# function drawCircle3D(r; n̂=[0 0 1], center=[0 0 0], ax=mplot3d.Axes3D(figure()), n=500, alter_axes=true, kwargs...)
#     Θ = range(0, stop=2π, length=n)
#     x = r.*cos.(Θ)
#     y = r.*sin.(Θ)
#     w = nullspace(n̂)

#     X, Y, Z, zorder = fill(NaN, (n, n)), fill(NaN, (n, n)), fill(NaN, (n, n)), fill(1, (n, n))
#     for i = 1:n, j = 1:n
#         points3D = planeTo3D([x[i] y[j] 0], w, center=center)
#         (X[i,j], Y[i,j], Z[i,j]) = points3D + center
#     end
#     plot_surface(X, Y, Z; kwargs...)

#     if alter_axes
#         maybeAxisLims("x", center[1]-r, center[1]+r)
#         maybeAxisLims("y", center[2]-r, center[2]+r)
#         maybeAxisLims("z", center[3]-r, center[3]+r)
#         ax.set_aspect("equal")
#     end
# end


function drawSquare(r::Real; ax=gca(), position=[0, 0], angle=0, color="blue", edgecolor="none", hatch="", linewidth=1, linestyle="-", label="", alpha=1, zorder=0)
    patch = pyimport("matplotlib.patches")
    linestyle = sanitizeLinestyle(linestyle)
    s = patch.Rectangle([position[1], position[2]], r, r, angle=rad2deg(angle), facecolor=color, edgecolor=edgecolor, linewidth=linewidth, label=label, linestyle=linestyle, alpha=alpha, zorder=zorder)
    ax.add_patch(s)

    wh = r*(abs(cos(angle)) + abs(sin(angle)))
    ww = r*(abs(sin(angle)) + abs(cos(angle)))
    xc = position[1] + r*cos(angle+0.25π)/√2
    yc = position[2] + r*sin(angle+0.25π)/√2

    maybeAxisLims("x", xc - 0.5ww, xc + 0.5ww)
    maybeAxisLims("y", yc - 0.5wh, yc + 0.5wh)
    # ax[:set_aspect]("equal")
end


function drawRectangle(w::Real, h::Real; ax=gca(), position=[0, 0], angle=0, color="blue", edgecolor="none", hatch="", linewidth=1, linestyle="-", label="", alpha=1, zorder=0)
    linestyle = sanitizeLinestyle(linestyle)
    patch = pyimport("matplotlib.patches")
    r = patch.Rectangle([position[1], position[2]], w, h, angle=rad2deg(angle), hatch=hatch, facecolor=color, edgecolor=edgecolor, linewidth=linewidth, label=label, linestyle=linestyle, alpha=alpha, zorder=zorder)
    ax.add_patch(r)

    wh = h*abs(cos(angle)) + w*abs(sin(angle))
    ww = h*abs(sin(angle)) + w*abs(cos(angle))
    l  = 0.5hypot(w, h)
    ϕ  = atan(h, w)
    xc = position[1] + l*cos(angle+ϕ)
    yc = position[2] + l*sin(angle+ϕ)

    maybeAxisLims("x", xc - 0.5ww, xc + 0.5ww)
    maybeAxisLims("y", yc - 0.5wh, yc + 0.5wh)
    # ax[:set_aspect]("equal")
end

function roundedRectangle(w::Real, h::Real; ax=gca(), position=[0,0], boxkind="round", pad=0.3, angle=0, color="blue", edgecolor="none", hatch="", linewidth=1, linestyle="-", label="", alpha=1, zorder=0)
    patch = pyimport("matplotlib.patches")
    h -= 2pad # fancyBox pad usually surrounds rectangle. Change w, h, position so that fancyBox is rectangle.
    w -= 2pad
    position += pad
    
    boxstyle = boxkind * ",pad=$(pad)"
    rr = patch.FancyBboxPatch([position[1], position[2]], w, h, boxstyle=boxstyle, hatch=hatch, facecolor=color, edgecolor=edgecolor, linewidth=linewidth, label=label, linestyle=linestyle, alpha=alpha, zorder=zorder)
    ax.add_patch(rr)
    
    wh = h*abs(cos(angle)) + w*abs(sin(angle))
    ww = h*abs(sin(angle)) + w*abs(cos(angle))
    l  = 0.5hypot(w, h)
    ϕ  = atan(h, w)
    xc = position[1] + l*cos(angle+ϕ)
    yc = position[2] + l*sin(angle+ϕ)

    maybeAxisLims("x", xc - 0.5ww, xc + 0.5ww)
    maybeAxisLims("y", yc - 0.5wh, yc + 0.5wh)
    # ax[:set_aspect]("equal")
end


function drawPolygon(n::Int64; r=1, angle=0, center=[0,0], ax=gca(), hatch="", color="blue", edgecolor="none", label="", linewidth=1, linestyle="-", alpha=1, zorder=0)
    patch = pyimport("matplotlib.patches")    
    linestyle = sanitizeLinestyle(linestyle)
    points = zeros(n, 2)
    for i = 1:n
        ϕ = 2π*i/n + angle
        points[i,:] = [r*cos(ϕ)+center[1] r*sin(ϕ)+center[2]]
    end    
    pol = patch.Polygon(points, closed=true, hatch=hatch, facecolor=color, edgecolor=edgecolor, linewidth=linewidth, linestyle=linestyle, label=label, alpha=alpha, zorder=zorder)
    ax.add_patch(pol)

    maybeAxisLims("x", minimum(points[:,1]), maximum(points[:,1]))
    maybeAxisLims("y", minimum(points[:,2]), maximum(points[:,2]))
    ax.set_aspect("equal")
end


# P is array of tuples, e.g. [0 0 0; 1 0 1; 1 1 0]
function draw3DPolygon(points; ax=mplot3d.Axes3D(figure()), kwargs...)
    mplot3d = pyimport("mpl_toolkits.mplot3d")
    poly = mplot3d.art3d.Poly3DCollection([points]; kwargs...)
    ax.add_collection3d(poly)
end

function drawCuboid(; C=zeros(3), L=ones(3), ax=axes(projection = "3d"), color="blue", hatch="", linewidth=0, edgecolor="none", alpha=1, zorder=0)
    xi, xf = C[1] - 0.5L[1], C[1] + 0.5L[1]
    yi, yf = C[2] - 0.5L[2], C[2] + 0.5L[2]
    zi, zf = C[3] - 0.5L[3], C[3] + 0.5L[3]

    P = [(xi, yi, zi),
         (xf, yi, zi),
         (xf, yf, zi),
         (xi, yf, zi),
         (xi, yi, zf),
         (xf, yi, zf),
         (xf, yf, zf),
         (xi, yf, zf)]

    face = [1 2 3 4; 5 6 7 8; 3 4 8 7; 1 2 6 5; 2 3 7 6; 1 4 8 5]
    for i = 1:6
        Q = P[face[i,:]]
        draw3DPolygon(Q, ax=ax, color=color, hatch=hatch, linewidth=linewidth, edgecolor=edgecolor, alpha=alpha, zorder=zorder)
    end
end


function drawAxes3D(; L=1, names=("x", "y", "z"), ax=subplot(projection="3d"))
    quiver(0, 0, 0, L, 0, 0)
    quiver(0, 0, 0, 0, L, 0)
    quiver(0, 0, 0, 0, 0, L)

    ax.text(1.05L, 0, 0, names[1], "x", ha="center", va="center")
    ax.text(0, 1.05L, 0, names[2], "y", ha="center", va="center")
    ax.text(0, 0, 1.05L, names[3], "z", ha="center", va="center")
end


function sanitizeLinestyle(ls)
    ls == "-"  && return "solid"
    ls == "--" && return "dashed"
    ls == "-." && return "dashdot"
    ls == ":"  && return "dotted"
                  return ls
end


function drawSegment(r; center=[0 0], direction=[0 1], halfAngle=0.25π, ax=gca(), n=100,  hatch="", color="blue", edgecolor="none", label="", linewidth=1, linestyle="-", alpha=1, zorder=0)
    patch = pyimport("matplotlib.patches")
    dir = atan(direction[2]-center[2], direction[1]-center[1])
    t = range(dir-halfAngle, stop=dir+halfAngle, length=n)
    points = ([r*cos.(t)+center[1] r*sin.(t)+center[2]; [center[1] center[2]]])
    points = [points; [center[1] center[2]]]

    linestyle = sanitizeLinestyle(linestyle)
    # @pyimport matplotlib.patches as patch
    seg = patch.Polygon(points, closed=true, hatch=hatch, facecolor=color, edgecolor=edgecolor, linewidth=linewidth, linestyle=linestyle, label=label, alpha=alpha, zorder=zorder)
    ax.add_patch(seg)

    maybeAxisLims("x", minimum(points[:,1]), maximum(points[:,1]))
    maybeAxisLims("y", minimum(points[:,2]), maximum(points[:,2]))
    ax.set_aspect("equal")
end


function drawGSEqProj(S::Real; ϕ=0, λ=0, halfAngle=π/4, r=R⨁, ax=gca(), n=1000, hatch="", color="blue", edgecolor="none", label="", linewidth=1, linestyle="-", alpha=1, zorder=0)
    patch = pyimport("matplotlib.patches")

    linestyle = sanitizeLinestyle(linestyle)
    ϕ = abs(ϕ)
    (halfAngle ≤ ϕ   ) && @error "Error: No equatorial coverage at latitudes greater than the view half-angle."
    (halfAngle ≥ 0.5π) && @error "Error: No equatorial coverage at latitudes greater than the view half-angle."

    a = (cos(ϕ)*tan(halfAngle))^2 - sin(ϕ)^2
    b = r*(tan(halfAngle)^2)*cos(ϕ)

    x₀ = b/a # Center of hyperbola

    ψ = atan(√a)
    arc = ψ - asin((x₀-r*cos(ϕ))*sin(ψ)/S)

    rMin = r*tan(halfAngle) / (tan(halfAngle)*cos(ϕ) - sin(ϕ)) # Minimum view distance.
    (rMin > S) && @error "Error: Minimum equator distance is less than maximum view distance S."

    # Get coordinates of straight line
    x = range(rMin*(1+eps()), stop=r*cos(ϕ) + S*cos(arc), length=n)
    y = .√((tan(halfAngle)*((2r*cos(ϕ)-x)*cos(ϕ) - r*cos(2*ϕ))).^2 - (sin(ϕ)*x).^2)

    # Get coordinates of curved part
    yMax = √((tan(halfAngle)*((r*cos(ϕ) - S*cos(arc))*cos(ϕ) - r*cos(2*ϕ))).^2 - (sin(ϕ)*(r*cos(ϕ) + S*cos(arc))).^2)
    ySeg = range(0, stop=yMax, length=n)
    xSeg = r*cos(ϕ) + .√(S^2 - ySeg.^2 - (r*sin(ϕ))^2)

    # Combine coordinates of straight and curved line segments.
    xTot, yTot = [x' reverse(xSeg)' xSeg' reverse(x)'], [y' reverse(ySeg)' -ySeg' reverse(-y)']

    # Rotate by 'λ' around origin.
    xTot2 = xTot*cos(λ) - yTot*sin(λ)
    yTot2 = xTot*sin(λ) + yTot*cos(λ)

    coords = [xTot2' yTot2']

    # Plot segment
    seg = patch.Polygon(coords, closed=true, hatch=hatch, facecolor=color, edgecolor=edgecolor, linewidth=linewidth, linestyle=linestyle, label=label, alpha=alpha, zorder=zorder)
    ax.add_patch(seg)

    maybeAxisLims("x", minimum(coords[:,1]), maximum(coords[:,1]))
    maybeAxisLims("y", minimum(coords[:,2]), maximum(coords[:,2]))
    ax.set_aspect("equal")
end


function verticalLine(x; ax=gca(), ymin=ylim()[1], ymax=ylim()[2], color="r", linewidth=1, linestyle="-", label="", zorder=0)
    # plot([x, x], [ymin, ymax], color=color, linewidth=linewidth, linestyle=linestyle, label=label, zorder=zorder)
    # maybeAxisLims("y", ymin, ymax)
    ymin2 = (ymin-ylim()[1])/(ylim()[2]-ylim()[1])
    ymax2 = (ymax-ylim()[1])/(ylim()[2]-ylim()[1])
    axvline(x=x, ymin=ymin2, ymax=ymax2, color=color, lw=linewidth, linestyle=linestyle, label=label, zorder=zorder)
end


function horizontalLine(y; ax=gca(), xmin=xlim()[1], xmax=xlim()[2], color="r", linewidth=1, linestyle="-", label="", zorder=0)
    # plot([xmin, xmax], [y, y], color=color, linewidth=linewidth, linestyle=linestyle, label=label, zorder=zorder)
    # maybeAxisLims("x", xmin, xmax)
    xmin2 = (xmin-xlim()[1])/(xlim()[2]-xlim()[1])
    xmax2 = (xmax-xlim()[1])/(xlim()[2]-xlim()[1])
    axhline(y=y, xmin=xmin2, xmax=xmax2, color=color, lw=linewidth, linestyle=linestyle, label=label, zorder=zorder)
end


function verticalBand(xmin, xmax; ax=gca(), ymin=ylim()[1], ymax=ylim()[2], color="g", alpha=1, hatch="", edgecolor="none", linewidth=1, linestyle="-", label="", zorder=0)
    # linestyle=sanitizeLinestyle(linestyle)
    axvspan(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, facecolor=color, hatch=hatch, alpha=alpha, linestyle=linestyle, linewidth=linewidth, edgecolor=edgecolor, label=label, zorder=zorder)
    maybeAxisLims("x", xmin, xmax)
    maybeAxisLims("y", ymin, ymax)
end


function horizontalBand(ymin, ymax; ax=gca(), xmin=xlim()[1], xmax=xlim()[2], color="g", alpha=1, hatch="", edgecolor="none", linewidth=1, linestyle="-", label="", zorder=0)
    # linestyle=sanitizeLinestyle(linestyle)
    axhspan(ymin=ymin, ymax=ymax, xmin=xmin, xmax=xmax, facecolor=color, hatch=hatch, alpha=alpha, linestyle=linestyle, linewidth=linewidth, edgecolor=edgecolor, label=label, zorder=zorder)
    maybeAxisLims("x", xmin, xmax)
    maybeAxisLims("y", ymin, ymax)
end


function maybeAxisLims(axis, axmin, axmax; ax=gca())
    axis == "x" && return ax.set_xlim(min(ax.get_xlim()[1], axmin), max(ax.get_xlim()[2], axmax))
    axis == "y" && return ax.set_ylim(min(ax.get_ylim()[1], axmin), max(ax.get_ylim()[2], axmax))
    axis == "z" && return ax.set_zlim(min(ax.get_zlim()[1], axmin), max(ax.get_zlim()[2], axmax))
    @error """`axis` must be "x", "y", or "z"."""
end


function drawCone(cotθ::Real, h₀::Real; n̂=[0.0 0.0 1.0], ax=mplot3d.Axes3D(figure()), origin=[0.0; 0.0; 0.0], n=100, rstride=1, cstride=1, antialiased=true, kwargs...)
    @assert norm(n̂) ≈ 1.0 "`n̂` must be a unit direction vector."
    r = h₀/cotθ
    q = range(-r, stop=r, length=n)
    
    if n̂[1] ≈ n̂[2] ≈ 0 # +Z or -Z
        k = (n̂[3] ≈ 1.0) ? 1 : -1
        R = k*[1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    else
        û = normalize([-n̂[2] n̂[1] 0])
        R = rotMat(acos(n̂[3]), û)
    end

    X, Y, Z = fill(NaN, (n, n)), fill(NaN, (n, n)), fill(NaN, (n, n))
    for i = 1:n, j = 1:n
        v⃗ = [q[i]; q[j]; hypot(q[i], q[j])*cotθ]
        v⃗ *= h₀/max(v⃗[3],h₀)
        (X[i,j], Y[i,j], Z[i,j]) = R*v⃗ + origin
    end
    plot_surface(X, Y, Z, rstride=rstride, cstride=cstride, antialiased=antialiased; kwargs...)
end


meshgrid(x, y) = repeat(x', length(y), 1), repeat(y, 1, length(x))

# Grom GMAT mathspec
# d is positive for clockwise ecef motion. d = sign(x*ẏ - ẋ*y)
function updateGroundTrack(λᵢ, λᵢ₋₁, ϕᵢ, ϕᵢ₋₁, dᵢ, dᵢ₋₁)
    mλ⁺ᵢ   = mod2pi(λᵢ    )
    mλ⁺ᵢ₋₁ = mod2pi(λᵢ₋₁  )
    mλ⁻ᵢ   = mod(λᵢ,   -2π)
    mλ⁻ᵢ₋₁ = mod(λᵢ₋₁, -2π)

    if     (dᵢ == dᵢ₋₁ ==  1) && (mλ⁺ᵢ₋₁ <  π && mλ⁺ᵢ   >  π) # Wraps off RHS
        m  = (ϕᵢ - ϕᵢ₋₁)/(mλ⁺ᵢ - mλ⁺ᵢ₋₁)
        ϕb = m*(π - mλ⁺ᵢ) + ϕᵢ
        plot([λᵢ₋₁, π], [ϕᵢ₋₁, ϕb], "b")
        plot([-π,  λᵢ], [ϕb,   ϕᵢ], "b")
    elseif (dᵢ == dᵢ₋₁ == -1) && (mλ⁻ᵢ   < -π && mλ⁻ᵢ₋₁ > -π) # Wraps off LHS
        m  = (ϕᵢ - ϕᵢ₋₁)/(mλ⁻ᵢ - mλ⁻ᵢ₋₁)
        ϕb = m*(-π - mλ⁻ᵢ) + ϕᵢ
        plot([λᵢ₋₁, -π], [ϕᵢ₋₁, ϕb], "b")
        plot([π,    λᵢ], [ϕb,   ϕᵢ], "b")
    else # No wrap
        plot([λᵢ₋₁,  λᵢ], [ϕᵢ₋₁, ϕᵢ], "b")
    end 
end

# Takes in ECEF footprint coords from footprint function.
# TODO: Make this function more elegant.
# TODO: Check "inverse split patch"
# DONE: Check "inverse patch"
# DONE: Stop plotting side lines during wraps.
function plot_footprint2D(coords_ECEF::Array{Float64}, T_ECEF::Array{Float64}; ax=gca(), color="none", hatch="", linewidth=1, linestyle="-", label="", alpha_fill=0, alpha_line=1, zorder=0, edgecolor="b", markTarget=false, marker="+", markerSize=40, markerColor="k") 
    patch = pyimport("matplotlib.patches")

    T = ecef2gdc(T_ECEF)[1:2]
    n = size(coords_ECEF)[1]
    coords = zeros(n, 2)
    for i = 1:n
        coords[i,:] = ecef2gdc(coords_ECEF[i,:])[1:2]
    end
        
    arr  = find(abs.(diff([coords[:,1]; coords[1,1]])) .> π)
    # if 0 breaks, check for inverse fix.
    # if 1 break, need to do edges and polar fixes.
    # if 2 breaks need to do edges and check for inverse fix., 
    lines   = isempty(arr) ? [coords] : Array{Float64, 2}[]
    patches = isempty(arr) ? [coords] : Array{Float64, 2}[]
        
    for i = 1:length(arr)
        len  = i ≠ length(arr) ? arr[i+1] : size(coords)[1] + arr[1]
        len -= arr[i]
        
        range = ((arr[i]-1 : arr[i] + len) .% n) + 1
        
        # Add edges to each end of string
        str = zeros(len+2, 2)
        str[2:end-1,:] = coords[range[2:end-1],:]
        str[1      ,:] = groundtrack_interpolate(coords[range[1  ],:], coords[range[2    ],:])
        str[  end  ,:] = groundtrack_interpolate(coords[range[end],:], coords[range[end-1],:])
          
        push!(lines, str)

        if length(arr) == 1 # Polar case, add edges
            indMin     = findmin([norm(coords[j,:] - T) for j = 1:n])[2]
            p          = (coords[indMin,2] > T[2]) ? -0.5 : 0.5
            
            str = [π*[sign(coords[range[2    ],1]) p]; str]
            str = [str; π*[sign(coords[range[end-1],1]) p]]
        end
        
        push!(patches, str)    
    end

    if length(arr) == 0 && !PIPoly(patches[1], T) # check for inverse fix 
        patches = [invertFootprint(patches[1])]
    end

    for i = 1:length(lines)
        plot(rad2deg.(lines[i][:,1]), rad2deg.(lines[i][:,2]), color=edgecolor, linewidth=linewidth, linestyle=linestyle, label=label, alpha=alpha_line, zorder=zorder) 
        shape = patch.Polygon(rad2deg.(patches[i]), closed=true, hatch=hatch, facecolor=color, edgecolor="none", linewidth=0, alpha=alpha_fill, zorder=zorder)
        ax.add_patch(shape)
    end
    markTarget && scatter(T[1], T[1], marker=marker, s=markerSize, color=markerColor, zorder=zorder+1)
end


function groundtrack_interpolate(p₁::Array{T}, p₂::Array{T}) where T <: Float64
    q  = sign(p₁[1])
    x₃ = p₂[1] + q*2π
    m  = (p₁[1] == x₃) ? 0 : abs(p₁[1]-q*π)/abs(p₁[1]-x₃)
    return [-q*π m*(p₂[2]-p₁[2]) + p₁[2]]
end

# Returns inverse footprint given a loop of GDC coordinates.
function invertFootprint(coords::Array{Float64,2})
    maxInd = findmax(coords[:,2])[2] # Index of the top of the loop
    coords2 = zeros(size(coords, 1)+7, 2)
    coords2[1:maxInd,:] = coords[1:maxInd,:]
    coords2[maxInd+1,:] = [coords[maxInd,1] 0.5π]
    coords2[maxInd+2,:] = [-π  0.5π]
    coords2[maxInd+3,:] = [-π -0.5π]
    coords2[maxInd+4,:] = [ π -0.5π]
    coords2[maxInd+5,:] = [ π  0.5π]
    coords2[maxInd+6,:] = [coords[maxInd,1] 0.5π]
    coords2[maxInd+7,:] = coords[maxInd,:]
    if maxInd ≠ size(coords, 1)
        coords2[maxInd+8:end,:] = coords[maxInd+1:end,:]
    end
    return coords2
end


"""
    flattenTransparentImage{T <: String}(infile::T, outfile::T; RGBA_bg=ones(4))

Export a partially transparent image with binary alpha channel (0 or 1) e.g. for gifs.
"""
function flattenTransparentImage(infile::T, outfile::T; RGBA_bg=ones(4)) where T <: Float64
    im1 = imread(infile)
    im2 = zeros(Float32, size(im1))
    
    for i = 1:size(im1)[1], j = 1:size(im1)[2]
        im2[i,j,:] = mimicAlpha(im1[i,j,:], dst=RGBA_bg)
    end
        
    imsave(outfile, im2)
end


"""
    mimicAlpha(src; dst=ones(4))

Returns a blend of a foreground color RGBA (`src`) with background RGBA (`dst`, default white).
Source here: https://en.wikipedia.org/wiki/Alpha_compositing#Alpha_blending
"""
function mimicAlpha(src; dst=ones(4))
    out      = zeros(4)
    out[4]   = src[4] + dst[4]*(1-src[4])
    out[1:3] = (out[4] == 0) ? zeros(3) : (src[1:3]*src[4] + dst[1:3]*dst[4]*(1-src[4]))/out[4]
    return out
end





#=
http://www.tapir.caltech.edu/~dtsang/python.html
function plotCoords(t, systemCoords, covData, bodies, points, vehicles, stations, covFrom, covTo, view, plot, plotType, animate, skipE)
    #clf();
    #fig = figure()

    if plotType == "position"
        inds = getPlotInds(view);

        if skipE != 0
            #drawEarth(fig);
        end

        if !animate
            if length(inds) == 3
                for k = 1:size(systemCoords)[3] # For each item
                    PyPlot.plot3D(systemCoords[:,1,k], systemCoords[:,2,k], systemCoords[:,3,k]);
                    hold(true);
                end
                axis("equal")
                hold(false)
            else
                for k = 1:size(systemCoords)[3] # For each item
                    PyPlot.plot(systemCoords[:,inds[1],k], systemCoords[:,inds[2],k]);
                    hold(true);
                end
                axis("equal")
                hold(false)
            end
            # save as png in plots/static/ folder
        else
            error("Have not implemented animations yet.");
            #
            @pyimport matplotlib.animation as anim
            # First set up the figure, the axis, and the plot element we want to animate
            fig = figure()
            #ax = plt.axes(xlim=(0, 2), ylim=(-2, 2))
            ax = plt.axes()
            global line = ax[:plot]([], [], lw = 2)[1]

            # call the animator.  blit=True means only re-draw the parts that have changed.
            myanim = anim.FuncAnimation(fig, animate, init, interval=20, blit=true);

            myanim[:save]("/plots/animations/" * filename * ".mp4", extra_args=["-vcodec", "libx264", "-pix_fmt", "yuv420p"]);
            #
            #
            # Using PyPlot in external GUI window
            pygui(false);
            PyPlot.ion();
            PyPlot.show();
            PyPlot.hold(false);
            for i = 1:size(systemCoords)[1] # For each timestep
                for j = 1:size(systemCoords)[2] # For each body
                    PyPlot.scatter(systemCoords[i,j,1], systemCoords[i,j,2]);
                end
                axis("equal");
                PyPlot.draw();
                IJulia.clear_output(wait=true)
                sleep(0.08);
            end
            #

            #
            #
            function myplot(data)
                Winston.plot(data[1,:], data[2,:], "b.")
            end

            # smallish system size
            N=30
            # initial configuration
            x = ljconfig(N, 0.1)
            p = zeros(size(x))
            # create an "input" object
            myinput = Input(x)
            reactive.lift(myplot, myinput)

            # lift(pyplot, myinput)
            #
            #using Plotly
            trace1 = [  "x"         => [],
                        "y"         => [],
                        "type"      => "scatter",
                        "mode"      => "markers",
                        "stream"    => ["token" => "xc05zyowlc"]
                        ];

            plot   = Plotly.plot([trace1], ["filename" => "plot", "world_readable" => false, "fileopt" => "overwrite"]);

            s1 = Plotly.Stream(["token" => "xc05zyowlc"]);
            s1.open();

            for i = 1:size(systemCoords)[1] # For each timestep
                s1.write([  "x" => [systemCoords[i,:,1]],
                            "y" => [systemCoords[i,:,2]]]);
                sleep(0.08);
            end

            s1.close();

            s = string("""<iframe height='750' id='igraph' scrolling='no' seamless='seamless' src='", plot["url"], "/700/700' width='750'></iframe>""");
            display("text/html", s);
            #
             Using Winstone & Reactive

            p = Input(0.0:0.0)
            function f(x)
                plot(sin(5*x), sin(2π*x))
            end
            lift(f, p)

            for i = 0:0.05:10
                sleep(1/60)
                push!(p, 0.0:0.01:i)
            end
            #
        #elseif animate == "combo"
        end
    elseif plotType == "coverage"
        println(size(covData)[2])
        for k = 1:size(covData)[2] # For each combination
            PyPlot.plot(Date(julian2datetime(t)), covData[:,k]);
            hold(true);
        end

    else
        error("""'plot' must be either "position" or "coverage".""");
    end
end
=#

 # initialization function: plot the background of each frame
#=
function init()
    global line
    line[:set_data]([], [])
    return (line,None)
end

# animation function.  This is called sequentially
function animate(i)
    x = systemCoords[i,:,1];
    y = systemCoords[i,:,2];
    global line
    line[:set_data](x, y)
    return (line,None)
end
=#