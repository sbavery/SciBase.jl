module AstroSciKit

# global dName = dirname(Base.source_path()) # used in vsop87d.jl, spacecraft.jl, rf.jl to load data
global azElGS # assigned in main.jl, used in frames.jl
global azElGSEcef # assigned in main.jl, used in frames.jl
global azElGSGDC # assigned in main.jl, used in frames.jl

using Polynomials, PyPlot, PyCall, HTTP, SpecialFunctions, Dates, JLD, HDF5, DelimitedFiles, LinearAlgebra, Printf
import Statistics.mean

files = [
		"astrodynamics",
		"attenuation",
		"utils", # Put first as other scripts need these functions
		"bodies",
		"calc",
		"constants",
		"frames",
		"geodesy",
		"geometry",
		"stations",
		"main",
		"maneuvers",
		"plot",
		"rf",
		"spacecraft",
		"stk",
		"time",
		"vsop87d"
		]

for f in files
    include(f * ".jl")
end

# !isfile(joinpath(datapath, "scLib.jld")) && updateSCLib()#dName * "/data/scLib.jld") # <- this for .jld library

# if isdir(dName * "/data/scLib.jld")
# 	rm(dName * "/data/scLib.jld")
# end

# if !isdir(dName * "/data/scLib")
# 	mkdir(dName * "/data/scLib")
# 	updateSCLib()
# end

# if isfile(dName * "/data/scLib.jld")
# 	rm(dName * "/data/scLib.jld")
# 	updateSCLib()
# end

!isdir("Plots") && mkdir("Plots")
!isdir("Plots/Orbit_Plots") && mkdir("Plots/Orbit_Plots")

end # module