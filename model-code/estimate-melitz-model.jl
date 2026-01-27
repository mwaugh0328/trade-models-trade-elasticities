include("gravity-tools.jl")
include("trade-environment.jl")
include("simmulate-trade-models.jl")
include("estimate-trade-models.jl")
using CSV
using DataFrames
using Plots
using MINPACK
using HypothesisTests
using Optimization
using OptimizationOptimJL
using OptimizationPRIMA
using LinearAlgebra
################################################################

# year = ["2004", "2011", "2017"]
year = ["2017"]
model = "melitz"
Nruns = 12
Nboots = 10
Ngoods = 50000

dirname = "./data/"

##############################################################################################################################
method = "over"
σ = 1.5

dfout = estimate_all(year, σ, model, method, Nruns, Nboots, dirname; Ngoods = Ngoods)

# CSV.write("./results/melitz-estimate-"*method*".csv", dfout; writeheader = true)

# ##############################################################################################################################

# method = "exact"

# dfout = estimate_all(year, model, method, Nruns, Nboots, dirname; Ngoods = Ngoods)

# CSV.write("./results/melitz-estimate-"*method*".csv", dfout; writeheader = true)