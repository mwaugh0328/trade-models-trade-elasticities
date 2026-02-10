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

year = [ "2017"]
#year = ["2017"]
model = "melitz-model2"
Nruns = 36
Nboots = 2
Ngoods = 50000

σ = 2.00

dirname = "./data/"
date = "12926-sigma-theta"

# ##############################################################################################################################
# method = "over"


# dfout = estimate_all(year, σ, model, method, Nruns, Nboots, dirname; Ngoods = Ngoods)

# CSV.write("./results/"*model*"-estimate-"*method*"-"*date*".csv", dfout; writeheader = true)

##############################################################################################################################

method = "exact"
    
dfout = estimate_all(year, σ, model, method, Nruns, Nboots, dirname; Ngoods = Ngoods)

# CSV.write("./results/"*model*"-estimate-"*method*"-"*date*".csv", dfout; writeheader = true)