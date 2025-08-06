include("gravity-tools.jl")
include("trade-environment.jl")
include("simmulate-eaton-kortum.jl")
include("estimate-models.jl")
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

year = ["2011"]

model = "ek"
method = "over"
Nruns = 36
Nboots = 20

dirname = "./data/"

dfout = estimate_all(year, model, method, Nruns, Nboots, dirname)

# CSV.write("./results/ek-estimate-"*method*".csv", dfout; writeheader = true)