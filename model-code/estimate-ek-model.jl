include("gravity-tools.jl")
include("trade-environment.jl")
include("simulate-trade-models.jl")
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
using Dates
################################################################


year = ["2004", "2011", "2017"]

#year = ["2017"]
model = "ek"
Nruns = 36
Nboots = 100
Ngoods = 100000
σ = 1.5

repo_root = normpath(joinpath(@__DIR__, ".."))
dirname = joinpath(repo_root, "data")
results_dir = joinpath(repo_root, "results")
date = Dates.format(today(), "yyyy-mm-dd")

# ##############################################################################################################################
method = "over"

dfout = estimate_all(year, σ, model, method, Nruns, Nboots, dirname; Ngoods = Ngoods)

CSV.write(joinpath(results_dir, model*"-estimate-"*method*"-"*date*".csv"), dfout; writeheader = true)

# ##############################################################################################################################

method = "exact"

dfout = estimate_all(year, σ, model, method, Nruns, Nboots, dirname; Ngoods = Ngoods)

CSV.write(joinpath(results_dir, model*"-estimate-"*method*"-"*date*".csv"), dfout; writeheader = true)