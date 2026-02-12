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

year = ["2004", "2011", "2017"]

model = "krugman-model2"
Nruns = 36
Nboots = 100
Ngoods = 1000
σ = 1.5

dirname = "./data/"
date = "021126"

################################################################################################################################
method = "over"

dfout = estimate_all(year, σ, model, method, Nruns, Nboots, dirname; Ngoods = Ngoods)

CSV.write("./results/"*model*"-estimate-"*method*"-"*date*".csv", dfout; writeheader = true)

################################################################################################################################

method = "exact"

dfout = estimate_all(year, σ, model, method, Nruns, Nboots, dirname; Ngoods = Ngoods)

CSV.write("./results/"*model*"-estimate-"*method*"-"*date*".csv", dfout; writeheader = true)