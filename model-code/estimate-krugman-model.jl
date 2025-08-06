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

year = ["2004", "2011", "2017"]
model = "krugman"
Nruns = 12
Nboots = 100
Ngoods = 10000

dirname = "./data/"

##############################################################################################################################
method = "over"

dfout = estimate_all(year, model, method, Nruns, Nboots, dirname; Ngoods = Ngoods)

CSV.write("./results/krugman-estimate-"*method*".csv", dfout; writeheader = true)

##############################################################################################################################

method = "exact"

dfout = estimate_all(year, model, method, Nruns, Nboots, dirname; Ngoods = Ngoods)

CSV.write("./results/krugman-estimate-"*method*".csv", dfout; writeheader = true)