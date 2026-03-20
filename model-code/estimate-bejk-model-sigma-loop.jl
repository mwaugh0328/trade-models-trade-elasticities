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
################################################################

year = ["2004", "2011", "2017"]
# year = ["2017"]
model = "bejk"

Nruns = 36
Nboots = 100
Ngoods = 100000

dirname = "./data/"

# Loop over different sigma values
for σ in 1.5:0.5:3.0
    
    println("Running estimation for σ = ", σ)
    
    # Create date string based on sigma value
    sigma_str = replace(string(σ), "." => "")
    date = "21126-sigma"*sigma_str
    
    ###############################################################################################################################
    method = "over"
    
    dfout = estimate_all(year, σ, model, method, Nruns, Nboots, dirname; Ngoods = Ngoods)
    
    CSV.write("./results/"*model*"-estimate-"*method*"-"*date*".csv", dfout; writeheader = true)
    
    ###############################################################################################################################
    
    method = "exact"
    
    dfout = estimate_all(year, σ, model, method, Nruns, Nboots, dirname; Ngoods = Ngoods)
    
    CSV.write("./results/"*model*"-estimate-"*method*"-"*date*".csv", dfout; writeheader = true)
    
    println("Completed estimation for σ = ", σ)
    println("="^80)
    
end