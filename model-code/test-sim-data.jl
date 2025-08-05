include("gravity-tools.jl")
include("trade-environment.jl")
include("simmulate-eaton-kortum.jl")
using CSV
using DataFrames
using Plots
using MINPACK
using HypothesisTests
using Optimization
using OptimizationOptimJL
using OptimizationPRIMA
using LinearAlgebra

df = DataFrame(CSV.File("./data/pricegap-df-2004.csv"))

grav_file = "./data/top30_gravity_data.csv"

dfgrav = DataFrame(CSV.File(grav_file))

rename!(df, Dict("exporter" => "iso_o", "importer" => "iso_d"))

df = innerjoin(df, dfgrav, on = ["iso_o", "iso_d"])

filter!(row -> ~(row.Xni ≈ 1.0), df);

# filter!(row -> ~(row.Xni ≈ 0.0), df);

# println( mean(df.logXni) / mean(df.dni) )

test = CorrelationTest(log.(df.dist), df.dni2)


################################################################

dftrade = DataFrame(CSV.File("./data/tradeshare-df-2004.csv"))

dftrade[!,"trade"] = log.(dftrade[!,"norm_tradeshare"] )

# removing the home trade flows
filter!(row -> ~(row.norm_tradeshare ≈ 1.0), dftrade);
dfcountryfix = deepcopy(dftrade)

# remove the zero trade flows
filter!(row -> ~(row.norm_tradeshare ≈ 0.0), dftrade);

Ncntry = 30

L = ones(Ncntry)

θ = 4.0
σ = 2.5

grv_params = gravity_params(Ncntry = Ncntry, θ = θ, L = L, dfcntryfix = dfcountryfix )

# dfsim, trd_prm, grvdata = generate_simmulated_data(4.0, 0.36, dftrade, grv_params; model = "ek", code = 100)

# estθ = boot_strap_simulation(4.0, 0.36, dftrade, grv_params; 
#     model = "ek", method = "exact", code = 100)


estimate_θ_dni(df, grv_params, trd_prm, 
grvdata; model = "bejk", method = "exact", Wmat = "optimal", display = true, Nruns = 8)

            
# p = SciMLBase.NullParameters()

# f(x, p) = estimate_θ_dni(x[1], dfsim.dni, dfsim.dni2, dfsim.dist,  grv_params, trd_prm,
#  grvdata; model = "ek", Wmat = "optimal", display = true, Nruns = 9)

# lb = [2.5,]
# ub = [10.0,]

# # Define the optimization problem
# prob = OptimizationProblem(f, [4.5], p, lb = lb, ub = ub)

# # Solve the problem using BOBYQA with options
# @time sol = Optimization.solve(prob, BOBYQA(); rhoend = 1e-4)

