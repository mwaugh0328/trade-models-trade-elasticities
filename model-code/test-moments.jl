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
################################################################
# builds the EK dataset

# dftrade, dfcntryfix, dflabor = make_ek_dataset()
# # this one has the country numbers which allows for the construction of the 
# # trade costs given the estimated fixed effects from the gravity regressiondf

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

# ################################################################
# # Run the Gravity regression

grvdata = gravity(dftrade, display = true);

# # ################################################################
# # # Recover the trade costs and technology parameters

d = zeros(Ncntry,Ncntry)
T = zeros(Ncntry)
W = ones(Ncntry)

make_trade_costs!(grvdata, d, grv_params)

make_technology!(grvdata, T, W, grv_params)

# # ################################################################
# # # Now simmulate the EK model

τ = zeros(Ncntry,Ncntry)

trd_prm = trade_params(θ = grv_params.θ, σ = σ, d = d, S = exp.(grvdata.S), Ncntry = grv_params.Ncntry, N = grv_params.L)

# # outreg = reg(pricegap_df, @formula(log(dist) ~ fe(importer) + fe(exporter) + dni), save = true, tol = 1e-10)

# outreg = reg(df, @formula(dni ~ fe(importer) + fe(exporter) + border + log(dist)), save = true, tol = 1e-10)


p = SciMLBase.NullParameters()

f(x, p) = estimate_θ_dni(x[1], df.dni, grv_params, trd_prm, grvdata; model = "ek", display = true, Nruns = 10)


lb = [2.5,]
ub = [10.0,]

# Define the optimization problem
prob = OptimizationProblem(f, [4.5], p, lb = lb, ub = ub)

# Solve the problem using BOBYQA with options
@time sol = Optimization.solve(prob, BOBYQA())