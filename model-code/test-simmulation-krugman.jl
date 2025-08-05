include("gravity-tools.jl")
include("trade-environment.jl")
include("simmulate-eaton-kortum.jl")
using CSV
using DataFrames
using Plots
using MINPACK
using Optimization
using OptimizationOptimJL
using OptimizationPRIMA
using LinearAlgebra

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


trd_prm = trade_params(θ = grv_params.θ, d = d, S = exp.(grvdata.S), Ncntry = grv_params.Ncntry, N = grv_params.L)

# @time πshares, foo = sim_trade_pattern_krugman(trd_prm.S, d, grv_params.θ);

@time πshares, foo = sim_trade_pattern_krugman(trd_prm);

dfmodel = make_trademodel_dataframe(πshares, grv_params)

# dfmodel_BEJK = plot_trade(πshares_BEKK, Ncntry);

plot(log.(dfmodel.tradeshare), log.(dfmodel.tradedata), seriestype = :scatter, alpha = 0.75,
    xlabel = "model",
    ylabel = "data",
    legend = false)

