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

τ = zeros(Ncntry,Ncntry)

trd_prm = trade_params(θ = grv_params.θ, d = d, S = exp.(grvdata.S), Ncntry = grv_params.Ncntry, N = grv_params.L)

# # @time πshares, foo = sim_trade_pattern_ek(exp.(grvdata.S), d, τ, grv_params.θ, 1.5);

# println("Simmulated trade pattern EK")

# πshares, foo = sim_trade_pattern_ek(trd_prm);

# beta_moment_model(foo, πshares)


# println("Simmulated trade pattern BEJK")

# πshares_bejk, foo = sim_trade_pattern_bejk(trd_prm);

# beta_moment_model(foo, πshares_bejk)


# test = estimate_θ(4.0, grv_params, trd_prm, grvdata)

# @time out2 = generate_moments(trd_prm, 100; code = 300, Nprices = 1000)


# p = SciMLBase.NullParameters()

# f(x, p) = norm(-5.60 - mean(estimate_θ(x[1], grv_params, trd_prm, grvdata)))

# lb = [2.0,]
# ub = [10.0,]

# # Define the optimization problem
# prob = OptimizationProblem(f, [4.5], p, lb = lb, ub = ub)

# # Solve the problem using BOBYQA with options
# sol = Optimization.solve(prob, BOBYQA())




# @time πshares_BEKK, foo = sim_trade_pattern_BEJK(exp.(grvdata.S), d, grv_params.θ, 1.5);

# πshares = average_trade_pattern(exp.(grvdata.S), d, grv_params.θ, 1.5, Nruns = 30);

# ################################################################
# # Check to see it they line up with the data


# dfmodel = make_trademodel_dataframe(πshares, grv_params)

# # dfmodel_BEJK = plot_trade(πshares_BEKK, Ncntry);

# plot(log.(dfmodel.tradeshare), log.(dfmodel.tradedata), seriestype = :scatter, alpha = 0.75,
#     xlabel = "model",
#     ylabel = "data",
#     legend = false)