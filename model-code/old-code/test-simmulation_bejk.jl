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

directory = "./data/"

yyy = "2011"

df = DataFrame(CSV.File(directory*"pricegap-df-"*yyy*".csv"))

grav_file = directory*"top30_gravity_data.csv"

dfgrav = DataFrame(CSV.File(grav_file))

rename!(df, Dict("exporter" => "iso_o", "importer" => "iso_d"))

df = innerjoin(df, dfgrav, on = ["iso_o", "iso_d"])

dftariffs = DataFrame(CSV.File(directory*"tariffs-"*yyy*".csv"))

rename!(dftariffs, Dict("exporter" => "iso_o", "importer" => "iso_d"))

df = innerjoin(df, dftariffs, on = ["iso_o", "iso_d"])

filter!(row -> ~(row.Xni ≈ 1.0), df);

# filter!(row -> ~(row.Xni ≈ 0.0), df);

# println( mean(df.logXni) / mean(df.dni) )

test = CorrelationTest(log.(df.dist), df.dni2)


################################################################

dftrade = DataFrame(CSV.File(directory*"tradeshare-df-"*yyy*".csv"))

dftrade[!,"trade"] = log.(dftrade[!,"norm_tradeshare"] )

# removing the home trade flows
filter!(row -> ~(row.norm_tradeshare ≈ 1.0), dftrade);
dfcountryfix = deepcopy(dftrade)

# remove the zero trade flows
filter!(row -> ~(row.norm_tradeshare ≈ 0.0), dftrade);

Ncntry = 30

θ = 5.0
σ = 1.5

grv_params = gravity_params(Ncntry = Ncntry, θ = θ, L = ones(Ncntry), dfcntryfix = dfcountryfix )

# ################################################################
# # Run the Gravity regression

grvdata = gravity(dftrade, display = false);

# # ################################################################
# # # Recover the trade costs and technology parameters

d = zeros(Ncntry,Ncntry)
T = zeros(Ncntry)
W = ones(Ncntry)

make_trade_costs!(grvdata, d, grv_params)

make_technology!(grvdata, T, W, grv_params)

trd_prm = trade_params(θ = grv_params.θ, σ = σ, d = d, S = exp.(grvdata.S), Ncntry = grv_params.Ncntry, N = grv_params.L)

# p = SciMLBase.NullParameters()

out = estimate_θ_σ_dni(df, grv_params, trd_prm,
            grvdata; model = "bejk", Wmat = [], display = true, Nruns = 12, Ngoods = 50000)


# prob = OptimizationProblem(g, [4.5, 2.5 ], p, lb = lb, ub = ub)

# sol = Optimization.solve(prob, BOBYQA(); rhoend = 1e-4)


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

# #Todo list:
# ```

# 1) First run the gravity regression
# 2) Take predicted trade flows and add error terms
# 3) Run the gravity regression on the simmulated trade data
# 4) Take paremters from 3) and simmulate the model (EK, BEJK, etc)
# 5) Sample prices (maybe with measurement error)
# 6) Compute moments
# 7) Construct weighting matrix and J-stats? (neet do fill this in)

# ```
