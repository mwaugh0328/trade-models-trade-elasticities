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

model = "ek"
method = "over"
Nruns = 36
Nboots = 100

dirname = "./data/"

dfout = estimate_all(year, model, method, Nruns, Nboots, dirname)

################################################################

# # dfout = DataFrame()

# for yyy in year

#     if yyy == "2004"

#         dfout = DataFrame()

#     end

#     df = DataFrame(CSV.File("./data/pricegap-df-"*yyy*".csv"))

#     grav_file = "./data/top30_gravity_data.csv"

#     dfgrav = DataFrame(CSV.File(grav_file))

#     rename!(df, Dict("exporter" => "iso_o", "importer" => "iso_d"))

#     df = innerjoin(df, dfgrav, on = ["iso_o", "iso_d"])

#     filter!(row -> ~(row.Xni ≈ 1.0), df);

#     # filter!(row -> ~(row.Xni ≈ 0.0), df);

#     # println( mean(df.logXni) / mean(df.dni) )

#     test = CorrelationTest(log.(df.dist), df.dni2)


#     ################################################################

#     dftrade = DataFrame(CSV.File("./data/tradeshare-df-"*yyy*".csv"))

#     dftrade[!,"trade"] = log.(dftrade[!,"norm_tradeshare"] )

#     # removing the home trade flows
#     filter!(row -> ~(row.norm_tradeshare ≈ 1.0), dftrade);
#     dfcountryfix = deepcopy(dftrade)

#     # remove the zero trade flows
#     filter!(row -> ~(row.norm_tradeshare ≈ 0.0), dftrade);

#     Ncntry = 30

#     θ = 4.0
#     σ = 2.5

#     grv_params = gravity_params(Ncntry = Ncntry, θ = θ, L = ones(Ncntry), dfcntryfix = dfcountryfix )

#     # ################################################################
#     # # Run the Gravity regression

#     grvdata = gravity(dftrade, display = true);

# # # ################################################################
# # # # Recover the trade costs and technology parameters

#     d = zeros(Ncntry,Ncntry)
#     T = zeros(Ncntry)
#     W = ones(Ncntry)

#     make_trade_costs!(grvdata, d, grv_params)

#     make_technology!(grvdata, T, W, grv_params)

#     trd_prm = trade_params(θ = grv_params.θ, σ = σ, d = d, S = exp.(grvdata.S), Ncntry = grv_params.Ncntry, N = grv_params.L)

# # # ################################################################
# # # # Now simmulate the EK model

#     θest, Jstat =  estimate_θ_dni(df, grv_params, trd_prm,  grvdata; 
#         model = model, method = method, Wmat = "optimal", display = true, Nruns = Nruns)

#     println("Estimated θ: ", θest, " J-stat: ", Jstat)
#     println(" ")

#     θinterval, Jinterval = boot_strap_simulation(θest, grvdata.σν, dftrade, grv_params; 
#         model = model, method = method, code = 269, Nboots = Nboots, Nruns = Nruns)

#     p90 = quantile(θinterval, 0.9)
#     p10 = quantile(θinterval, 0.10)

#     println("90% confidence interval: ", p10, " - ", p90)

#     Jpercentile = sum(x -> x <= Jstat, Jinterval) / length(Jinterval)

#     println("J-stat percentile: ", Jpercentile)

#     dfout = vcat(dfout, DataFrame(θ = θest, J = Jstat, 
#         p10 = p10, p90 = p90, Jpercentile = Jpercentile,
#         model = model, method = method, year = yyy))

# end


# CSV.write("./results/ek-estimate-"*method*".csv", dfout; writeheader = true)