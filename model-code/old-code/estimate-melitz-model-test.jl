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

#year = ["2004", "2011", "2017"]
year = ["2017"]
model = "melitz"
Nruns = 36
Nboots = 10
Ngoods = 50000
σ = 1.5

# dirname = "./data/"
# date = "12926-sigma-theta"

###############################################################################################################################
# method = "over"


# dfout = estimate_all(year, σ, model, method, Nruns, Nboots, dirname; Ngoods = Ngoods)

# CSV.write("./results/"*model*"-estimate-"*method*"-"*date*".csv", dfout; writeheader = true)

##############################################################################################################################

# method = "exact"
    
# dfout = estimate_all(year, σ, model, method, Nruns, Nboots, dirname; Ngoods = Ngoods)

# CSV.write("./results/"*model*"-estimate-"*method*"-"*date*".csv", dfout; writeheader = true)

directory = "./data/"

yyy = "2017"

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

#################################################################
#Run the Gravity regression

grvdata = gravity(dftrade, display = false);

################################################################
#Recover the trade costs and technology parameters

d = zeros(Ncntry,Ncntry)
T = zeros(Ncntry)
W = ones(Ncntry)

make_trade_costs!(grvdata, d, grv_params)

make_technology!(grvdata, T, W, grv_params)

trd_prm = trade_params(θ = grv_params.θ, σ = σ, d = d, S = exp.(grvdata.S), Ncntry = grv_params.Ncntry, N = grv_params.L)

model = "melitz-model2"
method = "exact"
Nprices = 60

sim_dni = generate_moments(trd_prm, Nruns; model = model, method = method, Nprices = Nprices, Ngoods = Ngoods) 


# sim_dni, sim_dni2 = generate_moments(foo, Nruns; model = model, method = "over", Nprices = Nprices, Ngoods = Ngoods) 