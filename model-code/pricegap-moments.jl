include("gravity-tools.jl")
include("trade-environment.jl")
include("simmulate-eaton-kortum.jl")
using CSV
using DataFrames
using Plots
using MINPACK

################################################################
# builds the EK dataset

# dftrade, dfcntryfix, dflabor = make_ek_dataset()
# # this one has the country numbers which allows for the construction of the 
# # trade costs given the estimated fixed effects from the gravity regressiondf

df = DataFrame(CSV.File("./data/pricegap-df-2017.csv"))

grav_file = "./data/top30_gravity_data.csv"

dfgrav = DataFrame(CSV.File(grav_file))

rename!(df, Dict("exporter" => "iso_o", "importer" => "iso_d"))

pricegap_df = innerjoin(df, dfgrav, on = ["iso_o", "iso_d"])

filter!(row -> ~(row.Xni ≈ 1.0), pricegap_df);

filter!(row -> ~(row.Xni ≈ 0.0), pricegap_df);

println( mean(pricegap_df.logXni) / mean(pricegap_df.dni) )

# outreg = reg(pricegap_df, @formula(log(dist) ~ fe(importer) + fe(exporter) + dni), save = true, tol = 1e-10)

outreg = reg(pricegap_df, @formula(τni ~ fe(importer) + fe(exporter) + border + log(dist)), save = true, tol = 1e-10)