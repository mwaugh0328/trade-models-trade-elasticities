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

year = "2017"

df = DataFrame(CSV.File("./data/pricegap-df-"*year*".csv"))

grav_file = "./data/top30_gravity_data.csv"

dfgrav = DataFrame(CSV.File(grav_file))

rename!(df, Dict("exporter" => "iso_o", "importer" => "iso_d"))

dftariffs = DataFrame(CSV.File("./data/tariffs-"*year*".csv"))

rename!(dftariffs, Dict("exporter" => "iso_o", "importer" => "iso_d"))

df = innerjoin(df, dftariffs, on = ["iso_o", "iso_d"])

pricegap_df = innerjoin(df, dfgrav, on = ["iso_o", "iso_d"])

filter!(row -> ~(row.Xni ≈ 1.0), pricegap_df);

filter!(row -> ~(row.Xni ≈ 0.0), pricegap_df);

println( mean(pricegap_df.logXni) / mean(pricegap_df.dni) )

################################################################################################################


test = CorrelationTest(log.(1.0 .+ 0.01*pricegap_df.tariff ), pricegap_df.dni)
println("Correlation of Price Gaps and Tariffs")
println(test)

test = CorrelationTest(log.(pricegap_df.dist), pricegap_df.dni)

println("Correlation of Price Gaps and Distance")
println(test)

test = CorrelationTest(pricegap_df.border, pricegap_df.dni)

println("Correlation of Price Gaps and Border")
println(test)

################################################################################################################

# outreg = reg(pricegap_df, @formula(dni ~ fe(importer) + fe(exporter) + border + log(dist)), save = true, tol = 1e-10)

# outreg = reg(pricegap_df, @formula(logXni ~ fe(importer) + fe(exporter) + border + log(dist) + log(1.0 + 0.01*tariff)), save = true, tol = 1e-10)

# print(outreg)

outreg = reg(pricegap_df, @formula((dni) ~ fe(importer) + fe(exporter) + border + log(dist) + log(1.0 + 0.01*tariff)), save = true, tol = 1e-10)

println(outreg)

outreg = reg(pricegap_df, @formula((dni) ~  border + log(dist) + log(1.0 + 0.01*tariff)), save = true, tol = 1e-10)

println(outreg)

# results = compute_triplet_ratios(pricegap_df)

# clean_results = filter(row -> !isnan(row.lhs) && !isnan(row.rhs) && 
#                             !isinf(row.lhs) && !isinf(row.rhs) && 
#                             row.lhs > 0 && row.rhs > 0, results)


# # # Run regression
# reg_model = reg(clean_results, @formula(log(lhs) ~ log(rhs)))

# println(reg_model)