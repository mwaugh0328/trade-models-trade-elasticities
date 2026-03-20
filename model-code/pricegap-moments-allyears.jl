include("gravity-tools.jl")
include("trade-environment.jl")
include("simulate-trade-models.jl")
using CSV
using DataFrames
using Plots
using MINPACK
 using HypothesisTests

################################################################

years = ["2004", "2011", "2017"]
all_dfs = DataFrame[]

for year in years
    df = DataFrame(CSV.File("./data/pricegap-df-" * year * ".csv"))
    dftariffs = DataFrame(CSV.File("./data/tariffs-" * year * ".csv"))
    grav_file = "./data/top30_gravity_data.csv"
    dfgrav = DataFrame(CSV.File(grav_file))

    rename!(df, Dict("exporter" => "iso_o", "importer" => "iso_d"))
    
    rename!(dftariffs, Dict("exporter" => "iso_o", "importer" => "iso_d"))

    df = innerjoin(df, dftariffs, on = ["iso_o", "iso_d"])

    pricegap_df = innerjoin(df, dfgrav, on = ["iso_o", "iso_d"])

    pricegap_df.year = fill(year, nrow(pricegap_df))  # Add year column

    push!(all_dfs, pricegap_df)
end

big_df = vcat(all_dfs...)

# Filter as before
filter!(row -> ~(row.Xni ≈ 1.0), big_df)
filter!(row -> ~(row.Xni ≈ 0.0), big_df)


################################################################

test = CorrelationTest(log.(big_df.dist), (big_df.dni))

println(" ")
println(" ")

println("Correlation of Price Gaps and Distance")
println(test)
ci_90 = confint(test, 0.90)
println("90% confidence interval: ", ci_90)


################################################################

test = CorrelationTest(big_df.border, (big_df.dni))

println(" ")
println(" ")
println("Correlation of Price Gaps and Border")
println(test)
ci_90 = confint(test, 0.90)
println("90% confidence interval: ", ci_90)

################################################################
test = CorrelationTest(log.(1.0 .+ 0.01*big_df.tariff ), (big_df.dni))
println(" ")
println(" ")
println("Correlation of Price Gaps and Tariffs")
println(test)
ci_90 = confint(test, 0.90)
println("90% confidence interval: ", ci_90)



outreg = reg(big_df, @formula(log(dni) ~  fe(year) + border + log(dist) + log(1.0 + 0.01*tariff)), save = true, tol = 1e-10)

# # Run regression with year fixed effect

println(outreg)


outreg = reg(big_df, @formula(log(dni) ~  fe(year) + fe(importer) + fe(exporter) +border + log(dist) + log(1.0 + 0.01*tariff)), save = true, tol = 1e-10)

# # Run regression with year fixed effect

println(outreg)