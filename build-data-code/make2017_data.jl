using DelimitedFiles
using XLSX
using DataFrames
using CSV
using MAT
include("functions_make_trade_share.jl")

grav_file = "./build-data-code/make-gravity-var/top30_gravity_data.csv"

df = DataFrame(CSV.File(grav_file))

###########################################################################################
# to make the 2017 dataset

gross_output_file = "./build-data-code/UNIDO-data/UNIDO-grossoutput-data-2017.csv"

gross_output = DataFrame(CSV.File(gross_output_file))

trade_file = "./build-data-code/baci-trade-data/top_30_trade_2017.csv"

tradedf = DataFrame(CSV.File(trade_file))

tradeshare, tradeshare_df = construct_tradeshare_baci(tradedf, gross_output)

###########################################################################################

rename!(tradeshare_df, Dict("exporter" => "iso_o", "importer" => "iso_d"))

trade_grav_df = innerjoin(tradeshare_df, df, on = ["iso_o", "iso_d"])

outfile = "tradeshare-df-2017.csv"

CSV.write(outfile, trade_grav_df, header = true)  # Write the DataFrame to a CSV file

###########################################################################################

price_file = "./build-data-code/2017-price-gap/2017_traded_prices.csv"

pricedf = DataFrame(CSV.File(price_file))

pricemat = Matrix(pricedf[:, 2:end])

df = build_dni(pricemat, gross_output, tradeshare)

outfile = "pricegap-df-2017.csv"

CSV.write(outfile, df, header = true)  # Write the DataFrame to a CSV file