using DelimitedFiles
using XLSX
using DataFrames
using CSV
using MAT
include("functions_make_trade_share.jl")
###########################################################################################

input_isic = readdlm("./make-tradeshare-pricegap/trade-data/isic_3digit_manuf.txt")

drop_30 = DataFrame(CSV.File("./make-tradeshare-pricegap/drop_30.csv"))

###########################################################################################

data_set = XLSX.readxlsx("./make-tradeshare-pricegap/price-data/clean_basic_headings_2011.xlsx")

# Extract the sheet (replace "Sheet1" with the actual sheet name if needed)
price_data = data_set["Sheet1"][1:136,1:160]

pricemat = adjust_price_data(price_data)

trade_flow_file = "./make-tradeshare-pricegap/trade-data/wtf_11mat.txt"

input_wtf = readdlm(trade_flow_file, ',', Float64)

# Construct the aggregate bilateral trade flow matrix

trademat = construct_tradematrix(input_wtf, input_isic)

new_trade_mat = aggregate_drop(trademat, drop_30[:, "new top 30"])

# grab the gross output data
gross_output_file = "./make-tradeshare-pricegap/UNIDO-data/UNIDO-grossoutput-data-2011.csv"

gross_output = DataFrame(CSV.File(gross_output_file))

tradeshare, tradeshare_df = construct_tradeshare(new_trade_mat, gross_output)

# df = build_dni(pricemat, gross_output, tradeshare)

# outfile = "pricegap-df-2011.csv"

# CSV.write(outfile, df, header = true)  # Write the DataFrame to a CSV file

###########################################################################################
###########################################################################################

# mat_data = matread("./old-code/estimation_mat_30_2014-version.mat")

# pricemat = mat_data["pmat_30"]

# istraded = mat_data["istraded"]

# ttt = istraded .≈ 0.0

# pricemat = pricemat[:, .!vec(ttt)]

# trade_flow_file = "./make-tradeshare-pricegap/trade-data/wtf_04mat.txt"

# input_wtf = readdlm(trade_flow_file, ',', Float64)

# # Construct the aggregate bilateral trade flow matrix

# trademat = construct_tradematrix(input_wtf, input_isic)

# new_trade_mat = aggregate_drop(trademat, drop_30[:, "new top 30"])

# # grab the gross output data
# gross_output_file = "./make-tradeshare-pricegap/UNIDO-data/UNIDO-grossoutput-data-2004.csv"

# gross_output = DataFrame(CSV.File(gross_output_file))

# tradeshare, tradeshare_df = construct_tradeshare(new_trade_mat, gross_output)

# df = build_dni(pricemat, gross_output, tradeshare)

# outfile = "pricegap-df-2004.csv"

# CSV.write(outfile, df, header = true)  # Write the DataFrame to a CSV file






