using DelimitedFiles
using XLSX
using DataFrames
using CSV
include("functions_make_trade_share.jl")

###########################################################################################

drop_30 = DataFrame(CSV.File("drop_30.csv"))

input_isic = readdlm("./trade-data/isic_3digit_manuf.txt")

grav_file = "../make-gravity-var/top30_gravity_data.csv"

df = DataFrame(CSV.File(grav_file))

years = ["2004", "2011"]

#year = "2004"

for year in years

    # Load the trade flow data from the text file
    # Assuming the file name is constructed as "wtf_<year>.mat.txt"

    trade_flow_file = "./trade-data/wtf_" * year[end-1:end] * "mat.txt"

    input_wtf = readdlm(trade_flow_file, ',', Float64)

    # Construct the aggregate bilateral trade flow matrix

    trademat = construct_tradematrix(input_wtf, input_isic)

    new_trade_mat = aggregate_drop(trademat, drop_30[:, "new top 30"])

    # grab the gross output data
    gross_output_file = "./UNIDO-data/UNIDO-grossoutput-data-" * year * ".csv"

    gross_output = DataFrame(CSV.File(gross_output_file))

    tradeshare, tradeshare_df = construct_tradeshare(new_trade_mat, gross_output)

    #now add in the gravity variables

    rename!(tradeshare_df, Dict("exporter" => "iso_o", "importer" => "iso_d"))

    trade_grav_df = innerjoin(tradeshare_df, df, on = ["iso_o", "iso_d"])

    outfile = "tradeshare-df-"*year*".csv"

    CSV.write(outfile, trade_grav_df, header = true)  # Write the DataFrame to a CSV file

end




###########################################################################################

# year = "2004"

# trade_flow_file = "./trade-data/wtf_" * year[end-1:end] * "mat.txt"

# # input_wtf = readdlm(trade_flow_file, ',', Float64)

# input_wtf = readdlm(trade_flow_file, ',', Any)  # Read as Any to handle mixed types
# input_wtf = replace(input_wtf, "" => NaN)      # Replace empty strings with NaN
# input_wtf = Float64.(input_wtf)                # Convert to Float64


# # Load the total output data from the Excel file
# sheet = XLSX.readxlsx("./old-code/make-tradeshare-2011/output_data_2011.xlsx")["Sheet1"]

# total_output = [x isa Float64 ? x : tryparse(Float64, x) === nothing ? NaN : tryparse(Float64, x) for x in sheet[2:135, 6]]

# drop_30 = sheet[2:135, 7]


# # Construct the trade share matrix
# # Assuming `construct_tradeshare` is implemented elsewhere
# tradeshare = construct_tradeshare(new_trade_mat, new_output)

# # Uncomment the following lines if gravity variables are needed
# # dist_mat = readdlm("dist_mat.txt")
# # d_mat, b_mat, e_code, i_code = construct_gravity_var(dist_mat, drop_30)

# # Save the trade share matrix and other variables if needed
# # save("trade_data_estimation.jld2", "tradeshare", tradeshare, "d_mat", d_mat, "b_mat", b_mat)

# # Uncomment to create gravity regression dataset
# # home_share = diagm(tradeshare)
# # grav_trade = tradeshare ./ permutedims(home_share, (2, 1))
# # grav_data_set = hcat(vec(i_code), vec(e_code), vec(grav_trade), vec(d_mat) ./ 1.6, vec(b_mat))
# # writedlm("grav_data.csv", grav_data_set, ',')

# println("Trade share matrix construction complete.")

# "./old-code/make-tradeshare-2011/output_data_2011.xlsx"

# df = DataFrame(XLSX.readtable("./old-code/make-tradeshare-2011/output_data_2011.xlsx", "Sheet1")

# selected_columns = df[:, ["New code", "Country code", "new top 30"]]

# top30 = selected_columns[selected_columns[:, "new top 30"] .== 1, :]