using LinearAlgebra  # For approximate comparisons
using Statistics


########################################################################################
########################################################################################

function construct_tradematrix(input_wtf, input_isic)
    # Adjust the trade data based on the input ISIC codes
    # and construct the bilateral trade flow matrix.

    cntries = 1:134  # Hardcoded country codes
    Ncntry = length(cntries)
    Ncode = length(input_isic)

    # Initialize the trade matrix
    trademat = zeros(Ncntry, Ncntry)

    for importer in 1:Ncntry
        # Filter trade flows for the current importer
        yyy = input_wtf[:, 5] .== cntries[importer]
        trade_flows = input_wtf[yyy, :]

        for exporter in 1:Ncntry
            # Filter trade flows for the current exporter
            www = trade_flows[:, 4] .== cntries[exporter]
            trade_flow_new = trade_flows[www, :]

            if !isempty(trade_flow_new)
                for isicode in 1:Ncode
                    # Filter trade flows for the current ISIC code
                    zzz = trade_flow_new[:, 1] .== input_isic[isicode]

                    if sum(zzz) == 0
                        # Skip if no trade flows exist for this ISIC code
                        continue
                    end

                    # Record the trade flows for the current ISIC code
                    trade_flow_newc = trade_flow_new[zzz, :]
                    trademat[exporter, importer] += sum(trade_flow_newc[:, 2])
                end
            end
        end
    end

    # Calculate and display the percentage of zeros in the trade matrix
    zeross = sum(trademat .== 0) - Ncntry
    println("Percent Zeros: ", zeross / (Ncntry^2 - Ncntry))

    return trademat
end

########################################################################################
########################################################################################



function aggregate_drop(trade_mat, drop_30)
    # It aggregates and then drops the countries for which we do not have data.

    new_trade_mat = copy(trade_mat)
    # new_gross_output = copy(total_gross_output)

    # Helper function for approximate equality
    function approx_equal(a, b; tol=1e-6)
        return abs(a - b) < tol
    end

    # Aggregate Belgium with Netherlands
    new_trade_mat[:, 9] .= trade_mat[:, 9] .+ trade_mat[:, 92]
    
    new_trade_mat[9, :] .= trade_mat[9, :] .+ trade_mat[92, :]
    
    # new_gross_output[9] += total_gross_output[92]

    # Aggregate Singapore with Malaysia
    new_trade_mat[:, 80] .= trade_mat[:, 80] .+ trade_mat[:, 109]
    
    new_trade_mat[80, :] .= trade_mat[80, :] .+ trade_mat[109, :]
    
    # new_gross_output[80] += total_gross_output[109]

    # Aggregate China with Hong Kong and Macau
    new_trade_mat[:, 27] .= trade_mat[:, 27] .+ trade_mat[:, 54] .+ trade_mat[:, 76]
    
    new_trade_mat[27, :] .= trade_mat[27, :] .+ trade_mat[54, :] .+ trade_mat[76, :]
    
    # # new_gross_output[27] += total_gross_output[54] + total_gross_output[76]

    # # # Define the countries to drop
    # drop = [20, 54, 73, 75, 76, 78, 80, 82, 87, 92, 107, 109, 112, 134]
    # # # Malaysia (80) and Singapore (109) are dropped anyways
    # # # saudia arabia (76) too 
    # # # drop south africa (112)


    # # Remove rows and columns corresponding to dropped countries
    # new_trade_mat = new_trade_mat[setdiff(1:end, drop), setdiff(1:end, drop)]
    # new_gross_output = new_gross_output[setdiff(1:end, drop)]

    # # Adjust drop_30 to exclude dropped countries
    # drop_30 = drop_30[setdiff(1:end, drop)]

    # Take only the top 30 countries using approximate comparison
    logical_indices = drop_30 .== 1

    new_trade_mat = new_trade_mat[logical_indices, logical_indices]
    # new_gross_output = new_gross_output[drop_30 .≈ 1]

    return new_trade_mat

end

########################################################################################
########################################################################################

function construct_tradeshare(trademat::Matrix{Float64}, gross_output)
    """
    Constructs the trade share matrix from the trade matrix and gross output data.
    Rows represent exporters, and columns represent importers.

    Arguments:
    - trademat: A matrix where rows are exporters and columns are importers.
    - gross_output: A vector of gross output values for each country.

    Returns:
    - tradeshare: A matrix of trade shares.
    """

    # Create an empty DataFrame
    df = DataFrame()

    df = DataFrame(exporter = String[], importer=String[], tradeshare=Float64[], norm_tradeshare=Float64[])

    Ncntry = nrow(gross_output)

    # Initialize matrices and vectors
    tradeshare = zeros(Float64, Ncntry, Ncntry)

    denominator = zeros(Float64, Ncntry)
    
    total_exports = zeros(Float64, Ncntry)
    
    total_imports = zeros(Float64, Ncntry)


    # Loop through each country
    for cntry = 1:Ncntry

        total_exports[cntry] = sum(trademat[cntry, :])

        total_imports[cntry] = sum(trademat[:, cntry])
        
        denominator[cntry] = gross_output[cntry, "Value"]

        tradeshare[:, cntry] .= trademat[:, cntry] ./ (denominator[cntry] + total_imports[cntry] - total_exports[cntry])
        
        tradeshare[cntry, cntry] = 1.0 - sum(tradeshare[:, cntry]) + tradeshare[cntry, cntry]
        # Adjust for positive trade flows with a country's self in the trademat matrix
        # not sure why this last addition is necessary

        # push!(df, (gross_output[!,"ISO_Code"], gross_output[cntry,"ISO_Code"], tradeshare[:, cntry]))

    end

    #####################################################################################

    for importer = 1:Ncntry
        for exporter = 1:Ncntry
            # Append the data to the DataFrame
            push!(df, (gross_output[exporter, "ISO_Code"], gross_output[importer, "ISO_Code"], tradeshare[exporter, importer],
            tradeshare[exporter, importer] / tradeshare[importer, importer]))
        end
    end

    return tradeshare, df
end

########################################################################################
########################################################################################


function adjust_price_data(price_data)
    # Extract the 'istraded' column
    istraded = price_data[2, 6:end]

    # Extract 'drop_30' and 'fx_rate'
    drop_30 = price_data[3:end, 1]
    fx_rate = price_data[3:end, 5]

    # Extract 'ppp_basic_headings'
    ppp_basic_headings = price_data[3:end, 6:end]

    # List of countries to drop for problematic reasons
    drop = [20, 54, 73, 75, 76, 78, 80, 82, 87, 92, 107, 109, 112, 134]

     # Remove rows corresponding to 'drop' using logical indexing
    keep_rows = setdiff(1:size(ppp_basic_headings, 1), drop)  # Rows to keep
    drop_30 = drop_30[keep_rows]
    ppp_basic_headings = ppp_basic_headings[keep_rows, :]
    fx_rate = fx_rate[keep_rows]

    # Filter rows where 'drop_30' is not zero
    valid_rows = drop_30 .!= 0
    drop_30 = drop_30[valid_rows]
    ppp_basic_headings = ppp_basic_headings[valid_rows, :]
    fx_rate = fx_rate[valid_rows]

    # Create 'fx_mat' by repeating 'fx_rate' across columns
    fx_mat = repeat(fx_rate, 1, size(istraded, 1))

    # Compute 'pricemat'
    pricemat = ppp_basic_headings ./ fx_mat

    ttt = istraded .≈ 0.0

    pricemat = pricemat[:, .!vec(ttt)]

    return pricemat
end

########################################################################################
########################################################################################

function build_dni(pricemat, gross_output, tradeshare)

        cntry = nrow(gross_output)

        df = DataFrame(exporter = String[], importer=String[], Xni = Float64[], logXni = Float64[], dni=Float64[], τni=Float64[])
    
        # Log of price matrix
        log_p = log.(pricemat)

        # Initialize matrices
        dni = zeros(cntry, cntry)
        τni = zeros(cntry, cntry)
        dni2 = zeros(cntry, cntry)
    
        # Compute price differences
        for importer in 1:cntry

            for exporter in 1:cntry
                # Compute the price difference
                pdiff = log_p[importer, :] .- log_p[exporter, :]
    
                # Sort the differences
                h = sortperm(pdiff)
    
                # Take the max and second max
                num = pdiff[h[end]]
                num2 = pdiff[h[end - 1]]
    
                # Compute the mean price difference
                den = mean(pdiff)
    
                # Compute proxies for aggregate price differences
                dni[exporter, importer] = num - den
                dni2[exporter, importer] = num2 - den
                τni[exporter, importer] = num

                push!(df, (gross_output[exporter, "ISO_Code"], gross_output[importer, "ISO_Code"], 
                tradeshare[exporter, importer] / tradeshare[exporter, exporter],
                log(tradeshare[exporter, importer] / tradeshare[exporter, exporter]),
                dni[exporter, importer],
                τni[exporter, importer]))
            
            end
        
        end
    
        # # Set up the normalized trade matrix trdx
        # trdx = trade_mat ./ reshape(diag(trade_mat), 1, cntry)
    
        # # Exclude diagonal entries and zeros
        # vvv = (trdx .== 0) .| (trdx .== 1)
    
    
        # Output
        return df
    end

#########################################################################################
#########################################################################################

function construct_tradeshare_baci(tradedf, gross_output)
    
#         # Create an empty DataFrame
         df = DataFrame()
    
         df = DataFrame(exporter = String[], importer=String[], tradeshare=Float64[], norm_tradeshare=Float64[])
    
         Ncntry = nrow(gross_output)
    
#         # Initialize matrices and vectors
         tradeshare = zeros(Float64, Ncntry, Ncntry)
    
         denominator = zeros(Float64, Ncntry)
        
         total_exports = zeros(Float64, Ncntry)
        
         total_imports = zeros(Float64, Ncntry)
    
#         # Loop through each country
         for (index, cntry) in enumerate(eachrow(gross_output))

            println("A: ", cntry.ISO_Code)

            total_exports[index] = sum(tradedf[tradedf.exporter .== cntry.ISO_Code, "value"])
    
            total_imports[index] = sum(tradedf[tradedf.importer .== cntry.ISO_Code, "value"])
            
            denominator[index] = cntry["Value"]

            # println("B: ", cntry.ISO_Code, " ", total_exports[index], " ", total_imports[index], " ", denominator[index])

            for (expr_index, expr) in enumerate(eachrow(gross_output))

                foo = tradedf[tradedf.importer .== cntry.ISO_Code, :]
                # grab all the stuff the country is importing

                if expr_index != index

                    # print(foo[foo.exporter .== expr.ISO_Code, "value"][1])
                    
                    tradeshare[expr_index, index] = foo[foo.exporter .== expr.ISO_Code, "value"][1] / (denominator[index] + total_imports[index] - total_exports[index])
                    # given the importer (index), this is how much they purchase from the exporter (expr_index)
                    # then it is divided by the denominator (gross output) of the importer (index) + total imports - total exports
                    # or how much they consume
                end

            end

            tradeshare[index, index] = 1.0 - sum(tradeshare[:, index])

        end
    
#         #####################################################################################
    
        for importer = 1:Ncntry
            for exporter = 1:Ncntry
                # Append the data to the DataFrame
                push!(df, (gross_output[exporter, "ISO_Code"], gross_output[importer, "ISO_Code"], tradeshare[exporter, importer],
                tradeshare[exporter, importer] / tradeshare[importer, importer]))
            end
        end
    
        return tradeshare, df
    end