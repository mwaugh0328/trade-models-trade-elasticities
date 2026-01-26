function estimate_all(years, model, method, Nruns, Nboots, directory; Ngoods = 100000)

    dfout = DataFrame()

    for yyy in years

        if yyy == "2004"

            dfout = DataFrame()

        end

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
        σ = 2.5

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

    # # ################################################################
    # # # Now simmulate the EK model

        @time θest, Jstat, model_moments, data_moments =  estimate_θ_dni(df, grv_params, trd_prm,  grvdata; 
            model = model, method = method, Wmat = "optimal", display = false, Nruns = Nruns, Ngoods = Ngoods)


        println("Year: ", yyy)
        println(" ")
        println("Estimated θ: ", θest, " J-stat: ", Jstat)
        println(" Model moments: ", model_moments, " Data moments: ", data_moments)
        println(" ")

        θinterval, Jinterval = boot_strap_simulation(θest, grvdata.σν, dftrade, trd_prm, grv_params; 
            model = model, method = method, code = 269, Nboots = Nboots, Nruns = Nruns, Ngoods = Ngoods)

        p90 = quantile(θinterval, 0.9)
        p10 = quantile(θinterval, 0.10)

        println("90% confidence interval: ", p10, " - ", p90)

        Jpercentile = sum(x -> x <= Jstat, Jinterval) / length(Jinterval)

        println("J-stat percentile: ", Jpercentile)

        # dfout = vcat(dfout, DataFrame(θ = θest, J = Jstat, 
        #     p10 = p10, p90 = p90, Jpercentile = Jpercentile,
        #     model = model, method = method, year = yyy))

        if method == "over"

            dfout = vcat(dfout, DataFrame(θ = θest, J = Jstat, 
            p10 = p10, p90 = p90, Jpercentile = Jpercentile,
            model_moment_1 = model_moments[1], data_moment_1 = data_moments[1],
            model_moment_2 = model_moments[2], data_moment_2 = data_moments[2],
            model_moment_3 = model_moments[3], data_moment_3 = data_moments[3], 
                model = model, method = method, year = yyy))

        else

            dfout = vcat(dfout, DataFrame(θ = θest, J = Jstat,
            p10 = p10, p90 = p90, Jpercentile = Jpercentile,
            model_moment = model_moments[1], data_moment = data_moments[1],
                model = model, method = method, year = yyy))

        end

        dfout = dfout

    end

    return dfout

end