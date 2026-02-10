function make_trade_data(directory, year)
    # A function to load and prepare the trade data for estimation
    # df, dftrade, dfcountryfix = make_trade_data(directory, year[1])

    df = DataFrame(CSV.File(directory*"pricegap-df-"*year*".csv"))

    grav_file = directory*"top30_gravity_data.csv"

    dfgrav = DataFrame(CSV.File(grav_file))

    rename!(df, Dict("exporter" => "iso_o", "importer" => "iso_d"))

    df = innerjoin(df, dfgrav, on = ["iso_o", "iso_d"])

    dftariffs = DataFrame(CSV.File(directory*"tariffs-"*year*".csv"))

    rename!(dftariffs, Dict("exporter" => "iso_o", "importer" => "iso_d"))

    df = innerjoin(df, dftariffs, on = ["iso_o", "iso_d"])

    filter!(row -> ~(row.Xni ≈ 1.0), df);

    # filter!(row -> ~(row.Xni ≈ 0.0), df);

    # println( mean(df.logXni) / mean(df.dni) )

    # test = CorrelationTest(log.(df.dist), df.dni2)


    ################################################################

    dftrade = DataFrame(CSV.File(directory*"tradeshare-df-"*year*".csv"))

    dftrade[!,"trade"] = log.(dftrade[!,"norm_tradeshare"] )

    # removing the home trade flows
    filter!(row -> ~(row.norm_tradeshare ≈ 1.0), dftrade);
    dfcountryfix = deepcopy(dftrade)

    # remove the zero trade flows
    filter!(row -> ~(row.norm_tradeshare ≈ 0.0), dftrade);

    return df, dftrade, dfcountryfix

end

###############################################################
###############################################################


function estimate_all(years, σ, model, method, Nruns, Nboots, directory; Ngoods = 100000)

    dfout = DataFrame()


    for yyy in years

        if yyy == "2004"

            dfout = DataFrame()

        end

        if yyy == "2004"

            Nprices = 62

        elseif yyy == "2017"

            Nprices = 64

        else # this is if 2011

            Nprices = 71

        end
        
        # Load the trade data for use in the estimation
        df, dftrade, dfcountryfix = make_trade_data(directory, yyy)

        Ncntry = 30

        θ = 5.0 # initial value for the trade elasticity
        σ = σ

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

        E = get_E_by_importer(df)

        trd_prm = trade_params(θ = grv_params.θ, σ = σ, d = d, S = exp.(grvdata.S), Ncntry = grv_params.Ncntry, N = grv_params.L, E = E)

    # # ################################################################
    # # # Now simmulate the EK model

        @time θest, Jstat, model_moments, data_moments =  estimate_θ_dni(df, grv_params, trd_prm,  grvdata; 
            model = model, method = method, Wmat = "optimal", display = false, Nruns = Nruns, Nprices = Nprices, Ngoods = Ngoods)


        println("Year: ", yyy)
        println(" ")
        println("Estimated θ: ", θest, " J-stat: ", Jstat)
        println(" Model moments: ", model_moments, " Data moments: ", data_moments)
        println(" ")

        θinterval, Jinterval = boot_strap_simulation(θest, grvdata.σν, dftrade, trd_prm, grv_params; 
            model = model, method = method, code = 269, Nboots = Nboots, Nruns = Nruns, Nprices = Nprices, Ngoods = Ngoods)

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



###############################################################
###############################################################

function dni_moment_model(pricemat, πshares)
        # Number of countries
    Ncntry = size(πshares, 1)

    # Log of price matrix
    log_p = log.(pricemat)

    dni = Array{Float64}(undef, Ncntry, Ncntry)
    dni2 = Array{Float64}(undef, Ncntry, Ncntry)
    dni3 = Array{Float64}(undef, Ncntry, Ncntry)
    Xni = Array{Float64}(undef, Ncntry, Ncntry)

    c = size(log_p, 2)

    s5p = Int(round.(0.75 .* c ))
    e5p = Int(round.(0.85 .* c ))

        # Compute price differences
    for importer in 1:Ncntry

        for exporter in 1:Ncntry
            # Compute the price difference
            pdiff = log_p[importer, :] .- log_p[exporter, :]

            # Sort the differences
            h = sortperm(pdiff)

            # Take the max
            num = pdiff[h[end]]
            num2 = pdiff[h[e5p]];
            num3 = pdiff[h[s5p]];

            # Compute the mean price difference
            den = mean(pdiff)

            # Compute proxies for aggregate price differences
            dni[exporter, importer] = num - den
            dni2[exporter, importer] = num2 - den
            dni3[exporter, importer] = num3 - den

            Xni[exporter, importer] = πshares[importer, exporter] / πshares[exporter, exporter]
            # IMPORTANT the shares are fliped where row is always the importer and the column is the exporter
            # so when we map into the Xni we need to flip the order of the indices


        end

    end
        # Exclude zeros and diagonal entries
    # notzeros = (Xni .≈ 0.0) .| (Xni .≈ 1.0)
    notzeros = (Xni .≈ 1.0)

    return log.(Xni[.!notzeros]), dni[.!notzeros], dni2[.!notzeros], dni3[.!notzeros] 

end

###############################################################
###############################################################

function beta_moment_model(pricemat, πshares)
    # Number of countries
    Ncntry = size(πshares, 1)

    # Log of price matrix
    log_p = log.(pricemat)

    # Initialize matrices
    dni = Array{Float64}(undef, Ncntry, Ncntry)
    Xni = Array{Float64}(undef, Ncntry, Ncntry)


    # Compute price differences
    for importer in 1:Ncntry

        for exporter in 1:Ncntry
            # Compute the price difference
            pdiff = log_p[importer, :] .- log_p[exporter, :]

            # Sort the differences
            h = sortperm(pdiff)

            # Take the max
            num = pdiff[h[end]]

            # Compute the mean price difference
            den = mean(pdiff)

            # Compute proxies for aggregate price differences
            dni[exporter, importer] = num - den

            Xni[exporter, importer] = πshares[importer, exporter] / πshares[exporter, exporter]
            # IMPORTANT the shares are fliped where row is always the importer and the column is the exporter
            # so when we map into the Xni we need to flip the order of the indices


        end

    end

    # Exclude zeros and diagonal entries
    notzeros = (Xni .≈ 0.0) .| (Xni .≈ 1.0)

    # println("trade", mean( log.(Xni[.!notzeros]) ), mean(πshares[.!notzeros]) )

    β = mean( log.(Xni[.!notzeros]) ) / mean( dni[.!notzeros] )

    return β
end

###############################################################
###############################################################

function estimate_θ(θ, gravity_parameters, trade_parameters, gravity_results; model = "ek", Nruns = 10, Nprices = 70, Ngoods = 100000)
    # A function to estimate the θ parameter using the gravity equation and the EK model

    # println(θ)

    @unpack Ncntry = gravity_parameters

    d = rescale_trade_cots(θ, gravity_parameters, gravity_results)

    foo = trade_params(θ = θ, d = d, trade_parameters)

    # println( mean(foo.d) )

    avgβ = mean( generate_moments(foo, Nruns; model = model, Nprices = Nprices, Ngoods = Ngoods) )

    return avgβ

end

##############################################################################################################################
##############################################################################################################################

function estimate_θ_σ_dni(dfdni, gravity_parameters, trade_parameters,
            gravity_results; model = "bejk", Wmat = "optimal", display = true, Nruns = 10, Ngoods = 100000)

    if model != "bejk"
        error("This function is only for the BEJK model")
    end

    lb = [2.5, 1.25]
    ub = [10.0, 3.0]

    p = SciMLBase.NullParameters()

    g(x, p) = estimate_θ_σ_dni(x[1], x[2], dfdni.dni, dfdni.dni2, dfdni.dist,  gravity_parameters, trade_parameters,
            gravity_results; model = model, Wmat = Wmat, display = display, return_moments = false, Nruns = Nruns, Ngoods = Ngoods)


                # Define the function to estimate θ using the over method
    prob = OptimizationProblem(g, [4.5, 1.35], p, lb = lb, ub = ub)
        

    sol = Optimization.solve(prob, BOBYQA(); rhoend = 1e-4)

    # Compatibility fix for different versions of Optimization.jl
    solution_value = hasfield(typeof(sol), :u) ? sol.u : sol.x


    ~, data_moments, model_moments = estimate_θ_σ_dni(solution_value[1], solution_value[2], dfdni.dni, dfdni.dni2, dfdni.dist,  gravity_parameters, trade_parameters,
             gravity_results; model = model, Wmat = Wmat, display = display, return_moments = true, Nruns = Nruns, Ngoods = Ngoods)


    return solution_value, length(dfdni.dni) * sol.objective, model_moments, data_moments
    #return solution_value
end

##############################################################################################################################
##############################################################################################################################

function estimate_θ_dni(dfdni, gravity_parameters, trade_parameters,
            gravity_results; model = "ek", method = "exact", Wmat = "optimal", display = true, Nruns = 10, Nprices = 70, Ngoods = 100000)

    lb = [1.5,]
    ub = [10.0,]

    p = SciMLBase.NullParameters()

    # println("value of σ: ", trade_parameters.σ)

    f(x, p) = estimate_θ_dni(x[1], dfdni.dni, gravity_parameters, trade_parameters,
            gravity_results; model = model, display = display, return_moments = false, Nruns = Nruns, Nprices = Nprices, Ngoods = Ngoods)

    g(x, p) = estimate_θ_dni(x[1], dfdni.dni, dfdni.dni2, dfdni.dist,  gravity_parameters, trade_parameters,
            gravity_results; model = model, Wmat = Wmat, display = display, return_moments = false, Nruns = Nruns, Nprices = Nprices, Ngoods = Ngoods)

    if method == "exact"

        # Define the function to estimate θ using the exact method

            prob = OptimizationProblem(f, [4.5], p, lb = lb, ub = ub)

    else
                # Define the function to estimate θ using the over method
            prob = OptimizationProblem(g, [4.5], p, lb = lb, ub = ub)
        
    end

    sol = Optimization.solve(prob, BOBYQA(); rhoend = 1e-4)

    # Compatibility fix for different versions of Optimization.jl
    solution_value = hasfield(typeof(sol), :u) ? sol.u[1] : sol.x[1]

    if method == "exact"

        ~, data_moments, model_moments = estimate_θ_dni(solution_value, dfdni.dni, gravity_parameters, trade_parameters,
            gravity_results; model = model, display = display, return_moments = true, Nruns = Nruns, Nprices = Nprices, Ngoods = Ngoods)

    else

        ~, data_moments, model_moments = estimate_θ_dni(solution_value, dfdni.dni, dfdni.dni2, dfdni.dist,  gravity_parameters, trade_parameters,
            gravity_results; model = model, Wmat = Wmat, display = display, return_moments = true, Nruns = Nruns, Nprices = Nprices, Ngoods = Ngoods)

            

    end



    return solution_value, length(dfdni.dni) * sol.objective, model_moments, data_moments

end


##############################################################################################################################
##############################################################################################################################

function estimate_θ_σ_dni(θ, σ, dni, dni2, dist, gravity_parameters, trade_parameters, gravity_results; model = "bejk", Wmat = "optimal",
     Nruns = 10, Nprices = 70, Ngoods = 100000, display = false, return_moments = false)
    # A function to estimate the θ and sigma paramter for the BEJK model

    # println(θ)

    @unpack Ncntry = gravity_parameters

    d = rescale_trade_cots(θ, gravity_parameters, gravity_results)

    foo = trade_params(θ = θ, σ = σ, d = d, trade_parameters)

    # println( mean(foo.d) )

    sim_dni, sim_dni2 = generate_moments(foo, Nruns; model = model, method = "over", Nprices = Nprices, Ngoods = Ngoods) 

    sim_dni = mean(sim_dni, dims = 2)

    sim_dni2 = mean(sim_dni2, dims = 2)

    sim_cov = (sim_dni .- mean(sim_dni, dims = 1)) .* (log.(dist) .- mean(log.(dist)))
    
    model_moments = [sim_dni sim_dni2 sim_cov]
    # println(mean(sim_cov))

    data_cov =  (dni .- mean(dni, dims = 1)) .* (log.(dist) .- mean(log.(dist))) 

    data_moments = [dni dni2 data_cov]
    # println(mean(data_cov))
    
    hθ = mean(  data_moments .- model_moments, dims = 1)
    # print(size(hθ))

    if Wmat == "optimal"
        # Compute the optimal weighting matrix
        W = cov(data_moments .- model_moments)

    else
        # Use the identity matrix as the weighting matrix
        W = I(3)
    end

    zero_fun = hθ*(W^-1)*hθ'

    # println(hθ)
    # println(size(W))

    if display == true

        println("Zero function: ", zero_fun, " Value of θ: ", θ, " Value of σ: ", σ)

    end

    if return_moments == true

        return zero_fun, mean(data_moments, dims = 1), mean(model_moments, dims = 1)

    else

        return zero_fun

    end

end

##############################################################################################################################
##############################################################################################################################

function estimate_θ_dni(θ, dni, dni2, dist, gravity_parameters, trade_parameters, gravity_results; model = "ek", Wmat = "optimal",
     Nruns = 10, Nprices = 70, Ngoods = 100000, display = false, return_moments = false)
    # A function to estimate the θ parameter using the gravity equation and the EK model

    # println(θ)

    @unpack Ncntry = gravity_parameters

    d = rescale_trade_cots(θ, gravity_parameters, gravity_results)

    foo = trade_params(θ = θ, d = d, trade_parameters)

    # println( mean(foo.d) )

    # println("Number of Prices: ", Nprices)

    sim_dni, sim_dni2 = generate_moments(foo, Nruns; model = model, method = "over", Nprices = Nprices, Ngoods = Ngoods) 

    sim_dni = mean(sim_dni, dims = 2)

    sim_dni2 = mean(sim_dni2, dims = 2)

    sim_cov = (sim_dni .- mean(sim_dni, dims = 1)) .* (log.(dist) .- mean(log.(dist)))
    
    model_moments = [sim_dni sim_dni2 sim_cov]
    # println(mean(sim_cov))


    data_cov =  (dni .- mean(dni, dims = 1)) .* (log.(dist) .- mean(log.(dist))) 

    data_moments = [dni dni2 data_cov]
    # println(mean(data_cov))
    
    hθ = mean(  data_moments .- model_moments, dims = 1)
    # print(size(hθ))

    if Wmat == "optimal"
        # Compute the optimal weighting matrix
        W = cov(data_moments .- model_moments)

    else
        # Use the identity matrix as the weighting matrix
        W = I(3)
    end

    zero_fun = hθ*(W^-1)*hθ'

    # println(hθ)
    # println(size(W))

    if display == true

        println("Zero function: ", zero_fun, " Value of θ: ", θ)

    end

    if return_moments == true

        return zero_fun, mean(data_moments, dims = 1), mean(model_moments, dims = 1)

    else

        return zero_fun

    end

end

##############################################################################################################################

function estimate_θ_dni(θ, dni, gravity_parameters, trade_parameters, gravity_results; model = "ek", 
    Nruns = 10, Nprices = 70, Ngoods = 100000, display = false, return_moments = false)
    # A function to estimate the θ parameter for the exactly identified case using the dni moment only

    # println(θ)

    @unpack Ncntry = gravity_parameters

    d = rescale_trade_cots(θ, gravity_parameters, gravity_results)

    foo = trade_params(θ = θ, d = d, trade_parameters)

    # println( mean(foo.d) )

    sim_dni = mean( generate_moments(foo, Nruns; model = model, method = "exact", Nprices = Nprices, Ngoods = Ngoods) , dims = 2)
    #average accross the different simmulations

    hθ = mean( dni .- sim_dni, dims = 1)
    #this is the average differnce accros n,i pairs

    zero_fun = hθ*hθ'

    if display == true

        println("Zero function: ", zero_fun, " Value of θ: ", θ)

    end

    if return_moments == true

        return zero_fun, mean(dni, dims = 1), mean(sim_dni, dims = 1)

    else

    return zero_fun

    end

end

##############################################################################################################################


function generate_moments(trade_parameters, Nruns; model = "ek", method = "exact", code = 1, Nprices = 70, Ngoods = 200000)
    # multiple dispatch version of the generate_moments function to generate a bunch of betas

    # β = Array{Float64}(undef, Nruns)

    if method == "exact"

        dni = Array{Float64}(undef, trade_parameters.Ncntry^2 - trade_parameters.Ncntry, Nruns)

        Threads.@threads for xxx = 1:Nruns

        dni[:, xxx] = generate_moments(trade_parameters; model = model, method = method, code = code + xxx, Nprices = Nprices, Ngoods = Ngoods)

        end
        
        return dni

    elseif method == "over"

        dni = Array{Float64}(undef, trade_parameters.Ncntry^2 - trade_parameters.Ncntry, Nruns)

        dni2 = Array{Float64}(undef, trade_parameters.Ncntry^2 - trade_parameters.Ncntry, Nruns)

        Threads.@threads for xxx = 1:Nruns

        dni[:, xxx], dni2[:, xxx] = generate_moments(trade_parameters; model = model, method = method, code = code + xxx, Nprices = Nprices, Ngoods = Ngoods)

        end
        
        return dni, dni2

    else
         error("Please specify estimation method, over or exact'.")

    end

end

##############################################################################################################################
##############################################################################################################################

function generate_moments(trade_parameters; model = "ek", method = "over", code = 1, Nprices = 70, Ngoods = 100000)
    # A function to simmulate a pattern of trade and then generate a 
    # random sample of final goods prices, then compute the moments

    @unpack Ncntry = trade_parameters

    πshares = Array{Float64}(undef, length(Ncntry), length(Ncntry))

    prices = Array{Float64}(undef, length(Ncntry), Ngoods)

    if model == "ek"
        # println("this is the EK model")
        πshares, prices = sim_trade_pattern_ek(trade_parameters; Ngoods = Ngoods, code = code)

    elseif model == "bejk"
        # print("this is the BEJK model")
        #println("BEJK model with σ = ", trade_parameters.σ)
        πshares, prices = sim_trade_pattern_bejk(trade_parameters; Ngoods = Ngoods, code = code)

    elseif model == "krugman"
        # print("this is the Krugman model")
        πshares, prices = sim_trade_pattern_krugman(trade_parameters; Ngoods = Ngoods, code = code)

    elseif model == "krugman-model2"

        πshares, prices = sim_trade_pattern_krugman_model2(trade_parameters; Ngoods = Ngoods, code = code)


    elseif model == "melitz"
        # print("this is the Melitz model")
        prices = Array{Float64}(undef, length(Ncntry), Ngoods*Ncntry)

        common_set = falses(Ncntry * Ngoods)

        πshares, prices, common_set = sim_trade_pattern_melitz_optimized(trade_parameters; Ngoods = Ngoods, code = code)

        num_prices = size(prices[:,common_set],2)

    elseif model == "melitz-model2"

        πshares, prices = sim_trade_pattern_melitz_model2(trade_parameters; Ngoods = Ngoods, code = code)

        num_prices = size(prices, 2)


    else
        error("Model not recognized. Use 'ek' or 'bejk' or 'krugman' or 'krugman-model2' or 'melitz'.")
        
    end

    # print(size(prices))

    if model != "melitz"
        
        sampled_prices= sample(MersenneTwister(09212013 + code), 1:Ngoods, Nprices; replace=false)

        pmat = prices[:, sampled_prices]

    elseif model == "melitz"

        prices = prices[:,common_set]

        sampled_prices= sample(MersenneTwister(09212013 + code), 1:num_prices, Nprices; replace=false)

        pmat = prices[:, sampled_prices]

    elseif model == "melitz-model2"
        
        sampled_prices= sample(MersenneTwister(09212013 + code), 1:num_prices, Nprices; replace=false)

        pmat = prices[:, sampled_prices]

    end

    # need to compute the moments of the sampled prices

    # β = beta_moment_model(pmat, πshares)

    if method == "exact"

        dni = dni_moment_model(pmat, πshares)[2]

        return dni

    elseif method == "over"

        dni, dni2 = dni_moment_model(pmat, πshares)[2:3]

        return dni, dni2

    else
        error("Please specify estimation method, over or exact'.")
    end

end
##############################################################################################################################
##############################################################################################################################

function rescale_trade_cots(θ, gravity_parameters, gravity_results)
    
    @unpack Ncntry = gravity_parameters
    
    d = Array{Float64}(undef, Ncntry, Ncntry) 
    
    foo = gravity_params(θ = θ, gravity_parameters)

    make_trade_costs!(gravity_results, d, foo)

    return d

end
    
##############################################################################################################################
##############################################################################################################################

function get_E_by_importer(df)
    # Get the importer column name (could be "importer" or "iso_d")
    imp_col = "iso_d" in names(df) ? :iso_d : :importer
    
    # Group by importer and take the first E_importer for each
    grp = groupby(df, imp_col)
    return [first(g.E_importer) for g in grp]
end


function generate_simmulated_data(θ, σν, tradedata, trade_parameters, gravity_parameters; model = "ek", code = 1, Nprices = 70, Ngoods = 100000)
    # A function to simmulate a pattern of trade and then generate a
    # random sample of final goods prices, then compute the moments

    @unpack Ncntry = gravity_parameters

    grvity_results = gravity(tradedata, σν; code = code, trade_cost_type = "ek", display = false)

    d = zeros(Ncntry,Ncntry)
    T = zeros(Ncntry)
    W = ones(Ncntry)

    foo_gravity_parameters = gravity_params(θ = θ, gravity_parameters)

    make_trade_costs!(grvity_results, d, foo_gravity_parameters)

    make_technology!(grvity_results, T, W, foo_gravity_parameters)

    foo_trade_parameters= trade_params(θ = θ, σ = trade_parameters.σ, d = d, S = exp.(grvity_results.S), Ncntry = foo_gravity_parameters.Ncntry, N = foo_gravity_parameters.L, E = trade_parameters.E)

    πshares = Array{Float64}(undef, length(Ncntry), length(Ncntry))

    prices = Array{Float64}(undef, length(Ncntry), Ngoods)

    if model == "ek"
        # println("this is the EK model")
        πshares, prices = sim_trade_pattern_ek(foo_trade_parameters; Ngoods = Ngoods, code = code)

    elseif model == "bejk"
        # print("this is the BEJK model")
        println("BEJK model with σ = ", foo_trade_parameters.σ)
        πshares, prices = sim_trade_pattern_bejk(foo_trade_parameters; Ngoods = Ngoods, code = code)

    elseif model == "krugman"
        # print("this is the Krugman model")
        πshares, prices = sim_trade_pattern_krugman(foo_trade_parameters; Ngoods = Ngoods, code = code)

    elseif model == "krugman-model2"
        # print("this is the Krugman model with variable markups")
        πshares, prices = sim_trade_pattern_krugman_model2(foo_trade_parameters; Ngoods = Ngoods, code = code)

    elseif model == "melitz"
        # print("this is the Melitz model")
        prices = Array{Float64}(undef, length(Ncntry), Ngoods*Ncntry)

        common_set = falses(Ncntry * Ngoods)

        πshares, prices, common_set = sim_trade_pattern_melitz_optimized(foo_trade_parameters; Ngoods = Ngoods, code = code)

        num_prices = size(prices[:,common_set],2)

    elseif model == "melitz-model2"

        # print("this is the Krugman model with variable markups")
        πshares, prices = sim_trade_pattern_melitz_model2(foo_trade_parameters; Ngoods = Ngoods, code = code)

        num_prices = size(prices, 2)

    else
        error("Model not recognized. Use 'ek' or 'bejk' or 'krugman'.")
    end

    # print(size(prices))

    if model != "melitz"
        
        sampled_prices= sample(MersenneTwister(09212013 + code), 1:Ngoods, Nprices; replace=false)

        pmat = prices[:, sampled_prices]

    elseif model == "melitz"

        prices = prices[:,common_set]

        sampled_prices= sample(MersenneTwister(09212013 + code), 1:num_prices, Nprices; replace=false)

        pmat = prices[:, sampled_prices]

    elseif model == "melitz-model2"

        sampled_prices= sample(MersenneTwister(09212013 + code), 1:num_prices, Nprices; replace=false)

        pmat = prices[:, sampled_prices]

    end

    dni, dni2 = dni_moment_model(pmat, πshares)[2:3]

    sim_df = DataFrame(dni = dni, dni2 = dni2, dist = gravity_parameters.dfcntryfix.dist)

    return sim_df, foo_trade_parameters, grvity_results


end

##############################################################################################################################
##############################################################################################################################

function boot_strap_simulation(θ, σν, tradedata, trade_parameters, gravity_parameters; 
    model = "ek", method = "exact", Wmat = "optimal", code = 1, Nprices = 70, Ngoods = 100000, Nboots = 100, Nruns  = 10)
    # A function to run a boot strap simulation to get the distribution of the θ estimator
    # the function returns the θ estimates and the value of the J-Stats
    
    lb = [2.5,]
    ub = [10.0,]

    p = SciMLBase.NullParameters()

    θval = Array{Float64}(undef, Nboots)
    Jval = Array{Float64}(undef, Nboots)

    for xxx = 1:Nboots

        sim_df, foo_trade_parameters, gravity_results = generate_simmulated_data(θ, σν, tradedata, trade_parameters, gravity_parameters; 
                model = model, code = code + xxx, Nprices = Nprices, Ngoods = Ngoods)

        f(x, p) = estimate_θ_dni(x[1], sim_df.dni, gravity_parameters, foo_trade_parameters,
            gravity_results; model = model, display = false, Ngoods = Ngoods, Nruns = Nruns, Nprices = Nprices )

        g(x, p) = estimate_θ_dni(x[1], sim_df.dni, sim_df.dni2, sim_df.dist,  gravity_parameters, foo_trade_parameters,
            gravity_results; model = model, Wmat = Wmat, display = false, Ngoods = Ngoods, Nruns = Nruns, Nprices = Nprices )

        if method == "exact"

            # Define the function to estimate θ using the exact method

            prob = OptimizationProblem(f, [4.5], p, lb = lb, ub = ub)

        else
                # Define the function to estimate θ using the over method
            prob = OptimizationProblem(g, [4.5], p, lb = lb, ub = ub)
        
        end

        @time sol = Optimization.solve(prob, BOBYQA(); rhoend = 1e-4)

        # Compatibility fix for different versions of Optimization.jl
        θval[xxx] = hasfield(typeof(sol), :u) ? sol.u[1] : sol.x[1]
        Jval[xxx] = length(sim_df.dni) * sol.objective

        println("Run: ", xxx, " θ: ", θval[xxx], " J-stat: ", Jval[xxx])

    end

    return θval, Jval

end