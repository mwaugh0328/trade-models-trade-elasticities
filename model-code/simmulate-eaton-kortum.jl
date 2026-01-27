using Random
using StatsBase
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
            gravity_results; model = "ek", method = "exact", Wmat = "optimal", display = true, Nruns = 10, Ngoods = 100000)

    lb = [2.5,]
    ub = [10.0,]

    p = SciMLBase.NullParameters()

    f(x, p) = estimate_θ_dni(x[1], dfdni.dni, gravity_parameters, trade_parameters,
            gravity_results; model = model, display = display, return_moments = false, Nruns = Nruns, Ngoods = Ngoods)

    g(x, p) = estimate_θ_dni(x[1], dfdni.dni, dfdni.dni2, dfdni.dist,  gravity_parameters, trade_parameters,
            gravity_results; model = model, Wmat = Wmat, display = display, return_moments = false, Nruns = Nruns, Ngoods = Ngoods)

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
            gravity_results; model = model, display = display, return_moments = true, Nruns = Nruns, Ngoods = Ngoods)

    else

        ~, data_moments, model_moments = estimate_θ_dni(solution_value, dfdni.dni, dfdni.dni2, dfdni.dist,  gravity_parameters, trade_parameters,
            gravity_results; model = model, Wmat = Wmat, display = display, return_moments = true, Nruns = Nruns, Ngoods = Ngoods)

            

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

    elseif model == "melitz"
        # print("this is the Melitz model")
        prices = Array{Float64}(undef, length(Ncntry), Ngoods*Ncntry)

        common_set = falses(Ncntry * Ngoods)

        πshares, prices, common_set = sim_trade_pattern_melitz_optimized(trade_parameters; Ngoods = Ngoods, code = code)

        num_prices = size(prices[:,common_set],2)

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

    foo_trade_parameters= trade_params(θ = θ, σ = trade_parameters.σ, d = d, S = exp.(grvity_results.S), Ncntry = foo_gravity_parameters.Ncntry, N = foo_gravity_parameters.L)

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
    elseif model == "melitz"
        # print("this is the Melitz model")
        prices = Array{Float64}(undef, length(Ncntry), Ngoods*Ncntry)

        common_set = falses(Ncntry * Ngoods)

        πshares, prices, common_set = sim_trade_pattern_melitz_optimized(foo_trade_parameters; Ngoods = Ngoods, code = code)

        num_prices = size(prices[:,common_set],2)

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
            gravity_results; model = model, display = false, Ngoods = Ngoods, Nruns = Nruns )

        g(x, p) = estimate_θ_dni(x[1], sim_df.dni, sim_df.dni2, sim_df.dist,  gravity_parameters, foo_trade_parameters,
            gravity_results; model = model, Wmat = Wmat, display = false, Ngoods = Ngoods, Nruns = Nruns )

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


###############################################################
###############################################################

# function average_trade_pattern(S, d, θ, σ; Ngoods = 100000, Nruns = 30)
#     # Computes the trade shares when avergaged over diffrent simmulaitons
#     # of the trade pattern. The function returns the average trade shares.

#     πshares = zeros(size(d))    

#     Threads.@threads for xxx = 1:Nruns

#         πshares = πshares .+ (1.0 / Nruns)*sim_trade_pattern_ek(S, d, θ, σ, Ngoods = Ngoods, code = xxx)[1]

#     end

#     return πshares

# end

# function average_trade_pattern(S, d, τ, θ, σ; Ngoods = 50000, Nruns = 5)
#     # Computes the trade shares when avergaged over diffrent simmulaitons
#     # of the trade pattern. The function returns the average trade shares.

#     πshares = Array{Float64}(undef, length(S), length(S), Nruns)
#     τ_revenue = Array{Float64}(undef, length(S), length(S), Nruns)
#     Φ = Array{Float64}(undef, length(S), Nruns)

#     avg_πshares = zeros(size(d))
#     avg_τ_revenue = zeros(size(S))
#     avg_Φ = zeros(size(S))

#     Threads.@threads for xxx = 1:Nruns

#        πshares[:,:, xxx], Φ[:,xxx], τ_revenue[:, :, xxx] = sim_trade_pattern_ek(S, d, θ, σ, Ngoods = Ngoods, code = xxx)

#     end

#     for xxx = 1:Nruns

#         avg_πshares = avg_πshares .+ (1.0 / Nruns)*πshares[:,:, xxx]

#         avg_τ_revenue = avg_τ_revenue .+ (1.0 / Nruns).*τ_revenue[:,:,xxx]

#         avg_Φ = avg_Φ .+ (1.0 / Nruns).*Φ[:,xxx]

#     end

#     return avg_πshares, avg_τ_revenue, avg_Φ

# end

###############################################################
###############################################################

function sim_trade_pattern_ek(trade_parameters; Ngoods = 100000, code = 1)
    # multiple dispatch version of the sim_trade_pattern_ek function
    # this allows me to pass the trade_parameters structure and it will work

    return sim_trade_pattern_ek(trade_parameters.S, trade_parameters.d,  trade_parameters.θ, 
        trade_parameters.σ; Ngoods = Ngoods, code = code)

end


function sim_trade_pattern_ek(S, d,  θ, σ; Ngoods = 100000, code = 1)
    # Constructs pattern of trade for the perfectly competitive model with Frechet
    # distributed productivity shocks. The function returns the trade shares and the
    # lowest price of each good in each country.
    #
    # S are parameters from gravity equation and are sufficient to simmulate marginal costs
    # d is the trade costs matrix with rows being imports, columns being exports
    # θ is the Frechet shape parameter
    # σ is the elasticity of substitution
    # options include number of goods and a code for the random number generator
    
    Ncntry = length(S)

    inv_Ngoods = one(σ) / Ngoods

    ###############################################################
    # Draw productivities and and compute unit costs to produce each good in
    # each country
    
    p = Array{Float64}(undef, Ncntry, Ngoods)

    u = Array{Float64}(undef, Ncntry, Ngoods)

    rand!(MersenneTwister(03281978 + code ), u)

    #println(code)

    @inbounds @views Threads.@threads for j in 1:Ncntry

        p[j, :] .= marginal_cost.(u[j,:], S[j], θ) 
        #not sure that this helped... may need to return back to the original code

    end

    ###############################################################

    # Loop to calculate the low price and country suppliers
    m = zeros(Ncntry, Ncntry) # need to be zero as I'm summing over stuff
    sum_price = zeros(Ncntry)

    rec_low_price = Array{Float64}(undef, Ncntry, Ngoods)

    @inbounds for gd in 1:Ngoods  # Loop over goods # threading here messes stuff up  

        @inbounds for im in 1:Ncntry  # Loop over importing countries

            low_price = p[im, gd]
            min_ex = im

            @inbounds for ex in 1:Ncntry

                cif_price = d[im, ex] * p[ex, gd] # price of exporter

                # low_price, min_ex = ifelse(cif_price < low_price, (cif_price, ex), (low_price, min_ex)) 

                if cif_price < low_price # if the price is lower than the current low price

                    low_price = cif_price # it is the low price
                    
                    min_ex = ex # and the exporter is the one with the lowest price
                end

            end

            ###############################################################
            # This is an alternative way to find the low cost exporter

            # cif_price = d[im, :] .* p[:, gd]

            # sorted_price = sort(cif_price)

            # low_price = sorted_price[1]

            # min_ex = findfirst(==(low_price), cif_price)

            # # This ==(low_price) creates an anonymous function that checks if an element is equal to low_price.

            # ###############################################################

            # Update trade matrix `m`

            # lp = low_price^(one(σ) - σ)

            m[im, min_ex] += low_price^(one(σ) - σ)
            # Update sum price and record lowest price

            sum_price[im] += low_price^(one(σ) - σ)

            rec_low_price[im, gd] = low_price
        end

    end

    # Loop to calculate aggregate price index and the trade shares.

    @inbounds for im in 1:Ncntry

        g_val = (sum_price[im] * inv_Ngoods)

        @inbounds for ex in 1:Ncntry

            m[im, ex] = inv_Ngoods*( m[im, ex] ) / g_val

        end

    end

    return m, rec_low_price

end

###############################################################
###############################################################

function marginal_cost(u, S, θ)
    # takes random number u, productivity S and frechet shape parameters
    # θ and returns the marginal cost of producing a good

    return ( log(u) / (-S) )^ ( one(θ) / θ )

    # (log.(u) ./ (-λ[j])) .^ (1/θ) 

end

function marginal_cost_second(u, first_price, S, θ)
    # takes random number u, productivity S and frechet shape parameters
    # θ and returns the marginal cost of producing a good

    return ( log(u) / (-S) + (first_price)^ (θ) )^ (one(θ) / θ)

    # (log(u2)./(-S(j)) + (p1const(:,j).^-1).^(-theta)).^(-1./theta);

    # Second draw, second best productivity, this comes from
    # 1-exp(-S*z_two^(-θ) + S*z_one^(-θ)) 

end

###############################################################
###############################################################

function sim_trade_pattern_bejk(trade_parameters; Ngoods = 100000, code = 1)
    # multiple dispatch versin of the sim_trade_pattern_ek function
    # this allows me to pass the trade_parameters structure and it will work

    return sim_trade_pattern_bejk(trade_parameters.S, trade_parameters.d,  trade_parameters.θ, 
        trade_parameters.σ; Ngoods = Ngoods, code = code)

end

###############################################################
###############################################################

function sim_trade_pattern_bejk(S, d, θ, σ; Ngoods = 100000, code = 1)
    # A function to simmulate a pattern of trade and then generate a trade
    # share matrix and a random sample of final goods prices given bertrand
    # pricing from BEJK (2003). 

    Ncntry = length(S)

    inv_Ngoods = 1 / Ngoods

    markup = σ / (σ - one(σ))

    

    ###########################################################
    # Draw productivities and compute unit costs to produce each good in 
    # each country

    p1const = Array{Float64}(undef, Ncntry, Ngoods)

    p2const = Array{Float64}(undef, Ncntry, Ngoods)

    u = Array{Float64}(undef, Ncntry, Ngoods)

    rand!(MersenneTwister(03281978 + code ), u)

    @inbounds @views Threads.@threads for j in 1:Ncntry

        p1const[j, :] .= marginal_cost.(u[j,:], S[j], θ) 
        # Invert to convert to unit costs. Here assume w[j] = 1 ?

    end

    rand!(MersenneTwister(02071983 + code ), u)
    
    @inbounds @views Threads.@threads for j in 1:Ncntry

        p2const[j, :] .= marginal_cost_second.(u[j,:], p1const[j, :], S[j], θ)

    end

    ###########################################################
    # Loop to calculate the low price and country suppliers

    m = zeros(Ncntry, Ncntry)

    sum_price = zeros(Ncntry)

    rec_low_price = Array{Float64}(undef, Ncntry, Ngoods)

    # rec_cost = Array{Float64}(undef, Ncntry, Ngoods)

    @inbounds for gd in 1:Ngoods # This is the good

        @inbounds for im in 1:Ncntry # This is the country importing the good

            # In BEJK there are two seperate issues:
            # 1. Needs to find who is the low cost producer AND second lowest cost producer.
            # 2. Need to find the price charged by the low cost producer which
            # is either the 2nd domestic low cost producer, the 2nd foreign low cost producer or the monopolist price.

            # this is how the old matlab code worked, the first julia BEJK version put together earlier was not correct (unsure why)
            # The first step is to find the two lowest international cost producers
        
            cif_price1 = d[im, 1] * p1const[1,gd]
            
            cif_price2 = d[im, 2] * p1const[2,gd]
           
            if cif_price1 < cif_price2 # if the first country is lower, its low cost

               low_cost = cif_price1

               low2_cost = cif_price2 # second country is second lowest
                
               min_ex = 1
            else # if the second country is lower, its low cost 1 is second lowest
               
                low_cost = cif_price2
               
                low2_cost = cif_price1 # first country is second lowest
               
                min_ex = 2
            end


            for ex in 3:Ncntry
                # now walk through remaining potential exporters

                cif_price = d[im, ex] * p1const[ex, gd]
                # this is the price of the exporter

                low2_cost = min(max(cif_price, low_cost), low2_cost)

                if cif_price < low_cost # if the exporter price is lower than the current low price

                    low_cost = cif_price # the low costs is that exporters price
                    
                    min_ex = ex # and the exporter is the one with the lowest price

                end

            end

            price_charged = min(min(d[im, min_ex] * p2const[min_ex, gd], low2_cost), markup * low_cost)

            m[im, min_ex] += (price_charged )^(one(σ) - σ)

            sum_price[im] += (price_charged )^(one(σ) - σ)

            rec_low_price[im, gd] = price_charged

            # rec_cost[im, gd] = low_cost

        end

    end

    ###########################################################
    # Loop to calculate aggregate price index and the trade shares

    for im in 1:Ncntry

        g_val = (sum_price[im] * inv_Ngoods)

        for ex in 1:Ncntry

            m[im, ex] = inv_Ngoods*( m[im, ex] ) / g_val

        end

    end

    return m, rec_low_price

end

##############################################################################################################################
##############################################################################################################################

function sim_trade_pattern_krugman(trade_parameters; Ngoods = 10000, code = 1)
    # multiple dispatch version of the sim_trade_pattern_krugman function

    return sim_trade_pattern_krugman(trade_parameters.S, trade_parameters.d,  trade_parameters.θ; Ngoods = Ngoods, code = code)

end

function sim_trade_pattern_krugman(S, d, θ; Ngoods = 10, code = 1)

    Ncntry = length(S)

    η = θ + 1
    markup = η / (η - 1)

    m = zeros(Ncntry, Ncntry) # need to be zero as I'm summing over stuff

    sum_price = zeros(Ncntry)

    price_matrix = zeros(Ncntry, Ngoods)

    final_price = zeros(Ncntry, Ncntry * Ngoods)

    u = fill(0.5, Ngoods)

    # Assign values in the price matrix
    @inbounds @views Threads.@threads for j in 1:Ncntry

        price_matrix[j, :] .= markup .* marginal_cost.(u, S[j], θ) 
        
    end

    # print(marginal_cost.(u, S[1], θ) )

    @inbounds for im in 1:Ncntry

        carry_prices = zeros(Ncntry, Ngoods)

        @inbounds for gd in 1:Ngoods

            # Loop over importing countries
            @inbounds for ex in 1:Ncntry

                carry_prices[ex, gd] = d[im, ex] * price_matrix[ex, gd]

                m[im, ex] += carry_prices[ex, gd]^(one(η) - η)

                sum_price[im] += carry_prices[ex, gd]^(one(η) - η)

            end

        end

        @views final_price[im, :] .= vec(carry_prices)
        # basicly this whole vector now reflects all the prices for each variety
    end

     @inbounds for im in 1:Ncntry

        g_val = sum_price[im]

        @inbounds for ex in 1:Ncntry

            m[im, ex] = ( m[im, ex] ) / g_val

        end

    end


    sampled_prices= sample(MersenneTwister(09111943 + code), 1:(Ncntry * Ngoods), Ngoods; replace=false)

    rec_low_price = final_price[:, sampled_prices]

    return m, rec_low_price

end

##############################################################################################################################
##############################################################################################################################


function sim_trade_pattern_melitz(trade_parameters; Ngoods = 10000, code = 1)
    # multiple dispatch version of the sim_trade_pattern_krugman function

    return sim_trade_pattern_melitz(trade_parameters.S, trade_parameters.d,  trade_parameters.θ; Ngoods = Ngoods, code = code)

end

function sim_trade_pattern_melitz(S, d, θ; Ngoods = 10000, code = 1)

    Ncntry = length(S)
    
    η = θ + 1

    markup = η / (η - 1)    

    # rng = MersenneTwister(03281978 + code)

    pi_nn = zeros(Ncntry)

    cost_cutoff = zeros(Ncntry)

    for j in 1:Ncntry
        
        phi_sum = sum(S .* d[j, :].^(-θ))

        pi_nn[j] = S[j] / phi_sum
        
        cost_cutoff[j] = (pi_nn[j] / S[j])^(1 / θ)
    end

    max_goods = maximum(pi_nn)

    m = zeros(Ncntry, Ncntry)

    sum_price = zeros(Ncntry)
    
    price_matrix = zeros(Ncntry, Ngoods)

    # Assign values in the price matrix
    for j in 1:Ncntry

        cgoods = round(Int, Ngoods * (pi_nn[j] / max_goods))

        u = Array{Float64}(undef, cgoods)

        rand!(MersenneTwister(03281978 + code + j), u)

        u = pi_nn[j] .* u

        price_matrix[ j , 1:cgoods] .= markup .* (u ./ S[j]).^(1 / θ)
        # Optionally shuffle: price_matrix[:, jjj] = price_matrix[shuffle(rng, 1:Ngoods), jjj]
    end

    final_price = Array{Float64}(undef, Ncntry, Ncntry * Ngoods)

    common_set = falses(Ncntry * Ngoods)

    # Compute exporting decisions and trade shares

    
   @inbounds for im in 1:Ncntry

        carry_prices = zeros(Ncntry, Ngoods)

        @inbounds for gd in 1:Ngoods

            @inbounds for ex in 1:Ncntry  
            
                if ex == im

                    carry_prices[ex, gd] = price_matrix[ex, gd]

                
                elseif d[im, ex] * price_matrix[ex, gd] <= markup * cost_cutoff[im]
            

                    carry_prices[ex, gd] = d[im, ex] * price_matrix[ex, gd]

                end

                if carry_prices[ex, gd] != 0.0
                   
                    m[im, ex] += carry_prices[ex, gd]^(one(η) - η)

                    sum_price[im] += carry_prices[ex, gd]^(one(η) - η)

                else
                    carry_prices[ex, gd] = NaN

                end

            end
        end

        @views final_price[im, :] .= vec(carry_prices)

    end

    for im in 1:Ncntry

        g_val = sum_price[im]

        for ex in 1:Ncntry

            m[im, ex] = ( m[im, ex] ) / g_val

        end

    end


    # Filter: keep only rows where all countries have a price (no NaN)
    common_set = [all(.!isnan.(final_price[:, gd])) for gd in 1:size(final_price, 2)]

    # final_price = final_price[:, common_set]

    # sampled_prices= sample(MersenneTwister(09111943 + code), 1:(Ncntry * Ngoods), Ngoods; replace=false)

    # rec_low_price = final_price[:, sampled_prices]
    # sample_common_set = common_set[sampled_prices]

    return m, final_price, common_set
end

# Multiple dispatch version for optimized Melitz function
function sim_trade_pattern_melitz_optimized(trade_parameters; Ngoods = 10000, code = 1)
    # multiple dispatch version of the sim_trade_pattern_melitz_optimized function

    return sim_trade_pattern_melitz_optimized(trade_parameters.S, trade_parameters.d,  trade_parameters.θ; Ngoods = Ngoods, code = code)

end

# # Optimized version of sim_trade_pattern_melitz
# function sim_trade_pattern_melitz_optimized(S, d, θ; Ngoods = 10000, code = 1)

#     Ncntry = length(S)
    
#     η = θ + 1
#     markup = η / (η - 1)    
#     inv_θ = 1 / θ
#     one_minus_η = one(η) - η

#     # Vectorized computation of pi_nn and cost_cutoff
#     phi_sum = [sum(S .* d[j, :].^(-θ)) for j in 1:Ncntry]

#     pi_nn = S ./ phi_sum
    
#     cost_cutoff = (pi_nn ./ S).^inv_θ
    
#     markup_cutoff = markup .* cost_cutoff

#     max_goods = maximum(pi_nn)

#     m = zeros(Ncntry, Ncntry)

#     sum_price = zeros(Ncntry)
    
#     # Pre-allocate price matrix
#     price_matrix = zeros(Ncntry, Ngoods)

#     # Vectorized price matrix assignment
#     for j in 1:Ncntry
#         cgoods = round(Int, Ngoods * (pi_nn[j] / max_goods))
#         if cgoods > 0
#             u = rand(MersenneTwister(03281978 + code + j), cgoods)
#             u .*= pi_nn[j]
#             price_matrix[j, 1:cgoods] .= markup .* (u ./ S[j]).^inv_θ
#         end
#     end

#     final_price = Array{Float64}(undef, Ncntry, Ncntry * Ngoods)

#     # Pre-allocate arrays
#     carry_prices = zeros(Ncntry, Ngoods)
    
#     # Optimized main computation loop
#     @inbounds for im in 1:Ncntry
#         fill!(carry_prices, 0.0)  # Reset instead of allocating new array
        
#         # Vectorized operations where possible
#         for gd in 1:Ngoods
#             # Home country always included
#             carry_prices[im, gd] = price_matrix[im, gd]
            
#             # Check export conditions for other countries
#             for ex in 1:Ncntry
#                 if ex != im
#                     cif_price = d[im, ex] * price_matrix[ex, gd]
#                     if cif_price <= markup_cutoff[im] && price_matrix[ex, gd] > 0
#                         carry_prices[ex, gd] = cif_price
#                     end
#                 end
                
#                 # Compute trade shares and price aggregates
#                 price = carry_prices[ex, gd]
#                 if price > 0
#                     price_power = price^one_minus_η
#                     m[im, ex] += price_power
#                     sum_price[im] += price_power
#                 else
#                     carry_prices[ex, gd] = NaN
#                 end
#             end
#         end

#         # Use view instead of copying
#         @views final_price[im, :] .= vec(carry_prices)
#     end

#     # Normalize trade shares efficiently
#     @inbounds for im in 1:Ncntry
#         g_val = sum_price[im]
#         if g_val > 0
#             for ex in 1:Ncntry
#                 m[im, ex] /= g_val
#             end
#         end
#     end

#     # Vectorized common_set computation
#     common_set = [all(!isnan, @view final_price[:, gd]) for gd in 1:(Ncntry * Ngoods)]

#     return m, final_price, common_set
# end

###############################################################################################
# NOTES on the MELITZ implementation in the function sim_trade_pattern_melitz_optimized

# Simulates bilateral trade patterns and micro-level prices from the Melitz model
# following the Simonovska-Waugh (SW) parameterization.

# # Model Setup (SW / Chaney 2008 Assumptions)

# This implementation assumes fixed export costs are proportional to origin country size:
#     fᵢⱼ = f ⋅ Lᵢ

# Under this assumption, the fraction of country i's firms selling in country j simplifies to:
#     Nᵢⱼ/Nᵢ = (dᵢⱼ⁻θ ⋅ Sⱼ) / Φᵢ

# where Φᵢ = Σₖ dᵢₖ⁻θ ⋅ Sₖ is the "market access" term. This convenient form arises
# because the Lᵢ terms cancel in the ratio of cutoffs.

# Key Theoretical Points

# 1. **Domestic share**: πⱼⱼ = Sⱼ/Φⱼ gives the fraction of varieties in j produced domestically.
#    This differs from EK where all goods are potentially produced everywhere.

# 2. **Export selection**: A firm from country `ex` with domestic price p can export to `im` if:
#        p ⋅ d[im,ex] ≤ p̄ᵢₘ
#    where p̄ᵢₘ is the destination's domestic cutoff price. The proportional fixed cost
#    assumption makes the destination's domestic cutoff the relevant threshold.

# 3. **Price distribution**: Firm productivity φ ~ Pareto(φ*, θ) conditional on entry.
#    Under SW normalization, φ* = Φʲ^(1/θ), so drawing u ~ U(0,1) and computing
#        p = markup ⋅ (u ⋅ πⱼⱼ / Sⱼ)^(1/θ)
#    yields the correct truncated Pareto distribution of prices.

# 4. **Common set**: Only the most productive firms clear all export cutoffs and sell in
#    every market. This set is smaller than in EK, which is why θ_Melitz < θ_EK for
#    the same data — selection compresses observed price gaps.

# # Arguments
# - `S::Vector`: Technology parameters (analogous to Tᵢ in EK) for each country
# - `d::Matrix`: Bilateral trade costs (iceberg), d[i,j] = cost for j to sell in i
# - `θ::Real`: Pareto shape parameter (governs productivity dispersion)
# - `Ngoods::Int=10000`: Number of potential goods to simulate per country
# - `code::Int=1`: Seed modifier for reproducibility across different runs

# # Returns
# - `m::Matrix`: Bilateral trade share matrix, m[i,j] = share of i's expenditure on j
# - `final_price::Matrix`: Price matrix, final_price[i, :] = prices faced by importer i
# - `common_set::Vector{Bool}`: Indicator for goods sold in all countries (for SW estimation)

# # References
# - Simonovska, I. and Waugh, M.E. "Trade Models, Trade Elasticities, and the Gains from Trade"
# - Chaney, T. (2008) "Distorted Gravity: The Intensive and Extensive Margins of International Trade"
# - Melitz, M.J. (2003) "The Impact of Trade on Intra-Industry Reallocations and Aggregate Industry Productivity"

function sim_trade_pattern_melitz_optimized(S, d, θ; Ngoods = 10000, code = 1)
    
    Ncntry = length(S)
   
    # =========================================================================
    # Model Parameters
    # =========================================================================
    # η = θ + 1 is the elasticity of substitution across varieties (CES parameter)
    # In Melitz-Chaney, θ (Pareto shape) and η are linked; this follows SW's notation
    η = θ + 1
    markup = η / (η - 1)        # CES markup: p = (η/(η-1)) ⋅ mc
    inv_θ = 1 / θ
    one_minus_η = one(η) - η    # Exponent for CES price aggregation: p^(1-η)
   
    # =========================================================================
    # Compute Domestic Shares and Cutoffs
    # =========================================================================
    # Φⱼ = Σᵢ Sᵢ ⋅ dⱼᵢ⁻θ is the "multilateral resistance" / market access term
    # Under SW's proportional fixed cost assumption, this fully characterizes
    # the competitive environment in each market
    phi_sum = [sum(S .* d[j, :].^(-θ)) for j in 1:Ncntry]
   
    # πⱼⱼ = Sⱼ/Φⱼ = domestic expenditure share = fraction of varieties from home
    # This is the Melitz analog of the EK domestic trade share
    pi_nn = S ./ phi_sum
   
    # Cost cutoff: firms with marginal cost c ≤ c* produce
    # Under SW normalization: c* = (πⱼⱼ/Sⱼ)^(1/θ) = Φⱼ^(-1/θ)
    # Equivalently, productivity cutoff φ* = Φⱼ^(1/θ)
    cost_cutoff = (pi_nn ./ S).^inv_θ
   
    # Maximum price at which a firm can sell domestically
    # Firms with p > markup ⋅ c* cannot cover fixed costs and exit
    markup_cutoff = markup .* cost_cutoff
   
    # =========================================================================
    # Initialize Storage
    # =========================================================================
    max_goods = maximum(pi_nn)
    m = zeros(Ncntry, Ncntry)           # Trade share matrix
    sum_price = zeros(Ncntry)            # Price index aggregator (Σ p^(1-η))
    price_matrix = zeros(Ncntry, Ngoods) # Domestic prices by origin country
   
    # =========================================================================
    # Generate Domestic Prices via Pareto Productivity Draws
    # =========================================================================
    # For each country j, we simulate the prices of firms that successfully enter.
    #
    # Theory: φ ~ Pareto(φ*, θ) ⟹ drawing u ~ U(0,1), set φ = φ* ⋅ u^(-1/θ)
    # Price: p = markup/φ = markup ⋅ u^(1/θ) / φ*
    #
    # Under SW normalization where φ* = Φⱼ^(1/θ) = (Sⱼ/πⱼⱼ)^(1/θ):
    #   p = markup ⋅ u^(1/θ) ⋅ (πⱼⱼ/Sⱼ)^(1/θ) = markup ⋅ (u ⋅ πⱼⱼ/Sⱼ)^(1/θ)
    #
    # We scale the number of goods by πⱼⱼ/max(π) so larger countries have more
    # varieties, matching the extensive margin predictions of Melitz.
   
    for j in 1:Ncntry
        # Number of goods produced in country j (scaled by relative domestic share)
        cgoods = round(Int, Ngoods * (pi_nn[j] / max_goods))
       
        if cgoods > 0
            # Seeded RNG for reproducibility across runs
            # Note: Creating MersenneTwister per country is intentional for reproducibility
            # but has some performance cost; acceptable for typical Ncntry (< 100)
            u = rand(MersenneTwister(03281978 + code + j), cgoods)
           
            # Scale uniform draws and transform to get Pareto-distributed prices
            u .*= pi_nn[j]
            price_matrix[j, 1:cgoods] .= markup .* (u ./ S[j]).^inv_θ
        end
        # Note: Goods cgoods+1 : Ngoods remain at zero for country j
        # This represents goods that don't exist because no firm in j had
        # sufficient productivity to enter. This is intentional and handled below.
    end
   
    # =========================================================================
    # Storage for Final Prices and Trade Computation
    # =========================================================================
    # final_price[im, :] stores all prices (domestic and imported) faced by importer im
    # Dimensions: Ncntry importers × (Ncntry × Ngoods) potential goods
    final_price = Array{Float64}(undef, Ncntry, Ncntry * Ngoods)
   
    # Temporary matrix for prices available in each destination market
    carry_prices = zeros(Ncntry, Ngoods)
   
    # =========================================================================
    # Main Loop: Compute Trade Shares and Delivered Prices
    # =========================================================================
    @inbounds for im in 1:Ncntry
        fill!(carry_prices, 0.0)  # Reset for each importing country
       
        for gd in 1:Ngoods
            # -----------------------------------------------------------------
            # Home country prices: always "available" if the good exists
            # -----------------------------------------------------------------
            # If price_matrix[im, gd] = 0, no domestic firm produces this good
            # (productivity draw was below cutoff). This is intentional — not
            # all goods exist in all countries under Melitz.
            carry_prices[im, gd] = price_matrix[im, gd]
           
            # -----------------------------------------------------------------
            # Check export conditions for foreign countries
            # -----------------------------------------------------------------
            for ex in 1:Ncntry
                if ex != im
                    # CIF price = FOB price × iceberg trade cost
                    cif_price = d[im, ex] * price_matrix[ex, gd]
                   
                    # Export selection condition (SW/Chaney):
                    # A firm exports from ex to im if its delivered price clears
                    # the destination's cutoff. Under proportional fixed costs,
                    # the relevant threshold is the importer's domestic cutoff.
                    #
                    # Also require price_matrix[ex, gd] > 0 (good must exist in origin)
                    if cif_price <= markup_cutoff[im] && price_matrix[ex, gd] > 0
                        carry_prices[ex, gd] = cif_price
                    end
                end
            end
           
            # -----------------------------------------------------------------
            # Aggregate into trade shares and price index
            # -----------------------------------------------------------------
            for ex in 1:Ncntry
                price = carry_prices[ex, gd]
               
                if price > 0
                    # CES expenditure share: proportional to p^(1-η)
                    price_power = price^one_minus_η
                    m[im, ex] += price_power
                    sum_price[im] += price_power
                else
                    # Mark as NaN: good not available from this origin in this market
                    # Used later to identify common set of goods sold everywhere
                    carry_prices[ex, gd] = NaN
                end
            end
        end
       
        # Store all prices for this importer (flattened across origins × goods)
        @views final_price[im, :] .= vec(carry_prices)
    end
   
    # =========================================================================
    # Normalize Trade Shares
    # =========================================================================
    # Convert expenditure-weighted sums to shares: m[im, ex] / Σₖ m[im, k]
    @inbounds for im in 1:Ncntry
        g_val = sum_price[im]
        if g_val > 0
            for ex in 1:Ncntry
                m[im, ex] /= g_val
            end
        end
    end
   
    # =========================================================================
    # Identify Common Set of Goods
    # =========================================================================
    # The "common set" contains goods sold in ALL countries — these are produced
    # by the most productive firms that clear every export cutoff.
    #
    # Key insight (SW): The common set is smaller in Melitz than EK because
    # selection is more stringent. Only top-productivity firms export everywhere.
    # This is why θ_Melitz < θ_EK when estimated on the same price data —
    # selection compresses the observed distribution of price gaps.
    #
    # For SW estimation, moments are computed on this common set to ensure
    # comparability of prices across countries.
    # For each good, check if it's sold in all countries (no NaN prices)
    common_set = Vector{Bool}(undef, Ncntry * Ngoods)

    for gd in 1:(Ncntry * Ngoods)
    
        prices_this_good = @view final_price[:, gd]
    
        sold_everywhere = all(p -> !isnan(p), prices_this_good)
    
        common_set[gd] = sold_everywhere
    end
   
    return m, final_price, common_set
end