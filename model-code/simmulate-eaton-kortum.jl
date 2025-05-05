using Random
using StatsBase
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

    # println("trade", mean( log.(Xni[.!notzeros]) ))

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



function generate_moments(trade_parameters, Nruns; model = "ek", code = 1, Nprices = 70, Ngoods = 100000)
    # multiple dispatch version of the generate_moments function to generate a bunch of betas

    β = Array{Float64}(undef, Nruns)

    Threads.@threads for xxx = 1:Nruns

        β[xxx] = generate_moments(trade_parameters; model = model, code = code + xxx, Nprices = Nprices, Ngoods = Ngoods)[2]

    end

    return β

end
    


function generate_moments(trade_parameters; model = "ek", code = 1, Nprices = 70, Ngoods = 100000)
    # A function to simmulate a pattern of trade and then generate a 
    # random sample of final goods prices, then compute the moments

    @unpack Ncntry = trade_parameters

    πshares = Array{Float64}(undef, length(Ncntry), length(Ncntry))

    prices = Array{Float64}(undef, length(Ncntry), Ngoods)

    if model == "ek"
        # this is the EK model
        πshares, prices = sim_trade_pattern_ek(trade_parameters; Ngoods = Ngoods, code = code)

    elseif model == "bejk"
        # this is the BEJK model
        πshares, prices = sim_trade_pattern_bejk(trade_parameters.S, trade_parameters.d, trade_parameters.θ, 
            trade_parameters.σ; Ngoods = Ngoods, code = code)
    end
        
    sampled_prices= sample(MersenneTwister(09212013 + code), 1:Ngoods, Nprices; replace=false)

    pmat = prices[:, sampled_prices]

    # need to compute the moments of the sampled prices

    β = beta_moment_model(pmat, πshares)

    return pmat, β

end



function rescale_trade_cots(θ, gravity_parameters, gravity_results)
    
    @unpack Ncntry = gravity_parameters
    
    d = Array{Float64}(undef, Ncntry, Ncntry) 
    
    foo = gravity_params(θ = θ, gravity_parameters)

    make_trade_costs!(gravity_results, d, foo)

    return d

end
    




###############################################################
###############################################################

function average_trade_pattern(S, d, θ, σ; Ngoods = 100000, Nruns = 30)
    # Computes the trade shares when avergaged over diffrent simmulaitons
    # of the trade pattern. The function returns the average trade shares.

    πshares = zeros(size(d))    

    Threads.@threads for xxx = 1:Nruns

        πshares = πshares .+ (1.0 / Nruns)*sim_trade_pattern_ek(S, d, θ, σ, Ngoods = Ngoods, code = xxx)[1]

    end

    return πshares

end

function average_trade_pattern(S, d, τ, θ, σ; Ngoods = 50000, Nruns = 5)
    # Computes the trade shares when avergaged over diffrent simmulaitons
    # of the trade pattern. The function returns the average trade shares.

    πshares = Array{Float64}(undef, length(S), length(S), Nruns)
    τ_revenue = Array{Float64}(undef, length(S), length(S), Nruns)
    Φ = Array{Float64}(undef, length(S), Nruns)

    avg_πshares = zeros(size(d))
    avg_τ_revenue = zeros(size(S))
    avg_Φ = zeros(size(S))

    Threads.@threads for xxx = 1:Nruns

       πshares[:,:, xxx], Φ[:,xxx], τ_revenue[:, :, xxx] = sim_trade_pattern_ek(S, d, θ, σ, Ngoods = Ngoods, code = xxx)

    end

    for xxx = 1:Nruns

        avg_πshares = avg_πshares .+ (1.0 / Nruns)*πshares[:,:, xxx]

        avg_τ_revenue = avg_τ_revenue .+ (1.0 / Nruns).*τ_revenue[:,:,xxx]

        avg_Φ = avg_Φ .+ (1.0 / Nruns).*Φ[:,xxx]

    end

    return avg_πshares, avg_τ_revenue, avg_Φ

end

###############################################################
###############################################################

function sim_trade_pattern_ek(trade_parameters; Ngoods = 200000, code = 1)
    # multiple dispatch version of the sim_trade_pattern_ek function
    # this allows me to pass the trade_parameters structure and it will work

    return sim_trade_pattern_ek(trade_parameters.S, trade_parameters.d,  trade_parameters.θ, 
        trade_parameters.σ, Ngoods = Ngoods, code = code)

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

    inv_Ngoods = 1.0 / Ngoods

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
            m[im, min_ex] += low_price^(1.0 - σ)  

            # Update sum price and record lowest price
            sum_price[im] += low_price^(1.0 - σ) 

            rec_low_price[im, gd] = low_price
        end

    end

    # Loop to calculate aggregate price index and the trade shares.

    for im in 1:Ncntry

        g_val = (sum_price[im] * inv_Ngoods)

        for ex in 1:Ncntry

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

    return ( log(u) / (-S) + (1 / first_price)^ (-θ) )^ (one(θ) / θ)

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

function sim_trade_pattern_bejk(S, d, θ, σ; Ngoods = 100000, code = 1)
    # A function to simmulate a pattern of trade and then generate a trade
    # share matrix and a random sample of final goods prices given bertrand
    # pricing from BEJK (2003). 

    Ncntry = length(S)

    inv_Ngoods = 1 / Ngoods

    markup = σ / (σ - 1)

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

    rec_cost = Array{Float64}(undef, Ncntry, Ngoods)

    @inbounds for gd in 1:Ngoods # This is the good

        @inbounds for im in 1:Ncntry # This is the country importing the good

            # In BEJK there are two seperate issues:
            # 1. Needs to find who is the low cost producer AND second lowest cost producer.
            # 2. Need to find the price charged by the low cost producer which
            # is either the 2nd domestic low cost producer, the 2nd foreign low cost producer or the monopolist price.

            low_cost = p1const[im, gd]

            low2_cost = p2const[im, gd]

            min_ex = im

            for ex in 1:Ncntry

                cif_price = d[im, ex] * p1const[ex, gd]

                if cif_price < low_cost # if the price is lower than the current low price

                    low2_cost = low_cost # low_cost is now second lowest

                    low_cost = cif_price # cif_price is the low price
                    
                    min_ex = ex # and the exporter is the one with the lowest price

                else
                    # if not one scnario to check is that the cif_price may be 
                    # the second lowest price

                    low2_cost = min(low2_cost, cif_price)

                end

            end

            price_charged = min(min(d[im, min_ex] * p2const[min_ex, gd], low2_cost), markup * low_cost)
           
            m[im, min_ex] += price_charged^(1 - σ)

            sum_price[im] += price_charged^(1 - σ)

            rec_low_price[im, gd] = price_charged

            rec_cost[im, gd] = low_cost

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

###############################################################
# this is the version with good-specific tariffs using muliple dispatch so when I call the function
# I can just pass the tariff matrix and it will work

# function sim_trade_pattern_ek(S, d, τ, θ, σ; Ngoods = 100000, code = 1)
#     # Constructs pattern of trade for the perfectly competitive model with Frechet
#     # distributed productivity shocks. The function returns the trade shares and the
#     # lowest price of each good in each country.
#     #
#     # S are parameters from gravity equation and are sufficient to simmulate marginal costs
#     # d is the trade costs matrix with rows being imports, columns being exports
#     # θ is the Frechet shape parameter
#     # σ is the elasticity of substitution
#     # options include number of goods and a code for the random number generator
    
#     Ncntry = length(S)

#     inv_Ngoods = 1.0 / Ngoods

#     ###############################################################
#     # Draw productivities and and compute unit costs to produce each good in
#     # each country
    
#     p = Array{Float64}(undef, Ncntry, Ngoods)

#     u = Array{Float64}(undef, Ncntry, Ngoods)

#     rand!(MersenneTwister(03281978 + code ), u)

#     #println(code)

#     @inbounds @views Threads.@threads for j in 1:Ncntry

#         p[j, :] .= marginal_cost.(u[j,:], S[j], θ) 
#         #not sure that this helped... may need to return back to the original code

#     end

#     ###############################################################

#     # Loop to calculate the low price and country suppliers
#     m = zeros(Ncntry, Ncntry) # need to be zero as I'm summing over stuff

#     τ_rev = zeros(Ncntry, Ncntry)

#     Φ = zeros(Ncntry)

#     sum_price = zeros(Ncntry)

#     rec_low_price = Array{Float64}(undef, Ncntry, Ngoods)

#     @inbounds for gd in 1:Ngoods  # Loop over goods # threading here messes stuff up  

#         @inbounds for im in 1:Ncntry  # Loop over importing countries

#             low_price = p[im, gd]
#             min_ex = im

#             @inbounds for ex in 1:Ncntry

#                 cif_price = (1.0 + τ[im, ex] ) * d[im, ex] * p[ex, gd] # price of exporter

#                 # low_price, min_ex = ifelse(cif_price < low_price, (cif_price, ex), (low_price, min_ex)) 

#                 if cif_price < low_price # if the price is lower than the current low price

#                     low_price = cif_price # it is the low price
                    
#                     min_ex = ex # and the exporter is the one with the lowest price
#                 end

#             end

#             ###############################################################
#             # This is an alternative way to find the low cost exporter

#             # cif_price = d[im, :] .* p[:, gd]

#             # sorted_price = sort(cif_price)

#             # low_price = sorted_price[1]

#             # min_ex = findfirst(==(low_price), cif_price)

#             # # This ==(low_price) creates an anonymous function that checks if an element is equal to low_price.

#             # ###############################################################

#             # Update trade matrix `m`
#             m[im, min_ex] += low_price^(1.0 - σ)
            
#             τ_rev[im, min_ex] += τ[im, min_ex] * ( low_price^(1.0 - σ) / (1.0 + τ[im , min_ex]) )
#             # this then computes tariff revenue with two caveats...
#             # need to divide through by the price index in the importing country, just like the m matrix
#             # this is on a per-unit basis, so need to multiply by total demand to get total tariff revenue

#             # Update sum price and record lowest price
#             sum_price[im] += low_price^(1.0 - σ) 

#             rec_low_price[im, gd] = low_price
#         end

#     end

#     # Loop to calculate aggregate price index and the trade shares.

#     for im in 1:Ncntry

#         Φ[im] = (sum_price[im] * inv_Ngoods)

#         for ex in 1:Ncntry

#             m[im, ex] = inv_Ngoods*( m[im, ex] ) / Φ[im]

#             τ_rev[im, ex] = inv_Ngoods*( τ_rev[im, ex] ) / Φ[im]

#         end

#     end

#     P = Φ.^(1.0 / (1.0 - σ)) # this is the CES price index

#     return m, P, τ_rev

# end



