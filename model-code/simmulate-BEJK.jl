using Random

###############################################################
###############################################################

function marginal_cost_second(u, first_price, S, θ)
    # takes random number u, productivity S and frechet shape parameters
    # θ and returns the marginal cost of producing a good

    return ( log(u) / (-S) + (1 / first_price)^ (-θ) )^ (one(θ) / θ)

    # Second draw, second best productivity, this comes from
    # 1-exp(-S*z_two^(-θ) + S*z_one^(-θ)) 

end


function sim_trade_pattern_BEJK(S, d, θ, σ; Ngoods = 100000, code = 1)
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

    for j in 1:Ncntry

        p1const[j, :] .= marginal_cost.(u[j,:], S[j], θ) 
        # Invert to convert to unit costs. Here assume w[j] = 1 ?

    end

    rand!(MersenneTwister(02071983 + code ), u)
    
    for j in 1:Ncntry

        p2const[j, :] .= marginal_cost_second.(u[j,:], p1const[j, :], S[j], θ)

    end

    ###########################################################
    # Loop to calculate the low price and country suppliers

    m = zeros(Ncntry, Ncntry)

    sum_price = zeros(Ncntry)

    rec_low_price = Array{Float64}(undef, Ncntry, Ngoods)

    rec_cost = Array{Float64}(undef, Ncntry, Ngoods)

    issecond = Array{Float64}(undef, Ncntry, Ngoods)

    istraded = Array{Float64}(undef, Ncntry, Ngoods)

    for gd in 1:Ngoods # This is the good

        for im in 1:Ncntry # This is the country importing the good

            # In BEJK there are two seperate issues:
            # First, one needs to findwho is the low cost producer and second lowest cost producer.
            # Second, one needs to find the price charged by the low cost producer which
            # is either the 2nd domestic low cost producer, the 2nd foreign low cost producer or the monopolist price.

            # Here try using 'sort' function directly. Need to check whether this is 
            # too slow and implement Mike's strategy instead.

            cif_price = d[im, :] .* p1const[:, gd]

            sorted_price = sort(cif_price)

            low_cost, low2_cost = sorted_price[1], sorted_price[2]

            price_charged = 0.0

            # First check if home guy is the low cost
            if low_cost == cif_price[im]

                istraded[im, gd] = 0

                price_charged = min(min(p2const[im, gd], low2_cost), markup * low_cost)
                # There is no τ on p2const because this is home/home situation

                m[im, im] += price_charged ^ (1 - σ)

            else

                istraded[im, gd] = 1

                ex = findfirst(==(low_cost), cif_price)

                price_charged = min(min(d[im, ex] * p2const[ex, gd], low2_cost), markup * low_cost)

                m[im, ex] += price_charged ^ (1 - σ)

            end

            sum_price[im] += price_charged ^ (1 - σ)

            rec_low_price[im, gd] = price_charged

            rec_cost[im, gd] = low_cost

            issecond[im, gd] = price_charged == low2_cost ? 1 : 0
            # where the ? is a shorthand for an if-else statement
            # so if price_charged == low2_cost then 1 else 0


        end

    end

    ###########################################################
    # Loop to calculate aggregate price index and the trade shares

    for im in 1:Ncntry

        g = (sum_price[im] * inv_Ngoods) 

        for ex in 1:Ncntry
        
            m[im, ex] = inv_Ngoods * m[im, ex] / g

        end

    end

    return m, rec_low_price, rec_cost, issecond, istraded

end


###############################################################

function sim_trade_pattern_BEJK(S, d, θ, σ; Ngoods = 100000, code = 1)
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

    return m, rec_low_price, rec_cost

end