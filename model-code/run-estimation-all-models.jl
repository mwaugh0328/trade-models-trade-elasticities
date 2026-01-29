

println("Running estimation for Eaton-Kortum model...")

include("estimate-ek-model.jl")

println("Running estimation for BEJK model...")

include("estimate-bejk-model.jl")

println("Running estimation for Krugman model...")

include("estimate-krugman-model.jl")

println("Running estimation for Melitz model...")

include("estimate-melitz-model.jl")