module OrrSommerfeld

    using LinearAlgebra
    using Optim
    using UnPack

    export Cheb
    export os_matrices, os_eigen, optimalgrowth

    
    include("chebyshev.jl")
    include("utilities.jl")
    include("baseflows.jl")
    include("osmatrices.jl")
    include("transientgrowth.jl")
    include("fields.jl")

end