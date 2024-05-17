module OrrSommerfeld

    using LinearAlgebra
    using UnPack

    export Cheb
    export compute_OS_matrices

    
    include("chebyshev.jl")
    include("baseflows.jl")
    include("osmatrices.jl")

end