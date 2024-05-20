module OrrSommerfeld

    using LinearAlgebra
    using Optim
    using UnPack

    export Cheb
    export os_matrices, os_eigen, optimal_growth, disturbance_energy
    export velocity_x, velocity_y, velocity_z, vorticity_x, vorticity_y, vorticity_z

    
    include("chebyshev.jl")
    include("utilities.jl")
    include("baseflows.jl")
    include("osmatrices.jl")
    include("transientgrowth.jl")
    include("fields.jl")

end