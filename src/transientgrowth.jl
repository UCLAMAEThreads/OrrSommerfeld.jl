struct OptimalGrowth{DT,GT}
    # OS matrix data
    d :: OSMatrix

    # Range of times over which the growth was inspected
    T :: Tuple

     # OS eigenvalues/vectors
     e :: OSEigen

    # Maximum energy growth and the associated time
    tformax :: Float64
    mgrowth :: Float64

    # Optimal energy growth t and amplitude
    t :: DT
    G :: GT

    # These are the optimal disturbance and response, expressed in
    # Chebyshev expansion coefficients 
    flow0 :: Vector{ComplexF64}
    flowt :: Vector{ComplexF64}

    # These are the coordinates of the optimal disturbance and response
    # in the OS eigenvector space
    κ0 :: Vector{ComplexF64}
    κt :: Vector{ComplexF64}
end

"""
    optimal_growth(d::OSMatrix,T::Tuple[; minval = -1.5, ntime = 100]) -> OptimalGrowth

Returns the optimal growth data for the system described by Orr-Sommerfeld matrix `d`, in a range of times `T`.
Only eigenvalues with imaginary parts above `minval` are used (default -1.5).
The results are returned in the `OptimalGrowth` type variable.
"""
function optimal_growth(d::OSMatrix,T; minval = -1.5, ntime = 100)
    @unpack M, ak2 = d

    eos = os_eigen(d,matrixpart=:full,ilims=(minval,1.0),normalize=:true)
    λu, Xu = eos.values, eos.vectors

    Qb, invF = qbmat(M,Xu,λu)

    maxer(t::Real) = -opnorm(exp(t*Qb))

    # Find the time and value of maximum energy 
    gcheck = maxer(0.01)^2 >= 1
    tformax, mgrowth = _max_growth(maxer,T,Val(gcheck))

    #= 
    Perform an SVD of the evolution matrix at the time of maximum growth
    to find the optimal disturbance and the corresponding growth rate
    =#
    evol = exp(tformax*Qb)
    Fevol = svd(evol)
    mgrowth = Fevol.S[1]^2

    # Reconstruct the optimal q0 and qt (in Chebyshev coeffients)
    κ0 = invF*Fevol.V[:,1]
    flow0 = sqrt(2*ak2)*Xu*κ0

    κt = invF*Fevol.U[:,1]
    flowt = sqrt(2*ak2)*Xu*κt
    
    # Compute the maximum growth rate curve
    tid, gg = _collect_growth(maxer,T,ntime)    

    return OptimalGrowth(d,T,eos,tformax,mgrowth,tid,gg,flow0,flowt,κ0,κt)

end

_max_growth(fcn,T,::Val{false}) = 0.0, 1.0

_max_growth(fcn,T::Real,::Val{true})  = T, fcn(T)^2

function _max_growth(fcn,T::Tuple{R,R},::Val{true}) where R <: Real
    ts, tf = T
    mindata = optimize(fcn,ts,tf)
    return mindata.minimizer, mindata.minimum^2
end

_collect_growth(fcn,T::Real,ntime::Int) = T, fcn(T)^2

function _collect_growth(fcn,T::Tuple{R,R},ntime::Int) where R <: Real
    tid = range(T...,length=ntime)
    gg = [fcn(t)^2 for t in tid]
    return tid, gg
end

"""
    disturbance_energy(κ0::Vector,d::OSMatrix,eos::OSEigen,t::Real) -> Real

Given a disturbance vector `κ0` (expressed as coefficients of OS eigenvectors) for an O-S system `d`
with associated eigenvalues/vectors `eos`, compute the disturbance energy
at time `t`. If `t` is a vector of times, then a vector of energies is returned.
"""
function disturbance_energy(κ0::Vector{ComplexF64},d::OSMatrix{N}, eos::OSEigen, t::Real) where N
    @unpack α, β, ak2, C, M = d

    λu, Xu = eos.values, eos.vectors
    κt = exp(-im*Diagonal(λu)*t)*κ0
    qhat = Xu*κt
    work = M*qhat

    return norm(work)^2

end

disturbance_energy(κ0::Vector{ComplexF64},d::OSMatrix{N}, eos::OSEigen, tr::AbstractVector) where N =
    [disturbance_energy(κ0,d,eos,t) for t in tr]

"""
    qbmat(M,X,λ) -> Matrix{ComplexF64}, Matrix{ComplexF64}

Computes the transient growth matrix ``Q = -i F \\lambda F^{-1}``,
where ``\\lambda`` are the eigenvalues of the OS matrix and ``F`` the factorization
matrix of the energy weight matrix `M`. It also returns ``F^{-1}``.
"""
function qbmat(M,X,λ)
    work = M*X
    M̃ = work'*work

    FM = svd(M̃)
    σ = sqrt.(FM.S)
    F = Diagonal(σ)*FM.U'
    invF = FM.U*Diagonal(1.0./σ)
    return -im*F*Diagonal(λ)*invF, invF

end

