struct OptimalGrowth
    d :: OSMatrix
    T :: Tuple{Float64,Float64}
    ak2 :: Real
    t :: Vector{Float64}
    G :: Vector{Float64}
    flow0 :: Vector{ComplexF64}
    flowt :: Vector{ComplexF64}
end

function optimal(d::OSMatrix,T::Tuple{Float64,Float64},M::Matrix,ak2::Real; minval = -1.5, ntime = 100)
    eOS = os_eigen(d,matrixpart=:full)

    λ = eOS.values
    X = copy(eOS.vectors)
    normalize_columns!(X,M)

    imin, imax = count_eigenvalues(λ,minval;maxval=1.0)
    Xu = view(X,:,imin:imax)
    λu = view(λ,imin:imax)
    
    Qb, invF = qbmat(M,Xu,λu)

    maxer(t::Real) = -opnorm(exp(t*Qb))

    gcheck = maxer(0.01)^2
    if gcheck < 1
        tformax, mgrowth = 0.0, 1.0
    else
        ts, tf = T
        mindata = optimize(maxer,ts,tf)
        tformax, mgrowth = mindata.minimizer, mindata.minimum^2
    end

    evol = exp(tformax*Qb)
    Fevol = svd(evol)
    mgrowth = Fevol.S[1]^2

    flowin = sqrt(2*ak2)*Xu*invF*Fevol.V[:,1]
    flowout = sqrt(2*ak2)*Xu*invF*Fevol.U[:,1]
    
    tid = range(ts,tf,length=ntime)
    gg = [maxer(t)^2 for t in tid]
    return OptimalGrowth(d,T,ak2,tid,gg,flowin,flowout)

end

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

