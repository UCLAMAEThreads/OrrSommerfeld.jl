struct Cheb{N}
    y :: Vector{Float64}
    D0 :: Matrix{Float64}
    D1 :: Matrix{Float64}
    D2 :: Matrix{Float64}
    D3 :: Matrix{Float64}
    D4 :: Matrix{Float64}
end

"""
    Cheb(N::Integer)

Returns a structure containing differentiation matrices of `N` points of Chebyshev-distributed
data in [-1,1]. If `C = Cheb(10)`, then
- `C.D0` contains the Chebyshev modes 0 through N in columns. Each column represents the mode evaluated on the Chebyshev points. Thus, a Chebyshev expansion of a function could be written as `C.D0*a`, where `a` is a vector of coefficients. 
- `C.D1` contains the first derivatives of the Chebyshev modes, in columns.
- `C.D2` contains the second derivatives of the Chebyshev modes, in columns
- `C.D3` contains the third derivatives of the Chebyshev modes, in columns
- `C.D4` contains the fourth derivatives of the Chebyshev modes, in columns
"""
function Cheb(N::Integer)
    vec = 0:N
    y = _cheb_points(N)
    D0 = cos.(π*vec*vec'/N)
    D1, D2, D3, D4 = zero(D0), zero(D0), zero(D0), zero(D0)
    
    D1[:,2] = D0[:,1]
    D1[:,3] = 4D0[:,2]
    D2[:,3] = 4D1[:,2]
    for n = 3:N
        D1[:,n+1] = 2*n*D0[:,n] + n/(n-2)*D1[:,n-1]
        D2[:,n+1] = 2*n*D1[:,n] + n/(n-2)*D2[:,n-1]
        D3[:,n+1] = 2*n*D2[:,n] + n/(n-2)*D3[:,n-1]
        D4[:,n+1] = 2*n*D3[:,n] + n/(n-2)*D4[:,n-1]
    end
    return Cheb{N}(y, D0, D1, D2, D3, D4)
end

@inline _cheb_points(n) = cos.(π*(0:n)/n)