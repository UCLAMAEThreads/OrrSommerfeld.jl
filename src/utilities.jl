"""
    normalize_columns!(X::AbstractMatrix,M::AbstractMatrix) -> AbstractMatrix

Returns X with columns normalized by the 2-norm of the corresponding column of M*X,
so that the resulting 2-norm of each column of M*X is equal to 1
"""
function normalize_columns!(X::AbstractMatrix{T},M::AbstractMatrix) where {T<:Union{Float64,ComplexF64}}
    Y = M*X
    for (x,y) in zip(eachcol(X),eachcol(Y))
        x .= x/norm(y)
    end
    return X
end

"""
    count_eigenvalues(evals::Vector,a::Real[;maxval=0.5])

For a set of eigenvalues `evals` ordered with decreasing
imaginary part, find the first and last indices such that 

    a <= imag(evals[i]) <= maxval
 
for ifirst <= i <= ilast
"""
function count_eigenvalues(evals::Vector,a::Real;maxval=0.5)
    ifirst = 1
    while imag(evals[ifirst]) > maxval
        ifirst += 1
    end
    ilast = ifirst
    while imag(evals[ilast]) > a
        ilast += 1
    end
    return ifirst, ilast - 1
end