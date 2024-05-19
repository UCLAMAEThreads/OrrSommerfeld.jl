struct OSMatrix{N}
    Aos :: Matrix{ComplexF64}
    Bos :: Matrix{ComplexF64}
    Asq :: Matrix{ComplexF64}
    Bsq :: Matrix{ComplexF64}
    Aoq :: Matrix{ComplexF64}
end


"""
    os_matrices(α,β,Re,C::Cheb[;baseflow=:poiseuille]) -> OSMatrix

Returns the O-S matrices for a given base flow (between walls at y = +1/-1), for
streamwise wavenumber `α`, spanwise wavenumber `β`, and Reynolds number `Re`.
The default base flow is Poiseuille flow, but this can be changed to `:couette` for
Couette flow.
"""
function os_matrices(α::Real,β::Real,Re::Real,C::Cheb{N};baseflow=:poiseuille) where N
    @unpack y, D0, D1, D2, D3, D4 = C


    ak2 = α^2 + β^2
    Nos = N+1
    Nsq = N+1

    # Base flow and its derivative
    U, dU, d2U  = eval(baseflow)(y)    

    Aos, Aoq, Asq = zeros(ComplexF64,Nos,Nos), zeros(ComplexF64,Nos,Nos), zeros(ComplexF64,Nos,Nos)
    Bos, Bsq = zeros(ComplexF64,Nos,Nos), zeros(ComplexF64,Nos,Nos)
    er = -200.0*im # This ensures that spurious eigenvalues from bc's are far away

    # Form the Orr-Sommerfeld operators Aos and Bos
    Bos .= D2 - ak2*D0
    Aos .= -(D4 - 2*ak2*D2 + ak2^2*D0)/(im*Re)
    Aos .+= α*Diagonal(U)*Bos - α*Diagonal(d2U)*D0
    Aos[1,:] = er*D0[1,:] # v(1) = 0
    Aos[2,:] = er*D1[1,:] # dv/dy(1) = 0
    Aos[Nos-1,:] = er*D1[Nos,:] # dv/dy(-1) = 0
    Aos[Nos,:] = er*D0[Nos,:] # v(-1) = 0

    Bos[1,:] = D0[1,:]
    Bos[2,:] = D1[1,:]
    Bos[Nos-1,:] = D1[Nos,:]
    Bos[Nos,:] = D0[Nos,:]

    # Form the Squire operators Asq and Bsq
    Asq .= α*Diagonal(U)*D0 - (D2-ak2*D0)/(im*Re)
    Bsq .= D0
    Asq[1,:] = er*D0[1,:]  # eta(1) = 0
    Asq[Nsq,:] = er*D0[Nsq,:] # eta(-1) = 0

    # Form the coupling operator from OS to Squire
    Aoq .= β*Diagonal(dU)*D0 
    Aoq[1,:] = zeros(1,Nos)
    Aoq[Nsq,:] = zeros(1,Nos)

    return OSMatrix{N}(Aos, Bos, Asq, Bsq, Aoq)
    
end

"""
    energy_matrix(Nos,Nsq,ak2) -> Matrix{Float64}

Computes the energy weight matrix for Couette and Poiseuille flow,
with `Nos` OS modes and `Nsq` Squire modes, and `ak2` equal to ``k^2 = α^2 + β^2``.

"""
function energy_matrix(Nos,Nsq,ak2)
    Cos = twonorm_weights(Nos)
    Dos = deven(Nos)

    Wos = Dos'*Cos*Dos + ak2*Cos
    Wsq = twonorm_weights(Nsq)

    Fos = svd(Wos)
    Mos = Diagonal(sqrt.(Fos.S))*Fos.U'

    Fsq = svd(Wsq)
    Msq = Diagonal(sqrt.(Fsq.S))*Fsq.U'

    return [Mos zeros(Nos,Nsq); zeros(Nsq,Nos) Msq]

end

"""
    os_eigen(d::OSMatrix[;matrixpart = :os])

Compute the eigenvalues and eigenvectors of the OS matrix `d`. The
optional argument `matrixpart` defaults to `:os`, which returns just
the OS modes. Instead, this can be set to `:sq` for Squire modes only,
or `:full` for all modes.
"""
os_eigen(d::OSMatrix; matrixpart = :os, reduce = (-Inf,Inf,I)) = _os_eigen_raw(d,reduce...,Val(matrixpart))

function _os_eigen_raw(d,minval,maxval,M,::Val{:os})
    F = eigen(d.Aos,d.Bos,sortby=(y -> -imag(y)))
    #_reduce_eigen(F,minval,maxval,M)
end

_os_eigen_raw(d,minval,maxval,M,::Val{:sq}) = eigen(d.Asq,d.Bsq,sortby=(y -> -imag(y)))

function _os_eigen_raw(d,minval,maxval,M,::Val{:full})
    nos, nsq = size(d.Aoq)
    A = [d.Aos zeros(nos,nsq); d.Aoq d.Asq]
    B = [d.Bos zeros(nos,nsq); zeros(nsq,nos) d.Bsq]
    eigen(A,B,sortby=(y -> -imag(y)))
end


#=
λ = F.values
    X = copy(F.vectors)
    normalize_columns!(X,M)

    imin, imax = count_eigenvalues(λ,minval;maxval=1.0)
    Xu = view(X,:,imin:imax)
    λu = view(λ,imin:imax)
    =#