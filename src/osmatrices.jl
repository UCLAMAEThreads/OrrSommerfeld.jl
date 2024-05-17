struct OSMatrix{N}
    Aos :: Matrix{ComplexF64}
    Bos :: Matrix{ComplexF64}
    Asq :: Matrix{ComplexF64}
    Bsq :: Matrix{ComplexF64}
    Aoq :: Matrix{ComplexF64}
end


"""
compute_OS_matrices(α,β,Re,C::Cheb[;baseflow=:poiseuille]) -> OSMatrix

Returns the O-S matrices for a given base flow (between walls at y = +1/-1), for
streamwise wavenumber `α`, spanwise wavenumber `β`, and Reynolds number `Re`.
The default base flow is Poiseuille flow, but this can be changed to `:couette` for
Couette flow.
"""
function compute_OS_matrices(α::Real,β::Real,Re::Real,C::Cheb{N};baseflow=:poiseuille) where N
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

