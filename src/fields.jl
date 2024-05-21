for fname in (:velocity_x, :velocity_y, :velocity_z, :vorticity_x, :vorticity_y, :vorticity_z)

    ufname = Symbol("_",fname)

    @eval function $fname(qhat::Vector{ComplexF64},d::OSMatrix{N}; glims = (0,4), Ng = 200, plane = :xy) where {N}
        @unpack α, β, ak2, C = d

        vhat_cheb, ηhat_cheb = qhat[1:N+1], qhat[N+2:2N+2]
        fhat = $ufname(vhat_cheb,ηhat_cheb,α,β,ak2,C)

        xg = range(glims...,length=Ng)
        wavenum = _get_wavenumber(α,β,Val(plane))

        return collect(xg), reverse(C.y), real(reverse(fhat)*exp.(im*(wavenum*xg')))
    end

    @eval function $fname(κ0::Vector{ComplexF64},d::OSMatrix{N}, eos::OSEigen, t::Real; glims = (0,4), Ng = 200, plane = :xy) where {N}
        @unpack α, β, ak2, C = d

        λu, Xu = eos.values, eos.vectors
        κt = exp(-im*Diagonal(λu)*t)*κ0
        qhat = Xu*κt

        vhat_cheb, ηhat_cheb = qhat[1:N+1], qhat[N+2:2N+2]
        fhat = $ufname(vhat_cheb,ηhat_cheb,α,β,ak2,C)

        xg = range(glims...,length=Ng)
        wavenum = _get_wavenumber(α,β,Val(plane))

        return collect(xg), reverse(C.y), real(reverse(fhat)*exp.(im*(wavenum*xg')))
    end

end

_velocity_x(vhat,ηhat,α,β,ak2,C) = im/ak2*(α*C.D1*vhat - β*C.D0*ηhat)
_velocity_y(vhat,ηhat,α,β,ak2,C) = C.D0*vhat
_velocity_z(vhat,ηhat,α,β,ak2,C) = im/ak2*(β*C.D1*vhat + α*C.D0*ηhat)
_vorticity_x(vhat,ηhat,α,β,ak2,C) = -im*β/ak2*(ak2*C.D0 .- C.D2)*vhat + im*α/ak2*C.D1*ηhat
_vorticity_y(vhat,ηhat,α,β,ak2,C) = C.D0*ηhat
_vorticity_z(vhat,ηhat,α,β,ak2,C) = im*α/ak2*(ak2*C.D0 .- C.D2)*vhat + im*β/ak2*C.D1*ηhat


_get_wavenumber(α,β,::Val{:xy}) = α
_get_wavenumber(α,β,::Val{:yz}) = β