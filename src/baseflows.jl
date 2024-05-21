# Base flow functions

const TANHFACTOR=7.0

for bflow in (:poiseuille, :couette, :sharpmixing, :tanhmixing)
    Ufcn = Symbol("_U_",bflow)
    dUfcn = Symbol("_dU_",bflow)
    d2Ufcn = Symbol("_d2U_",bflow)
    @eval function $bflow(y::AbstractVector)
        U = $Ufcn.(y)
        dU = $dUfcn.(y)
        d2U = $d2Ufcn.(y)
        return U, dU, d2U
    end

    @eval function $bflow(y::Float64)
        U = $Ufcn(y)
        dU = $dUfcn(y)
        d2U = $d2Ufcn(y)
        return U, dU, d2U
    end
end

_U_poiseuille(y) = 1 - y^2
_dU_poiseuille(y) = -2*y
_d2U_poiseuille(y) = -2


_U_couette(y) = y
_dU_couette(y) = 1
_d2U_couette(y) = 0

_U_sharpmixing(y) = y >= 0 ? 1 : -1
_dU_sharpmixing(y) = 0
_d2U_sharpmixing(y) = 0

_U_tanhmixing(y) = 0.5*(1+tanh(TANHFACTOR*y))
_dU_tanhmixing(y) = 0.5*TANHFACTOR*sech(TANHFACTOR*y)^2
_d2U_tanhmixing(y) = -TANHFACTOR^2*tanh(TANHFACTOR*y)*sech(TANHFACTOR*y)^2

