# Base flow functions

for bflow in (:poiseuille, :couette)
    Ufcn = Symbol("_U_",bflow)
    dUfcn = Symbol("_dU_",bflow)
    d2Ufcn = Symbol("_d2U_",bflow)
    @eval function $bflow(y::AbstractVector)
        U = $Ufcn.(y)
        dU = $dUfcn.(y)
        d2U = $d2Ufcn.(y)
        return U, dU, d2U
    end
end

_U_poiseuille(y) = 1 - y^2
_dU_poiseuille(y) = -2*y
_d2U_poiseuille(y) = -2


_U_couette(y) = y
_dU_couette(y) = 1
_d2U_couette(y) = 0