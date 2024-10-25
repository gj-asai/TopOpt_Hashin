struct StressLimits
    Xt::Real
    Xc::Real
    Yt::Real
    Yc::Real
    Sln::Real
    Stn::Real
end
StressLimits() = StressLimits(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

struct Material{dim}
    C::SymmetricTensor{4,dim}
    stress_limits::StressLimits
end

function Isotropic2D(; E, ν, limits=StressLimits())
    S = [1/E -ν/E 0;
        -ν/E 1/E 0;
        0 0 2(1+ν)/E]
    C = fromvoigt(SymmetricTensor{4,2}, inv(S))
    return Material(C, limits)
end

function Orthotropic2D(; El, Et, νlt, Glt, limits=StressLimits())
    S = [1/El -νlt/El 0;
        -νlt/El 1/Et 0;
        0 0 1/Glt]
    C = fromvoigt(SymmetricTensor{4,2}, inv(S))
    return Material(C, limits)
end
