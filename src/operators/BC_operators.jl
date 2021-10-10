
#Cylinder impedance BC
function Z(t , θ::Float64 , p::Param)

    Z0 = p.Z0
    ϵ  = p.ϵ
    Ω  = p.Ω
    
    return Z0*(1 + ϵ*cos(θ - Ω*t)^2)
end

function dZ(t,θ,p::Param)
    
    Z0 = p.Z0
    ϵ  = p.ϵ
    Ω  = p.Ω

    return 2 * Z0 * Ω * ϵ * cos(θ - Ω * t) * sin(θ - Ω * t)
end