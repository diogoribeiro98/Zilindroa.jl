#Create Initially traveling gaussian wavepacket
function get_gaussian_pulse(p::Param , Opp::Opps)

    #Get Variables
    Nr  = p.rnodes
    Nθ  = p.thnodes

    A0   = p.A0 
    σ   = p.σ
    ω   = p.ω
    r0  = p.r0
    m   = p.m

    #Get operators
    r   = Opp.rcoord
    θ   = Opp.tcoord

    Dr  = Opp.Dr

    #Create initial state
    ψ_matrix = zeros( Nr , Nθ )
    dψ_matrix = zeros( Nr , Nθ )

    for i in 1:Nr , j in 1:Nθ
    
        ψ_matrix[i,j] = A0*cos( m*θ[j] ) * sin( ω*r[i] ) * exp( - 0.5* ( (r[i]-r0)/ σ)^2 )
        dψ_matrix[i,j] = A0*cos( m*θ[j] ) * exp( - 0.5* ( (r[i]-r0)/ σ)^2 ) * ( ω*cos(ω*r[i]) - (r[i]-r0)*sin(ω*r[i])/σ^2  )

    end
    
    ψ  = reshape(ψ_matrix  , Nr*Nθ)
    dψ = reshape(dψ_matrix , Nr*Nθ)

    #Make boundaries zeros
    for i in 1:Nθ
    
        #Inner BC 
        idx = (i-1)*Nr+1 
        ψ[idx] = 0
        dψ[idx] = 0
    
        #Outer BC 
        idx = i*Nr 
        ψ[idx] = 0
        dψ[idx] = 0
    end



    #Return initial state
    return ψ , dψ

end