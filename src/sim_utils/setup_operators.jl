function set_up_operator(p::Param)

    # Variables
    rmin        = p.rmin
    rmax        = p.rmax
    rnodes      = p.rnodes

    thnodes     = p.thnodes

    #Create position vectors
    dr = (rmax - rmin)/(rnodes-1)
    r = [rmin + (i-1)*dr for i in 1:rnodes]
    
    dth = 2*pi/(thnodes)
    theta = [ (i-1)*dth for i in 1:thnodes]
    
    #Create Differential Operators 
    Dr  =   Diff_Operator_2D(1 ,1, [dr,dth] , [rnodes,thnodes]) 
    Drr =   Diff_Operator_2D(2 ,1, [dr,dth] , [rnodes,thnodes])

    Dp  =  Diff_Operator_2D_Periodic(1 ,2, [dr,dth] , [rnodes,thnodes]) 
    Dpp =  Diff_Operator_2D_Periodic(2 ,2, [dr,dth] , [rnodes,thnodes])
    
    one_over_r2 = spzeros(rnodes*thnodes,rnodes*thnodes)

    for i in 1:rnodes*thnodes
        one_over_r2[i,i] = (1/r[1 + (i-1)%rnodes])^2
    end

    RHS_Opp = Drr + one_over_r2*Dpp + (0.25)*one_over_r2
    #sparse(Drr + (0.25)*one_over_r2*I  + one_over_r2*Dpp);

    av1 = zeros(rnodes*thnodes)
    av2 = zeros(rnodes*thnodes)
    
    #Operator structure
    return Opps(r, theta , av1 , av2  , Dr, Dp, RHS_Opp)
    
end