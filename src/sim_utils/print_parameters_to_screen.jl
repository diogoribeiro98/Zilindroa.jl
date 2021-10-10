function print_parameters_to_screen(p::Param)

    #Grid information
    ti          = p.t_sim_init
    tf          = p.t_sim_final

    t_save      = p.out_every_t
    deltat      = p.deltat
    
    rmin        = p.rmin
    rmax        = p.rmax
    Nr          = p.rnodes

    Nt          = p.tnodes

    dr      = (xmax - xmin)/(Nr-1)
    dt      = 2*pi/Nt

    #Print information
    println("---------------------------------------------")
    println("Parameter Reader:")
    println("---------------------------------------------")
    println("Grid Information:")
    println("r      = [" , tmin , "," , tmax ,"]" )
    println("t      = [" , 0    , "," , 2*pi ,"]" )
    println("dr     = "  , dr) 
    println("dtheta = "  , dt)
    println("---------------------------------------------")
    println("Time Information:")
    println("t      = [" , ti , "," , tf ,"]" )
    println("dt     = "  , deltat)
    println("---------------------------------------------")
    println("Initial Conditions:")
    println("A0     = " , p.A0)
    println("σ      = " , p.σ)
    println("r0     = " , p.r0)
    println("ω      = " , p.ω)
    println("m      = " , p.m)
    println("---------------------------------------------")
    println("Cylinder data:")
    println("Z           = " , p.Z0)
    println("epsilon     = " , p.ϵ )
    println("Omega       = " , p.Ω )
    println("---------------------------------------------")
    println("saving to file  = " ,  string(p.folder ,"/", p.fname , ".h5"))
    println("---------------------------------------------")

end