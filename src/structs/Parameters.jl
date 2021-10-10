@with_kw struct Param

    #Temporal variables for simulation
    t_sim_init          :: Float64   = 0.0
    t_sim_final         :: Float64   = 0.0   
    out_every_t         :: Float64   = 1.0
    deltat              :: Float64   = 0.05

    #Max possible sim
    max_runtime         ::String    = "70:00:00"
    
   #Spacial variables
   rmin            :: Real          =  1.0
   rmax            :: Float64       =  30.0
   rnodes          :: Integer       =  128

   thnodes         :: Integer       =  128
    
    #Initial Conditions
    A0              ::Float64       = 1.0
    σ               ::Float64       = 2.0
    r0              ::Float64       = 5.0
    ω               ::Float64       = 5.0
    m               ::Integer       = 0.0

    #Cylinder Variables
    ϵ               ::Real      = 0.1
    Ω               ::Real      = 0.3
    Z0              ::Real      = 0.001

    #Data variables
    folder          :: String    = "data/"
    fname           :: String    = "default_file_name"

end