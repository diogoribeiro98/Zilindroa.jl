#Wave equation Differential Equation
function Zilndroa_Diffeq!(df,f,param,t)

    #Get Operators and variables
    Opps        = param[1]
    vars        = param[2]

    r = Opps.rcoord
    θ = Opps.tcoord

    RHS = Opps.DiffO_1
    Dr  = Opps.Dr
    Dθ  = Opps.Dθ

    dr_ψ  = Opps.aux_vex_1
    dθ_ψ  = Opps.aux_vex_2

    #Fields of our problem (Julia changes the rhs if you change the lhs)
    ψ   = f.x[1]
    ϕ   = f.x[2]

    dψ  = df.x[1]
    dϕ  = df.x[2]

    mul!(dr_ψ , Dr , ψ)
    mul!(dθ_ψ , Dθ , ψ)

    #Differential equation
    mul!(dψ , I ,  ϕ)
    mul!(dϕ , RHS ,  ψ)

    Nr      = vars.rnodes
    Ntheta  = vars.thnodes
    Ω       = vars.Ω
    #Boundary conditions
    
    for i in 1:Ntheta
    
        #Inner BC (Impedance)   
        idx = (i-1)*Nr+1 
        #dψ[idx] = Z(t , θ[i] ,  vars ) * ( dr_ψ[idx] - 0.5 * ψ[idx] / r[1] ) - Ω * dθ_ψ[idx]
        dψ[idx] = 0
        
        #Outer BC (reflective)
        idx = i*Nr 
        dψ[idx] = 0
    end

end

function solve_zilindroa_system(p_in::Param)

    #
    # Setup Routine
    #

    #Welcome message
    welcome_message();

    #Start time
    init_time = time()
    max_runtime = time_string_to_seconds(p_in.max_runtime)

    #Setup routine setup check
    if( !setup_folder(p_in) ) 
        print_exit_message(  time() - init_time )
        return false
    end

    #Read parameters
    p = p_in
    print_parameters_to_screen(p)  
    println("OK")

    #Storage file
    println("Creating HDF5 file...")

    fdata   = "./data"
    fname   = fdata * p.folder *  "$(p.fname).h5"
    fid = h5open(fname, "w")

    save_parameters_to_file( p , fid )
    println("OK")

    #Create operators structures
    println("Setting up Operators...")
    Operators = set_up_operator(p)
    println("OK")

    #Time variables
    tmin = p.t_sim_init
    tmax = p.t_sim_final
    tspan       = ( tmin , tmax )

    iter        = 0
    iter_save   = 0
    every       = floor(Int , p.out_every_t/p.deltat + 1.0)
    save_time   = 0

    #Initial Conditions
    println("Setting up inicial configuration...")
    ψ , dψ = get_gaussian_pulse( p , Operators)

    U = ArrayPartition( ψ , dψ )
    println("OK")

    my_params = ( Operators , p)

    #Problem ODE and integrator
    prob = ODEProblem(  Zilndroa_Diffeq! , U , tspan , my_params )

    integrator = init( prob , RK4() , save_everystep=false , dt=p.deltat , adaptive=false )

    #Create a group for the functions
    g1 =  create_group(fid , "ψ");
    g2 =  create_group(fid , "dψ");
    
    #Save first configuration of field

    dset_name = string("iteration_",iter_save)
    
    data1 = integrator.u.x[1]  
    dset1, dtype1  = create_dataset(g1, dset_name , data1 )
    write_dataset(dset1, dtype1 , data1)

    data2 = integrator.u.x[2]
    dset2, dtype2  = create_dataset(g2, dset_name , data2 )
    write_dataset(dset2, dtype2 , data2)

    #Add HDF5.attributes time to both datasets
    HDF5.attributes(dset1)["time"] =  integrator.t
    HDF5.attributes(dset2)["time"] =  integrator.t

    #Progress bar and time evolution
    println("Starting up simulation...")
    prog1 = Progress( Int(floor(p.t_sim_final-p.t_sim_init))*100 )
    
    for (u,t) in tuples(integrator)
    
        #Increase iteration counter
        iter += 1
    
        #Step integrator
        DifferentialEquations.step!(integrator)
        
        #Update bar
        ProgressMeter.update!(prog1, Int(floor(integrator.t*100)))

        #If save iteration
        if iter % every == 0
            
            iter_save += 1
            save_time = integrator.t
            
            #data set names
            dset_name = string("iteration_",iter_save)
    
            #data set for both displacement and derivative
            data1 = integrator.u.x[1]  
            dset1, dtype1  = create_dataset(g1, dset_name , data1 )
            write_dataset(dset1, dtype1 , data1)

            data2 = integrator.u.x[2]
            dset2, dtype2  = create_dataset(g2, dset_name , data2 )
            write_dataset(dset2, dtype2 , data2)

            #Add HDF5.attributes time to both datasets
            HDF5.attributes(dset1)["time"] =  integrator.t
            HDF5.attributes(dset2)["time"] =  integrator.t

        end
  
        #Check if time has exceeded
        if( time_since(init_time) > max_runtime )
            println("\nMaximum allowed runtime reached , Ending Simulation!")
            break
        end


    end

    #Add time information to file
    HDF5.attributes(fid)["max_iter"]    = iter_save
    HDF5.attributes(fid)["t_min"]       = p.t_sim_init
    HDF5.attributes(fid)["t_max"]       = save_time
    HDF5.attributes(fid)["out_every_t"] = p.out_every_t
    HDF5.attributes(fid)["dt"]          = p.deltat

    close(fid)
    println("Saving Files...")
   
    print_exit_message(  time() - init_time )

    return true

end