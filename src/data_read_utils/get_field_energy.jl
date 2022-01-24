function simpson_integration(x_vec,y_vec)
    
    sum = 0.
    dx = x_vec[2]-x_vec[1]

    #First point
    sum += 3*y_vec[1]/8 + 3*y_vec[end]/8
    sum += 7*y_vec[2]/6 + 7*y_vec[end-1]/6
    sum += 23*y_vec[3]/24 + 23*y_vec[end-2]/24

    for idx in 4:(length(x_vec)-3)
        sum += y_vec[idx]
    end
    return sum*dx
end

function simpson_integration_2D(x_vec,y_vec,z_vec)
    
    sum_vector = zeros(length(y_vec))

    for idx in 1:length(sum_vector)
        sum_vector[idx] = simpson_integration(x_vec,z_vec[:,idx])
    end

    sum = simpson_integration(y_vec , sum_vector)
     
    return sum
end

function get_field_energy( fname::String, every::Int , tlimit::Float64)

    #Get Variables  
    fid         = h5open( fname , "r")

    rmin = read_attribute(fid, "rmin")
    rmax = read_attribute(fid, "rmax")
    rnodes   = read_attribute(fid, "rnodes")

    thnodes   = read_attribute(fid, "thnodes")

    #Discretize space
    dr = (rmax - rmin)/(rnodes-1)
    r = [rmin + (i-1)*dr for i in 1:rnodes]
    
    dth = 2*pi/(thnodes)
    theta = [ (i-1)*dth for i in 1:thnodes]

    #Create Differential Operators 
    Dr  =   Diff_Operator_2D(1 ,1, [dr,dth] , [rnodes,thnodes]) 
    Drr =   Diff_Operator_2D(2 ,1, [dr,dth] , [rnodes,thnodes])

    Dp  =  Diff_Operator_2D_Periodic(1 ,2, [dr,dth] , [rnodes,thnodes]) 
    Dpp =  Diff_Operator_2D_Periodic(2 ,2, [dr,dth] , [rnodes,thnodes])
    
    one_over_r  = spzeros(rnodes*thnodes,rnodes*thnodes)
    r2_matrix = spzeros(rnodes*thnodes,rnodes*thnodes)

    for i in 1:rnodes*thnodes
        one_over_r[i,i] = (1/r[1 + (i-1)%rnodes])
        r2_matrix[i,i] = (r[1 + (i-1)%rnodes])^2
    end
    
    #Get max time iteration 
    _ , tmax , max_sim_iter  = get_time_variables(fname)    
    max_plot_iter = max_sim_iter

    if(tlimit <=  tmax)
        max_plot_iter = floor(Int , max_sim_iter * tlimit/tmax )
    else
        println("ATTENTION: max time input is larger than simulation time, using tmax")
    end
    
    # Total volume of space
    Volume = ((ymax-ymin)*(xmax-xmin))
    dV = dx*dy
    I_factor = dV / Volume

    #Create storage matrix
    Energy_matrix = zeros( rnodes , thnodes )

    #Storage vectors
    E_vector = Float64[]
    t_vector = Float64[]

    t = 0
    E = 0
    
    #Loop to get energy
    for i in 0:every:max_plot_iter
        
        #Time and fields
        t, _ , _ , ψ , dψ = get_fields(fname,i)
        
        #Gradient term
        Psi  = reshape(ψ, rnodes*thnodes)

        vr = reshape( Drr * Psi ,  (rnodes,thnodes))
        vp = reshape( one_over_r * Dp * Psi , (rnodes,thnodes))

        Energy_matrix = (vr.^2 .+ vy.^2 .+ dψ.^2) * r2_matrix
    
        E = simpson_integration_2D(x , y ,Energy_matrix)

        #Add to storage vectors  and normalize
        push!(t_vector, t)
        push!(E_vector, I_factor*E  )

    end

    return t_vector , E_vector

end