function run_Zilindroa(;
        
    #Time variables
    tspan::Array{Float64 , 1}   =   [0. , 0.]   ,
    dt::Float64                 =   0.05        ,
    out_every::Float64          =   1.0         ,
    max_run_time::String        =   "24:00:00"  ,

    #Spatial variables
    Nodes::Array{Int64 , 1}   =   [128 , 128]   ,
    Rcylinder::Float64        =   1.       ,
    Rcavity::Float64          =   30.       ,

    #Cylinder parameters
    Z::Float64              =   1.0     ,
    Omega::Float64          =   0.5     ,
    epsilon::Float64        =   0.0    ,
            
    #Initial configuration
    Gaussian_pulse:: Array{Float64, 1} = [ 3.5 , 3.7 , 20.  , 0.1 , 2 ],
            
    output_folder              = "" 

    )

if(tspan[2] - tspan[1] < 0)
   println("ERROR: t_final must be larger than t_initial")
    return
end
   
#filename
folder_name = output_folder

if(output_folder == "")
    folder_name = string("/", Dates.today() , "/")    
end   
    
filename = string(
            "Zilindroa"             , 
            "_R1_"      , Rcylinder , 
            "_R2_"      , Rcavity   , 
            "_Z_"       , Z         , 
            "_eps_"     , epsilon   ,
            "_Omega_"   , Omega     , 
            "_Nr_"      , Nodes[1]  ,
            "_Nt_"      , Nodes[2]  ,
            "_ti_"      , tspan[1]  ,
             "_tf_"     , tspan[2] )

#Setup simulation
p0 = Param(

#Temporal variables
t_sim_init      = tspan[1]      ,
t_sim_final     = tspan[2]      ,
out_every_t     = out_every     ,
deltat          = dt            ,

max_runtime     = max_run_time  ,
#Spacial variables
rmin           =    Rcylinder   ,
rmax           =    Rcavity     ,
xnodes         =    Nnodes[1]   ,

#Spacial variables
ymin           =    -Rcavity    ,
ymax           =    Rcavity     ,
ynodes         =    Nnodes      ,

#Initial Conditions
A0              =   Gaussian_pulse[1] , #3.5         ,
σ               =   Gaussian_pulse[2] , #4.5         ,
r0              =   Gaussian_pulse[3] , #40.0        ,
ω               =   Gaussian_pulse[4] , #0.5         ,
m               =   Gaussian_pulse[5] , #2           ,

Z0              = Z         ,
ϵ               = epsilon   , 
Ω               = Omega     ,

#Data variables
folder          = folder_name   ,
fname           = filename
    
 )                        

solve_zilindroa_system(p0)

return 1
end