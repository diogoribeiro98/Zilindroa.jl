# Get displacement at a given iteration
function get_fields( fname::String , itt::Int)

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

    #Check it time iteration is OK
    if( itt > max_iter)
        println("ERROR: out of bounds (max iteration is" , max_iter,")")
        return
    end

    #Get Displacement
    g1      = open_group(fid, "ψ")

    dname   = string("iteration_" , itt)
    dset    = open_dataset(g1, dname )
    ttt     = read_attribute(dset,"time")
    
    ψ = reshape(Array(dset) , rnodes,thnodes )
   
    #Get Derivative
    g1      = open_group(fid, "dψ")

    dname   = string("iteration_" , itt)
    dset    = open_dataset(g1, dname )
    
    
    dψ = reshape(Array(dset),rnodes,thnodes)

    close(fid)

    #Return all data
    return ttt , x , y , ψ , dψ

end