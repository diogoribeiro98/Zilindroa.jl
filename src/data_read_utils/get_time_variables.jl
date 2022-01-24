function get_time_variables( fname::String )

    #Get Variables  
    fid         = h5open( fname , "r")

    #Max iteration and time
    max_iter = read_attribute(fid, "max_iter")
    t_min = read_attribute(fid, "t_min")
    t_max = read_attribute(fid, "t_max")

    close(fid)
    
    return t_min , t_max , max_iter 

end