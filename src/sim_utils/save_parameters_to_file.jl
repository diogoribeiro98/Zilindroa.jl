function save_parameters_to_file( p::Param , fid::HDF5.File  )

    #Time variables 
    # -> saved at end of simulation

    #Spatial variables
    HDF5.attributes(fid)["rmin"]      = p.rmin
    HDF5.attributes(fid)["rmax"]      = p.rmax
    HDF5.attributes(fid)["rnodes"]    = p.rnodes

    HDF5.attributes(fid)["tnodes"]    = p.tnodes

    #Initial condition
    HDF5.attributes(fid)["A0"]      = p.A0
    HDF5.attributes(fid)["sigma"]   = p.σ
    HDF5.attributes(fid)["r0"]      = p.r0
    HDF5.attributes(fid)["omega"]   = p.ω
    HDF5.attributes(fid)["m_mode"]  = p.m
    
     #Cylinder data
     HDF5.attributes(fid)["epsilon"]    = p.ϵ
     HDF5.attributes(fid)["Omega"]      = p.Ω
     HDF5.attributes(fid)["Z"]          = p.Z0
    
end