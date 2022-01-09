struct Opps

    #Coordinates
    rcoord      :: Array{Float64,1}
    tcoord      :: Array{Float64,1}

    aux_vex_1     :: Array{Float64,1}
    aux_vex_2     :: Array{Float64,1}
    
    #Operators
    Dr          :: SparseMatrixCSC{Float64,Int64}
    Dθ         :: SparseMatrixCSC{Float64,Int64}
    
    #RHS Differential Operator
    DiffO_1     :: SparseMatrixCSC{Float64,Int64}

   
end