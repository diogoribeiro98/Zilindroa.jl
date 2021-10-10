#
#   Differential Operators 1D of order 2 (Q = 4)
#

function Diff_Operator_1D(  derivative_order::Int , dx::T , N::Int ) where T

    #Check dimension of array
    if(N < 4)
        print("ERROR: Dimension of operator must be greater than 3")
        return
    end

    #Create operator
    Opp = zeros(T , N,N)

    #1st order derivative
    if( derivative_order == 1)
       
        d1_1 = Array{T}([ -11.0/6.0  , 	3.0	, -3.0/2.0 , 	1.0/3.0	]/dx)
        d1_2 = Array{T}([ -1.0/4.0  , -5.0/6.0	, 3.0/2.0 , 	-1.0/2.0	 , 1.0/12.0	]/dx)
        d1_c = Array{T}([ 1.0/12.0  , 	-2.0/3.0	, 0 , 	2.0/3.0	 , -1.0/12.0	]/dx)
        
        #Fill first line  and last
        for i in 1:length(d1_1)

            #First row (left to right)
            Opp[ 1 , i ]        = d1_1[i]

            #Last row (right to left)
            Opp[ N , N-(i-1) ]  = -d1_1[i]
        end

        #Fill second line and second lasts
        for i in 1:length(d1_2)

            #Second row (left to right)
            Opp[ 2 , i ]        = d1_2[i]

            #Last row (right to left)
            Opp[ N-1 , N-(i-1) ]  = -d1_2[i]    
        end
    
        #Fill remaining lines
        for i in 3:(N-2)
            for j in 1:length(d1_c)
                Opp[i , i+j-3 ] = d1_c[j]
            end
        end

        #Return operator
        return Opp 

    end

    #2nd order derivative
    if( derivative_order == 2)

        blah = 1/ (dx^2)
        d2_1 = Array{T}([ 35.0 / 12.0 , 	-26.0/3.0	, 19.0/2.0 , 	-14.0/3.0 , 11.0/12.0	] * blah) 
        d2_2 = Array{T}([ 5.0 / 6.0 , 	-5.0/4.0	, -1.0/3.0 , 	7.0/6.0 , -1.0/2.0 , 1.0/12.0	] * blah )
        d2_c = Array{T}([ -1.0/12.0  , 	4.0/3.0	, -5.0/2.0 , 	4.0/3.0	 , -1.0/12.0	]* blah )



        #Fill first and last line
        for i in 1:length(d2_1)

            #First row (left to right)
            Opp[ 1 , i ]        = d2_1[i]

            #Last row (right to left)
            Opp[ N , N-(i-1) ]  = d2_1[i]
        end

        #Fill second and second last line
        for i in 1:length(d2_2)

            #First row (left to right)
            Opp[ 2 , i ]        = d2_2[i]

            #Last row (right to left)
            Opp[ N-1 , N-(i-1) ]  = d2_2[i]
        end

        #Fill remaining lines
        for i in 3:(N-2)
            for j in 1:length(d2_c)
                Opp[i , i+j-3 ] = d2_c[j]
            end
        end

        #Return operator
        return Opp

    end

    #Higher derivatives not avaliable
    if(derivative_order > 2 )
        print("ERROR: Derivative Order still not implementd")
        return
    end

end


function Diff_Operator_1D_Periodic(  derivative_order::Int , dx::T , N::Int )  where T

    #Check dimension of array
    if(N < 4)
        print("ERROR: Dimension of operator must be greater than 3")
        return
    end

    #Create operator
    Opp = spzeros(T , N,N)

    #1st order derivative
    if( derivative_order == 1)

        #Higher Order Derivative ( Q = 16)
        for i in 1:N , j in 1:N
            i-j==1   && (Opp[i,j] =  -2.0 / 3.0)
            j-i==1   && (Opp[i,j] =   2.0 / 3.0)
            j-i==N-1 && (Opp[i,j] = - 2.0 / 3.0)
            i-j==N-1 && (Opp[i,j] =   2.0 / 3.0)
            
            i-j==2   && (Opp[i,j] =   1.0 / 12.0)
            j-i==2   && (Opp[i,j] =  -1.0 / 12.0)
            j-i==N-2 && (Opp[i,j] =   1.0 / 12.0)
            i-j==N-2 && (Opp[i,j] =  -1.0 / 12.0)

        end

        return Opp/dx
    end

    #2nd order derivative
    if( derivative_order == 2)

        #Higher order derivative ( Q = 16)
        for i in 1:N , j in 1:N
            abs(i-j)==2   && (Opp[i,j] = -1.0 / 12.0 )
            abs(i-j)==1   && (Opp[i,j] =  4.0 / 3.0  )
            i == j        && (Opp[i,j] = -5.0 / 2.0  )
            abs(i-j)==N-1 && (Opp[i,j] =  4.0 / 3.0  )
            abs(i-j)==N-2 && (Opp[i,j] = -1.0 / 12.0)

        end

        
        return Opp/dx^2
    end

    #Higher derivatives not avaliable
    if(derivative_order > 2 )
        print("ERROR: Derivative Order still not implementd")
        return
    end

end

#
#   Differential Operators 2D 
#

function Diff_Operator_2D( derivative_order::Int , index::Int , delta , dim)

    if index == 1 
        A = Diff_Operator_1D(derivative_order , delta[1] , dim[1] )
        return kron(sparse( I, dim[2], dim[2] ) , A)
    end
    
    if index == 2
        A = Diff_Operator_1D(derivative_order , delta[2] , dim[2] )
        return kron(A, sparse( I, dim[1], dim[1] ))
    end

    #Derivative order 3
    if(index > 2 )
        print("ERROR: Out of bounds")
        return
    end
end


function Diff_Operator_2D_Periodic( derivative_order::Int , index::Int , delta , dim)

    if index == 1 
        A = Diff_Operator_1D_Periodic(derivative_order , delta[1] , dim[1] )
        return kron(sparse( I, dim[2], dim[2] ) , A)
    end
    
    if index == 2
        A = Diff_Operator_1D_Periodic(derivative_order , delta[2] , dim[2] )
        return kron(A, sparse( I, dim[1], dim[1] ))
    end

    #Derivative order 3
    if(index > 2 )
        print("ERROR: Out of bounds")
        return
    end
end