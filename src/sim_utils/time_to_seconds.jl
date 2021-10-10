#Get current time since t0
time_since(t0) = time() - t0

#String to time
function time_string_to_seconds( time::String)

    #Split time string by ':'
    v = split(time , ":" , limit = 3)

    #Calculate time in seconds
    t_seconds = 0
    for i in 1:3
        t_seconds += parse(Float64 , v[i])*60^(3-i)
    end

    return t_seconds
end

#String to time
function time_seconds_to_string( seconds::Float64)

    #
    hour = round( Int, seconds / (60*60))
    seconds %= (60*60)
    minutes = round(Int , seconds / 60, ) 
    seconds %= 60
    
    return string(  lpad(hour,2,"0") , ":" ,lpad(minutes,2,"0") , ":" , lpad(round(Int , seconds),2,"0") )

end