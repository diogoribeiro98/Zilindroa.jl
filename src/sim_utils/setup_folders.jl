function setup_folder(p::Param)

    println("Checking folder structure...")
    setup = false

    fdata   = "./data"
    folder  = fdata * p.folder
    fname   = fdata * p.folder *  "$(p.fname).h5"

     # if no name specified error out
     if folder == ""
        println("ERROR: Output folder cannot be empty string.")
        return setup
    end

    #Create data folder
   if !isdir(fdata)
        println(" './data' folder did not exist, creating it...")
        mkdir(fdata)
    end

    #Create specified folder inside data folder if it does not exist
    if !isdir(   folder )
        println("Output folder did not exist, creating it...")
        mkdir(folder)
    end

    if isfile( fname )

        println("File " * fname * " already exists" )
        return setup
        
    end

    setup = true
    println("OK")

    return setup
end