#' Compiling TMB C++ templates
#'
#' Compiles the  TMB templates into a shared object
#' file. Must be done a single time following installation or
#' updating of the package.
#'
#' @export
compile.mmre <- function(){
    wd <- getwd()
    dir <- paste(system.file(package = "mmre"), "/src", sep = "")
    setwd(dir)
    if (!dir.exists("../bin")){
        dir.create("../bin")
    }
    files <- strsplit(list.files(), "[.]")
    base <- sapply(files, function(x) x[1])
    ext <- sapply(files, function(x) x[2]) 
    for (i in base[ext == "cpp"]){
        compile(paste(i, ".cpp", sep = ""))
        unlink(paste(i, ".o", sep = ""))
        file.rename(paste(i, ".so", sep = ""),
                    paste("../bin/", i, ".so", sep = ""))
    }
    setwd(wd)
}

#' load DLLs for C++ templtes
#'
#' Loads required DLLs for \code{fit.mmre}
#'
#' @export
dll.mmre <- function(){
    dll.dir <- paste(system.file(package = "mmre"), "/bin/", sep = "")
    for (i in paste(dll.dir, list.files(dll.dir), sep = "")){
        dyn.load(i)
    }
}
