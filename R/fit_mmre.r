#' Funtions requred to fit two state Markov model with or without individual effects
#' @details
#' Must compile required TMB templates to use

              
#' The mmre class
#' @section Slots:
#' \describe{
#' \item{data}{a data frame which must contain the following named elements ID,
#' state, time to be used with \code{fit.mmre}}
#' \item{summary}{a summary of data, automatiacally created when calling \code{fit.mmre},
#' filled using \code{mmre.summary}. Includes raw transition probabilities.}
#' \item{fit_data}{a list of the data required to use \code{fit.mmre} (automatically filled)}
#' \item{continuous}{logical, should a continuous model be used (default)}
#' \item{fitted}{a character of the chosen model based on arguments given to \code{fit.mmre},
#' chosen using \code{mmre.mod}}
#' \item{cov_names}{a charater string or vector of covariate names (elements in \code{data})
#' (only for multiple individuals)}
#' \item{parameters}{a named list of starting values for the model, assumed model is chosen based
#' on the contents of this list (i.e., to fit a correlated random effect model \code{parameters} must contain
#' the elements \code{u}, \code{log_sig_u} and \code{cov_par})}
#' \item{fit}{contains the model output}
#' \item{sdreport}{summary information of the model}
#' }
#' @name mmre
#' @rdname mmre
#' @aliases mmre-class
#' @exportClass mmre
mmre <- setClass("mmre",
                 slots = c(data = "data.frame",
                           summary = "list",
                           fit_data = "list",
                           fitted = "character",
                           cov_names = "character",
                           decay = "logical",
                           parameters = "list",
                           fit = "list",
                           sdreport = "matrix"),
                 validity = function(object){
                     if(sum(c("ID", "state","time") %in% names(object@data))!=3) {
                         return("Data must contain elements named ID, state, and time. Try ?fit.mmre")
                     }
                     return(TRUE)
                 })


#' Funtion to get summary information from a \code{mmre} object
#' @param data a \code{mmre} object
mmre.summary <- function(data){
    n_wh <- length(table(data$ID))
    wh <- names(table(data$ID))
    if(n_wh > 1){
        individual <- sapply(lapply(split(data$state,data$ID),function(x) matrix(x,nrow = 1)),
                                   function(x) trans.matrix(x,prob = TRUE))
        trans <- individual
        }else{
            all <- matrix(data$state,nrow = 1)
            trans <- trans.matrix(all, prob = TRUE)
        }
    names(trans) <- wh
    sum <- list(n_wh, wh, trans)
    names(sum) <- c("Number of animals", "IDs", "Raw transition propabilities")
    return(sum)
}

#' Function to combine argumets and data into a \code{mmre} object for using \code{fit.mmre}
#' @param data a \code{mmre} object
#' @param parameters a named list of starting values for the model, assumed model is chosen based on the contents of
#' this list (i.e., to fit a correlated random effect model \code{parameters} must contain the elements \code{u},
#' \code{log_sig_u} and \code{cov_par})
#' @param cov.names a charater string or vector of covariate names (elements in \code{data}) (only for multiple individuals)
#' @param decay logical decay formulation to be used or not
get.mmre.data <- function(data, parameters, cov.names = "none",decay = FALSE, truncation = NULL){
    obj <- mmre(data = data)
    obj@summary <-  mmre.summary(obj@data)
    obj@parameters <- parameters
    obj@cov_names <- cov.names
    obj@decay <- decay
    n <- obj@summary$Number
    if(n > 1) {
        states <- split(obj@data$state,obj@data$ID)
        times <- split(obj@data$time,obj@data$ID)
        ID <- as.factor(levels(obj@data$ID))
        obj@fit_data$response <- list(ID = ID,states = states, times = times)
        
    }else{
        states <- obj@data$state
        times <- obj@data$time
        obj@fit_data$response <- list(states = states, times = times)
    }
    if(!("none"%in%obj@cov_names)){
        covs <- cbind(obj@data[cov.names])
        if(n > 1) {
            covs <- lapply(split(covs,obj@data$ID),as.matrix,ncol = ncol(covs))
        }
        obj@fit_data$response$covariates <- covs
    }
    if(!is.null(truncation)){obj@fit_data$response$truncation <- truncation}
    return(obj)
}

#' Function that choses C++ template to use for model
#' @param mmre.data a \code{mmre} object
mmre.mod <- function(mmre.data = NULL){
    par_names <- names(mmre.data@parameters)
    RE <- ifelse("u"%in% par_names,TRUE, FALSE)
    cov <- ifelse(!("none"%in%mmre.data@cov_names),TRUE,FALSE)
    betmat <- ifelse("betas_matrix"%in%par_names,TRUE,FALSE)
    fixed <- !mmre.data@decay
    trunc <- !is.null(mmre.data@fit_data$response$truncation)
    if(!cov & betmat){stop("Please include covariate names to match parameter matrix")}
    if(cov & !betmat){stop("Please include parameter matrix of covariate statring values")}
    if(!cov & !betmat & trunc){stop("Cannot truncate a missing covariate")}
    if(trunc & !fixed){stop("Please shoose either decay or truncated")}
    if(!RE & !cov  & !betmat & !trunc){
        mod_chosen <- "multiple_continuous"
    }
    if(RE & !cov & !betmat & !trunc){
        mod_chosen <- "multiple_continuous_independentRE"
    }
    if(RE & cov & betmat & !fixed & !trunc){
        mod_chosen <- "multiple_continuous_offdiag_decay_independentRE"
    }
    if(!RE & cov & betmat & !fixed & !trunc){
        mod_chosen <- "multiple_continuous_offdiag_decay"
    }
    if(!RE & cov & betmat & fixed & !trunc){
        mod_chosen <- "multiple_continuous_offdiag"
    }
    if(RE & cov & betmat & fixed & !trunc){
        mod_chosen <- "multiple_continuous_offdiag_independentRE"
    }
    if(!RE & cov & betmat & fixed & trunc){
        mod_chosen <- "multiple_continuous_offdiag_truncated"
    }
    if(RE & cov & betmat & fixed & trunc){
        mod_chosen <- "multiple_continuous_offdiag_truncated_independentRE"
    }
    return(mod_chosen)
}

#' Function to fit two state Markov model
#' @inheritParams get.mmre.data
#' @export
fit.mmre <- function(data = NULL, parameters, cov.names = "none", decay = FALSE, truncation = NULL){
    get_mmre <- get.mmre.data(data = data, parameters = parameters, cov.names = cov.names, decay = decay, truncation = truncation)
    mmre_mod <- mmre.mod(get_mmre)
    print(mmre_mod)
    get_mmre@fitted <- mmre_mod
    data <- get_mmre@fit_data$response 
    params <- get_mmre@parameters
    obj <- MakeADFun(data, params, DLL = get_mmre@fitted)
    opt <- optim(obj$par, obj$fn, gr = obj$gr)
    res <- sdreport(obj)
    get_mmre@fit <- obj
    get_mmre@sdreport <- summary(res)
    return(get_mmre)
}
   





