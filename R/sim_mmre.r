#' Function to simulate Matrkov transition model where transition intensities are
#' given by beta0 + beta1*exp(beta2*cov) + individual.random.effect
#' Number simulated will correspond to length of given times and covariates.
#' @param par.sim a 2 by three matrix of parameters. Each row coresponds to the off diagonal transition intensities.
#' Each column corresponds to the covariates.
#' @param times numeric vector of observed times
#' @param ID a character vector of whale Ids for zero meaned mvn independent random effects
#' @param start.state numeric initial state (can only be 1 or 2)
#' @param cov a vector of covariate values (i.e., days since exposure)
#' @param random logical, should individual random effects be included.
#' @param fit a mmre fit object, either supplied or all other armugents
#' @export
setGeneric("sim.mmre.decay",
           function(par.sim , cov, times, ID, start.state, random ,fit){
               standardGeneric("sim.mmre.decay")
           })

setMethod("sim.mmre.decay",
          c("matrix","numeric","numeric","character", "numeric","logical","missing"),
          function(par.sim , cov, times, ID, start.state, random, fit){
              n <- length(table(ID))
              ids <- names(table(ID))
              if(random){
                  mvn <- mvtnorm::rmvnorm(n,mean = rep(0,2), sigma = diag(rep(1,2)))
                  ## zero mean mvn with fixed sd
                  }
              states <- list()
              for(i in 1:n){
                  tms <- times[ID==ids[i]]
                  covs <- cov[ID==ids[i]]
                  start.state <- start.state ## at the moment same start state for each whale
                  state <- numeric(length(tms)) ## initialise state vector
                  state[1] <- start.state
                  for(j in 2:length(state)){
                      t <- tms[j] - tms[j-1]
                      c <- covs[j]
                      if(random){
                          ## off diag q s1 to s2
                          q_1.2 <- exp(par.sim[1,1] +
                                       par.sim[1,2]*exp(-exp(par.sim[1,3])*c) + mvn[i,1])
                          ## off diag q s2 to s1
                          q_2.1 <-  exp(par.sim[2,1] +
                                        par.sim[2,2]*exp(-exp(par.sim[2,3])*c) + mvn[i,2])
                      }else{
                          ## off diag q s1 to s2
                          q_1.2 <- exp(par.sim[1,1] + par.sim[1,2]*exp(-exp(par.sim[1,3])*c) )
                          ## off diag q s2 to s1
                          q_2.1 <-  exp(par.sim[2,1] + par.sim[2,2]*exp(-exp(par.sim[2,3])*c) )
                      }
                      ## Q matrix
                      Q <- matrix(c(-q_1.2,q_1.2,q_2.1,-q_2.1),nrow = 2,byrow = TRUE)
                      Qt <- Q*t
                      P <- expm::expm(Qt)
                      state[j] <- sample(c(1,2),1, prob = P[state[j - 1],] )
                  }
                  states[[i]] <- state
              }
              states <- unlist(states)
              res <- data.frame(times = times ,state = states, ID = ID,cov = cov)
              names(res) <- c("time", "state","ID","cov")
              if(random){
                  return(list(sim = res, mvn = mvn))
              }else{
                  return(list(sim = res))
              }
              
          }
          )
setMethod("sim.mmre.decay",
         c("missing","missing","missing","missing","missing","missing","mmre"),
          function(par.sim , cov, times, ID, start.state, fit){
              random <- "u"%in%names(fit@parameters)
              sim <- sim.mmre.decay(par.sim = get.params(fit,FALSE),
                                    cov = fit@data[,fit@cov_names], times = fit@data$time,
                                    ID = as.character(fit@data$ID) , start.state = fit@data$state[1],
                                    random = random)
              return(sim)
          }
          )
#' Function to simulate a two state Markov transition log-linear model
#' Number simulated will correspond to length of given times and covariate.
#' @inheritParams sim.mmre.decay
#' @export
setGeneric("sim.mmre.offdiag",
           function(par.sim , cov, times, ID, start.state, random ,fit){
               standardGeneric("sim.mmre.offdiag")
           })

setMethod("sim.mmre.offdiag",
          c("matrix","numeric","numeric","character", "numeric","logical","missing"),
          function(par.sim , cov, times, ID, start.state, random, fit){
              n <- length(table(ID))
              ids <- names(table(ID))
              if(random){
                  mvn <- mvtnorm::rmvnorm(n,mean = rep(0,2), sigma = diag(rep(1,2)))
                  ## zero mean mvn with fixed sd
                  }
              states <- list()
              for(i in 1:n){
                  tms <- times[ID==ids[i]]
                  covs <- cov[ID==ids[i]]
                  start.state <- start.state ## at the moment same start state for each whale
                  state <- numeric(length(tms)) ## initialise state vector
                  state[1] <- start.state
                  n.c <- ncol(par.sim)
                  for(j in 2:length(state)){
                      t <- tms[j] - tms[j-1]
                      c <- covs[j]
                      if(random){
                          ## off diag q s1 to s2
                          q_1.2 <- exp(par.sim[1,1] +
                                        c*par.sim[1,(2:n.c)] + mvn[i,1])
                          ## off diag q s2 to s1
                          q_2.1 <-  exp(par.sim[2,1] +
                                        c*par.sim[2,(2:n.c)] + mvn[i,2])
                      }else{
                          ## off diag q s1 to s2
                          q_1.2 <- exp(par.sim[1,1] + c*par.sim[1,(2:n.c)])
                          ## off diag q s2 to s1
                          q_2.1 <- exp(par.sim[2,1] + c*par.sim[2,(2:n.c)])
                      }
                      ## Q matrix
                      Q <- matrix(c(-q_1.2,q_1.2,q_2.1,-q_2.1),nrow = 2,byrow = TRUE)
                      Qt <- Q*t
                      P <- expm::expm(Qt)
                      state[j] <- sample(c(1,2),1, prob = P[state[j - 1],] )
                  }
                  states[[i]] <- state
              }
              states <- unlist(states)
              res <- data.frame(times = times ,state = states, ID = ID,cov = cov)
              names(res) <- c("time", "state","ID","cov")
              if(random){
                  return(list(sim = res, mvn = mvn))
              }else{
                  return(list(sim = res))
              }
              
          }
          )
setMethod("sim.mmre.offdiag",
         c("missing","missing","missing","missing","missing","missing","mmre"),
          function(par.sim , cov, times, ID, start.state,random, fit){
              random <- "u"%in%names(fit@parameters)
              sim <- sim.mmre.offdiag(par.sim = get.params(fit,FALSE),
                                    cov = fit@data[,fit@cov_names], times = fit@data$time,
                                    ID = as.character(fit@data$ID) , start.state = fit@data$state[1],
                                    random = random)
              
              return(sim)
          }
          )

#' Function to simulate a two state Markov transition baseline model
#' Number simulated will correspond to length of given times and covariate.
#' @inheritParams sim.mmre.decay
#' @export
setGeneric("sim.mmre",
           function(par.sim, times, ID, start.state, random ,fit){
               standardGeneric("sim.mmre")
           })

setMethod("sim.mmre",
          c("matrix","numeric","character", "numeric","logical","missing"),
          function(par.sim, times, ID, start.state, random, fit){
              n <- length(table(ID))
              ids <- names(table(ID))
              if(random){
                  mvn <- mvtnorm::rmvnorm(n,mean = rep(0,2), sigma = diag(rep(1,2)))
                  ## zero mean mvn with fixed sd
                  }
              states <- list()
              for(i in 1:n){
                  tms <- times[ID==ids[i]]
                  start.state <- start.state ## at the moment same start state for each whale
                  state <- numeric(length(tms)) ## initialise state vector
                  state[1] <- start.state
                  for(j in 2:length(state)){
                      t <- tms[j] - tms[j-1]
                      if(random){
                          ## off diag q s1 to s2
                          q_1.2 <- exp(par.sim[1,1] + mvn[i,1])
                          ## off diag q s2 to s1
                          q_2.1 <-  exp(par.sim[2,1] + mvn[i,2])
                      }else{
                          ## off diag q s1 to s2
                          q_1.2 <- exp(par.sim[1,1])
                          ## off diag q s2 to s1
                          q_2.1 <- exp(par.sim[2,1])
                      }
                      ## Q matrix
                      Q <- matrix(c(-q_1.2,q_1.2,q_2.1,-q_2.1),nrow = 2,byrow = TRUE)
                      Qt <- Q*t
                      P <- expm::expm(Qt)
                      state[j] <- sample(c(1,2),1,prob = P[state[j - 1],] )
                  }
                  states[[i]] <- state
              }
              states <- unlist(states)
              res <- data.frame(times = times ,state = states, ID = ID)
              names(res) <- c("time", "state","ID")
              if(random){
                  return(list(sim = res, mvn = mvn))
              }else{
                  return(list(sim = res))
              }
              
          }
          )
setMethod("sim.mmre",
         c("missing","missing","missing","missing","missing","mmre"),
          function(par.sim , times, ID, start.state,random, fit){
              random <- "u"%in%names(fit@parameters)
              sim <- sim.mmre(par.sim = as.matrix(fit@sdreport[1:2,1]),
                                    times = fit@data$time,
                                    ID = as.character(fit@data$ID) , start.state = fit@data$state[1],
                                    random = random)
              
              return(sim)
          }
          )

#' Function to simulate a two state Markov transition truncated log-linear model
#' Number simulated will correspond to length of given times and covariate.
#' @inheritParams sim.mmre.decay
#' @param trunc numeric truncation at which covariate has no effect
#' @export
setGeneric("sim.mmre.trunc",
           function(par.sim , cov, times, ID, start.state, random, trunc,fit){
               standardGeneric("sim.mmre.trunc")
           })

setMethod("sim.mmre.trunc",
          c("matrix","numeric","numeric","character", "numeric","logical","numeric","missing"),
          function(par.sim , cov, times, ID, start.state, random,trunc, fit){
              n <- length(table(ID))
              ids <- names(table(ID))
              if(random){
                  mvn <- mvtnorm::rmvnorm(n,mean = rep(0,2), sigma = diag(rep(1,2)))
                  ## zero mean mvn with fixed sd
                  }
              states <- list()
              for(i in 1:n){
                  tms <- times[ID==ids[i]]
                  covs <- cov[ID==ids[i]]
                  start.state <- start.state ## at the moment same start state for each whale
                  state <- numeric(length(tms)) ## initialise state vector
                  state[1] <- start.state
                  n.c <- ncol(par.sim)
                  for(j in 2:length(state)){
                      t <- tms[j] - tms[j-1]
                      c <- covs[j]
                      if(random){
                          if(c == 0){
                              ## off diag q s1 to s2
                              q_1.2 <- exp(par.sim[1,1] + par.sim[1,2] + mvn[i,1])
                              ## off diag q s2 to s1
                              q_2.1 <-  exp(par.sim[2,1] +
                                            par.sim[2,2] + mvn[i,2])
                          }else{
                              if(c >= trunc){
                                  ## off diag q s1 to s2
                                  q_1.2 <- exp(par.sim[1,1] + mvn[i,1])
                                  ## off diag q s2 to s1
                                  q_2.1 <-  exp(par.sim[2,1] +
                                                mvn[i,2])
                              }else{
                                  ## off diag q s1 to s2
                                  q_1.2 <- exp(par.sim[1,1] + par.sim[1,3]*c + mvn[i,1])
                                  ## off diag q s2 to s1
                                  q_2.1 <-  exp(par.sim[2,1] +
                                                par.sim[2,3]*c + mvn[i,2])
                              }
                          }
                      }else{
                          if(c == 0){
                              ## off diag q s1 to s2
                              q_1.2 <- exp(par.sim[1,1] + par.sim[1,2] )
                              ## off diag q s2 to s1
                              q_2.1 <-  exp(par.sim[2,1] +
                                            par.sim[2,2] )
                          }else{
                              if(c >= trunc){
                                  ## off diag q s1 to s2
                                  q_1.2 <- exp(par.sim[1,1])
                                  ## off diag q s2 to s1
                                  q_2.1 <-  exp(par.sim[2,1])
                              }else{
                                  ## off diag q s1 to s2
                                  q_1.2 <- exp(par.sim[1,1] + par.sim[1,3]*c )
                                  ## off diag q s2 to s1
                                  q_2.1 <-  exp(par.sim[2,1] +
                                                par.sim[2,3]*c )
                              }
                          }
                      }
                      ## Q matrix
                      Q <- matrix(c(-q_1.2,q_1.2,q_2.1,-q_2.1),nrow = 2,byrow = TRUE)
                      Qt <- Q*t
                      P <- expm::expm(Qt)
                      state[j] <- sample(c(1,2),1, prob = P[state[j - 1],] )
                  }
                  states[[i]] <- state
              }
              states <- unlist(states)
              res <- data.frame(times = times ,state = states, ID = ID,cov = cov)
              names(res) <- c("time", "state","ID","cov")
              if(random){
                  return(list(sim = res, mvn = mvn))
              }else{
                  return(list(sim = res))
              }
              
          }
          )
setMethod("sim.mmre.trunc",
         c("missing","missing","missing","missing","missing","missing","missing","mmre"),
          function(par.sim , cov, times, ID, start.state,trunc, fit){
              random <- "u"%in%names(fit@parameters)
              sim <- sim.mmre.trunc(par.sim = get.params(fit,FALSE),
                                    cov = fit@data[,fit@cov_names], times = fit@data$time,
                                    ID = as.character(fit@data$ID) , start.state = fit@data$state[1],
                                    random = random, trunc = fit@fit_data$response$truncation)
              
              return(sim)
          }
          )
