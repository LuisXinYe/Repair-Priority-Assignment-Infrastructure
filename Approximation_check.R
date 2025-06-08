library(parallel)

cur_addr <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(cur_addr)

epsilon <- 0.001

NumScenarios <- 200
NumReplications <- 5000
NumBaseStation <- 10
NumIteration <- 50000
#Simulation window parameters
NumState <- 0
for (i in 0:NumBaseStation) {
  NumState <- NumState+(NumBaseStation-i+1)*choose(NumBaseStation,i)
}
# Definition of states
State <- list("vector", length=NumState)
for (i in 0:NumBaseStation) {
  temp_vec <- rep(0, NumBaseStation+1)
  
  temp_vec[1] <- i
  if (i == 0) {
    remaining_combs <- expand.grid(rep(list(c(0, 1)), NumBaseStation))
    row_idle <- nrow(remaining_combs)
    for (j in 1:row_idle) {
      tmp_df <- as.vector(t(remaining_combs[j, ]))
      State[[j]] <- c(0, tmp_df)
    }
  }
  else {
    remaining_combs <- expand.grid(rep(list(c(0, 1)), NumBaseStation-1))
    row_serve <- nrow(remaining_combs)
    for (j in 1:row_serve) {
      tmp_df <- as.vector(t(remaining_combs[j, ]))
      if (i == 1) State[[row_idle+j]] <- c(i, 1, tmp_df)
      else if (i == NumBaseStation) State[[row_idle+(i-1)*row_serve+j]] <- c(i, tmp_df, 1)
      else State[[row_idle+(i-1)*row_serve+j]] <- c(i, tmp_df[1:(i-1)], 1, tmp_df[i:length(tmp_df)])
    }
  }
}

numberOfCores <- parallel::detectCores()
cat("numberOfCores: ", numberOfCores, "\n")
cl <- makeCluster(numberOfCores-1)

clusterExport(cl, c("NumBaseStation", "NumScenarios", "NumReplications", "State", "NumIteration", "NumState", "epsilon"))
clusterEvalQ(cl, library(Rsolnp))

allResult <- parSapply(cl, 1:NumScenarios, function(j){
  
  xMin=0
  xMax=12
  yMin=0
  yMax=12
  xDelta=xMax-xMin
  yDelta=yMax-yMin #rectangle dimensions
  
  #Point process parameters
  
  #Simulate Poisson point process
  numbPoints=NumBaseStation # Poisson number of points
  xx=xDelta*runif(numbPoints)+xMin # x coordinates of Poisson points
  yy=yDelta*runif(numbPoints)+yMin # y coordinates of Poisson points
  depot_xx=6
  depot_yy=6
  
  #travel_speed <- 2/3 # unit: km/min
  travel_speed <- 1/2 # unit: km/min
  on_scene_service_time <- 75 # unit: min
  
  euclidean_distance <- vector("numeric", length = NumBaseStation)
  for (i in 1:NumBaseStation) {
    p1 <- c(xx[i],yy[i])
    depot <- c(depot_xx, depot_yy)
    euclidean_distance[i] <- sqrt((depot[1] - p1[1])^2 + (depot[2] - p1[2])^2) # unit: km
  }
  travel_time <- euclidean_distance/travel_speed # unit: min
  service_time <- travel_time + on_scene_service_time # unit: min
  mu <- 1/service_time*60 # service rate; unit: 1/hr
  
  # optimal: 1
  # heuristic: 2
  lambda <- runif(NumBaseStation,0.004,0.02) # arrival rate; unit: 1/hr
  cost <- runif(NumBaseStation,50,500) # unit: CNY per hour
  Lambda <- sum(lambda+mu)
  rho <- 0.0239 # real interest rate
  alpha <- 0.0045
  initial_state <- c(0, rep(0, NumBaseStation)) # initial state
  
  ValueFuncTmp <- vector("list", length = NumIteration)
  for (i in 1:NumIteration) {
    
    if (i > 1) ValueFunc <- ValueFuncTmp[[i-1]]
    else ValueFunc <- rep(0, NumState) # initial values
    
    for (j in 1:NumState) {
      tmp_value <- 1/(Lambda+alpha)*sum(cost*State[[j]][-1])
      lambdaSum <- sum(lambda[which(State[[j]][-1] == 0)])
      #cat(lambdaSum,"")
      
      ActionIndex <- State[[j]][1]
      
      #cat(ActionIndex,"")
      
      if (ActionIndex > 0) {
        # the fourth term
        tmp_value <- tmp_value+1/(Lambda+alpha)*(Lambda-mu[ActionIndex]-lambdaSum)*ValueFunc[j]
        
        if (length(which(State[[j]][-1]==0))>0) {
          
          # find the index of the second term
          NewState2 <- State[[j]]
          NewState2[1] <- 0
          NewState2[ActionIndex+1] <- 0
          if (length(which(NewState2[-1] == 1)) == 0) Index2 <- 1
          else {
            max_index <- max(which(NewState2[-1] == 1))
            for (m in (1+2^(max_index-1)):2^NumBaseStation) {
              if (identical(State[[m]], NewState2)) {
                Index2 <- m
                break
              }
            }
          }
          tmp_value <- tmp_value+1/(Lambda+alpha)*mu[ActionIndex]*ValueFunc[Index2]
          
          # find the index of the third term
          Index3 <- rep(0,length(which(State[[j]][-1]==0)))
          FailureIndex <- which(State[[j]][-1]==0)+1
          NewState3 <- list("numeric", length(which(State[[j]][-1]==0)))
          for (m in 1:length(which(State[[j]][-1]==0))) {
            NewState3[[m]] <- State[[j]]
            NewState3[[m]][FailureIndex[m]] <- 1
            
            max_index <- max(which(NewState3[[m]][-c(1,State[[j]][1]+1)] == 1))
            for (n in ((NewState3[[m]][1]+1)*2^(NumBaseStation-1)+2^(max_index-1)+1):((NewState3[[m]][1]+2)*2^(NumBaseStation-1))) {
              if (identical(State[[n]], NewState3[[m]])) {
                Index3[m] <- n
                break
              }
            }
            tmp_value <- tmp_value+1/(Lambda+alpha)*lambda[FailureIndex[m]-1]*ValueFunc[Index3[m]]
          }
        }
        else {
          # find the index of the second term
          NewState2 <- State[[j]]
          NewState2[1] <- 0
          NewState2[ActionIndex+1] <- 0
          
          max_index <- max(which(NewState2[-1] == 1))
          for (m in (1+2^(max_index-1)):2^NumBaseStation) {
            if (identical(State[[m]], NewState2)) {
              Index2 <- m
              break
            }
          }
          tmp_value <- tmp_value+1/(Lambda+alpha)*mu[ActionIndex]*ValueFunc[Index2]
        }
      }
      else {
        AvailableAction <- which(State[[j]][-1] == 1)
        AvailableAction <- append(AvailableAction,0)
        TmpTerm <- rep(0,length(AvailableAction))
        #cat("Number of available actions:",length(AvailableAction))
        for (k in 1:length(AvailableAction)) {
          if (AvailableAction[k]>0) {
            # the second term
            NewState2 <- State[[j]]
            NewState2[AvailableAction[k]+1] <- 0
            if (length(which(NewState2[-1] == 1)) == 0) Index2 <- 1
            else {
              max_index <- max(which(NewState2[-1] == 1))
              for (m in (1+2^(max_index-1)):2^NumBaseStation) {
                if (identical(State[[m]], NewState2)) {
                  Index2 <- m
                  break
                }
              }
            }
            TmpTerm[k] <- 1/(Lambda+alpha)*mu[AvailableAction[k]]*ValueFunc[Index2]
            
            # the third term
            if (length(which(State[[j]][-1]==0))>0) {
              Index3 <- rep(0,length(which(State[[j]][-1]==0)))
              FailureIndex3 <- which(State[[j]][-1]==0)+1
              NewState3 <- list("numeric", length(which(State[[j]][-1]==0)))
              for (m in 1:length(which(State[[j]][-1]==0))) {
                NewState3[[m]] <- State[[j]]
                NewState3[[m]][1] <- AvailableAction[k]
                NewState3[[m]][FailureIndex3[m]] <- 1
                
                max_index <- max(which(NewState3[[m]][-c(1,AvailableAction[k]+1)] == 1))
                for (n in ((NewState3[[m]][1]+1)*2^(NumBaseStation-1)+2^(max_index-1)+1):((NewState3[[m]][1]+2)*2^(NumBaseStation-1))) {
                  if (identical(State[[n]], NewState3[[m]])) {
                    Index3[m] <- n
                    break
                  }
                }
                TmpTerm[k] <- TmpTerm[k]+1/(Lambda+alpha)*lambda[FailureIndex3[m]-1]*ValueFunc[Index3[m]]
              }
            }
            
            # the fourth term
            NewState4 <- State[[j]]
            NewState4[1] <- AvailableAction[k]
            
            if (length(which(NewState4[-c(1,AvailableAction[k]+1)] == 1)) == 0) Index4 <- (NewState4[1]+1)*2^(NumBaseStation-1)+1
            else {
              max_index <- max(which(NewState4[-c(1,AvailableAction[k]+1)] == 1))
              for (n in ((NewState4[1]+1)*2^(NumBaseStation-1)+2^(max_index-1)+1):((NewState4[1]+2)*2^(NumBaseStation-1))) {
                if (identical(State[[n]], NewState4)) {
                  Index4 <- n
                  break
                }
              }
            }
            TmpTerm[k] <- TmpTerm[k]+1/(Lambda+alpha)*(Lambda-mu[AvailableAction[k]]-lambdaSum)*ValueFunc[Index4]
          }
          else {
            # the third term
            if (length(which(State[[j]][-1]==0))>0) {
              Index3 <- rep(0,length(which(State[[j]][-1]==0)))
              FailureIndex3 <- which(State[[j]][-1]==0)+1
              NewState3 <- list("numeric", length(which(State[[j]][-1]==0)))
              for (m in 1:length(which(State[[j]][-1]==0))) {
                NewState3[[m]] <- State[[j]]
                NewState3[[m]][1] <- 0
                NewState3[[m]][FailureIndex3[m]] <- 1
                
                max_index <- max(which(NewState3[[m]][-1] == 1))
                for (n in (1+2^(max_index-1)):2^NumBaseStation) {
                  if (identical(State[[n]], NewState3[[m]])) {
                    Index3[m] <- n
                    break
                  }
                }
                TmpTerm[k] <- TmpTerm[k]+1/(Lambda+alpha)*lambda[FailureIndex3[m]-1]*ValueFunc[Index3[m]]
              }
            }
            # the fourth term
            TmpTerm[k] <- TmpTerm[k]+1/(Lambda+alpha)*(Lambda-lambdaSum)*ValueFunc[j]
          }
        }
        tmp_value <- tmp_value+min(TmpTerm)
      }
      ValueFuncTmp[[i]][j] <- tmp_value
    }
    if (i ==1) next
    else {
      if (max(ValueFuncTmp[[i]]-ValueFuncTmp[[i-1]]) < epsilon) break
      else next
    }
  }
  
  
  result_policy <-c(0,0)
  for (l in 1:2) {
    index_policy <- l # indicate the applied policy
    total_cost <- 0  # total interruption cost
    for (h in 1:NumReplications) {
      time <- 0       # time
      state <- initial_state  # system state
      
      ObservationWindow <- 24*365
      
      #fail_time <- rep(ObservationWindow, length = NumBaseStation) # record the failure time
      fail_time <- runif(NumBaseStation, 0, ObservationWindow) # record the failure time
      
      # simulate the process
      while (time < ObservationWindow) {  # set stop time of the simulation
        failed_set <- which(state[-1] == 1) # which base station is failed
        normal_set <- which(state[-1] == 0) # which base station is under normal operation
        next_event_time <- ObservationWindow # initialization of the next event time
        
        next_failure_time <- vector("numeric", length = length(normal_set))
        if (length(failed_set) == 0) {
          # next failure time
          for (i in 1:NumBaseStation) {
            next_failure_time[i] <- rexp(1, lambda[i])
          }
          next_event_time <- min(next_failure_time)
          fail_time[which.min(next_failure_time)] <- time + next_event_time
          state[which.min(next_failure_time)+1] <- 1
          total_cost <- total_cost
          time <- time + next_event_time # update time
        }
        else if (length(normal_set) == 0) {
          # next service completion time
          if (state[1] != 0) {
            next_completion_time <- rexp(1, mu[state[1]])
            next_event_time <- next_completion_time
            state[state[1]+1] <- 0
            fail_time[state[1]] <- ObservationWindow
            state[1] <- 0
            total_cost <- total_cost + sum(cost/alpha*(exp(-alpha*time)-exp(-alpha*(time+next_event_time))))
            time <- time + next_event_time # update time
          }
          else {
            if (index_policy == 1) { # optimal policy
              lambdaSum <- sum(lambda[which(state[-1] == 0)])
              AvailableAction <- which(state[-1] == 1)
              TmpTerm <- rep(0,length(AvailableAction))
              #cat("Number of available actions:",length(AvailableAction))
              for (k in 1:length(AvailableAction)) {
                # the second term
                NewState2 <- state
                NewState2[AvailableAction[k]+1] <- 0
                if (length(which(NewState2[-1] == 1)) == 0) Index2 <- 1
                else {
                  max_index <- max(which(NewState2[-1] == 1))
                  for (m in (1+2^(max_index-1)):2^NumBaseStation) {
                    if (identical(State[[m]], NewState2)) {
                      Index2 <- m
                      break
                    }
                  }
                }
                TmpTerm[k] <- 1/(Lambda+alpha)*mu[AvailableAction[k]]*ValueFunc[Index2]
                
                # the fourth term
                NewState4 <- state
                NewState4[1] <- AvailableAction[k]
                
                if (length(which(NewState4[-c(1,AvailableAction[k]+1)] == 1)) == 0) Index4 <- (NewState4[1]+1)*2^(NumBaseStation-1)+1
                else {
                  max_index <- max(which(NewState4[-c(1,AvailableAction[k]+1)] == 1))
                  for (n in ((NewState4[1]+1)*2^(NumBaseStation-1)+2^(max_index-1)+1):((NewState4[1]+2)*2^(NumBaseStation-1))) {
                    if (identical(State[[n]], NewState4)) {
                      Index4 <- n
                      break
                    }
                  }
                }
                TmpTerm[k] <- TmpTerm[k]+1/(Lambda+alpha)*(Lambda-mu[AvailableAction[k]]-lambdaSum)*ValueFunc[Index4]
              }
              optim_action <- which.min(TmpTerm)
              next_completion_time <- rexp(1,mu[AvailableAction[optim_action]])
              next_event_time <- next_completion_time
              total_cost <- total_cost + sum(cost*state[-1]/alpha*(exp(-alpha*time)-exp(-alpha*(time+next_event_time))))
              state[AvailableAction[optim_action]+1] <- 0
              state[1] <- 0
              time <- time + next_event_time # update time}]
            }
            if (index_policy == 2) { # heuristic
              if (which.max(cost*mu/lambda) == which.max(cost*mu)) {
                next_event_time <- rexp(1, mu[which.max(cost*mu/lambda)])
                state[1] <- 0
                state[which.max(cost*mu)+1] <- 0
                total_cost <- total_cost + sum(cost/alpha*(exp(-alpha*time)-exp(-alpha*(time+next_event_time))))
                time <- time + next_event_time # update time
              }
              else {
                if ((cost*mu)[which.max(cost*mu/lambda)]>(1-(mu[which.max(cost*mu)]-mu[which.max(cost*mu/lambda)])/(Lambda+alpha))*max(cost*mu)) {
                  next_event_time <- rexp(1, mu[which.max(cost*mu/lambda)])
                  state[1] <- 0
                  state[which.max(cost*mu/lambda)+1] <- 0
                  total_cost <- total_cost + sum(cost/alpha*(exp(-alpha*time)-exp(-alpha*(time+next_event_time))))
                  time <- time + next_event_time # update time
                }
                else {
                  next_event_time <- rexp(1, mu[which.max(cost*mu)])
                  state[1] <- 0
                  state[which.max(cost*mu)+1] <- 0
                  total_cost <- total_cost + sum(cost/alpha*(exp(-alpha*time)-exp(-alpha*(time+next_event_time))))
                  time <- time + next_event_time # update time
                }
              }
            }
          }
        }
        else {
          for (i in 1:length(normal_set)) {
            next_failure_time[i] <- rexp(1, lambda[normal_set[i]])
          }
          if (state[1] != 0) {
            next_completion_time <- rexp(1, mu[state[1]])
            if (next_completion_time < min(next_failure_time)) {
              next_event_time <- min(next_failure_time, next_completion_time)
              total_cost <- total_cost + sum(cost*state[-1]/alpha*(exp(-alpha*time)-exp(-alpha*(time+next_event_time))))
              state[state[1]+1] <- 0
              fail_time[state[1]] <- ObservationWindow
              state[1] <- 0
              time <- time + next_event_time # update time}
            }
            else {
              next_event_time <- min(next_failure_time)
              fail_time[normal_set[which.min(next_failure_time)]] <- time + next_event_time
              total_cost <- total_cost + sum(cost*state[-1]/alpha*(exp(-alpha*time)-exp(-alpha*(time+next_event_time))))
              state[normal_set[which.min(next_failure_time)]+1] <- 1
              time <- time + next_event_time # update time}
            }
          }
          else {
            if (index_policy == 1) { # optimal policy
              lambdaSum <- sum(lambda[which(state[-1] == 0)])
              
              AvailableAction <- which(state[-1] == 1)
              AvailableAction <- append(AvailableAction,0)
              TmpTerm <- rep(0,length(AvailableAction))
              #cat("Number of available actions:",length(AvailableAction))
              for (k in 1:length(AvailableAction)) {
                if (AvailableAction[k]>0) {
                  # the second term
                  NewState2 <- state
                  NewState2[AvailableAction[k]+1] <- 0
                  if (length(which(NewState2[-1] == 1)) == 0) Index2 <- 1
                  else {
                    max_index <- max(which(NewState2[-1] == 1))
                    for (m in (1+2^(max_index-1)):2^NumBaseStation) {
                      if (identical(State[[m]], NewState2)) {
                        Index2 <- m
                        break
                      }
                    }
                  }
                  TmpTerm[k] <- 1/(Lambda+alpha)*mu[AvailableAction[k]]*ValueFunc[Index2]
                  
                  # the third term
                  if (length(which(state[-1]==0))>0) {
                    Index3 <- rep(0,length(which(state[-1]==0)))
                    FailureIndex3 <- which(state[-1]==0)+1
                    NewState3 <- list("numeric", length(which(state[-1]==0)))
                    for (m in 1:length(which(state[-1]==0))) {
                      NewState3[[m]] <- state
                      NewState3[[m]][1] <- AvailableAction[k]
                      NewState3[[m]][FailureIndex3[m]] <- 1
                      
                      max_index <- max(which(NewState3[[m]][-c(1,AvailableAction[k]+1)] == 1))
                      for (n in ((NewState3[[m]][1]+1)*2^(NumBaseStation-1)+2^(max_index-1)+1):((NewState3[[m]][1]+2)*2^(NumBaseStation-1))) {
                        if (identical(State[[n]], NewState3[[m]])) {
                          Index3[m] <- n
                          break
                        }
                      }
                      TmpTerm[k] <- TmpTerm[k]+1/(Lambda+alpha)*lambda[FailureIndex3[m]-1]*ValueFunc[Index3[m]]
                    }
                  }
                  
                  # the fourth term
                  NewState4 <- state
                  NewState4[1] <- AvailableAction[k]
                  
                  if (length(which(NewState4[-c(1,AvailableAction[k]+1)] == 1)) == 0) Index4 <- (NewState4[1]+1)*2^(NumBaseStation-1)+1
                  else {
                    max_index <- max(which(NewState4[-c(1,AvailableAction[k]+1)] == 1))
                    for (n in ((NewState4[1]+1)*2^(NumBaseStation-1)+2^(max_index-1)+1):((NewState4[1]+2)*2^(NumBaseStation-1))) {
                      if (identical(State[[n]], NewState4)) {
                        Index4 <- n
                        break
                      }
                    }
                  }
                  TmpTerm[k] <- TmpTerm[k]+1/(Lambda+alpha)*(Lambda-mu[AvailableAction[k]]-lambdaSum)*ValueFunc[Index4]
                }
                else {
                  # the third term
                  if (length(which(state[-1]==0))>0) {
                    Index3 <- rep(0,length(which(state[-1]==0)))
                    FailureIndex3 <- which(state[-1]==0)+1
                    NewState3 <- list("numeric", length(which(state[-1]==0)))
                    for (m in 1:length(which(state[-1]==0))) {
                      NewState3[[m]] <- state
                      NewState3[[m]][1] <- 0
                      NewState3[[m]][FailureIndex3[m]] <- 1
                      
                      max_index <- max(which(NewState3[[m]][-1] == 1))
                      for (n in (1+2^(max_index-1)):2^NumBaseStation) {
                        if (identical(State[[n]], NewState3[[m]])) {
                          Index3[m] <- n
                          break
                        }
                      }
                      TmpTerm[k] <- TmpTerm[k]+1/(Lambda+alpha)*lambda[FailureIndex3[m]-1]*ValueFunc[Index3[m]]
                    }
                  }
                  # the fourth term
                  TmpTerm[k] <- TmpTerm[k]+1/(Lambda+alpha)*(Lambda-lambdaSum)*ValueFunc[j]
                }
              }
              optim_action <- which.min(TmpTerm)
              if (optim_action == 0) {
                next_event_time <- min(next_failure_time)
                state[normal_set[which.min(next_failure_time)]+1] <- 1
                total_cost <- total_cost + sum(cost*state[-1]/alpha*(exp(-alpha*time)-exp(-alpha*(time+next_event_time))))
                time <- time + next_event_time # update time}
              } else {
                next_completion_time <- rexp(1, mu[AvailableAction[optim_action]])
                next_event_time <- min(next_failure_time, next_completion_time)
                total_cost <- total_cost + sum(cost*state[-1]/alpha*(exp(-alpha*time)-exp(-alpha*(time+next_event_time))))
                if (next_completion_time < min(next_failure_time)) { # completion
                  state[AvailableAction[optim_action]+1] <- 0
                  state[1] <- 0
                  time <- time + next_event_time # update time}
                }
                else { # an new failure arrival
                  state[1] <- AvailableAction[optim_action]
                  state[normal_set[which.min(next_failure_time)]+1] <- 1
                  time <- time + next_event_time # update time}
                }
              }
            }
            if (index_policy == 2) { # heuristic
              if (which.max((cost*mu/lambda)[failed_set]) == which.max((cost*mu)[failed_set])) {
                next_completion_time <- rexp(1, mu[failed_set[which.max((cost*mu)[failed_set])]])
                if (next_completion_time < min(next_failure_time)) { # completion
                  next_event_time <- next_completion_time
                  total_cost <- total_cost + sum(cost*state[-1]/alpha*(exp(-alpha*time)-exp(-alpha*(time+next_event_time))))
                  state[failed_set[which.max((cost*mu)[failed_set])]+1] <- 0
                  time <- time + next_event_time # update time}
                }
                else { # an new failure arrival
                  next_event_time <- min(next_failure_time)
                  total_cost <- total_cost + sum(cost*state[-1]/alpha*(exp(-alpha*time)-exp(-alpha*(time+next_event_time))))
                  state[1] <- failed_set[which.max((cost*mu)[failed_set])]
                  state[normal_set[which.min(next_failure_time)]+1] <- 1
                  time <- time + next_event_time # update time}
                }
              }
              else {
                if ((cost*mu)[which.max(cost*mu/lambda)]>(1-(mu[which.max(cost*mu)]-mu[which.max(cost*mu/lambda)])/(Lambda+alpha))*max(cost*mu)) {
                  next_completion_time <- rexp(1, mu[failed_set[which.max((cost*mu/lambda)[failed_set])]])
                  if (next_completion_time < min(next_failure_time)) { # completion
                    next_event_time <- next_completion_time
                    total_cost <- total_cost + sum(cost*state[-1]/alpha*(exp(-alpha*time)-exp(-alpha*(time+next_event_time))))
                    state[failed_set[which.max((cost*mu/lambda)[failed_set])]+1] <- 0
                    time <- time + next_event_time # update time}
                  }
                  else { # an new failure arrival
                    next_event_time <- min(next_failure_time)
                    total_cost <- total_cost + sum(cost*state[-1]/alpha*(exp(-alpha*time)-exp(-alpha*(time+next_event_time))))
                    state[1] <- failed_set[which.max((cost*mu/lambda)[failed_set])]
                    state[normal_set[which.min(next_failure_time)]+1] <- 1
                    time <- time + next_event_time # update time}
                  }
                }
                else {
                  next_completion_time <- rexp(1, mu[failed_set[which.max((cost*mu)[failed_set])]])
                  if (next_completion_time < min(next_failure_time)) { # completion
                    next_event_time <- next_completion_time
                    total_cost <- total_cost + sum(cost*state[-1]/alpha*(exp(-alpha*time)-exp(-alpha*(time+next_event_time))))
                    state[failed_set[which.max((cost*mu)[failed_set])]+1] <- 0
                    time <- time + next_event_time # update time}
                  }
                  else { # an new failure arrival
                    next_event_time <- min(next_failure_time)
                    total_cost <- total_cost + sum(cost*state[-1]/alpha*(exp(-alpha*time)-exp(-alpha*(time+next_event_time))))
                    state[1] <- failed_set[which.max((cost*mu)[failed_set])]
                    state[normal_set[which.min(next_failure_time)]+1] <- 1
                    time <- time + next_event_time # update time}
                  }
                }
              }
            }
          }
        }
      }
    }
    result_policy[l] <- total_cost/NumReplications
  }
  result_policy
})

save(list = c("allResult"), file = paste0("Approx_", NumBaseStation, ".Rdata"))

stopCluster(cl)


