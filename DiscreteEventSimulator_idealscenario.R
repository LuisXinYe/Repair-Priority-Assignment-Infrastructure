require(parallel)
cur_addr <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(cur_addr)

NumScenarios <- 200
NumReplications <- 5000
NumBaseStation <- 5 #{5,10,15,20,25,30}

numberOfCores <- parallel::detectCores()
cat("numberOfCores: ", numberOfCores, "\n")
cl <- makeCluster(numberOfCores-1)

clusterExport(cl, c("NumBaseStation", "NumScenarios", "NumReplications"))
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
  
  # cmu-lambda-rule: 1
  # cmu-rule: 2
  # FCFS: 3
  # nearest-first: 4
  # least-mu: 5
  # heuristic: 6
  lambda <- runif(NumBaseStation,0.004,0.02) # arrival rate; unit: 1/hr
  cost <- runif(NumBaseStation,50,500) # unit: CNY per hour
  Lambda <- sum(lambda+mu)
  rho <- 0.0239 # real interest rate
  alpha <- 0.0045
  initial_state <- c(0, rep(0, NumBaseStation)) # initial state
  
  result_policy <-c(0,0,0,0,0,0)
  for (l in 1:6) {
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
        next_completion_time <- vector("numeric", length = length(failed_set))
        
        if (length(normal_set) > 0) {
          for (i in 1:length(normal_set)) {
            next_failure_time[i] <- rexp(1, lambda[normal_set[i]])
          }
        } else next_failure_time <- 0
        
        if (length(failed_set) > 0) {
          for (i in 1:length(failed_set)) {
            next_completion_time[i] <- rexp(1, mu[failed_set[i]])
          }
        } else next_completion_time <- 0
        
        if (length(failed_set) == 0) {
          next_event_time <- min(next_failure_time)
          fail_time[which.min(next_failure_time)] <- time + next_event_time
          state[which.min(next_failure_time)+1] <- 1
          total_cost <- total_cost
          time <- time + next_event_time # update time
        }
        else if (length(normal_set) == 0) {
          # next service completion time
          if (state[1] != 0) {
            next_completion_time <- next_completion_time[state[1]]
            next_event_time <- next_completion_time
            state[state[1]+1] <- 0
            fail_time[state[1]] <- ObservationWindow
            state[1] <- 0
            total_cost <- total_cost + sum(cost/alpha*(exp(-alpha*time)-exp(-alpha*(time+next_event_time))))
            time <- time + next_event_time # update time
          }
          else {
            if (index_policy == 1) { # cmu-lambda-rule
              next_completion_time <- next_completion_time[which.max(cost*mu/lambda)]
              next_event_time <- next_completion_time
              state[1] <- 0
              state[which.max(cost*mu/lambda)+1] <- 0
              total_cost <- total_cost + sum(cost/alpha*(exp(-alpha*time)-exp(-alpha*(time+next_event_time))))
              time <- time + next_event_time # update time
            }
            if (index_policy == 2) { # cmu-rule
              next_completion_time <- next_completion_time[which.max(cost*mu)]
              next_event_time <- next_completion_time
              state[1] <- 0
              state[which.max(cost*mu)+1] <- 0
              total_cost <- total_cost + sum(cost/alpha*(exp(-alpha*time)-exp(-alpha*(time+next_event_time))))
              time <- time + next_event_time # update time
            }
            if (index_policy == 3) { # FCFS
              next_completion_time <- next_completion_time[which.min(fail_time)]
              next_event_time <- next_completion_time
              state[1] <- 0
              state[which.min(fail_time)+1] <- 0
              fail_time[which.min(fail_time)] <- ObservationWindow
              total_cost <- total_cost + sum(cost/alpha*(exp(-alpha*time)-exp(-alpha*(time+next_event_time))))
              time <- time + next_event_time # update time
            }
            if (index_policy == 4) { # nearest-first
              next_completion_time <- next_completion_time[which.min(euclidean_distance)]
              next_event_time <- next_completion_time
              state[1] <- 0
              state[which.min(euclidean_distance)+1] <- 0
              total_cost <- total_cost + sum(cost/alpha*(exp(-alpha*time)-exp(-alpha*(time+next_event_time))))
              time <- time + next_event_time # update time
            }
            if (index_policy == 5) { # least-mu
              next_completion_time <- rexp(1, mu[which.min(mu)])
              next_event_time <- next_completion_time
              state[1] <- 0
              state[which.min(mu)+1] <- 0
              total_cost <- total_cost + sum(cost/alpha*(exp(-alpha*time)-exp(-alpha*(time+next_event_time))))
              time <- time + next_event_time # update time
            }
            if (index_policy == 6) { # heuristic
              if (which.max(cost*mu/lambda) == which.max(cost*mu)) {
                next_event_time <- next_completion_time[which.max(cost*mu/lambda)]
                state[1] <- 0
                state[which.max(cost*mu)+1] <- 0
                total_cost <- total_cost + sum(cost/alpha*(exp(-alpha*time)-exp(-alpha*(time+next_event_time))))
                time <- time + next_event_time # update time
              }
              else {
                if (mu[which.max(cost*mu)] >= mu[which.max(cost*mu/lambda)]) {
                  next_event_time <- next_completion_time[which.max(cost*mu/lambda)]
                  state[1] <- 0
                  state[which.max(cost*mu/lambda)+1] <- 0
                  total_cost <- total_cost + sum(cost/alpha*(exp(-alpha*time)-exp(-alpha*(time+next_event_time))))
                  time <- time + next_event_time # update time
                }
                else {
                  if ((cost*mu)[which.max(cost*mu/lambda)]>(1-(mu[which.max(cost*mu)]-mu[which.max(cost*mu/lambda)])/(Lambda+alpha))*max(cost*mu)) {
                    next_event_time <- next_completion_time[which.max(cost*mu/lambda)]
                    state[1] <- 0
                    state[which.max(cost*mu/lambda)+1] <- 0
                    total_cost <- total_cost + sum(cost/alpha*(exp(-alpha*time)-exp(-alpha*(time+next_event_time))))
                    time <- time + next_event_time # update time
                  }
                  else {
                    next_event_time <- next_completion_time[which.max(cost*mu)]
                    state[1] <- 0
                    state[which.max(cost*mu)+1] <- 0
                    total_cost <- total_cost + sum(cost/alpha*(exp(-alpha*time)-exp(-alpha*(time+next_event_time))))
                    time <- time + next_event_time # update time
                  }
                }
              }
            }
          }
        }
        else {
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
            if (index_policy == 1) { # cmu-lambda-rule
              next_completion_time <- rexp(1, mu[failed_set[which.max(cost[failed_set]*mu[failed_set]/lambda[failed_set])]])
              if (next_completion_time < min(next_failure_time)) { # completion
                next_event_time <- next_completion_time
                total_cost <- total_cost + sum(cost*state[-1]/alpha*(exp(-alpha*time)-exp(-alpha*(time+next_event_time))))
                state[failed_set[which.max(cost[failed_set]*mu[failed_set]/lambda[failed_set])]+1] <- 0
                time <- time + next_event_time # update time}
              }
              else { # an new failure arrival
                next_event_time <- min(next_failure_time)
                total_cost <- total_cost + sum(cost*state[-1]/alpha*(exp(-alpha*time)-exp(-alpha*(time+next_event_time))))
                state[1] <- failed_set[which.max(cost[failed_set]*mu[failed_set]/lambda[failed_set])]
                state[normal_set[which.min(next_failure_time)]+1] <- 1
                time <- time + next_event_time # update time}
              }
            }
            if (index_policy == 2) { # cmu-rule
              next_completion_time <- rexp(1, mu[failed_set[which.max(cost[failed_set]*mu[failed_set])]])
              if (next_completion_time < min(next_failure_time)) { # completion
                next_event_time <- next_completion_time
                total_cost <- total_cost + sum(cost*state[-1]/alpha*(exp(-alpha*time)-exp(-alpha*(time+next_event_time))))
                state[failed_set[which.max(cost[failed_set]*mu[failed_set])]+1] <- 0
                time <- time + next_event_time # update time}
              }
              else { # an new failure arrival
                next_event_time <- min(next_failure_time)
                total_cost <- total_cost + sum(cost*state[-1]/alpha*(exp(-alpha*time)-exp(-alpha*(time+next_event_time))))
                state[1] <- failed_set[which.max(cost[failed_set]*mu[failed_set])]
                state[normal_set[which.min(next_failure_time)]+1] <- 1
                time <- time + next_event_time # update time}
              }
            }
            if (index_policy == 3) { # FCFS
              next_completion_time <- rexp(1, mu[which.min(fail_time)])
              if (next_completion_time < min(next_failure_time)) { # completion
                next_event_time <- next_completion_time
                total_cost <- total_cost + sum(cost*state[-1]/alpha*(exp(-alpha*time)-exp(-alpha*(time+next_event_time))))
                state[which.min(fail_time)+1] <- 0
                fail_time[which.min(fail_time)] <- ObservationWindow
                time <- time + next_event_time # update time}
              }
              else { # an new failure arrival
                next_event_time <- min(next_failure_time)
                total_cost <- total_cost + sum(cost*state[-1]/alpha*(exp(-alpha*time)-exp(-alpha*(time+next_event_time))))
                state[1] <- which.min(fail_time)
                state[normal_set[which.min(next_failure_time)]+1] <- 1
                fail_time[normal_set[which.min(next_failure_time)]] <- time + next_event_time
                time <- time + next_event_time # update time}
              }
            }
            if (index_policy == 4) { # nearest-first
              next_completion_time <- rexp(1, mu[failed_set[which.min(euclidean_distance[failed_set])]])
              if (next_completion_time < min(next_failure_time)) { # completion
                next_event_time <- next_completion_time
                total_cost <- total_cost + sum(cost*state[-1]/alpha*(exp(-alpha*time)-exp(-alpha*(time+next_event_time))))
                state[failed_set[which.min(euclidean_distance[failed_set])]+1] <- 0
                time <- time + next_event_time # update time}
              }
              else { # an new failure arrival
                next_event_time <- min(next_failure_time)
                total_cost <- total_cost + sum(cost*state[-1]/alpha*(exp(-alpha*time)-exp(-alpha*(time+next_event_time))))
                state[normal_set[which.min(next_failure_time)]+1] <- 1
                state[1] <- failed_set[which.min(euclidean_distance[failed_set])]
                time <- time + next_event_time # update time}
              }
            }
            if (index_policy == 5) { # least-mu
              next_completion_time <- rexp(1, mu[failed_set[which.min(mu[failed_set])]])
              if (next_completion_time < min(next_failure_time)) { # completion
                next_event_time <- next_completion_time
                total_cost <- total_cost + sum(cost*state[-1]/alpha*(exp(-alpha*time)-exp(-alpha*(time+next_event_time))))
                state[failed_set[which.min(mu[failed_set])]+1] <- 0
                time <- time + next_event_time # update time}
              }
              else { # an new failure arrival
                next_event_time <- min(next_failure_time)
                total_cost <- total_cost + sum(cost*state[-1]/alpha*(exp(-alpha*time)-exp(-alpha*(time+next_event_time))))
                state[1] <- failed_set[which.min(mu[failed_set])]
                state[normal_set[which.min(next_failure_time)]+1] <- 1
                time <- time + next_event_time # update time}
              }
            }
            if (index_policy == 6) { # heuristic
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
                if (mu[failed_set[which.max((cost*mu)[failed_set])]] >= mu[failed_set[which.max((cost*mu/lambda)[failed_set])]]) {
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
                  if ((cost*mu)[failed_set[which.max((cost*mu/lambda)[failed_set])]]>(1-(mu[failed_set[which.max((cost*mu)[failed_set])]]-mu[failed_set[which.max((cost*mu/lambda)[failed_set])]])/(Lambda+alpha))*max((cost*mu)[failed_set])) {
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
    }
    result_policy[l] <- total_cost/NumReplications
  }
  
  result_policy
  
})

#save(list = c("allResult"), file = paste0("NumBaseStation", NumBaseStation, ".Rdata"))

stopCluster(cl)


