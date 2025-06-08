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
clusterEvalQ(cl, c(library(statmod)))

allResult <- parSapply(cl, 1:NumScenarios, function(j){
  #Simulation window parameters
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
  
  cv <- 1 # cv is 1
  index_dist <- 1 # indicate the assumed travel time distribution (1: lognormal, 2: inverse-gaussian)
  
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
        if (length(failed_set) == 0) {
          # next failure time
          for (i in 1:NumBaseStation) {
            next_failure_time[i] <- rexp(1, lambda[i])*60
          }
          if (state[1] == 0) {
            next_event_time <- min(next_failure_time)
            state[which.min(next_failure_time)+1] <- 1
            fail_time[which.min(next_failure_time)] <- time + next_event_time
            total_cost <- total_cost
            time <- time + next_event_time # update time
          } else {
            state_tmp <- -state[1]
            if (index_dist == 1) {
              cdf <- plnorm(time-depart_time, meanlog = log(travel_time[state_tmp]^2/sqrt((travel_time[state_tmp]*cv)^2 + travel_time[state_tmp]^2)), 
                            sdlog = sqrt(log(1 + ((travel_time[state_tmp]*cv)^2/travel_time[state_tmp]^2))))
              u <- runif(1)
              tmp <- qlnorm(u * (1 - cdf) + cdf, meanlog = log(travel_time[state_tmp]^2/sqrt((travel_time[state_tmp]*cv)^2 + travel_time[state_tmp]^2)), 
                            sdlog = sqrt(log(1 + ((travel_time[state_tmp]*cv)^2/travel_time[state_tmp]^2))))
            }
            else {
              cdf <- pinvgauss(time-depart_time, mean = travel_time[state_tmp], shape = travel_time[state_tmp]/cv^2)
              u <- runif(1)
              tmp <- qinvgauss(u * (1 - cdf) + cdf, mean = travel_time[state_tmp], shape = travel_time[state_tmp]/cv^2)
            }
            next_completion_time <- tmp-(time-depart_time)
            next_event_time <- min(next_completion_time, min(next_failure_time))
            if (next_completion_time < min(next_failure_time)) {
              state[1] <- 0
              total_cost <- total_cost
              time <- time + next_event_time # update time
              depart_time <- time
            }
            else {
              state[which.min(next_failure_time)+1] <- 1
              fail_time[which.min(next_failure_time)] <- time + next_event_time
              total_cost <- total_cost
              time <- time + next_event_time # update time
            }
          }
        }
        else if (length(normal_set) == 0) {
          # next service completion time
          if (state[1] > 0) {
            #next_completion_time <- rexp(1, mu[state[1]])
            if (index_dist == 1) {
              if (expected_arriv_time >= time-depart_time) {
                cdf <- plnorm(time-depart_time, meanlog = log(travel_time[state[1]]^2/sqrt((travel_time[state[1]]*cv)^2 + travel_time[state[1]]^2)), 
                              sdlog = sqrt(log(1 + ((travel_time[state[1]]*cv)^2/travel_time[state[1]]^2))))
                u <- runif(1)
                tmp <- qlnorm(u * (1 - cdf) + cdf, meanlog = log(travel_time[state[1]]^2/sqrt((travel_time[state[1]]*cv)^2 + travel_time[state[1]]^2)), 
                              sdlog = sqrt(log(1 + ((travel_time[state[1]]*cv)^2/travel_time[state[1]]^2))))
                next_completion_time <- tmp+on_scene_service_time
              }
              else {
                next_completion_time <- on_scene_service_time-(time-depart_time-expected_arriv_time)
              }
            }
            else {
              if (expected_arriv_time >= time-depart_time) {
                cdf <- pinvgauss(time-depart_time, mean = travel_time[state[1]], shape = travel_time[state[1]]/cv^2)
                u <- runif(1)
                tmp <- qinvgauss(u * (1 - cdf) + cdf, mean = travel_time[state[1]], shape = travel_time[state[1]]/cv^2)
                next_completion_time <- tmp+on_scene_service_time
              }
              else {
                next_completion_time <- on_scene_service_time-(time-depart_time-expected_arriv_time)
              }
            }
            next_event_time <- next_completion_time
            state[1] <- -state[1]
            state[-state[1]+1] <- 0
            fail_time[-state[1]] <- ObservationWindow
            total_cost <- total_cost + sum(cost/alpha*(exp(-alpha*time/60)-exp(-alpha*(time+next_event_time)/60)))
            time <- time + next_event_time # update time
            depar_time <- time
          }
          else if (state[1] < 0) {
            state_tmp <- -state[1]
            #next_completion_time <- rexp(1, mu[state[1]])
            if (index_dist == 1) {
              cdf <- plnorm(time-depart_time, meanlog = log(travel_time[state_tmp]^2/sqrt((travel_time[state_tmp]*cv)^2 + travel_time[state_tmp]^2)), 
                            sdlog = sqrt(log(1 + ((travel_time[state_tmp]*cv)^2/travel_time[state_tmp]^2))))
              u <- runif(1)
              tmp <- qlnorm(u * (1 - cdf) + cdf, meanlog = log(travel_time[state_tmp]^2/sqrt((travel_time[state_tmp]*cv)^2 + travel_time[state_tmp]^2)), 
                            sdlog = sqrt(log(1 + ((travel_time[state_tmp]*cv)^2/travel_time[state_tmp]^2))))
              next_completion_time <- tmp-(time-depart_time)
            }
            else {
              cdf <- pinvgauss(time-depart_time, mean = travel_time[state_tmp], shape = travel_time[state_tmp]/cv^2)
              u <- runif(1)
              tmp <- qinvgauss(u * (1 - cdf) + cdf, mean = travel_time[state_tmp], shape = travel_time[state_tmp]/cv^2)
              next_completion_time <- tmp-(time-depart_time)
            }
            next_event_time <- next_completion_time
            state[1] <- 0
            total_cost <- total_cost + sum(cost/alpha*(exp(-alpha*time/60)-exp(-alpha*(time+next_event_time)/60)))
            time <- time + next_event_time # update time
            depart_time <- time
          }
          else {
            if (index_policy == 1) { # cmu-lambda-rule
              if (index_dist == 1) {
                tmp <- rlnorm(1, meanlog = log(travel_time[which.max(cost*mu/lambda)]^2/sqrt((travel_time[which.max(cost*mu/lambda)]*cv)^2 + travel_time[which.max(cost*mu/lambda)]^2)), 
                              sdlog = sqrt(log(1 + ((travel_time[which.max(cost*mu/lambda)]*cv)^2/travel_time[which.max(cost*mu/lambda)]^2))))
                next_completion_time <- tmp+on_scene_service_time
              }
              if (index_dist == 2) {
                tmp <- rinvgauss(1, mean = travel_time[which.max(cost*mu/lambda)], shape = travel_time[which.max(cost*mu/lambda)]/cv^2)
                next_completion_time <- tmp+on_scene_service_time
              }
              next_event_time <- next_completion_time
              state[1] <- -which.max(cost*mu/lambda)
              state[which.max(cost*mu/lambda)+1] <- 0
              fail_time[which.max(cost*mu/lambda)] <- ObservationWindow
              total_cost <- total_cost + sum(cost/alpha*(exp(-alpha*time/60)-exp(-alpha*(time+next_event_time)/60)))
              time <- time + next_event_time # update time
              depart_time <- time
            }
            if (index_policy == 2) { # cmu-rule
              if (index_dist == 1) {
                tmp <- rlnorm(1, meanlog = log(travel_time[which.max(cost*mu)]^2/sqrt((travel_time[which.max(cost*mu)]*cv)^2 + travel_time[which.max(cost*mu)]^2)), 
                              sdlog = sqrt(log(1 + ((travel_time[which.max(cost*mu)]*cv)^2/travel_time[which.max(cost*mu)]^2))))
                next_completion_time <- tmp+on_scene_service_time
              }
              if (index_dist == 2) {
                tmp <- rinvgauss(1, mean = travel_time[which.max(cost*mu)], shape = travel_time[which.max(cost*mu)]/cv^2)
                next_completion_time <- tmp+on_scene_service_time
              }
              next_event_time <- next_completion_time
              state[1] <- -which.max(cost*mu)
              state[which.max(cost*mu)+1] <- 0
              fail_time[which.max(cost*mu)] <- ObservationWindow
              total_cost <- total_cost + sum(cost/alpha*(exp(-alpha*time/60)-exp(-alpha*(time+next_event_time)/60)))
              time <- time + next_event_time # update time
              depart_time <- time
            }
            if (index_policy == 3) { # FCFS
              #next_completion_time <- rexp(1, mu[which.min(fail_time)])
              if (index_dist == 1) {
                tmp <- rlnorm(1, meanlog = log(travel_time[which.min(fail_time)]^2/sqrt((travel_time[which.min(fail_time)]*cv)^2 + travel_time[which.min(fail_time)]^2)), 
                              sdlog = sqrt(log(1 + ((travel_time[which.min(fail_time)]*cv)^2/travel_time[which.min(fail_time)]^2))))
                expected_arriv_time <- tmp
                next_completion_time <- tmp+on_scene_service_time
              }
              if (index_dist == 2) {
                tmp <- rinvgauss(1, mean = travel_time[which.min(fail_time)], shape = travel_time[which.min(fail_time)]/cv^2)
                expected_arriv_time <- tmp
                next_completion_time <- tmp+on_scene_service_time
              }
              next_event_time <- next_completion_time
              state[1] <- -which.min(fail_time)
              state[which.min(fail_time)+1] <- 0
              fail_time[which.min(fail_time)] <- ObservationWindow
              total_cost <- total_cost + sum(cost/alpha*(exp(-alpha*time/60)-exp(-alpha*(time+next_event_time)/60)))
              time <- time + next_event_time # update time
              depart_time <- time
            }
            if (index_policy == 4) { # nearest-first
              if (index_dist == 1) {
                tmp <- rlnorm(1, meanlog = log(travel_time[which.min(euclidean_distance)]^2/sqrt((travel_time[which.min(euclidean_distance)]*cv)^2 + travel_time[which.min(euclidean_distance)]^2)), 
                              sdlog = sqrt(log(1 + ((travel_time[which.min(euclidean_distance)]*cv)^2/travel_time[which.min(euclidean_distance)]^2))))
                next_completion_time <- tmp+on_scene_service_time
              }
              if (index_dist == 2) {
                tmp <- rinvgauss(1, mean = travel_time[which.min(euclidean_distance)], shape = travel_time[which.min(euclidean_distance)]/cv^2)
                next_completion_time <- tmp+on_scene_service_time
              }
              next_event_time <- next_completion_time
              state[1] <- -which.min(euclidean_distance)
              state[which.min(euclidean_distance)+1] <- 0
              fail_time[which.min(euclidean_distance)] <- ObservationWindow
              total_cost <- total_cost + sum(cost/alpha*(exp(-alpha*time/60)-exp(-alpha*(time+next_event_time)/60)))
              time <- time + next_event_time # update time
              depart_time <- time
            }
            if (index_policy == 5) { # least-mu rule
              if (index_dist == 1) {
                tmp <- rlnorm(1, meanlog = log(travel_time[which.min(mu)]^2/sqrt((travel_time[which.min(mu)]*cv)^2 + travel_time[which.min(mu)]^2)), 
                              sdlog = sqrt(log(1 + ((travel_time[which.min(mu)]*cv)^2/travel_time[which.min(mu)]^2))))
                next_completion_time <- tmp+on_scene_service_time
              }
              if (index_dist == 2) {
                tmp <- rinvgauss(1, mean = travel_time[which.min(mu)], shape = travel_time[which.min(mu)]/cv^2)
                next_completion_time <- tmp+on_scene_service_time
              }
              next_event_time <- next_completion_time
              state[1] <- -which.min(mu)
              state[which.min(mu)+1] <- 0
              fail_time[which.min(mu)] <- ObservationWindow
              total_cost <- total_cost + sum(cost/alpha*(exp(-alpha*time/60)-exp(-alpha*(time+next_event_time)/60)))
              time <- time + next_event_time # update time
              depart_time <- time
            }
            if (index_policy == 6) { # heuristic
              if (which.max(cost*mu/lambda) == which.max(cost*mu)) {
                if (index_dist == 1) {
                  tmp <- rlnorm(1, meanlog = log(travel_time[which.max(cost*mu)]^2/sqrt((travel_time[which.max(cost*mu)]*cv)^2 + travel_time[which.max(cost*mu)]^2)), 
                                sdlog = sqrt(log(1 + ((travel_time[which.max(cost*mu)]*cv)^2/travel_time[which.max(cost*mu)]^2))))
                  next_completion_time <- tmp+on_scene_service_time
                }
                if (index_dist == 2) {
                  tmp <- rinvgauss(1, mean = travel_time[which.max(cost*mu)], shape = travel_time[which.max(cost*mu)]/cv^2)
                  next_completion_time <- tmp+on_scene_service_time
                }
                
                next_event_time <- next_completion_time
                state[1] <- -which.max(cost*mu)
                state[which.max(cost*mu)+1] <- 0
                fail_time[which.max(cost*mu)] <- ObservationWindow
                total_cost <- total_cost + sum(cost/alpha*(exp(-alpha*time/60)-exp(-alpha*(time+next_event_time)/60)))
                time <- time + next_event_time # update time
                depart_time <- time
              }
              else {
                if (mu[which.max(cost*mu)] >= mu[which.max(cost*mu/lambda)]) {
                  if (index_dist == 1) {
                    tmp <- rlnorm(1, meanlog = log(travel_time[which.max(cost*mu/lambda)]^2/sqrt((travel_time[which.max(cost*mu/lambda)]*cv)^2 + travel_time[which.max(cost*mu/lambda)]^2)), 
                                  sdlog = sqrt(log(1 + ((travel_time[which.max(cost*mu/lambda)]*cv)^2/travel_time[which.max(cost*mu/lambda)]^2))))
                    next_completion_time <- tmp+on_scene_service_time
                  }
                  if (index_dist == 2) {
                    tmp <- rinvgauss(1, mean = travel_time[which.max(cost*mu/lambda)], shape = travel_time[which.max(cost*mu/lambda)]/cv^2)
                    next_completion_time <- tmp+on_scene_service_time
                  }
                  next_event_time <- next_completion_time
                  state[1] <- -which.max(cost*mu/lambda)
                  state[which.max(cost*mu/lambda)+1] <- 0
                  fail_time[which.max(cost*mu/lambda)] <- ObservationWindow
                  total_cost <- total_cost + sum(cost/alpha*(exp(-alpha*time/60)-exp(-alpha*(time+next_event_time)/60)))
                  time <- time + next_event_time # update time
                  depart_time <- time
                }
                else {
                  if ((cost*mu)[which.max(cost*mu/lambda)]>(1-(mu[which.max(cost*mu)]-mu[which.max(cost*mu/lambda)])/(Lambda+alpha))*max(cost*mu)) {
                    if (index_dist == 1) {
                      tmp <- rlnorm(1, meanlog = log(travel_time[which.max(cost*mu/lambda)]^2/sqrt((travel_time[which.max(cost*mu/lambda)]*cv)^2 + travel_time[which.max(cost*mu/lambda)]^2)), 
                                    sdlog = sqrt(log(1 + ((travel_time[which.max(cost*mu/lambda)]*cv)^2/travel_time[which.max(cost*mu/lambda)]^2))))
                      next_completion_time <- tmp+on_scene_service_time
                    }
                    if (index_dist == 2) {
                      tmp <- rinvgauss(1, mean = travel_time[which.max(cost*mu/lambda)], shape = travel_time[which.max(cost*mu/lambda)]/cv^2)
                      next_completion_time <- tmp+on_scene_service_time
                    }
                    next_event_time <- next_completion_time
                    state[1] <- -which.max(cost*mu/lambda)
                    state[which.max(cost*mu/lambda)+1] <- 0
                    fail_time[which.max(cost*mu/lambda)] <- ObservationWindow
                    total_cost <- total_cost + sum(cost/alpha*(exp(-alpha*time/60)-exp(-alpha*(time+next_event_time)/60)))
                    time <- time + next_event_time # update time
                    depart_time <- time
                  }
                  else {
                    if (index_dist == 1) {
                      tmp <- rlnorm(1, meanlog = log(travel_time[which.max(cost*mu)]^2/sqrt((travel_time[which.max(cost*mu)]*cv)^2 + travel_time[which.max(cost*mu)]^2)), 
                                    sdlog = sqrt(log(1 + ((travel_time[which.max(cost*mu)]*cv)^2/travel_time[which.max(cost*mu)]^2))))
                      next_completion_time <- tmp+on_scene_service_time
                    }
                    if (index_dist == 2) {
                      tmp <- rinvgauss(1, mean = travel_time[which.max(cost*mu)], shape = travel_time[which.max(cost*mu)]/cv^2)
                      next_completion_time <- tmp+on_scene_service_time
                    }
                    next_event_time <- next_completion_time
                    state[1] <- -which.max(cost*mu)
                    state[which.max(cost*mu)+1] <- 0
                    fail_time[which.max(cost*mu)] <- ObservationWindow
                    total_cost <- total_cost + sum(cost/alpha*(exp(-alpha*time/60)-exp(-alpha*(time+next_event_time)/60)))
                    time <- time + next_event_time # update time
                    depart_time <- time
                  }
                }
              }
            }
          }
        }
        else {
          for (i in 1:length(normal_set)) {
            next_failure_time[i] <- rexp(1, lambda[normal_set[i]])*60
          }
          if (state[1] > 0) { # departing to site
            #next_completion_time <- rexp(1, mu[state[1]])
            if (index_dist == 1) {
              if (expected_arriv_time >= time-depart_time) {
                cdf <- plnorm(time-depart_time, meanlog = log(travel_time[state[1]]^2/sqrt((travel_time[state[1]]*cv)^2 + travel_time[state[1]]^2)), 
                              sdlog = sqrt(log(1 + ((travel_time[state[1]]*cv)^2/travel_time[state[1]]^2))))
                u <- runif(1)
                tmp <- qlnorm(u * (1 - cdf) + cdf, meanlog = log(travel_time[state[1]]^2/sqrt((travel_time[state[1]]*cv)^2 + travel_time[state[1]]^2)), 
                              sdlog = sqrt(log(1 + ((travel_time[state[1]]*cv)^2/travel_time[state[1]]^2))))
                expected_arriv_time <- tmp-(time-depart_time)
                next_completion_time <- expected_arriv_time+on_scene_service_time
              }
              else {
                next_completion_time <- on_scene_service_time-(time-depart_time-expected_arriv_time)
              }
            }
            else {
              if (expected_arriv_time >= time-depart_time) {
                cdf <- pinvgauss(time-depart_time, mean = travel_time[state[1]], shape = travel_time[state[1]]/cv^2)
                u <- runif(1)
                tmp <- qinvgauss(u * (1 - cdf) + cdf, mean = travel_time[state[1]], shape = travel_time[state[1]]/cv^2)
                expected_arriv_time <- tmp-(time-depart_time)
                next_completion_time <- expected_arriv_time+on_scene_service_time
              }
              else {
                next_completion_time <- on_scene_service_time-(time-depart_time-expected_arriv_time)
              }
            }
            next_event_time <- min(next_failure_time, next_completion_time)
            total_cost <- total_cost + sum(cost*state[-1]/alpha*(exp(-alpha*time/60)-exp(-alpha*(time+next_event_time)/60)))
            if (next_completion_time < min(next_failure_time)) {
              state[state[1]+1] <- 0
              fail_time[state[1]] <- ObservationWindow
              state[1] <- -state[1]
              time <- time + next_event_time # update time
              depart_time <- time # update departing time
            }
            else {
              state[normal_set[which.min(next_failure_time)]+1] <- 1
              fail_time[normal_set[which.min(next_failure_time)]] <- time + next_event_time
              time <- time + next_event_time # update time
            }
          }
          else if (state[1] < 0) { # departing to site
            #next_completion_time <- rexp(1, mu[state[1]])
            state_tmp <- -state[1]
            if (index_dist == 1) {
              cdf <- plnorm(time-depart_time, meanlog = log(travel_time[state_tmp]^2/sqrt((travel_time[state_tmp]*cv)^2 + travel_time[state_tmp]^2)), 
                            sdlog = sqrt(log(1 + ((travel_time[state_tmp]*cv)^2/travel_time[state_tmp]^2))))
              u <- runif(1)
              tmp <- qlnorm(u * (1 - cdf) + cdf, meanlog = log(travel_time[state_tmp]^2/sqrt((travel_time[state_tmp]*cv)^2 + travel_time[state_tmp]^2)), 
                            sdlog = sqrt(log(1 + ((travel_time[state_tmp]*cv)^2/travel_time[state_tmp]^2))))
              next_completion_time <- tmp-(time-depart_time)
            }
            else {
              cdf <- pinvgauss(time-depart_time, mean = travel_time[state_tmp], shape = travel_time[state_tmp]/cv^2)
              u <- runif(1)
              tmp <- qinvgauss(u * (1 - cdf) + cdf, mean = travel_time[state_tmp], shape = travel_time[state_tmp]/cv^2)
              next_completion_time <- tmp-(time-depart_time)
            }
            next_event_time <- min(next_failure_time, next_completion_time)
            total_cost <- total_cost + sum(cost*state[-1]/alpha*(exp(-alpha*time/60)-exp(-alpha*(time+next_event_time)/60)))
            if (next_completion_time < min(next_failure_time)) {
              state[1] <- 0
              time <- time + next_event_time # update time
            }
            else {
              state[normal_set[which.min(next_failure_time)]+1] <- 1
              fail_time[normal_set[which.min(next_failure_time)]] <- time + next_event_time
              time <- time + next_event_time # update time
            }
          }
          else {
            if (index_policy == 1) { # cmu-lambda-rule
              if (index_dist == 1) {
                tmp <- rlnorm(1, meanlog = log(travel_time[failed_set[which.max(cost[failed_set]*mu[failed_set]/lambda[failed_set])]]^2/sqrt((travel_time[failed_set[which.max(cost[failed_set]*mu[failed_set]/lambda[failed_set])]]*cv)^2 + travel_time[failed_set[which.max(cost[failed_set]*mu[failed_set]/lambda[failed_set])]]^2)), 
                              sdlog = sqrt(log(1 + ((travel_time[failed_set[which.max(cost[failed_set]*mu[failed_set]/lambda[failed_set])]]*cv)^2/travel_time[failed_set[which.max(cost[failed_set]*mu[failed_set]/lambda[failed_set])]]^2))))
                expected_arriv_time <- tmp
                next_completion_time <- tmp+on_scene_service_time
              }
              else {
                tmp <- rinvgauss(1, mean = travel_time[failed_set[which.max(cost[failed_set]*mu[failed_set]/lambda[failed_set])]], shape = travel_time[failed_set[which.max(cost[failed_set]*mu[failed_set]/lambda[failed_set])]]/cv^2)
                expected_arriv_time <- tmp
                next_completion_time <- tmp+on_scene_service_time
              }
              next_event_time <- min(next_failure_time, next_completion_time)
              total_cost <- total_cost + sum(cost*state[-1]/alpha*(exp(-alpha*time/60)-exp(-alpha*(time+next_event_time)/60)))
              if (next_completion_time < min(next_failure_time)) { # completion
                state[failed_set[which.max(cost[failed_set]*mu[failed_set]/lambda[failed_set])]+1] <- 0
                state[1] <- -which.min(fail_time)
                fail_time[which.min(fail_time)] <- ObservationWindow
                time <- time + next_event_time # update time}
                depart_time <- time
              }
              else { # an new failure arrival
                state[1] <- failed_set[which.max(cost[failed_set]*mu[failed_set]/lambda[failed_set])]
                state[normal_set[which.min(next_failure_time)]+1] <- 1
                fail_time[normal_set[which.min(next_failure_time)]] <- time + next_event_time
                depart_time <- time
                time <- time + next_event_time # update time}
              }
            }
            if (index_policy == 2) { # cmu-rule
              if (index_dist == 1) {
                tmp <- rlnorm(1, meanlog = log(travel_time[failed_set[which.max(cost[failed_set]*mu[failed_set])]]^2/sqrt((travel_time[failed_set[which.max(cost[failed_set]*mu[failed_set])]]*cv)^2 + travel_time[failed_set[which.max(cost[failed_set]*mu[failed_set])]]^2)), 
                              sdlog = sqrt(log(1 + ((travel_time[failed_set[which.max(cost[failed_set]*mu[failed_set])]]*cv)^2/travel_time[failed_set[which.max(cost[failed_set]*mu[failed_set])]]^2))))
                expected_arriv_time <- tmp
                next_completion_time <- tmp+on_scene_service_time
              }
              else {
                tmp <- rinvgauss(1, mean = travel_time[failed_set[which.max(cost[failed_set]*mu[failed_set])]], shape = travel_time[failed_set[which.max(cost[failed_set]*mu[failed_set])]]/cv^2)
                expected_arriv_time <- tmp
                next_completion_time <- tmp+on_scene_service_time
              }
              next_event_time <- min(next_failure_time, next_completion_time)
              total_cost <- total_cost + sum(cost*state[-1]/alpha*(exp(-alpha*time/60)-exp(-alpha*(time+next_event_time)/60)))
              if (next_completion_time < min(next_failure_time)) { # completion
                state[failed_set[which.max(cost[failed_set]*mu[failed_set])]+1] <- 0
                state[1] <- -which.min(fail_time)
                fail_time[which.min(fail_time)] <- ObservationWindow
                time <- time + next_event_time # update time}
                depart_time <- time
              }
              else { # an new failure arrival
                state[1] <- failed_set[which.max(cost[failed_set]*mu[failed_set])]
                state[normal_set[which.min(next_failure_time)]+1] <- 1
                fail_time[normal_set[which.min(next_failure_time)]] <- time + next_event_time
                depart_time <- time
                time <- time + next_event_time # update time}
              }
            }
            if (index_policy == 3) { # FCFS
              #next_completion_time <- rexp(1, mu[which.min(fail_time)])
              if (index_dist == 1) {
                tmp <- rlnorm(1, meanlog = log(travel_time[which.min(fail_time)]^2/sqrt((travel_time[which.min(fail_time)]*cv)^2 + travel_time[which.min(fail_time)]^2)), 
                              sdlog = sqrt(log(1 + ((travel_time[which.min(fail_time)]*cv)^2/travel_time[which.min(fail_time)]^2))))
                expected_arriv_time <- tmp
                next_completion_time <- tmp+on_scene_service_time
              }
              else {
                tmp <- rinvgauss(1, mean = travel_time[which.min(fail_time)], shape = travel_time[which.min(fail_time)]/cv^2)
                expected_arriv_time <- tmp
                next_completion_time <- tmp+on_scene_service_time
              }
              next_event_time <- min(next_failure_time, next_completion_time)
              total_cost <- total_cost + sum(cost*state[-1]/alpha*(exp(-alpha*time/60)-exp(-alpha*(time+next_event_time)/60)))
              if (next_completion_time < min(next_failure_time)) { # completion
                state[which.min(fail_time)+1] <- 0
                state[1] <- -which.min(fail_time)
                fail_time[which.min(fail_time)] <- ObservationWindow
                time <- time + next_event_time # update time}
                depart_time <- time
              }
              else { # an new failure arrival
                state[1] <- which.min(fail_time)
                state[normal_set[which.min(next_failure_time)]+1] <- 1
                fail_time[normal_set[which.min(next_failure_time)]] <- time + next_event_time
                depart_time <- time
                time <- time + next_event_time # update time}
              }
            }
            if (index_policy == 4) { # nearest-first
              #next_completion <- rexp(1, mu[failed_set[which.min(euclidean_distance[failed_set])]])
              if (index_dist == 1) {
                tmp <- rlnorm(1, meanlog = log(travel_time[failed_set[which.min(euclidean_distance[failed_set])]]^2/sqrt((travel_time[failed_set[which.min(euclidean_distance[failed_set])]]*cv)^2 + travel_time[failed_set[which.min(euclidean_distance[failed_set])]]^2)), 
                              sdlog = sqrt(log(1 + ((travel_time[failed_set[which.min(euclidean_distance[failed_set])]]*cv)^2/travel_time[failed_set[which.min(euclidean_distance[failed_set])]]^2))))
                expected_arriv_time <- tmp
                next_completion_time <- tmp+on_scene_service_time
              }
              else {
                tmp <- rinvgauss(1, mean = travel_time[failed_set[which.min(euclidean_distance[failed_set])]], shape = travel_time[failed_set[which.min(euclidean_distance[failed_set])]]/cv^2)
                expected_arriv_time <- tmp
                next_completion_time <- tmp+on_scene_service_time
              }
              next_event_time <- min(next_failure_time, next_completion_time)
              total_cost <- total_cost + sum(cost*state[-1]/alpha*(exp(-alpha*time/60)-exp(-alpha*(time+next_event_time)/60)))
              if (next_completion_time < min(next_failure_time)) { # completion
                state[failed_set[which.min(euclidean_distance[failed_set])]+1] <- 0
                state[1] <- -failed_set[which.min(euclidean_distance[failed_set])]
                fail_time[which.min(fail_time)] <- ObservationWindow
                time <- time + next_event_time # update time}
                depart_time <- time
              }
              else { # an new failure arrival
                state[1] <- failed_set[which.min(euclidean_distance[failed_set])]
                state[normal_set[which.min(next_failure_time[failed_set])]+1] <- 1
                fail_time[normal_set[which.min(next_failure_time)]] <- time + next_event_time
                depart_time <- time
                time <- time + next_event_time # update time}
              }
            }
            if (index_policy == 5) { # least-mu rule
              if (index_dist == 1) {
                tmp <- rlnorm(1, meanlog = log(travel_time[failed_set[which.min(mu[failed_set])]]^2/sqrt((travel_time[failed_set[which.min(mu[failed_set])]]*cv)^2 + travel_time[failed_set[which.min(mu[failed_set])]]^2)), 
                              sdlog = sqrt(log(1 + ((travel_time[failed_set[which.min(mu[failed_set])]]*cv)^2/travel_time[failed_set[which.min(mu[failed_set])]]^2))))
                expected_arriv_time <- tmp
                next_completion_time <- tmp+on_scene_service_time
              }
              else {
                tmp <- rinvgauss(1, mean = travel_time[failed_set[which.min(mu[failed_set])]], shape = travel_time[failed_set[which.min(mu[failed_set])]]/cv^2)
                expected_arriv_time <- tmp
                next_completion_time <- tmp+on_scene_service_time
              }
              next_event_time <- min(next_failure_time, next_completion_time)
              total_cost <- total_cost + sum(cost*state[-1]/alpha*(exp(-alpha*time/60)-exp(-alpha*(time+next_event_time)/60)))
              if (next_completion_time < min(next_failure_time)) { # completion
                state[failed_set[which.min(mu[failed_set])]+1] <- 0
                state[1] <- -which.min(fail_time)
                fail_time[which.min(fail_time)] <- ObservationWindow
                time <- time + next_event_time # update time}
                depart_time <- time
              }
              else { # an new failure arrival
                state[1] <- failed_set[which.min(mu[failed_set])]
                state[normal_set[which.min(next_failure_time)]+1] <- 1
                fail_time[normal_set[which.min(next_failure_time)]] <- time + next_event_time
                depart_time <- time
                time <- time + next_event_time # update time}
              }
            }
            if (index_policy == 6) { # heuristic
              if (which.max((cost*mu/lambda)[failed_set]) == which.max((cost*mu)[failed_set])) {
                if (index_dist == 1) {
                  tmp <- rlnorm(1, meanlog = log(travel_time[failed_set[which.max(cost[failed_set]*mu[failed_set])]]^2/sqrt((travel_time[failed_set[which.max(cost[failed_set]*mu[failed_set])]]*cv)^2 + travel_time[failed_set[which.max(cost[failed_set]*mu[failed_set])]]^2)), 
                                sdlog = sqrt(log(1 + ((travel_time[failed_set[which.max(cost[failed_set]*mu[failed_set])]]*cv)^2/travel_time[failed_set[which.max(cost[failed_set]*mu[failed_set])]]^2))))
                  expected_arriv_time <- tmp
                  next_completion_time <- tmp+on_scene_service_time
                }
                else {
                  tmp <- rinvgauss(1, mean = travel_time[failed_set[which.max(cost[failed_set]*mu[failed_set])]], shape = travel_time[failed_set[which.max(cost[failed_set]*mu[failed_set])]]/cv^2)
                  expected_arriv_time <- tmp
                  next_completion_time <- tmp+on_scene_service_time
                }
                next_event_time <- min(next_failure_time, next_completion_time)
                total_cost <- total_cost + sum(cost*state[-1]/alpha*(exp(-alpha*time/60)-exp(-alpha*(time+next_event_time)/60)))
                if (next_completion_time < min(next_failure_time)) { # completion
                  state[failed_set[which.max(cost[failed_set]*mu[failed_set])]+1] <- 0
                  state[1] <- -which.min(fail_time)
                  fail_time[which.min(fail_time)] <- ObservationWindow
                  time <- time + next_event_time # update time}
                  depart_time <- time
                }
                else { # an new failure arrival
                  state[1] <- failed_set[which.max(cost[failed_set]*mu[failed_set])]
                  state[normal_set[which.min(next_failure_time)]+1] <- 1
                  fail_time[normal_set[which.min(next_failure_time)]] <- time + next_event_time
                  depart_time <- time
                  time <- time + next_event_time # update time}
                }
              }
              else {
                if (mu[failed_set[which.max((cost*mu)[failed_set])]] >= mu[failed_set[which.max((cost*mu/lambda)[failed_set])]]) {
                  if (index_dist == 1) {
                    tmp <- rlnorm(1, meanlog = log(travel_time[failed_set[which.max(cost[failed_set]*mu[failed_set]/lambda[failed_set])]]^2/sqrt((travel_time[failed_set[which.max(cost[failed_set]*mu[failed_set]/lambda[failed_set])]]*cv)^2 + travel_time[failed_set[which.max(cost[failed_set]*mu[failed_set]/lambda[failed_set])]]^2)), 
                                  sdlog = sqrt(log(1 + ((travel_time[failed_set[which.max(cost[failed_set]*mu[failed_set]/lambda[failed_set])]]*cv)^2/travel_time[failed_set[which.max(cost[failed_set]*mu[failed_set]/lambda[failed_set])]]^2))))
                    expected_arriv_time <- tmp
                    next_completion_time <- tmp+on_scene_service_time
                  }
                  else {
                    tmp <- rinvgauss(1, mean = travel_time[failed_set[which.max(cost[failed_set]*mu[failed_set]/lambda[failed_set])]], shape = travel_time[failed_set[which.max(cost[failed_set]*mu[failed_set]/lambda[failed_set])]]/cv^2)
                    expected_arriv_time <- tmp
                    next_completion_time <- tmp+on_scene_service_time
                  }
                  next_event_time <- min(next_failure_time, next_completion_time)
                  total_cost <- total_cost + sum(cost*state[-1]/alpha*(exp(-alpha*time/60)-exp(-alpha*(time+next_event_time)/60)))
                  if (next_completion_time < min(next_failure_time)) { # completion
                    state[failed_set[which.max(cost[failed_set]*mu[failed_set]/lambda[failed_set])]+1] <- 0
                    state[1] <- -which.min(fail_time)
                    fail_time[which.min(fail_time)] <- ObservationWindow
                    time <- time + next_event_time # update time}
                    depart_time <- time
                  }
                  else { # an new failure arrival
                    state[1] <- failed_set[which.max(cost[failed_set]*mu[failed_set]/lambda[failed_set])]
                    state[normal_set[which.min(next_failure_time)]+1] <- 1
                    fail_time[normal_set[which.min(next_failure_time)]] <- time + next_event_time
                    depart_time <- time
                    time <- time + next_event_time # update time}
                  }
                }
                else {
                  if ((cost*mu)[failed_set[which.max((cost*mu/lambda)[failed_set])]]>(1-(mu[failed_set[which.max((cost*mu)[failed_set])]]-mu[failed_set[which.max((cost*mu/lambda)[failed_set])]])/(Lambda+alpha))*max((cost*mu)[failed_set])) {
                    if (index_dist == 1) {
                      tmp <- rlnorm(1, meanlog = log(travel_time[failed_set[which.max(cost[failed_set]*mu[failed_set]/lambda[failed_set])]]^2/sqrt((travel_time[failed_set[which.max(cost[failed_set]*mu[failed_set]/lambda[failed_set])]]*cv)^2 + travel_time[failed_set[which.max(cost[failed_set]*mu[failed_set]/lambda[failed_set])]]^2)), 
                                    sdlog = sqrt(log(1 + ((travel_time[failed_set[which.max(cost[failed_set]*mu[failed_set]/lambda[failed_set])]]*cv)^2/travel_time[failed_set[which.max(cost[failed_set]*mu[failed_set]/lambda[failed_set])]]^2))))
                      expected_arriv_time <- tmp
                      next_completion_time <- tmp+on_scene_service_time
                    }
                    else {
                      tmp <- rinvgauss(1, mean = travel_time[failed_set[which.max(cost[failed_set]*mu[failed_set]/lambda[failed_set])]], shape = travel_time[failed_set[which.max(cost[failed_set]*mu[failed_set]/lambda[failed_set])]]/cv^2)
                      expected_arriv_time <- tmp
                      next_completion_time <- tmp+on_scene_service_time
                    }
                    next_event_time <- min(next_failure_time, next_completion_time)
                    total_cost <- total_cost + sum(cost*state[-1]/alpha*(exp(-alpha*time/60)-exp(-alpha*(time+next_event_time)/60)))
                    if (next_completion_time < min(next_failure_time)) { # completion
                      state[failed_set[which.max(cost[failed_set]*mu[failed_set]/lambda[failed_set])]+1] <- 0
                      state[1] <- -which.min(fail_time)
                      fail_time[which.min(fail_time)] <- ObservationWindow
                      time <- time + next_event_time # update time}
                      depart_time <- time
                    }
                    else { # an new failure arrival
                      state[1] <- failed_set[which.max(cost[failed_set]*mu[failed_set]/lambda[failed_set])]
                      state[normal_set[which.min(next_failure_time)]+1] <- 1
                      fail_time[normal_set[which.min(next_failure_time)]] <- time + next_event_time
                      depart_time <- time
                      time <- time + next_event_time # update time}
                    }
                  }
                  else {
                    if (index_dist == 1) {
                      tmp <- rlnorm(1, meanlog = log(travel_time[failed_set[which.max(cost[failed_set]*mu[failed_set])]]^2/sqrt((travel_time[failed_set[which.max(cost[failed_set]*mu[failed_set])]]*cv)^2 + travel_time[failed_set[which.max(cost[failed_set]*mu[failed_set])]]^2)), 
                                    sdlog = sqrt(log(1 + ((travel_time[failed_set[which.max(cost[failed_set]*mu[failed_set])]]*cv)^2/travel_time[failed_set[which.max(cost[failed_set]*mu[failed_set])]]^2))))
                      expected_arriv_time <- tmp
                      next_completion_time <- tmp+on_scene_service_time
                    }
                    else {
                      tmp <- rinvgauss(1, mean = travel_time[failed_set[which.max(cost[failed_set]*mu[failed_set])]], shape = travel_time[failed_set[which.max(cost[failed_set]*mu[failed_set])]]/cv^2)
                      expected_arriv_time <- tmp
                      next_completion_time <- tmp+on_scene_service_time
                    }
                    next_event_time <- min(next_failure_time, next_completion_time)
                    total_cost <- total_cost + sum(cost*state[-1]/alpha*(exp(-alpha*time/60)-exp(-alpha*(time+next_event_time)/60)))
                    if (next_completion_time < min(next_failure_time)) { # completion
                      state[failed_set[which.max(cost[failed_set]*mu[failed_set])]+1] <- 0
                      state[1] <- -which.min(fail_time)
                      fail_time[which.min(fail_time)] <- ObservationWindow
                      time <- time + next_event_time # update time}
                      depart_time <- time
                    }
                    else { # an new failure arrival
                      state[1] <- failed_set[which.max(cost[failed_set]*mu[failed_set])]
                      state[normal_set[which.min(next_failure_time)]+1] <- 1
                      fail_time[normal_set[which.min(next_failure_time)]] <- time + next_event_time
                      depart_time <- time
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
save(list = c("allResult"), file = paste0("NumBaseStation", NumBaseStation, "relaxation1", ".Rdata"))

stopCluster(cl)
