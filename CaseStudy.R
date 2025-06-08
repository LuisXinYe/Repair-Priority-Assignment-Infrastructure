library(readxl)
library(statmod)
library(geosphere)
cur_addr <- dirname(rstudioapi::getActiveDocumentContext()$path)

# cmulambda-rule: 1
# cmu-rule: 2
# FCFS: 3
# nearest-first: 4
# least-mu: 5
# heuristic: 6

NumReplications <- 100
area <- 1 # 1:Urban 2:Suburban
if (area == 1) {
  NumBaseStation <- 30
  
  travel <- read_xlsx(paste0(cur_addr, '/CaseStudy.xlsx', sep = ""), sheet = 1)[["行驶时间"]][1:30]/3600
  on_scene_service_time <- 1.25 # unit: hour
  mu <- 1/(travel+on_scene_service_time) # service rate; unit: 1/hr
  cv <- 1 # cv is 1
  
  lambda <- read_xlsx(paste0(cur_addr, '/CaseStudy.xlsx', sep = ""), sheet = 1)[["累积失效次数"]][1:30]/240/24 # arrival rate; unit: 1/hr
  cost <- read_xlsx(paste0(cur_addr, '/CaseStudy.xlsx', sep = ""), sheet = 1)[["失效成本率"]][1:30]  # unit: CNY per hour
  longitude <- read_xlsx(paste0(cur_addr, '/CaseStudy.xlsx', sep = ""), sheet = 1)[["经度"]][1:30]  # longitude
  latitude <- read_xlsx(paste0(cur_addr, '/CaseStudy.xlsx', sep = ""), sheet = 1)[["纬度"]][1:30]  # latitude
  Lambda <- sum(lambda+mu)
  rho <- 0.0239 # real interest rate
  alpha <- 0.0045 #discount rate
  
  initial_state <- c(0, rep(0, NumBaseStation)) # initial state
  #initial_state <- c(0, rep(1, NumBaseStation)) # initial state
  
  haversine_distance <- vector("numeric", NumBaseStation)
  for (i in 1:NumBaseStation) { # Urban area
    p1 <- c(longitude[i],latitude[i])
    depot <- c(123.971168, 47.332965)
    haversine_distance[i] <- distHaversine(p1,depot)
  }
} else {
  NumBaseStation <- 16
  
  travel <- read_xlsx(paste0(cur_addr, '/CaseStudy.xlsx', sep = ""), sheet = 2)[["行驶时间"]][1:16]/3600
  on_scene_service_time <- 1.25 # unit: hour
  mu <- 1/(travel+on_scene_service_time) # service rate; unit: 1/hr
  cv <- 1 # cv is 1
  
  lambda <- read_xlsx(paste0(cur_addr, '/CaseStudy.xlsx', sep = ""), sheet = 2)[["累积失效次数"]][1:16]/240/24 # arrival rate; unit: 1/hr
  cost <- read_xlsx(paste0(cur_addr, '/CaseStudy.xlsx', sep = ""), sheet = 2)[["失效成本率"]][1:16]  # unit: CNY per hour
  longitude <- read_xlsx(paste0(cur_addr, '/CaseStudy.xlsx', sep = ""), sheet = 2)[["经度"]][1:16]  # longitude
  latitude <- read_xlsx(paste0(cur_addr, '/CaseStudy.xlsx', sep = ""), sheet = 2)[["纬度"]][1:16]  # latitude
  Lambda <- sum(lambda+mu)
  rho <- 0.0239 # real interest rate
  #alpha <- Lambda*rho #discount rate
  alpha <- rho/365/24 #discount rate
  initial_state <- c(0, rep(0, NumBaseStation)) # initial state
  #initial_state <- c(0, rep(1, NumBaseStation)) # initial state
  
  haversine_distance <- vector("numeric", NumBaseStation)
  for (i in 1:NumBaseStation) { # Suburban area
    p1 <- c(longitude[i],latitude[i])
    depot <- c(123.613877, 47.216194)
    haversine_distance[i] <- distHaversine(p1,depot)
  }
}

total_cost <- vector("list", length = 6)
for (l in 1:6) {
  index_policy <- l # indicate the applied policy
  total_cost[[l]] <- vector("numeric", length = NumReplications)  # total interruption cost
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
        if (state[1] == 0) {
          next_event_time <- min(next_failure_time)
          state[which.min(next_failure_time)+1] <- 1
          fail_time[which.min(next_failure_time)] <- time + next_event_time
          total_cost[[l]][h] <- total_cost[[l]][h]
          time <- time + next_event_time # update time
        } else {
          state_tmp <- -state[1]
          cdf <- pinvgauss(time-depart_time, mean = travel[state_tmp], shape = travel[state_tmp]/cv^2)
          u <- runif(1)
          tmp <- qinvgauss(u * (1 - cdf) + cdf, mean = travel[state_tmp], shape = travel[state_tmp]/cv^2)
          next_completion_time <- tmp-(time-depart_time)
          next_event_time <- min(next_completion_time, min(next_failure_time))
          if (next_completion_time < min(next_failure_time)) {
            state[1] <- 0
            total_cost[[l]][h] <- total_cost[[l]][h]
            time <- time + next_event_time # update time
            depart_time <- time
          }
          else {
            state[which.min(next_failure_time)+1] <- 1
            fail_time[which.min(next_failure_time)] <- time + next_event_time
            total_cost[[l]][h] <- total_cost[[l]][h]
            time <- time + next_event_time # update time
          }
        }
      }
      else if (length(normal_set) == 0) {
        # next service completion time
        if (state[1] > 0) {
          #next_completion_time <- rexp(1, mu[state[1]])
          if (expected_arriv_time >= time-depart_time) {
            cdf <- pinvgauss(time-depart_time, mean = travel[state[1]], shape = travel[state[1]]/cv^2)
            u <- runif(1)
            tmp <- qinvgauss(u * (1 - cdf) + cdf, mean = travel[state[1]], shape = travel[state[1]]/cv^2)
            expected_arriv_time <- tmp-(time-depart_time)
            next_completion_time <- expected_arriv_time+on_scene_service_time
          }
          else {
            next_completion_time <- on_scene_service_time-(time-depart_time-expected_arriv_time)
          }
          next_event_time <- next_completion_time
          state[1] <- -state[1]
          state[-state[1]+1] <- 0
          fail_time[-state[1]] <- ObservationWindow
          total_cost[[l]][h] <- total_cost[[l]][h] + sum(cost*state[-1]/alpha*(exp(-alpha*time)-exp(-alpha*(time+next_event_time))))
          time <- time + next_event_time # update time
          depar_time <- time
        }
        else if (state[1] < 0) {
          state_tmp <- -state[1]
          #next_completion_time <- rexp(1, mu[state[1]])
          cdf <- pinvgauss(time-depart_time, mean = travel[state_tmp], shape = travel[state_tmp]/cv^2)
          u <- runif(1)
          tmp <- qinvgauss(u * (1 - cdf) + cdf, mean = travel[state_tmp], shape = travel[state_tmp]/cv^2)
          next_completion_time <- tmp-(time-depart_time)
          next_event_time <- next_completion_time
          state[1] <- 0
          total_cost[[l]][h] <- total_cost[[l]][h] + sum(cost/alpha*(exp(-alpha*time)-exp(-alpha*(time+next_event_time))))
          time <- time + next_event_time # update time
          depart_time <- time
        }
        else {
          if (index_policy == 1) { # cmulambda-rule
            #next_completion_time <- rexp(1, mu[which.min(cost*mu)])
            tmp <- rinvgauss(1, mean = travel[which.max(cost*mu/lambda)], shape = travel[which.max(cost*mu/lambda)]/cv^2)
            expected_arriv_time <- tmp
            next_completion_time <- tmp+on_scene_service_time
            next_event_time <- next_completion_time
            state[1] <- -which.max(cost*mu/lambda)
            state[which.max(cost*mu/lambda)+1] <- 0
            total_cost[[l]][h] <- total_cost[[l]][h] + sum(cost/alpha*(exp(-alpha*time)-exp(-alpha*(time+next_event_time))))
            time <- time + next_event_time # update time
            depart_time <- time
          }
          if (index_policy == 2) { # cmu-rule
            #next_completion_time <- rexp(1, mu[which.min(cost*mu)])
            tmp <- rinvgauss(1, mean = travel[which.max(cost*mu)], shape = travel[which.max(cost*mu)]/cv^2)
            expected_arriv_time <- tmp
            next_completion_time <- tmp+on_scene_service_time
            next_event_time <- next_completion_time
            state[1] <- -which.max(cost*mu)
            state[which.max(cost*mu)+1] <- 0
            total_cost[[l]][h] <- total_cost[[l]][h] + sum(cost/alpha*(exp(-alpha*time)-exp(-alpha*(time+next_event_time))))
            time <- time + next_event_time # update time
            depart_time <- time
          }
          if (index_policy == 3) { # FCFS
            #next_completion_time <- rexp(1, mu[which.min(fail_time)])
            tmp <- rinvgauss(1, mean = travel[which.min(fail_time)], shape = travel[which.min(fail_time)]/cv^2)
            expected_arriv_time <- tmp
            next_completion_time <- tmp+on_scene_service_time
            next_event_time <- next_completion_time
            state[1] <- -which.min(fail_time)
            state[which.min(fail_time)+1] <- 0
            fail_time[which.min(fail_time)] <- ObservationWindow
            total_cost[[l]][h] <- total_cost[[l]][h] + sum(cost/alpha*(exp(-alpha*time)-exp(-alpha*(time+next_event_time))))
            time <- time + next_event_time # update time
            depart_time <- time
          }
          if (index_policy == 4) { # nearest-first
            #next_completion_time <- rexp(1, mu[which.min(haversine_distance)])
            tmp <- rinvgauss(1, mean = travel[which.min(haversine_distance)], shape = travel[which.min(haversine_distance)]/cv^2)
            expected_arriv_time <- tmp
            next_completion_time <- tmp+on_scene_service_time
            next_event_time <- next_completion_time
            state[1] <- -which.min(haversine_distance)
            state[which.min(haversine_distance)+1] <- 0
            total_cost[[l]][h] <- total_cost[[l]][h] + sum(cost/alpha*(exp(-alpha*time)-exp(-alpha*(time+next_event_time))))
            time <- time + next_event_time # update time
            depart_time <- time
          }
          if (index_policy == 5) { # least-mu
            #next_completion_time <- rexp(1, mu[which.min(mu)])
            tmp <- rinvgauss(1, mean = travel[which.min(mu)], shape = travel[which.min(mu)]/cv^2)
            expected_arriv_time <- tmp
            next_completion_time <- tmp+on_scene_service_time
            next_event_time <- next_completion_time
            state[1] <- -which.min(mu)
            state[which.min(mu)+1] <- 0
            total_cost[[l]][h] <- total_cost[[l]][h] + sum(cost/alpha*(exp(-alpha*time)-exp(-alpha*(time+next_event_time))))
            time <- time + next_event_time # update time
            depart_time <- time
          }
          if (index_policy == 6) { # heuristic
            AvailableAction <- which(state[-1] == 1)
            
            if (which.max(cost*mu/lambda) == which.max(cost*mu)) optim_action <- which.max(cost*mu/lambda)
            else {
              if (mu[which.max(cost*mu)] >= mu[which.max(cost*mu/lambda)]) {
                optim_action <- which.max(cost*mu/lambda)
              }
              else {
                if ((cost*mu)[which.max(cost*mu/lambda)]>(1-(mu[which.max(cost*mu)]-mu[which.max(cost*mu/lambda)])/(Lambda+alpha))*max(cost*mu)) {
                  optim_action <- which.max(cost*mu/lambda)
                }
                else optim_action <- which.max(cost*mu)
              }
            }
            
            tmp <- rinvgauss(1, mean = travel[AvailableAction[optim_action]], shape = travel[AvailableAction[optim_action]]/cv^2)
            expected_arriv_time <- tmp
            next_completion_time <- tmp+on_scene_service_time
            next_event_time <- next_completion_time
            total_cost[[l]][h] <- total_cost[[l]][h] + sum(cost*state[-1]/alpha*(exp(-alpha*time)-exp(-alpha*(time+next_event_time))))
            state[AvailableAction[optim_action]+1] <- 0
            state[1] <- -AvailableAction[optim_action]
            time <- time + next_event_time # update time}
            depart_time <- time
          }
        }
      }
      else {
        for (i in 1:length(normal_set)) {
          next_failure_time[i] <- rexp(1, lambda[normal_set[i]])
        }
        if (state[1] > 0) { # departing to base station
          #next_completion_time <- rexp(1, mu[state[1]])
          if (expected_arriv_time >= time-depart_time) {
            cdf <- pinvgauss(time-depart_time, mean = travel[state[1]], shape = travel[state[1]]/cv^2)
            u <- runif(1)
            tmp <- qinvgauss(u * (1 - cdf) + cdf, mean = travel[state[1]], shape = travel[state[1]]/cv^2)
            expected_arriv_time <- tmp-(time-depart_time)
            next_completion_time <- expected_arriv_time+on_scene_service_time
          }
          else {
            next_completion_time <- on_scene_service_time-(time-depart_time-expected_arriv_time)
          }
          next_event_time <- min(next_failure_time, next_completion_time)
          total_cost[[l]][h] <- total_cost[[l]][h] + sum(cost*state[-1]/alpha*(exp(-alpha*time)-exp(-alpha*(time+next_event_time))))
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
        else if (state[1] < 0) { # back to site
          #next_completion_time <- rexp(1, mu[state[1]])
          state_tmp <- -state[1]
          cdf <- pinvgauss(time-depart_time, mean = travel[state_tmp], shape = travel[state_tmp]/cv^2)
          u <- runif(1)
          tmp <- qinvgauss(u * (1 - cdf) + cdf, mean = travel[state_tmp], shape = travel[state_tmp]/cv^2)
          next_completion_time <- tmp-(time-depart_time)
          next_event_time <- min(next_failure_time, next_completion_time)
          total_cost[[l]][h] <- total_cost[[l]][h] + sum(cost*state[-1]/alpha*(exp(-alpha*time)-exp(-alpha*(time+next_event_time))))
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
          if (index_policy == 1) { # cmulambda-rule
            #next_completion <- rexp(1, mu[failed_set[which.min(cost[failed_set]*mu[failed_set])]])
            tmp <- rinvgauss(1, mean = travel[failed_set[which.max(cost[failed_set]*mu[failed_set]/lambda[failed_set])]], shape = travel[failed_set[which.max(cost[failed_set]*mu[failed_set]/lambda[failed_set])]]/cv^2)
            expected_arriv_time <- tmp
            next_completion_time <- tmp+on_scene_service_time
            next_event_time <- min(next_failure_time, next_completion_time)
            total_cost[[l]][h] <- total_cost[[l]][h] + sum(cost*state[-1]/alpha*(exp(-alpha*time)-exp(-alpha*(time+next_event_time))))
            if (next_completion_time < min(next_failure_time)) { # completion
              state[failed_set[which.max(cost[failed_set]*mu[failed_set]/lambda[failed_set])]+1] <- 0
              state[1] <- -failed_set[which.max(cost[failed_set]*mu[failed_set]/lambda[failed_set])]
              time <- time + next_event_time # update time}
              depart_time <- time
            }
            else { # an new failure arrival
              state[1] <- failed_set[which.max(cost[failed_set]*mu[failed_set]/lambda[failed_set])]
              state[normal_set[which.min(next_failure_time)]+1] <- 1
              depart_time <- time
              time <- time + next_event_time # update time}
            }
          }
          if (index_policy == 2) { # cmu-rule
            #next_completion <- rexp(1, mu[failed_set[which.min(cost[failed_set]*mu[failed_set])]])
            tmp <- rinvgauss(1, mean = travel[failed_set[which.max(cost[failed_set]*mu[failed_set])]], shape = travel[failed_set[which.max(cost[failed_set]*mu[failed_set])]]/cv^2)
            expected_arriv_time <- tmp
            next_completion_time <- tmp+on_scene_service_time
            next_event_time <- min(next_failure_time, next_completion_time)
            total_cost[[l]][h] <- total_cost[[l]][h] + sum(cost*state[-1]/alpha*(exp(-alpha*time)-exp(-alpha*(time+next_event_time))))
            if (next_completion_time < min(next_failure_time)) { # completion
              state[failed_set[which.max(cost[failed_set]*mu[failed_set])]+1] <- 0
              state[1] <- -failed_set[which.max(cost[failed_set]*mu[failed_set])]
              time <- time + next_event_time # update time}
              depart_time <- time
            }
            else { # an new failure arrival
              state[1] <- failed_set[which.max(cost[failed_set]*mu[failed_set])]
              state[normal_set[which.min(next_failure_time)]+1] <- 1
              depart_time <- time
              time <- time + next_event_time # update time}
            }
          }
          if (index_policy == 3) { # FCFS
            #next_completion_time <- rexp(1, mu[which.min(fail_time)])
            tmp <- rinvgauss(1, mean = travel[which.min(fail_time)], shape = travel[which.min(fail_time)]/cv^2)
            expected_arriv_time <- tmp
            next_completion_time <- tmp+on_scene_service_time
            next_event_time <- min(next_failure_time, next_completion_time)
            total_cost[[l]][h] <- total_cost[[l]][h] + sum(cost*state[-1]/alpha*(exp(-alpha*time)-exp(-alpha*(time+next_event_time))))
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
            tmp <- rinvgauss(1, mean = travel[failed_set[which.min(haversine_distance[failed_set])]], shape = travel[failed_set[which.min(haversine_distance[failed_set])]]/cv^2)
            expected_arriv_time <- tmp
            next_completion_time <- tmp+on_scene_service_time
            next_event_time <- min(next_failure_time, next_completion_time)
            total_cost[[l]][h] <- total_cost[[l]][h] + sum(cost*state[-1]/alpha*(exp(-alpha*time)-exp(-alpha*(time+next_event_time))))
            if (next_completion_time < min(next_failure_time)) { # completion
              state[failed_set[which.min(haversine_distance[failed_set])]+1] <- 0
              state[1] <- -failed_set[which.min(haversine_distance[failed_set])]
              time <- time + next_event_time # update time}
              depart_time <- time
            }
            else { # an new failure arrival
              state[1] <- failed_set[which.min(haversine_distance[failed_set])]
              state[normal_set[which.min(next_failure_time)]+1] <- 1
              depart_time <- time
              time <- time + next_event_time # update time}
            }
          }
          if (index_policy == 5) { # least-mu
            tmp <- rinvgauss(1, mean = travel[failed_set[which.min(mu[failed_set])]], shape = travel[failed_set[which.min(mu[failed_set])]]/cv^2)
            expected_arriv_time <- tmp
            next_completion_time <- tmp+on_scene_service_time
            next_event_time <- min(next_failure_time, next_completion_time)
            total_cost[[l]][h] <- total_cost[[l]][h] + sum(cost*state[-1]/alpha*(exp(-alpha*time)-exp(-alpha*(time+next_event_time))))
            if (next_completion_time < min(next_failure_time)) { # completion
              state[failed_set[which.min(mu[failed_set])]+1] <- 0
              state[1] <- -failed_set[which.min(mu[failed_set])]
              time <- time + next_event_time # update time}
              depart_time <- time
            }
            else { # an new failure arrival
              state[1] <- failed_set[which.min(mu[failed_set])]
              state[normal_set[which.min(next_failure_time)]+1] <- 1
              depart_time <- time
              time <- time + next_event_time # update time}
            }
          }
          if (index_policy == 6) { # heuristic
            AvailableAction <- which(state[-1] == 1)
            TmpTerm_1 <- rep(0,length(AvailableAction))
            TmpTerm_2 <- rep(0,length(AvailableAction))
            for (k in 1:length(AvailableAction)) {
              TmpTerm_1[k] <- (cost*mu/lambda)[AvailableAction[k]]
              TmpTerm_2[k] <- (cost*mu)[AvailableAction[k]]
            }
            
            if (which.max(TmpTerm_1) == which.max(TmpTerm_2)) optim_action <- which.max(TmpTerm_1)
            else {
              if (mu[which.max(cost*mu)] >= mu[which.max(cost*mu/lambda)]) {
                optim_action <- which.max(TmpTerm_1)
              }
              else {
                if ((cost*mu)[which.max(TmpTerm_1)]>(1-(mu[which.max(TmpTerm_2)]-mu[which.max(TmpTerm_1)])/(Lambda+alpha))*max(TmpTerm_2)) {
                  optim_action <- which.max(TmpTerm_1)
                }
                else optim_action <- which.max(TmpTerm_2)
              }
            }
             
            tmp <- rinvgauss(1, mean = travel[AvailableAction[optim_action]], shape = travel[AvailableAction[optim_action]]/cv^2)
            expected_arriv_time <- tmp
            next_completion_time <- tmp+on_scene_service_time
            next_event_time <- min(next_failure_time, next_completion_time)
            total_cost[[l]][h] <- total_cost[[l]][h] + sum(cost*state[-1]/alpha*(exp(-alpha*time)-exp(-alpha*(time+next_event_time))))
            if (next_completion_time < min(next_failure_time)) { # completion
              state[AvailableAction[optim_action]+1] <- 0
              state[1] <- -AvailableAction[optim_action]
              time <- time + next_event_time # update time}
              depart_time <- time
            }
            else { # an new failure arrival
              state[1] <- AvailableAction[optim_action]
              state[normal_set[which.min(next_failure_time)]+1] <- 1
              depart_time <- time
              time <- time + next_event_time # update time}
            }
          }
        }
      }
    }
  }
}
mean(total_cost[[1]])
mean(total_cost[[2]])
mean(total_cost[[3]])
mean(total_cost[[4]])
mean(total_cost[[5]])
mean(total_cost[[6]])
#save(list = c("allResult"), file = paste0("NumBaseStation", NumBaseStation, ".Rdata"))


