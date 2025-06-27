
# model 1: cultural traditions---------------------------

# parameters:
# N, number of individuals per population
# p, number of populations (so total number of individuals is pN)
# a, prob of social learning each trait of random target
# mu, expected number of traits to be invented
# m, expected number of individuals migrating out of population
# t_max, maximum number of timesteps

# outcomes:
# S, total number of unique traits in the population
# s, mean cultural similarity between groups (when p>1)


kempe_model1 <- function(N = 100, p = 5, a = 0.9, mu = 0.1, m = 0, t_max = 5000) {
  
  # trait matrix: each column is an agent, each row is a trait. 1=know that trait
  traits <- matrix(0, nrow = 100, ncol = p*N) 
  
  # record total number of traits in each population
  num_traits <- matrix(NA, nrow = t_max, ncol = p)
  
  # record final traits known in each population
  final_traits <- matrix(NA, nrow = 100, ncol = p)
  
  for (t in 1:t_max) {
    
    for (pop in 1:p) {
      
      # 1. pick random individual and replace with naive
      new_agent <- (pop-1)*N + sample(N, 1)
      traits[,new_agent] <- 0
      
      # 2. naive picks random other agent and copies their traits with prob a per trait
      demonstrator <- sample(setdiff(seq((pop-1)*N + 1, (pop-1)*N + N), new_agent), 1)
      demonstrator_traits <- which(traits[,demonstrator] == 1)
      copy_prob <- runif(nrow(traits)) < a
      for (i in demonstrator_traits) if (copy_prob[i]) traits[i,new_agent] <- 1
      
      # 3. naive innovates new traits with expectation mu, using poisson distribution
      num_innovate <- rpois(1, mu)
      if (sum(num_innovate) > 0) {
        for (i in 1:num_innovate) {
          traits[which(rowSums(traits[,((pop-1)*N+1):((pop-1)*N+N)]) == 0)[1], new_agent] <- 1
        }
      }
      
      # 4. naive migrates with prob m/2 to random other population
      if (p > 1 & runif(1) < (m/2)) {
        if (p == 2) destination <- (1:p)[-pop]
        if (p > 2) destination <- sample((1:p)[-pop], 1)
        return_migrant <- sample(((destination-1)*N+1):((destination-1)*N+N), 1)
        return_migrant_traits <- traits[,return_migrant]
        traits[,return_migrant] <- traits[,new_agent]
        traits[,new_agent] <- return_migrant_traits
      }
      
      # store number of traits for this population
      num_traits[t,pop] <- sum(rowSums(traits[,((pop-1)*N+1):((pop-1)*N+N)]) > 0)
      
      # store final traits if t=t_max
      if (t == t_max) final_traits[,pop] <- rowSums(traits[,((pop-1)*N+1):((pop-1)*N+N)]) > 0
    }
    
  }
  
  # expected number of traits, S, from Strimling et al. (2009) equation 3
  expected_S <- (mu*N/a)*log((1/(1-a))) + mu / (1-a)
  
  # observed mean number of traits in last 50% of timesteps
  observed_S <- mean(tail(num_traits, t_max/2))
  
  # observed final number of traits at last timestep
  S <- mean(tail(num_traits, 1))
  
  # mean similarity between populations, s
  # number of traits in either population / number of traits in both populations
  # averaged across all possible pairs of populations
  s <- 0
  if (p > 1) {
    
    for (i in 1:(p-1)) {
      for (j in (i+1):p) {
        # avoid divide by zero when neither group knows any traits
        if (sum(final_traits[,i] + final_traits[,j] > 0) > 0) {
          s <- s + sum(final_traits[,i] + final_traits[,j] > 1) / sum(final_traits[,i] + final_traits[,j] > 0)
        }
      }
    }
    
    s <- s / choose(p,2)
  }
  
  # output
  list(num_traits = num_traits,
       final_traits = final_traits,
       expected_S = expected_S,
       observed_S = observed_S,
       S = S,
       s = s)
}

## run model 1-------------------------------

kempe_model1_output <- kempe_model1()

## fig1a, time series of S-----------------------------------

kempe_fig1a <- function(kempe_model1_output) {
  
  # plot number of traits over time
  cols <- c("#E64B35","#4DBBD5","#00A087","#3C5488","#F39B7F","#8491B4","#91D1C2","#DC0000","#7E6148","#B09C85")
  
  num_traits <- kempe_model1_output$num_traits
  
  plot(num_traits[,1], type = 'l', col = cols[1],
       ylab = expression(Number~of~cultural~traits~italic(bar(S))),
       xlab = "Timestep")
  if (ncol(num_traits) > 1) {
    for (pop in 2:ncol(num_traits)) lines(num_traits[,pop], type = 'l', col = cols[pop])
    lines(rowMeans(num_traits), type = 'l', lwd = 2)
  }
  abline(h = kempe_model1_output$expected_S, lty = 2)
  
}

kempe_fig1a(kempe_model1_output)

data.frame(expected_S = kempe_model1_output$expected_S, 
           observed_S = kempe_model1_output$observed_S)

# save plot
png(file = "fig1a.png", width = 7, height = 5, units = "in", res = 600)
kempe_fig1a(kempe_model1_output)
dev.off()

## fig 1b, example trait profile for p groups--------------------

kempe_fig1b <- function(kempe_model1_output) {
  
  # highest known trait
  h <- max(which(rowSums(kempe_model1_output$final_traits) > 0))
  
  # truncate data to traits known by at least one population
  profiles <- kempe_model1_output$final_traits[1:h,]
  p <- ncol(profiles)
  
  plot(1, type="n", xlab="Trait", ylab="Population", 
       xlim=c(0, h), ylim=c(0, p),
       yaxt = "n", xaxt = "n",
       xaxs = "i", yaxs = "i")
  axis(1, at = 0.5:(h-0.5),
       labels = 1:h,
       tick = F, las = 1)
  axis(2, at = 0.5:(p-0.5),
       labels = 1:p,
       tick = F, las = 1)
  
  for (i in 1:h) {
    for (j in 1:p) {
      if (profiles[i,j]) rect(i-1,j-1,i,j, col = "grey", lwd = 0)
    }
  }
  
  grid(nx = h, ny = p, lty = 1)
  
}

kempe_fig1b(kempe_model1_output)

# save plot
png(file = "fig1b.png", width = 7, height = 4, units = "in", res = 600)
kempe_fig1b(kempe_model1_output)
dev.off()


## figure 2a, s against N------------------------

fig2a_data <- data.frame(N = seq(10,100,10),
                         s = rep(0, 10))
r_max <- 100  # independent runs per data point
counter <- 0

for (i in 1:nrow(fig2a_data)) {
  for (r in 1:r_max) {
    
    counter <- counter + 1
    cat("Run", counter, "of", nrow(fig2a_data)*r_max, fill=T)
    
    fig2a_data$s[i] <- fig2a_data$s[i] + kempe_model1(a = 0.9, p = 5, mu = 0.1, 
                                                      m = 0, t_max = 2000, 
                                                      N = fig2a_data$N[i])$s
    
  }
}
fig2a_data$s <- fig2a_data$s / r_max

# save data
save(fig2a_data, file = "fig2a_data.Rda")

# load data
#load(file = "fig2a_data.Rda")

kempe_fig2a <- function(fig2a_data) {
  
  plot(y = fig2a_data$s, x = fig2a_data$N, type = 'b', 
       ylab = expression(Mean~similarity~between~populations~italic(bar(s))),
       xlab = expression(Population~size~italic(N)))
  
}

kempe_fig2a(fig2a_data)

# save plot
png(file = "fig2a.png", width = 7, height = 5, units = "in", res = 600)
kempe_fig2a(fig2a_data)
dev.off()


## figure 2b, s against a, no migration---------------------------------

fig2b_data <- data.frame(a = seq(0.1,0.9,0.1),
                         s = rep(0, 9))
r_max <- 1500  # independent runs per data point
counter <- 0

for (i in 1:nrow(fig2b_data)) {
  for (r in 1:r_max) {
    
    counter <- counter + 1
    cat("Run", counter, "of", nrow(fig2b_data)*r_max, fill=T)
    
    fig2b_data$s[i] <- fig2b_data$s[i] + kempe_model1(N = 50, p = 5, mu = 0.1, 
                                                      t_max = 2000, m = 0,
                                                      a = fig2b_data$a[i])$s
    
  }
}
fig2b_data$s <- fig2b_data$s / r_max

# save data
save(fig2b_data, file = "fig2b_data.Rda")

# load data
#load(file = "fig2b_data.Rda")

kempe_fig2b <- function(fig2b_data) {
  
  plot(y = fig2b_data$s, x = fig2b_data$a, type = 'b', 
       ylab = expression(Mean~similarity~between~populations~italic(bar(s))),
       xlab = expression(Accuracy~of~social~learning~italic(a)))
  
}

kempe_fig2b(fig2b_data)

# save plot
png(file = "fig2b.png", width = 7, height = 5, units = "in", res = 600)
kempe_fig2b(fig2b_data)
dev.off()
  

## fig 2c, contour plot for s, a and N----------------------------

N_values <- seq(10,100,10)
a_values <- seq(0.1,0.9,0.1)
r_max <- 100

# cols are N values, rows are a values 
s_matrix <- matrix(0, ncol = length(N_values), nrow = length(a_values))

counter <- 0

for (N in N_values) {
  for (a in a_values) {
    for (r in 1:r_max) {
      
      counter <- counter + 1
      cat("Run", counter, "of", length(s_matrix)*r_max, fill=T)
      
      s_matrix[which(a_values == a), which(N_values == N)] <- 
        s_matrix[which(a_values == a), which(N_values == N)] +
        kempe_model1(t_max=2000, p=5, mu=0.1, m=0, N=N, a=a)$s
      
    }
  }
}
s_matrix <- s_matrix / r_max

# save data
save(s_matrix, file = "s_matrix.Rda")

# load data
#load(file = "s_matrix.Rda")

kempe_fig2c <- function(a_values, N_values, s_matrix) {
  
  contour(x = a_values, y = N_values, z = s_matrix,
          levels = seq(0.3,0.65,0.1),
          xlab = expression(Accuracy~of~social~learning~italic(a)),
          ylab = expression(Population~size~italic(N)),
          xaxs = "i", yaxs = "i")
  
}

kempe_fig2c(a_values, N_values, s_matrix)

# save plot
png(file = "fig2c.png", width = 7, height = 5, units = "in", res = 600)
kempe_fig2c(a_values, N_values, s_matrix)
dev.off()

## fig 3a, m and S------------------------

fig3_data <- data.frame(m = seq(0,0.5,0.05),
                         S = rep(0, 11),
                         s = rep(0, 11))

r_max <- 100  # independent runs per data point
counter <- 0

for (i in 1:nrow(fig3_data)) {
  for (r in 1:r_max) {
    
    counter <- counter + 1
    cat("Run", counter, "of", nrow(fig3_data)*r_max, fill=T)
    
    d <- kempe_model1(a = 0.9, p = 5, mu = 0.1, 
                      t_max = 2000, N = 50,
                      m = fig3_data$m[i])
    fig3_data$S[i] <- fig3_data$S[i] + d$S
    fig3_data$s[i] <- fig3_data$s[i] + d$s
    
  }
}

fig3_data$S <- fig3_data$S / r_max
fig3_data$s <- fig3_data$s / r_max

# save data
save(fig3_data, file = "fig3_data.Rda")

# load data
#load(file = "fig3_data.Rda")

kempe_fig3a <- function(fig3_data) {
  
  plot(y = fig3_data$S, x = fig3_data$m, type = 'b', 
       ylab = expression(Mean~number~of~traits~italic(bar(S))),
       xlab = expression(Migration~rate~italic(m)))
  
}

kempe_fig3a(fig3_data)

# save plot
png(file = "fig3a.png", width = 7, height = 5, units = "in", res = 600)
kempe_fig3a(fig3_data)
dev.off()


## fig 3b, m and s------------------------

kempe_fig3b <- function(fig3_data) {
  
  plot(y = fig3_data$s, x = fig3_data$m, type = 'b', 
       ylab = expression(Mean~similarity~between~populations~italic(bar(s))),
       xlab = expression(Migration~rate~italic(m)))
  
}

kempe_fig3b(fig3_data)

# save plot
png(file = "fig3b.png", width = 7, height = 5, units = "in", res = 600)
kempe_fig3b(fig3_data)
dev.off()


# model 2: cumulative culture---------------------

# parameters:
# N, number of individuals
# a, prob of social learning each trait level of n random targets
# n, number of demonstrators
# mu, expected number of traits to be invented
# t_max, maximum number of timesteps

# output
# l, mean trait level per timestep
# final traits at end of sim

kempe_model2 <- function(N = 100, a = 0.7, mu = 0.1, n = 3, t_max = 10000) {
  
  # trait vector: trait level of each agent
  traits <- rep(0,N) 
  
  # record mean trait level at each timestep
  mean_level <- rep(NA, t_max)
  
  for (t in 1:t_max) {
      
      # 1. pick random individual and replace with naive
      new_agent <- sample(N, 1)
      traits[new_agent] <- 0
      
      # 2. naive picks n random other agents and copies their traits with prob a per trait
      dems <- sample((1:N)[-new_agent], n)
      # for each demonstrator i
      for (i in dems) {
        # if dem has higher level than new agent
        if (traits[i] > traits[new_agent]) {
          # for each higher level j
          for (j in (traits[new_agent]+1):traits[i]) {
            # with probability a
            if (runif(1) < a) {
              # copy trait
              traits[new_agent] <- j
            } else {
              # otherwise stop and move to next dem
              break
            }
          }
        }
      }
      
      # 3. naive innovates new trait with probability mu
      if (runif(1) < mu) traits[new_agent] <- traits[new_agent] + 1
      
      # store mean trait level l
      mean_level[t] <- mean(traits)
    
  }
  
  # output
  list(mean_level = mean_level,
       final_traits = traits)
  
}

## run model 2-----------------------

kempe_model2_output <- kempe_model2()

## fig 4a, time series for l---------------

kempe_fig4a <- function(kempe_model2_output) {
  
  plot(kempe_model2_output$mean_level, type = 'l',
       ylab = expression(Average~level~~italic(bar(l))),
       xlab = "Timestep")
  
}

kempe_fig4a(kempe_model2_output)

# save plot
png(file = "fig4a.png", width = 7, height = 5, units = "in", res = 600)
kempe_fig4a(kempe_model2_output)
dev.off()


## fig 4b, trait distribution---------------

kempe_fig4b <- function(kempe_model2_output) {
  
  barplot(table(kempe_model2_output$final_traits) / length(kempe_model2_output$final_traits),
          ylab = "Proportion of population",
          xlab = "Level")
  
}

kempe_fig4b(kempe_model2_output)

# save plot
png(file = "fig4b.png", width = 6, height = 5, units = "in", res = 600)
kempe_fig4b(kempe_model2_output)
dev.off()


## fig 5a, l and n---------------

fig5a_data <- data.frame(n = 1:10,
                         l = rep(0, 10))
r_max <- 20  # independent runs per data point
counter <- 0

for (i in 1:nrow(fig5a_data)) {
  for (r in 1:r_max) {
    
    counter <- counter + 1
    cat("Run", counter, "of", nrow(fig5a_data)*r_max, fill=T)
    
    fig5a_data$l[i] <- fig5a_data$l[i] + tail(kempe_model2(a = 0.9, mu = 0.1, 
                                                           t_max = 10000, N = 100,
                                                           n = fig5a_data$n[i])$mean_level, 1)
    
  }
}
fig5a_data$l <- fig5a_data$l / r_max

# save data
save(fig5a_data, file = "fig5a_data.Rda")

# load data
#load(file = "fig5a_data.Rda")

kempe_fig5a <- function(fig5a_data) {
  
  plot(y = fig5a_data$l, x = fig5a_data$n, type = 'b', 
       ylab = expression(Average~level~italic(bar(l))),
       xlab = expression(Number~of~cultural~models~italic(n)))
  
}

kempe_fig5a(fig5a_data)

# save plot
png(file = "fig5a.png", width = 7, height = 5, units = "in", res = 600)
kempe_fig5a(fig5a_data)
dev.off()


## fig 5b, l and a---------------

fig5b_data <- data.frame(a = seq(0,0.95,0.05),
                         l = rep(0, length(seq(0,0.95,0.05))))
r_max <- 20  # independent runs per data point
counter <- 0

for (i in 1:nrow(fig5b_data)) {
  for (r in 1:r_max) {
    
    counter <- counter + 1
    cat("Run", counter, "of", nrow(fig5b_data)*r_max, fill=T)
    
    fig5b_data$l[i] <- fig5b_data$l[i] + tail(kempe_model2(n = 3, mu = 0.1, 
                                                           t_max = 10000, N = 100,
                                                           a = fig5b_data$a[i])$mean_level, 1)
    
  }
}
fig5b_data$l <- fig5b_data$l / r_max

# save data
save(fig5b_data, file = "fig5b_data.Rda")

# load data
#load(file = "fig5b_data.Rda")

kempe_fig5b <- function(fig5b_data) {
  
  plot(y = fig5b_data$l, x = fig5b_data$a, type = 'b', 
       ylab = expression(Average~level~italic(bar(l))),
       xlab = expression(Accuracy~of~social~learning~italic(a)))
  
}

kempe_fig5b(fig5b_data)

# save plot
png(file = "fig5b.png", width = 7, height = 5, units = "in", res = 600)
kempe_fig5b(fig5b_data)
dev.off()


## fig 5c, contour of l, a and n----------------------

n_values <- 1:10
a_values <- seq(0,0.95,0.05)
r_max <- 1

# cols are a values, rows are n values 
l_matrix <- matrix(0, ncol = length(a_values), nrow = length(n_values))

counter <- 0

for (a in a_values) {
  for (n in n_values) {
    for (r in 1:r_max) {
      
      counter <- counter + 1
      cat("Run", counter, "of", length(l_matrix)*r_max, fill=T)
      
      l_matrix[which(n_values == n), which(a_values == a)] <- 
        l_matrix[which(n_values == n), which(a_values == a)] +
        tail(kempe_model2(n = n, a = a, mu = 0.1, 
                          t_max = 10000, N = 100)$mean_level, 1)
      
    }
  }
}
l_matrix <- l_matrix / r_max

# save data
save(l_matrix, file = "l_matrix.Rda")

# load data
#load(file = "l_matrix.Rda")

kempe_fig5c <- function(n_values, a_values, l_matrix) {
  
  contour(x = n_values, y = a_values, z = l_matrix,
          levels = c(1,2,3,4,5,10,20,50),
          xlab = expression(Number~of~cultural~models~italic(n)),
          ylab = expression(Accuracy~of~social~learning~italic(a)),
          xaxs = "i", yaxs = "i")
  
}

kempe_fig5c(n_values, a_values, l_matrix)

# save plot
png(file = "fig5c.png", width = 7, height = 5, units = "in", res = 600)
kempe_fig5c(n_values, a_values, l_matrix)
dev.off()


## fig 5d (extra), l and mu---------------

mu_values <- seq(0.1,1,0.1)

fig5d_data <- data.frame(mu = mu_values,
                         l = rep(0, length(mu_values)))
r_max <- 20  # independent runs per data point
counter <- 0

for (i in 1:nrow(fig5d_data)) {
  for (r in 1:r_max) {
    
    counter <- counter + 1
    cat("Run", counter, "of", nrow(fig5d_data)*r_max, fill=T)
    
    fig5d_data$l[i] <- fig5d_data$l[i] + tail(kempe_model2(n = 3, a = 0.7, 
                                                           t_max = 10000, N = 100,
                                                           mu = fig5d_data$mu[i])$mean_level, 1)
    
  }
}
fig5d_data$l <- fig5d_data$l / r_max

# save data
save(fig5d_data, file = "fig5d_data.Rda")

# load data
#load(file = "fig5d_data.Rda")

kempe_fig5d <- function(fig5d_data) {
  
  plot(y = fig5d_data$l, x = fig5d_data$mu, type = 'b', 
       ylab = expression(Average~level~italic(bar(l))),
       xlab = expression(Probability~of~innovation~italic(mu)),
       ylim = c(0,max(fig5d_data$l)*1.5))
  
}

kempe_fig5d(fig5d_data)

# save plot
png(file = "fig5d.png", width = 7, height = 5, units = "in", res = 600)
kempe_fig5d(fig5d_data)
dev.off()

