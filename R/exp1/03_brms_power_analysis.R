library(tidyverse)
library(brms)
library(broom)

theme_set(theme_classic())

# first, load data from pilot, compute correlation between conditions using this value
#cc <- df_all %>% 
#  group_by(t_contr, participant) %>% 
#  summarize(nCorrect = mean(nCorrect, na.rm = T)) %>% 
#  spread(t_contr,nCorrect) %>% 
# rename(t1 = `1`, t2 = `2`, t3 = `3`, t4 = `4`) %>% 
#  select(-participant) %>% 
#  as.matrix() %>% 
#  Hmisc::rcorr()
#cc <- cc$r
#mean(cc[upper.tri(cc)])
# average correlation is 0.7519376
# correlation matrix for random generation

sim_d <- function(seed, n) {
  
  
  set.seed(seed)
  
  d <-
    tibble(id = rep(1:n, times = 3),
           group     = rep(c("c1", "c2","c3"), each = n)) %>% 
    mutate(treatment1 = ifelse(group == "c2", 1, 0),
           treatment2 = ifelse(group == "c3", 1, 0),
           y          = NA)
  dd <- MASS::mvrnorm(n=n, mu=c(mu_c1, mu_c2, mu_c3), Sigma=cm, empirical=TRUE)
  
  d$y[d$group == "c1"] <- dd[1,]
  d$y[d$group == "c2"] <- dd[2,]
  d$y[d$group == "c3"] <- dd[3,]
  
  
  d
}

# define the means
mu_c1 <- 0
mu_c2 <- 0.12
mu_c3 <- 0.24 #
cm <- matrix(c(1.0000000, 0.7186702, 0.6176024, 
               0.7186702, 1.0000000, 0.8548623, 
               0.6176024, 0.8548623, 1.0000000), nrow = 3) 


n <- 50

set.seed(1)

d <- sim_d(1, 50)

ez::ezANOVA(data = d,
            dv = y,
            wid = id,
            within = group)

# first define means of groups, so traditional anova would result in medium effect (eta_p = 0.06)

n <- 100
n_sim <- 1000
eta_ps <- tibble(isim=1:n_sim,eta_p = NA)

for (i in 1:n_sim) {
  d <- sim_d(i,n)
  aov1 <- ez::ezANOVA(data = d,
              dv = y,
              wid = id,
              within = group,detailed = T)$ANOVA
  eta_ps$eta_p[eta_ps$isim == i] <- aov1$SSn[2]/(aov1$SSn[2]+aov1$SSd[2])
}
eta_ps$eta_p %>% mean()
sd(eta_ps$eta_p)/n_sim
quantile(eta_ps$eta_p, c(0.025,0.975))

get_prior(data = d,
          family = gaussian,
          y ~ 0 + intercept + treatment1 + treatment2)

fit <-
  brm(data = d,
      family = gaussian,
      y ~ 0 + intercept + treatment1 + treatment2,
      prior = c(prior(normal(0, 2), class = b),
                prior(student_t(3, 1, 1), class = sigma)),
      seed = 1)

plot(fit)

print(fit)

tidy(fit, prob = .95)



# how many simulations would you like?
n_sim <- 100

# this will help us track time
t1 <- Sys.time()

# here's the main event!
s <-
  tibble(seed = 1:n_sim) %>% 
  mutate(d    = map(seed, sim_d, n = 100)) %>% 
  mutate(fit  = map2(d, seed, ~update(fit, newdata = .x, seed = .y)))

t2 <- Sys.time()

s %>% 
  mutate(treatment = map(fit, tidy, prob = .95)) %>% 
  unnest(treatment) %>% 
  filter(term == "b_treatment1") %>% 
  
  ggplot(aes(x = seed, y = estimate, ymin = lower, ymax = upper)) +
  geom_hline(yintercept = c(0, .5), color = "white") +
  geom_pointrange(fatten = 1/2) +
  labs(x = "seed (i.e., simulation index)",
       y = expression(beta[1]))

s %>% 
  mutate(treatment = map(fit, tidy, prob = .95)) %>% 
  unnest(treatment) %>% 
  filter(term == "b_treatment1") %>% 
  mutate(check = ifelse(lower > 0, 1, 0)) %>% 
  summarise(power = mean(check))