# fit data
library(quickpsy)
predict_weibull <- function(x, guess = 0.5, lapse = 1-acc[5]) {
  guess+(1-guess - lapse)*weibull_fun(x,c(fit$par$par[1],fit$par$par[2]))
}

acc <- c(2,2.28,2.74,3.01,3.05)/4

acc
yess <- c(rbinom(4000,1,acc[1]),
          rbinom(4000,1,acc[2]),
          rbinom(4000,1,acc[3]),
          rbinom(4000,1,acc[4]),
          rbinom(4000,1,acc[5]))

cs <- c(0,0.1, 0.16, 0.21, 0.25)

cs2 <- rep(cs, each =4000)
df <- tibble(cs = cs2, yess = yess)

df2 <- tibble(cs = cs, acc = acc)
fit <- quickpsy::quickpsy(df2, cs, acc,thresholds = F,fun = weibull_fun, guess = 0.5, lapses = 1-acc[5]) 
plot(fit)

target_p <- c(0.6,0.65,0.7)

# values find by hand
c1 <- 0.11679
c2 <- 0.14049
c3 <- 0.16675
new_cs <- c(c1,c2,c3)
stopifnot(abs(target_p-predict_weibull(new_cs)) < tol)


