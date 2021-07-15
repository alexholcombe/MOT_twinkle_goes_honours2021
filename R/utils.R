
gabor_patch <- function(szDeg = 1, contr = 1, lambda = 6, bandwidth = 1, ppd = 50, theta = 45, phase = 0.25, trim = 0.005) {
  #   Detailed explanation goes here
  imSize <- szDeg * ppd;                   # image size: n X n
  

  # sigmadeg = 0.1;                        # gaussian standard deviation in degrees
  sigmadeg <- (sqrt(log(4))/(2*pi*6))*((2^bandwidth + 1)/(2^bandwidth - 1))
  sigma    <- sigmadeg * ppd               # gaussian standard deviation in pixels
  

  X <- 1:imSize;                           # X is a vector from 1 to imageSize
  X0 <- (X / imSize) - 0.5;                 # rescale X -> -.5 to .5
  
  freq <- lambda;                          # compute frequency from wavelength
  phaseRad <- (phase * 2* pi);             #convert to radians: 0 -> 2*pi
  
  XXYY <- pracma::meshgrid(X0, X0);        # 2D matrices
  Xm <- XXYY$X
  Ym <- XXYY$Y 
  
  thetaRad <- (theta / 360) * 2 * pi;      # convert theta (orientation) to radians
  
  Xt <- Xm * cos(thetaRad);                # compute proportion of Xm for given orientation
  Yt <- Ym * sin(thetaRad);                # compute proportion of Ym for given orientation
  XYt <-  Xt + Yt;                         # sum X and Y components
  XYf <- XYt * freq * 2*pi;                # convert to radians and scale by frequency
  grating <- contr * sin( XYf + phaseRad); # make 2D sinewave
  
  
  # gaussian envelope
  
  s <- sigma / imSize;                     # gaussian width as fraction of imageSize
  gauss <- 
    exp( -(((Xm^2)+(Ym^2)) / (2* s^2)) );  # formula for 2D gaussian 
  
  gauss[gauss < trim] <- 0;                 # trim around edges (for 8-bit colour displays)
  gabor <- grating * gauss;                 # dot-product
  
  return(gabor)
  
  
}

pink_noise_2d <- function (x = 750, y = 750, degreesPerImage = 15) {
  n <- x * y
  noise <- matrix(data = runif(n)-0.5, nrow = x, ncol = y)
  aa <- fft(noise)
  fft_noise <- mrbsizeR::fftshift(aa)
  u <- freqSpace(x);
  v <- freqSpace(y);
  u <- u/degreesPerImage;
  v <- v/degreesPerImage;
  mesh_data <- pracma::meshgrid(v,u);
  meshv <- mesh_data$X;
  meshu <- mesh_data$Y;
  
  browner <- sqrt(meshu^2 + meshv^2)^(-1);
  browner[is.infinite(browner)] <- 0; # remove DC
  brown_noise <- fft_noise*browner;
  
  bb <- mrbsizeR::ifftshift(brown_noise)
  
  unclipped_noise <- fft(bb, inverse = T) / length(bb)
  stopifnot(max(Im(unclipped_noise)) < 0.0001)
  unclipped_noise <- Re(unclipped_noise)
  # We clip the noise at 2 std above and below.
  std_UN <- sd(as.vector(unclipped_noise));
  clipped_noise <- unclipped_noise;
  clipped_noise[unclipped_noise >= 2*std_UN] <- 2*std_UN
  clipped_noise[unclipped_noise <= -2*std_UN] <- -2*std_UN
  
  #output
  return(unclipped_noise) 
  
}

freqSpace <- function(nx) {
# Output frequency space for matlab discrete fourier transforms
  negativeHalf = ceiling((nx-1)/2);
  positiveHalf = nx-negativeHalf-1;
  return (-negativeHalf:positiveHalf)/(2*negativeHalf);
}

rescale_noise <- function(im, newrms){
  # rescale_noise Rescales noise by different factor
  m     <- mean(im)
  sd    <- sd(im)
  newsd <- 128 * newrms
  im_new <- 128 + (im - m) * (newsd / sd)
  return(im_new)
}

rescale_target <- function(tgt, new_contr){
  tgt <- tgt * new_cotnr
  return(im_new)
}

create_protocol_detection <- function(prot_id,t_contr, ntrials_per_level = 50, nNoises = 500) {
  p <- expand.grid(prot_id=prot_id, t_contr = t_contr, target_present = c(0,1)) %>% 
    as_tibble() %>% 
    mutate(t_contr = if_else(target_present == 1, t_contr, NA_integer_)) %>% 
    mutate(correctAns = if_else(target_present == 1, "left", "right"))
  p <- purrr::map_dfr(seq_len(ntrials_per_level), ~p) %>% 
    mutate(noise_id = sample(500, n(), replace = F))
  return(p)
}

create_protocol_mot_pilot <- function(prot_id, t_contr_id, ntrials_per_level = 40, start_time = 2) {
  block_order <- sample(t_contr_id)

  p <- tibble(t_contr = factor(rep(block_order, each = ntrials_per_level), levels = block_order)) %>% 
    mutate(prot_id=prot_id, trajectory_id = rep(1:ntrials_per_level, times = length(t_contr_id))) %>% 
    mutate(start_time = runif(n(),0,2) %>% round(2)) %>% 
    group_by(t_contr) %>% 
    sample_frac(1) %>% 
    ungroup() %>% 
    mutate(t_contr = as.character(t_contr) %>% as.numeric(), trial_id = 1:n(), noise_id = 1:n()) %>% 
    select(prot_id, trial_id, noise_id, trajectory_id, t_contr, start_time)
  return(p)
}

create_protocol_mot_exp1 <- function(prot_id, t_contr_id, mark_nomark, ntrials_per_level = 30) {
  # first noiseMOT than MOT
  
  block_order <- c(t_contr_id,0)
  block_order2 <- rep(block_order, each = 2)
  mark_nomark2 <- c(mark_nomark,"mark","noMark")
  df <- tibble(id = rep(1:4, each = 2), t_contr = block_order2)
  df2 <- df %>% 
    group_by(id) %>% 
    mutate(mark = if_else(row_number() <= n()/2,"first", "second")) %>% 
    ungroup() %>% 
    mutate(mark_val = mark_nomark2) %>% 
    select(-id)
  
  p <- tibble(t_contr = factor(rep(block_order2, each = ntrials_per_level), levels = block_order)) %>% 
    mutate(prot_id=prot_id, trajectory_id = 1:n(),noise_id = 1:n()) %>% 
    group_by(t_contr) %>% 
    sample_frac(1) %>% 
    mutate(mark = if_else(row_number() <= n()/2,"first", "second")) %>% 
    ungroup() %>% 
    mutate(t_contr = as.character(t_contr) %>% as.numeric(), trial_id = 1:n()) %>% 
    left_join(df2,by = c("t_contr", "mark")) %>% 
    rename(mark_type = mark_val) %>% 
    select(prot_id, trial_id, noise_id, trajectory_id, t_contr, mark_type)
    
    
  p
    
  return(p)
}

create_protocol_exp2_detection_range_estimation <- function() {
  contr_levels <- seq(0.12, 0.4, length.out = 8 )
  p <- tibble(contr_level = rep(contr_levels,times = 5))
  
  p
  
  return(p)
}

create_protocol_exp2_detection_calibration <- function(x = used_contrasts) {
  
  p <- tibble(contr_level = rep(used_contrasts,times = 5))
  
  p
  
  return(p)
}


perm <- function(v) {
  n <- length(v)
  if (n == 1) v
  else {
    X <- NULL
    for (i in 1:n) X <- rbind(X, cbind(v[i], perm(v[-i])))
    X
  }
}

plot_fit <- function(fit_par, df,subject_id = NULL) {
  if(is.null(subject_id)) {
    subject_id <- "NA"
  }
  df_summ <- df %>% group_by(contr_level) %>% summarize(m = mean(correct))  
  x <- seq(0.001, 1, by = 0.0001)
  y <- p_corr_fun(x, fit_par[1], fit_par[2])
  df_pred <- tibble(x,y)
  df_pred %>% 
    ggplot(aes(x,y)) + 
    geom_line() + scale_x_log10() +
    geom_point(data = df_summ, aes(x= contr_level, y = m)) + 
    theme(aspect.ratio = 1) +
    ggtitle(sprintf("Subject %d",subject_id)) + 
    ylim(0.4,1)
  
}


p_corr_fun <- function(c,c_T,beta) {
  pnorm(0.5 * (c/c_T)^beta)
}

inv_p_corr_fun <- function(p,c_T,beta) {
  return((qnorm(p)/0.5)^(1/beta) * c_T)
}

dprime_fun <- function(d, c_T,beta) {
  return(c_T * d^(1 / beta))
}


nll <- function(params,contrasts, correct) {
  c_T <- params[1]
  beta <- params[2]
  like <- 0
  if (c_T > 0 & c_T <= 1 & beta > 0) {
    like <- like + sum(log(p_corr_fun(contrasts[correct == 1], c_T, beta)))
    like <- like + sum(log(1 - p_corr_fun(contrasts[correct == 0], c_T, beta)))
    return(-like)
  } else {
    return(Inf)
  }
}

iterate_nll <- function(...) {
  current <- tibble(...)
  
  # return
  current %>%
    mutate(nll = nll(c(current$cT[1],current$beta[1]),current$data$contr_level, current$data$correct))
}

f <- function(...) {
  current <- tibble(...)
  # do cool stuff and access content from current row with
  plot_fit(fit_par = c(current$cT[1], current$beta[1]), current$data, subject_id = current$subject_id) %>% print()
  # no return
}

mle_SDT <- function(df) {
  contrasts <- df$contr_level
  correct   <- df$correct
  
  init_cT <- 0.12
  init_beta <- 1.6
  
  op <- optim(c(init_cT,init_beta), fn = nll, contrasts = contrasts, correct = correct)
  return(tibble(cT=op$pa[1],beta=op$par[2]))
}

mle_SDT_based_on_range_estim_results <- function(df) {
  # this function is intended to run from Rscript, thus we do not use tidyverse
  contrasts <- df$contr_level
  correct   <- df$correct
  
  init_cT <- 0.06
  init_beta <- 0.9
  
  op <- optim(c(init_cT,init_beta), fn = nll, contrasts = contrasts, correct = correct)
  return(data.frame(cT=op$pa[1],beta=op$par[2]))
}

change_contrasts_protocols <- function(px, contrasts_for_mot) {
  stopifnot(length(contrasts_for_mot)==unique(px$t_contr) %>% length())
  contr_tbl <- tibble(t_contr = px$t_contr %>% unique() %>% sort(), t_contr_val = contrasts_for_mot)
  px <- px %>% left_join(contr_tbl, by = "t_contr")
  return(px)
}

blend_images <- function(noise, gabor) {
  
  return(motrack::add_image_additive(noise, gabor, nrow(noise) / 2, ncol(noise) / 2))
}

cha<-function(x,y){
  chull(x,y)->i
  return(splancs::areapl(cbind(x[i],y[i])))
}

describe_trajectory <- function(traj1) {
  
  traj1 %>% group_by(t,target) %>% mutate(target_center_x = mean(valuex),
                                          target_center_y = mean(valuey),
                                          dist_from_center = sqrt(valuex^2+valuey^2), 
                                          dist_from_tgt_center = sqrt((valuex - target_center_x)^2 + (valuey - target_center_y)^2),
                                          chull = cha(valuex,valuey)) %>%
    arrange(trajectory_id,t) %>% 
    group_by(target) %>% 
    summarize(dist_from_center = mean(dist_from_center),
              dist_from_tgt_center_mean = mean(dist_from_tgt_center),
              dist_from_tgt_center_sd   = sd(dist_from_tgt_center),
              dist_from_tgt_center_min   = min(dist_from_tgt_center),
              dist_from_tgt_center_max   = max(dist_from_tgt_center),
              chull_mean =  mean(chull),
              chull_sd =  sd(chull))
  
}
ch_inter <- function(df1) {
  i <- geometry::intersectn(cbind(df1$valuex[df1$target],df1$valuey[df1$target]),cbind(df1$valuex[!df1$target],df1$valuey[!df1$target])) 
  if(is.null(i$ch$area)) {
    return (0)
  } else {
    i$ch$area
  }
  
}

describe_trajectory_tgtdist <- function(traj1) {
  
  tt <- traj1 %>% group_by(trajectory_id, t) %>% do(ch = ch_inter(.)) %>% mutate(ch = unlist(ch))
  
  tt <- tt %>% mutate(ch = unlist(ch))
 
  traj1_pw <- traj1 %>% left_join(traj1, by = c("trajectory_id","t")) %>% mutate(d = sqrt((valuex.x - valuex.y)^2 + (valuey.x - valuey.y)^2)) %>% filter(x.x < x.y) %>% 
    left_join(tt, by = c("trajectory_id","t"))
  traj1_pw %>% filter(target.x != target.y) %>% 
    group_by(trajectory_id) %>% 
    summarize(min_tgtdist = min(d),
              mean_tgtdist = mean(d),
              lower_quar_tgtdist = quantile(d,0.25),
              min_ch_inter = min(ch),
              mean_ch_inter = mean(ch),
              med_ch_inter = median(ch))
  
}

reshape_trajectories <- function(df) {
  
  df_x <- df %>% select(-starts_with("Y")) %>% gather("x","valuex",starts_with("X"))
  df_y <- df %>% select(-starts_with("X")) %>%  gather("y","valuey",starts_with("Y"))
  df1 <- df_x %>% left_join(df_y, by = c("trajectory_id", "t")) %>% 
    mutate(x = str_remove(x, "X"),
           y = str_remove(y, "y")) %>% 
    filter(x == y) %>% 
    mutate(target = x < 5) %>% select(-y)
  df1
}

download_files <- function(df, local_data_pth, should_overwrite = T) {
  # we need to set correct class as the current version of osfr does not works with dplyr properly
  class(df) <- c("osf_tbl_file","osf_tbl", class(df)) 
  df %>% 
    rowwise() %>% 
    do(osf_retrieve_file(.$id) %>% 
         osf_download(path = file.path(local_data_pth, .$name), 
                      overwrite = should_overwrite))
}

generate_noise_exp2 <- function(seed_id, out_pth,n_noises) {
  tm <- FDhelpers::create.time.measure(n_noises)
  
  for(i in 1:n_noises) {
    noise <- pink_noise_2d()
    noise <- rescale_noise(noise, 0.1)
    noise[1,1] <- 255
    noise[1,2] <- 0
    noise2 <- round(noise)
    png::writePNG(noise2, file.path(out_pth, sprintf("noise%03d.png",i)))
    
    tm <- FDhelpers::update.tm(tm)
    FDhelpers::print.tm(tm)
  }
  
  
}
generate_gabor_exp2 <- function(out_dir, file_name) {
  
  png::writePNG(gabor_patch(), file.path(out_dir, file_name))
  
}