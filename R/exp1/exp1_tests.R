library(testthat)
library(tidyverse)

test_that("verify protocol", {
          t_contr_all <- perm(t_contr_id_mot)
          mark_nomark_raw <-  expand.grid(c("mark","noMark"), c("mark","noMark"),c("mark","noMark")) %>% as.matrix()
          mark_nomark_all <- matrix(NA, ncol = 6, nrow = 8)
          for (i in 1:nrow(mark_nomark_raw)) {
            mnm <- c(mark_nomark_raw[i,])
            mark_nomark_raw_complement <- if_else(mnm == "mark", "noMark","mark")
            mark_nomark_all[i,c(1,3,5)] <- mnm
            mark_nomark_all[i,c(2,4,6)] <- mark_nomark_raw_complement
          }
          prot_id <- 10
          p <- create_protocol_mot_exp1(prot_id, t_contr_all[1,],mark_nomark_all[1,])
          expect_true(all(p$prot_id == prot_id))
          expect_true(all(sort(unique(p$t_contr)) == 0:3))
          expect_true(all(sort(p$noise_id) == 1:nrow(p)))
          expect_true(all(sort(p$trajectory_id) == 1:nrow(p)))
          expect_equal(tail(unique(p$t_contr),1),0)
          expect_true(all(table(p$mark_type,p$t_contr) == 30))
          expect_true(all(colSums(table(p$mark_type,p$trial_id)) == 1))
          expect_true(all(p$trial_id == 1:nrow(p)))
          expect_equal(p$noise_id,p$trajectory_id)
          }
          )
