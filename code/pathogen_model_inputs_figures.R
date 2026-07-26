# need to run parts of pathogen_model.R
# should separate out that script later ...
source("code/init.R")

# get onto times later....
surveillance_effort_points = alpha("blue", 1)
surv_cex = 5
# surv_pch = 0
# surv_pch = 15
surv_pch = 16
date_cex = 1
weight_column = "likelihood_weight"

purps = colorRampPalette(brewer.pal(9,"Purples"))(100)

# harsher palette of purples for surveillance effort:
extreme_purps = colorRampPalette(purps[c(seq(1,30,5),30:100)])(100)

library(lubridate)

dat <- read_rds("output/dat.rds")
# dat$df$x <- jitter(dat$df$x, amount = 0.5)
# dat$df$y <- jitter(dat$df$y, amount = 0.25)

# filter(dat$df, host_type == "mosquito") %>%
#   dplyr::select(x, y) %>%
#   plot()
# 
# filter(dat$df, host_type == "mosquito") %>%
#   mutate(x = jitter(x, amount = 0.2),
#          y = jitter(y, amount = 0.2)) %>%
#   dplyr::select(x, y) %>%
#   plot()

# jitter positions after grouping locations
jittered_locs <- dat$df %>%
  filter(!(host_type == "human" & cases == 0)) %>%
  dplyr::select(x, y) %>%
  unique() %>%
  mutate(xj = jitter(x, amount = 0.7),
         yj = jitter(y, amount = 0.4))
  
plot(st_geometry(aus_simple))
points(jittered_locs[,c("x", "y")])
points(jittered_locs[,c("xj", "yj")], col = "red")

source("secure_code/nearestLand.R")
idx <- which(is.na(raster::extract(hpop_blur, jittered_locs[,c("xj", "yj")])))
length(idx)

# this takes a while ...
landed <- nearestLand(jittered_locs[idx, c("x", "y")], 
                               hpop_blur, 
                               max_distance = 100000)

dat$df <- dat$df %>%
  left_join(jittered_locs, by = join_by(x==x, y==y)) %>%
  mutate(x = xj, y = yj)

# check everything was joined
filter(dat$df, is.na(xj)) %>% dplyr::select(host_type, virus_detection) %>% ftable()
# ^ (those are in fact -ves)
filter(dat$df, !is.na(xj)) %>% dplyr::select(host_type, virus_detection) %>% ftable()

############################################################################
# all date-handling comes from 2023_report_figures and the demogorgon script ...

all_humans_dated <- read.csv("data/secure/all_humans_with_semisampled_dates.csv") %>%
  mutate(diagnostic_detail = "PCR",
         virus_detection = "JEV positive")

all_dates = rbind(all_mozzies[,c("diagnostic_type", "diagnostic_detail", 
                                 "host_type", "virus_detection","exposure_weight")],
                  all_ferals[,c("diagnostic_type", "diagnostic_detail", 
                                "host_type", "virus_detection","exposure_weight")],
                  all_dompigs[,c("diagnostic_type", "diagnostic_detail", 
                                 "host_type", "virus_detection","exposure_weight")],
                  all_chickens[,c("diagnostic_type", "diagnostic_detail", 
                                  "host_type", "virus_detection","exposure_weight")],
                  all_humans_dated[,c("diagnostic_type", "diagnostic_detail", 
                                      "host_type", "virus_detection","exposure_weight")]) %>%
  left_join(diagnostic_weights_table,
            by=c("diagnostic_type", "diagnostic_detail"="diagnostic_endpoint", 
                 "host_type", "virus_detection")) %>%
  mutate(date=c(as.Date(all_mozzies$date),
                all_ferals$date,
                all_dompigs$date,
                all_chickens$date,
                all_humans_dated$date)) %>%
  mutate(ypos = case_when(host_type == "human" ~ 5,
                          host_type == "mosquito" ~ 4,
                          host_type == "domestic_pig" ~ 3,
                          host_type == "feral_pig" ~ 2,
                          host_type == "chicken" ~ 1),
         likelihood_weight = diagnostic_weight * exposure_weight) %>%
  drop_na(likelihood_weight) #%>%
  #filter(date > dmy("01/01/2021")) # there was a 2002 in the pigs at least

range(all_dates$date)
# noice .. just don't ask me how I attatched dates to the people ... 
# perhaps is better to leave no date? Idk? Would lose a bunch

###############################################################################

month_breaks <- seq(ymd("2020-02-01"),
                    max(as.Date(all_dates$date), na.rm=TRUE),
                    by="month")
quarter_breaks <- month_breaks[seq(3, length(month_breaks), 3)]

jit_fac = 2

# pal = brewer.pal(11, "RdYlBu")
# 
# neg_pal = colorRampPalette(c(pal[6],"blue"))(100)
# pos_pal = colorRampPalette(c(pal[6],"red"))(100)

neg_pal = alpha("blue", seq(0.01,0.5,length.out=100))
pos_pal = alpha("red", seq(0.01,0.5,length.out=100))


{png("figures/pathogen_inputs_points_only.png",
     width = 900*reso/72, 
     height = 1200*reso/72, 
     res = reso)
  #type = 'cairo')
  
  par(mfrow = c(3, 2),
      mar = c(0.1, 0.1, 2.1, 0.1), bty="n",
      oma=c(0,0,0,2))
  
  # humans
  plot(st_geometry(aus_simple), cex.main=1.6,
       main="Probable exposure sites of confirmed JEV in humans", 
       border='grey20')
  # points(dat$human_points, cex=surv_cex, col="red", pch=surv_pch)
  points(dat$df[dat$df$host_type == "human" & dat$df$cases > 0, c("x","y")], 
         col=pos_pal[dat$df[dat$df$host_type == "human" & dat$df$cases > 0, 
                            weight_column]*100], 
         pch=surv_pch, cex=surv_cex,)
  subfigure_label(par("usr"), 0.1, 0.9, "(a)", 1.6)
  
  # mosquitoes
  plot(st_geometry(aus_simple), lwd=1.2, border='grey20', cex.main=1.6,
       main="JEV positive/negative samples in mosquitoes")
  points(dat$df[dat$df$host_type == "mosquito" & dat$df$cases == 0, c("x","y")], 
         col=neg_pal[dat$df[dat$df$host_type == "mosquito" & dat$df$cases == 0, 
                            weight_column]*100], 
         pch=surv_pch)
  points(dat$df[dat$df$host_type == "mosquito" & dat$df$cases > 0, c("x","y")], 
         col=pos_pal[dat$df[dat$df$host_type == "mosquito" & dat$df$cases > 0, 
                            weight_column]*100], 
         pch=surv_pch)  
  subfigure_label(par("usr"), 0.1, 0.9, "(b)", 1.6)
  
  # domestic pigs
  plot(st_geometry(aus_simple), lwd=1.2, border='grey20', cex.main=1.6,
       main="JEV positive/negative samples in domestic pigs")
  points(dat$df[dat$df$host_type == "domestic_pig" & dat$df$cases == 0, c("x","y")], 
         col=neg_pal[dat$df[dat$df$host_type == "domestic_pig" & dat$df$cases == 0, 
                            weight_column]*100], 
         pch=surv_pch)
  points(dat$df[dat$df$host_type == "domestic_pig" & dat$df$cases > 0, c("x","y")], 
         col=pos_pal[dat$df[dat$df$host_type == "domestic_pig" & dat$df$cases > 0, 
                            weight_column]*100], 
         pch=surv_pch)  
  subfigure_label(par("usr"), 0.1, 0.9, "(c)", 1.6)
  
  # feral pigs
  plot(st_geometry(aus_simple), lwd=1.2, border='grey20', cex.main=1.6,
       main="JEV positive/negative samples in feral pigs")
  points(dat$df[dat$df$host_type == "feral_pig" & dat$df$cases == 0, c("x","y")], 
         col=neg_pal[dat$df[dat$df$host_type == "feral_pig" & dat$df$cases == 0, 
                            weight_column]*100], 
         pch=surv_pch)
  points(dat$df[dat$df$host_type == "feral_pig" & dat$df$cases > 0, c("x","y")], 
         col=pos_pal[dat$df[dat$df$host_type == "feral_pig" & dat$df$cases > 0, 
                            weight_column]*100], 
         pch=surv_pch)  
  subfigure_label(par("usr"), 0.1, 0.9, "(d)", 1.6)
  
  # chickens
  plot(st_geometry(aus_simple), lwd=1.2, border='grey20', cex.main=1.6,
       main="JEV positive/negative samples in sentinel chickens")
  points(dat$df[dat$df$host_type == "chicken" & dat$df$cases == 0, c("x","y")], 
         col=neg_pal[dat$df[dat$df$host_type == "chicken" & dat$df$cases == 0, 
                            weight_column]*100], 
         pch=surv_pch)
  points(dat$df[dat$df$host_type == "chicken" & dat$df$cases > 0, c("x","y")], 
         col=pos_pal[dat$df[dat$df$host_type == "chicken" & dat$df$cases > 0, 
                            weight_column]*100], 
         pch=surv_pch)  
  subfigure_label(par("usr"), 0.1, 0.9, "(e)", 1.6)
  
  par(mar=c(15.1,9.1,6.1,1.1), bty="o")
  plot(0, xlim = c(as.Date("2021-01-01"), as.Date("2023-12-04")), ylim=c(0.5,5.5), xlab="",
       main="JEV surveillance in time by host type", cex.main=1.6, ylab="",
       yaxt="n", xaxt="n", type="n", cex.lab=1.2)
  abline(v=quarter_breaks, col="grey90")
  abline(v=c(as.Date("2021-01-01"),as.Date("2022-01-01"),as.Date("2023-01-01")), col="grey60")
  # abline(v=as.Date("2022-03-04"), col="grey50", lty=2, lwd=2) # CDIoNS
  # plot all negatives
  points(as.Date(all_dates$date[all_dates$virus_detection == "JEV negative"]),
         jitter(all_dates$ypos[all_dates$virus_detection == "JEV negative"],
                factor=jit_fac), pch=surv_pch, cex=date_cex,
         col=neg_pal[round(all_dates$likelihood_weight[all_dates$virus_detection == "JEV negative"]*100)])
  # plot all positives
  points(as.Date(all_dates$date[all_dates$virus_detection == "JEV positive"]),
         jitter(all_dates$ypos[all_dates$virus_detection == "JEV positive"],
                factor=jit_fac), pch=surv_pch, cex=date_cex,
         col=pos_pal[round(all_dates$likelihood_weight[all_dates$virus_detection == "JEV positive"]*100)])
  # not sure why the humans disappear?
  # points(all_dates$date[all_dates$ypos == 5], cex=0.8,
  #        jitter(rep(5, length(all_dates$date[all_dates$ypos == 5])),
  #               factor=jit_fac*2), pch=surv_pch,
  #        col=pos_pal[round(all_dates$likelihood_weight[all_dates$virus_detection == "JEV positive"]*100)])

  axis(2, at=5:1, c("Human","Mosquito","Domestic pig","Feral pig","Chicken"), cex.axis=1.2, las=2)
  axis(1, at=quarter_breaks, format(quarter_breaks, format="%b %y"), las=2, cex.axis=1.2)
  mtext("Date", 1, line=4.5)
  mtext("(f)", 3, 2.5, adj=0, cex=1.1)
  
  
  par(new=TRUE, mar=c(0,0,0,0))
  plot(0, axes=FALSE, type="n", xlim=c(0,1), ylim=c(0,1), xlab="", ylab="")
  gradient.rect(0.38, 0.1,
                0.58, 0.12,
                col=rev(neg_pal),
                border='grey20',
                gradient='x')
  gradient.rect(0.58, 0.1,
                0.78, 0.12,
                col=(pos_pal),
                border='grey20',
                gradient='x')
  axis(1, at=seq(0.38,0.78,length.out=5), labels=seq(-1,1,0.5), line=-5.5)
  text(c(0.48, 0.68),c(0.03, 0.03), c("JEV Negative", "JEV Positive"), cex=1.2)
  text(0.58,0.15,"Likelihood weight", cex=1.4)
  
  # For legend over empty sixth panel:
  # par(mar=c(0.1,0.1,0.1,5.1))
  # plot(0, axes=FALSE, type="n", xlim=c(0,1), ylim=c(0,1), xlab="", ylab="")
  # gradient.rect(0.35, 0.5,
  #               0.58, 0.54,
  #               col=rev(neg_pal),
  #               border='grey20',
  #               gradient='x')
  # gradient.rect(0.58, 0.5,
  #               0.81, 0.54,
  #               col=(pos_pal),
  #               border='grey20',
  #               gradient='x')
  # axis(1, at=seq(0.38,0.78,length.out=5), labels=seq(-1,1,0.5), line=-20.5, cex.axis=1.6)
  # text(x = c(0.46, 0.7), y = c(0.41, 0.41), c("JEV Negative", "JEV Positive"), cex=1.6)
  # text(0.58,0.57,"Likelihood weight", cex=1.6)
  
  dev.off()}


# Just the surveillance density surfaces ...
{png("figures/pathogen_inputs_surface_only.png",
     width = 1000*reso/72, 
     height = 1100*reso/72, 
     res = reso)
  
  par(mfrow = c(3, 2),
      mar = c(0.1, 0.1, 3.1, 1), bty="n",
      oma=c(0,0,0,2))
  
  # humans
  plot(sqrt(hpop_blur/100), col=extreme_purps, legend.width=1.4,
       legend.args = list(text='Human population\nper square km', 
                          line=0.9, side=2, cex=1.1), 
       axis.args = list(at=c(sqrt(10),10,sqrt(250),sqrt(500),30), 
                        labels=c(10,100,250,500,900), cex.axis=1.2),
       box = FALSE, axes = FALSE, cex.main=1.6,
       main="Distance-weighted human population density")
  plot(st_geometry(aus_simple), add=TRUE, lwd=1.2, border='grey20')
  subfigure_label(par("usr"), 0.1, 0.9, "(a)", 1.6)
  
  # mosquitoes
  plot(#nat_mozz_density/100, # decemberist
    dat$densities$mosquito_sampling/100,
    box = FALSE, axes=FALSE, col=extreme_purps, 
    legend.width=1.4,
    legend.args = list(text='Mosquito catches screened\nper square km', 
                       line=0.9, cex=1.1, side=2), 
    axis.args = list(cex.axis=1.2),
    cex.main=1.6,
    main="Mosquito screening effort")
  plot(st_geometry(aus_simple), lwd=1.2, border='grey20', add=TRUE)
  subfigure_label(par("usr"), 0.1, 0.9, "(b)", 1.6)
  
  # domestic pigs
  plot(dat$densities$domestic_pig_sampling/100, box = FALSE, axes=FALSE, 
       col=extreme_purps,
       legend.args = list(text='Sampled domestic pigs\nper square km', 
                          cex=1.1, line=0.9, side=2), 
       axis.args = list(at=c(0,2e-4,4e-4,6e-4,8e-4), 
                        labels=c("0","2e-4","4e-4","6e-4","8e-4"),
                        cex.axis=1.2),
       legend.width=1.4, cex.main=1.6,
       main="Domestic pig surveillance density")
  plot(st_geometry(aus_simple), lwd=1.2, border='grey20', add=TRUE)
  subfigure_label(par("usr"), 0.1, 0.9, "(c)", 1.6)
  
  # feral pigs
  plot(dat$densities$feral_pig_sampling/100, box = FALSE, axes=FALSE, 
       col=extreme_purps,
       legend.args = list(text='Sampled feral pigs\nper square km', cex=1.1, 
                          line=0.9, side=2),
       axis.args = list(at=c(0,2e-3,4e-3,6e-3,8e-3),
                        labels=c("0","2e-3","4e-3","6e-3","8e-3"),
                        cex.axis=1.2),
       legend.width=1.4, cex.main=1.6,
       main="Feral pig screening effort")
  plot(st_geometry(aus_simple), lwd=1.2, border='grey20', add=TRUE) 
  subfigure_label(par("usr"), 0.1, 0.9, "(d)", 1.6)
  
  plot(dat$densities$sentinel_chicken_sampling/100, box=FALSE, axes=FALSE,
       col=extreme_purps,
       legend.args = list(text='Sampled sentinel chickens\nper square km', cex=1.1, 
                          line=0.9, side=2),
       axis.args = list(at=c(0,5e-4,1e-3,1.5e-3,2e-3), 
                        labels=c("0","5e-4","1e-3","1.5e-3","2e-3"),
                        cex.axis=1.2),
       legend.width=1.4, cex.main=1.6,
       main="Sentinel chicken surveillance density")
  plot(st_geometry(aus_simple), lwd=1.2, border='grey20', add=TRUE)
  subfigure_label(par("usr"), 0.1, 0.9, "(e)", 1.6)
  
  dev.off()}


# just making some adjustments (to label sizes) for thesis
# yep this is way better
sub_cex_main <- 2
{png("figures/pathogen_inputs_points_only_thesis.png",
     width = 900*reso/72, 
     height = 1200*reso/72, 
     res = reso)
  #type = 'cairo')
  
  par(mfrow = c(3, 2),
      mar = c(0.1, 0.1, 2.1, 0.1), bty="n",
      oma=c(0,0,0,2))
  
  # humans
  plot(st_geometry(aus_simple), cex.main=sub_cex_main,
       main="Probable exposure sites of confirmed JEV in humans", 
       border='grey20')
  # points(dat$human_points, cex=surv_cex, col="red", pch=surv_pch)
  points(dat$df[dat$df$host_type == "human" & dat$df$cases > 0, c("x","y")], 
         col=pos_pal[dat$df[dat$df$host_type == "human" & dat$df$cases > 0, 
                            weight_column]*100], 
         cex=surv_cex,
         pch=surv_pch)
  subfigure_label(par("usr"), 0.1, 0.9, "(a)", sub_cex_main)
  
  # mosquitoes
  plot(st_geometry(aus_simple), lwd=1.2, border='grey20', cex.main=sub_cex_main,
       main="JEV positive/negative samples in mosquitoes")
  points(dat$df[dat$df$host_type == "mosquito" & dat$df$cases == 0, c("x","y")], 
         col=neg_pal[dat$df[dat$df$host_type == "mosquito" & dat$df$cases == 0, 
                            weight_column]*100], 
         cex=surv_cex,
         pch=surv_pch)
  points(dat$df[dat$df$host_type == "mosquito" & dat$df$cases > 0, c("x","y")], 
         col=pos_pal[dat$df[dat$df$host_type == "mosquito" & dat$df$cases > 0, 
                            weight_column]*100], 
         cex=surv_cex,
         pch=surv_pch)  
  subfigure_label(par("usr"), 0.1, 0.9, "(b)", sub_cex_main)
  
  # domestic pigs
  plot(st_geometry(aus_simple), lwd=1.2, border='grey20', cex.main=sub_cex_main,
       main="JEV positive/negative samples in domestic pigs")
  points(dat$df[dat$df$host_type == "domestic_pig" & dat$df$cases == 0, c("x","y")], 
         col=neg_pal[dat$df[dat$df$host_type == "domestic_pig" & dat$df$cases == 0, 
                            weight_column]*100], 
         cex=surv_cex,
         pch=surv_pch)
  points(dat$df[dat$df$host_type == "domestic_pig" & dat$df$cases > 0, c("x","y")], 
         col=pos_pal[dat$df[dat$df$host_type == "domestic_pig" & dat$df$cases > 0, 
                            weight_column]*100], 
         cex=surv_cex,
         pch=surv_pch)  
  subfigure_label(par("usr"), 0.1, 0.9, "(c)", sub_cex_main)
  
  # feral pigs
  plot(st_geometry(aus_simple), lwd=1.2, border='grey20', cex.main=sub_cex_main,
       main="JEV positive/negative samples in feral pigs")
  points(dat$df[dat$df$host_type == "feral_pig" & dat$df$cases == 0, c("x","y")], 
         col=neg_pal[dat$df[dat$df$host_type == "feral_pig" & dat$df$cases == 0, 
                            weight_column]*100], 
         cex=surv_cex,
         pch=surv_pch)
  points(dat$df[dat$df$host_type == "feral_pig" & dat$df$cases > 0, c("x","y")], 
         col=pos_pal[dat$df[dat$df$host_type == "feral_pig" & dat$df$cases > 0, 
                            weight_column]*100], 
         cex=surv_cex,
         pch=surv_pch)  
  subfigure_label(par("usr"), 0.1, 0.9, "(d)", sub_cex_main)
  
  # chickens
  plot(st_geometry(aus_simple), lwd=1.2, border='grey20', cex.main=sub_cex_main,
       main="JEV positive/negative samples in sentinel chickens")
  points(dat$df[dat$df$host_type == "chicken" & dat$df$cases == 0, c("x","y")], 
         col=neg_pal[dat$df[dat$df$host_type == "chicken" & dat$df$cases == 0, 
                            weight_column]*100], 
         cex=surv_cex,
         pch=surv_pch)
  points(dat$df[dat$df$host_type == "chicken" & dat$df$cases > 0, c("x","y")], 
         col=pos_pal[dat$df[dat$df$host_type == "chicken" & dat$df$cases > 0, 
                            weight_column]*100], 
         cex=surv_cex,
         pch=surv_pch)  
  subfigure_label(par("usr"), 0.1, 0.9, "(e)", sub_cex_main)
  
  par(mar=c(15.1,9.1,6.1,1.1), bty="o")
  plot(0, xlim = c(as.Date("2021-01-01"), as.Date("2023-12-04")), ylim=c(0.5,5.5), xlab="",
       main="JEV surveillance in time by host type", cex.main=sub_cex_main, ylab="",
       yaxt="n", xaxt="n", type="n", cex.lab=1.2)
  abline(v=quarter_breaks, col="grey90")
  abline(v=c(as.Date("2021-01-01"),as.Date("2022-01-01"),as.Date("2023-01-01")), col="grey60")
  # abline(v=as.Date("2022-03-04"), col="grey50", lty=2, lwd=2) # CDIoNS
  # plot all negatives
  points(as.Date(all_dates$date[all_dates$virus_detection == "JEV negative"]),
         jitter(all_dates$ypos[all_dates$virus_detection == "JEV negative"],
                factor=jit_fac), pch=surv_pch, cex=date_cex,
         col=neg_pal[round(all_dates$likelihood_weight[all_dates$virus_detection == "JEV negative"]*100)])
  # plot all positives
  points(as.Date(all_dates$date[all_dates$virus_detection == "JEV positive"]),
         jitter(all_dates$ypos[all_dates$virus_detection == "JEV positive"],
                factor=jit_fac), pch=surv_pch, cex=date_cex,
         col=pos_pal[round(all_dates$likelihood_weight[all_dates$virus_detection == "JEV positive"]*100)])
  # not sure why the humans disappear?
  # points(all_dates$date[all_dates$ypos == 5], cex=0.8,
  #        jitter(rep(5, length(all_dates$date[all_dates$ypos == 5])),
  #               factor=jit_fac*2), pch=surv_pch,
  #        col=pos_pal[round(all_dates$likelihood_weight[all_dates$virus_detection == "JEV positive"]*100)])
  
  axis(2, at=5:1, c("Human","Mosquito","Domestic pig","Feral pig","Chicken"), cex.axis=1.7, las=2)
  axis(1, at=quarter_breaks, format(quarter_breaks, format="%b %y"), las=2, cex.axis=1.4)
  mtext("Date", 1, line=5, cex=1.2)
  mtext("(f)", 3, 2.5, adj=0, cex=1.4)
  
  
  par(new=TRUE, mar=c(0,0,0,0))
  plot(0, axes=FALSE, type="n", xlim=c(0,1), ylim=c(0,1), xlab="", ylab="")
  gradient.rect(0.38, 0.1,
                0.58, 0.12,
                col=rev(neg_pal),
                border='grey20',
                gradient='x')
  gradient.rect(0.58, 0.1,
                0.78, 0.12,
                col=(pos_pal),
                border='grey20',
                gradient='x')
  axis(1, at=seq(0.38,0.78,length.out=5), labels=seq(-1,1,0.5), line=-5.5, cex=1.5)
  text(c(0.48, 0.68),c(0.02, 0.02), c("JEV Negative", "JEV Positive"), cex=1.6)
  text(0.58,0.15,"Likelihood weight", cex=1.8)
  
  dev.off()}



sub_cex_main <- 2
leg_cex_main <- 1.3
leg_cex_axis <- 1.4
{png("figures/pathogen_inputs_surface_only_thesis.png",
     width = 1000*reso/72, 
     height = 1100*reso/72, 
     res = reso)
  
  par(mfrow = c(3, 2),
      mar = c(0.1, 0.1, 3.1, 1), bty="n",
      oma=c(0,0,0,2))
  
  # humans
  plot(sqrt(hpop_blur/100), col=extreme_purps, legend.width=1.4,
       legend.args = list(text='Human population\nper square km', 
                          line=0.9, side=2, cex=leg_cex_main), 
       axis.args = list(at=c(sqrt(10),10,sqrt(250),sqrt(500),30), 
                        labels=c(10,100,250,500,900), cex.axis=leg_cex_axis),
       box = FALSE, axes = FALSE, cex.main=sub_cex_main,
       main="Distance-weighted human population density")
  plot(st_geometry(aus_simple), add=TRUE, lwd=1.2, border='grey20')
  subfigure_label(par("usr"), 0.1, 0.9, "(a)", sub_cex_main)
  
  # mosquitoes
  plot(#nat_mozz_density/100, # decemberist
    dat$densities$mosquito_sampling/100,
    box = FALSE, axes=FALSE, col=extreme_purps, 
    legend.width=1.4,
    legend.args = list(text='Mosquito catches screened\nper square km', 
                       line=0.9, cex=leg_cex_main, side=2), 
    axis.args = list(cex.axis=leg_cex_axis),
    cex.main=sub_cex_main,
    main="Mosquito screening effort")
  plot(st_geometry(aus_simple), lwd=1.2, border='grey20', add=TRUE)
  subfigure_label(par("usr"), 0.1, 0.9, "(b)", sub_cex_main)
  
  # domestic pigs
  plot(dat$densities$domestic_pig_sampling/100, box = FALSE, axes=FALSE, 
       col=extreme_purps,
       legend.args = list(text='Sampled domestic pigs\nper square km', 
                          cex=leg_cex_main, line=0.9, side=2), 
       axis.args = list(at=seq(0,2e-3, length.out=5), 
                        #labels=c("0","2e-4","4e-4","6e-4","8e-4"),
                        cex.axis=leg_cex_axis),
       legend.width=1.4, cex.main=sub_cex_main,
       main="Domestic pig surveillance density")
  plot(st_geometry(aus_simple), lwd=1.2, border='grey20', add=TRUE)
  subfigure_label(par("usr"), 0.1, 0.9, "(c)", sub_cex_main)
  
  # feral pigs
  plot(dat$densities$feral_pig_sampling/100, box = FALSE, axes=FALSE, 
       col=extreme_purps,
       legend.args = list(text='Sampled feral pigs\nper square km', cex=leg_cex_main, 
                          line=0.9, side=2),
       axis.args = list(at=c(0,2e-3,4e-3,6e-3,8e-3),
                        labels=c("0","2e-3","4e-3","6e-3","8e-3"),
                        cex.axis=leg_cex_axis),
       legend.width=1.4, cex.main=sub_cex_main,
       main="Feral pig screening effort")
  plot(st_geometry(aus_simple), lwd=1.2, border='grey20', add=TRUE) 
  subfigure_label(par("usr"), 0.1, 0.9, "(d)", sub_cex_main)
  
  plot(dat$densities$sentinel_chicken_sampling/100, box=FALSE, axes=FALSE,
       col=extreme_purps,
       legend.args = list(text='Sampled sentinel chickens\nper square km', cex=leg_cex_main, 
                          line=0.9, side=2),
       axis.args = list(#at=c(0,5e-4,1e-3,1.5e-3,2e-3), 
                        at = seq(0,3e-3, length.out=4),
                        #labels=c("0","5e-4","1e-3","1.5e-3","2e-3"),
                        cex.axis=leg_cex_axis),
       legend.width=1.4, cex.main=sub_cex_main,
       main="Sentinel chicken surveillance density")
  plot(st_geometry(aus_simple), lwd=1.2, border='grey20', add=TRUE)
  subfigure_label(par("usr"), 0.1, 0.9, "(e)", sub_cex_main)
  
  dev.off()}




# just making some adjustments (to label sizes) for thesis
# yep this is way better
month_breaks <- seq(ymd("2019-11-01"),
                    max(as.Date(all_dates$date), na.rm=TRUE),
                    by="month")
quarter_breaks <- month_breaks[seq(3, length(month_breaks), 3)]
sub_cex_main <- 2
{png("figures/pathogen_inputs_pos_only.png",
     width = 900*reso/72, 
     height = 1200*reso/72, 
     res = reso)
  #type = 'cairo')
  
  par(mfrow = c(3, 2),
      mar = c(0.1, 0.1, 2.1, 0.1), bty="n",
      oma=c(0,0,0,2))
  
  # humans
  plot(st_geometry(aus_simple), cex.main=sub_cex_main,
       main="Probable exposure sites of confirmed JEV in humans", 
       border='grey20')
  # points(dat$human_points, cex=surv_cex, col="red", pch=surv_pch)
  points(dat$df[dat$df$host_type == "human" & dat$df$cases > 0, c("x","y")], 
         col=pos_pal[dat$df[dat$df$host_type == "human" & dat$df$cases > 0, 
                            weight_column]*100], 
         cex=surv_cex,
         pch=surv_pch)
  subfigure_label(par("usr"), 0.1, 0.9, "(a)", sub_cex_main)
  
  # mosquitoes
  plot(st_geometry(aus_simple), lwd=1.2, border='grey20', cex.main=sub_cex_main,
       main="JEV positive samples in mosquitoes")
  points(dat$df[dat$df$host_type == "mosquito" & dat$df$cases > 0, c("x","y")], 
         col=pos_pal[dat$df[dat$df$host_type == "mosquito" & dat$df$cases > 0, 
                            weight_column]*100], 
         pch=surv_pch, 
         cex=surv_cex)  
  subfigure_label(par("usr"), 0.1, 0.9, "(b)", sub_cex_main)
  
  # domestic pigs
  plot(st_geometry(aus_simple), lwd=1.2, border='grey20', cex.main=sub_cex_main,
       main="JEV positive samples in domestic pigs")
  points(dat$df[dat$df$host_type == "domestic_pig" & dat$df$cases > 0, c("x","y")], 
         col=pos_pal[dat$df[dat$df$host_type == "domestic_pig" & dat$df$cases > 0, 
                            weight_column]*100], 
         pch=surv_pch, cex=surv_cex)  
  subfigure_label(par("usr"), 0.1, 0.9, "(c)", sub_cex_main)
  
  # feral pigs
  plot(st_geometry(aus_simple), lwd=1.2, border='grey20', cex.main=sub_cex_main,
       main="JEV positive samples in feral pigs")
  points(dat$df[dat$df$host_type == "feral_pig" & dat$df$cases > 0, c("x","y")], 
         col=pos_pal[dat$df[dat$df$host_type == "feral_pig" & dat$df$cases > 0, 
                            weight_column]*100], 
         pch=surv_pch, cex=surv_cex)  
  subfigure_label(par("usr"), 0.1, 0.9, "(d)", sub_cex_main)
  
  # chickens
  plot(st_geometry(aus_simple), lwd=1.2, border='grey20', cex.main=sub_cex_main,
       main="JEV positive samples in sentinel chickens")
  points(dat$df[dat$df$host_type == "chicken" & dat$df$cases > 0, c("x","y")], 
         col=pos_pal[dat$df[dat$df$host_type == "chicken" & dat$df$cases > 0, 
                            weight_column]*100], 
         pch=surv_pch, cex=surv_cex)  
  subfigure_label(par("usr"), 0.1, 0.9, "(e)", sub_cex_main)
  
  par(mar=c(15.1,9.1,6.1,1.1), bty="o")
  plot(0, xlim = c(as.Date("2021-01-01"), as.Date("2023-12-04")), ylim=c(0.7,5.3), xlab="",
       main="JEV detections in time by host type", cex.main=sub_cex_main, ylab="",
       yaxt="n", xaxt="n", type="n", cex.lab=1.2)
  abline(h = seq(1.5,4.5), col="grey90")
  abline(v=quarter_breaks, col="grey90")
  abline(v=c(as.Date("2020-01-01"),as.Date("2021-01-01"),as.Date("2022-01-01"),as.Date("2023-01-01")), col="grey60")
  
  # abline(v=as.Date("2022-03-04"), col="grey50", lty=2, lwd=2) # CDIoNS
  # plot all negatives
  # plot all positives
  points(as.Date(all_dates$date[all_dates$virus_detection == "JEV positive"]),
         jitter(all_dates$ypos[all_dates$virus_detection == "JEV positive"],
                factor=jit_fac), pch=surv_pch, cex=date_cex,
         col=pos_pal[round(all_dates$likelihood_weight[all_dates$virus_detection == "JEV positive"]*100)])
  # not sure why the humans disappear?
  # points(all_dates$date[all_dates$ypos == 5], cex=0.8,
  #        jitter(rep(5, length(all_dates$date[all_dates$ypos == 5])),
  #               factor=jit_fac*2), pch=surv_pch,
  #        col=pos_pal[round(all_dates$likelihood_weight[all_dates$virus_detection == "JEV positive"]*100)])
  
  axis(2, at=5:1, c("Human","Mosquito","Domestic pig","Feral pig","Chicken"), cex.axis=1.7, las=2)
  axis(1, at=quarter_breaks, format(quarter_breaks, format="%b %y"), las=2, cex.axis=1.4)
  mtext("Date", 1, line=5, cex=1.2)
  mtext("(f)", 3, 2.5, adj=0, cex=1.4)
  
  
  par(new=TRUE, mar=c(0,0,0,0))
  plot(0, axes=FALSE, type="n", xlim=c(0,1), ylim=c(0,1), xlab="", ylab="")
  # gradient.rect(0.38, 0.1,
  #               0.58, 0.12,
  #               col=rev(neg_pal),
  #               border='grey20',
  #               gradient='x')
  # gradient.rect(0.58, 0.1,
  #               0.78, 0.12,
  #               col=(pos_pal),
  #               border='grey20',
  #               gradient='x')
  # axis(1, at=seq(0.38,0.78,length.out=5), labels=seq(-1,1,0.5), line=-5.5, cex=1.5)
  # text(c(0.48, 0.68),c(0.02, 0.02), c("JEV Negative", "JEV Positive"), cex=1.6)
  # text(0.58,0.15,"Likelihood weight", cex=1.8)
  
  gradient.rect(0.38, 0.1,
                0.78, 0.12,
                col=(pos_pal),
                border='grey20',
                gradient='x')
  axis(1, at=seq(0.38,0.78,length.out=5), labels=seq(0,1,0.25), line=-5.5, cex=1.5)
  #text(0.58,0.02, "JEV Positive", cex=1.6)
  text(0.58,0.15,"Likelihood weight - JEV positive", cex=1.8)
  
  dev.off()}

