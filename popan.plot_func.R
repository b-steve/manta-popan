popan.plot <- function(fit.ma, year.start = 2009, year.end = 2019) {

  ## Extracting parameter values from the model fit
  # ENs female + male
  EN.mf <- summary(fit.ma, par.fun = function(x, ...) x$ENs[, 1] + x$ENs[, 2])
  # ENs female
  EN.f <- summary(fit.ma, pars = "ENs", groups = 1)
  # ENs male
  EN.m <- summary(fit.ma, pars = "ENs", groups = 2)
  # ps female
  p.f <- summary(fit.ma, pars = "ps", groups = 1)
  # ps male
  p.m <- summary(fit.ma, pars = "ps", groups = 2)
  # phis female
  phi.f <- summary(fit.ma, pars = "phis", groups = 1)
  # phis male
  phi.m <- summary(fit.ma, pars = "phis", groups = 2)
  # rhos female
  rho.f <- summary(fit.ma, pars = "rhos", groups = 1)
  # rhos male
  rho.m <- summary(fit.ma, pars = "rhos", groups = 2)
  
  n.start <- (year.start-2009)+1
  n.end <- (year.end-2009)+1
  max.seq <- (n.end-n.start)+1
  
  ## PLOTTING ####
  par(mfrow=c(3,2)) # dividing plots into 6 panels (row, col)
  par(mar = c(3,3,2,1)) # adjust graph margin: below, left, top, right
  par(mgp = c(1.8,0.4,0))
  
  ## PANEL A = ENs for total population (female+male) and for females only ####
  # total population (female+male)
  plot(EN.mf[n.start:n.end, 1],
       col = "orangered2",
       ylim = c(0, max(EN.mf[, 4])),
       type = "b", 
       pch = 15,
       xlab = "Year",
       ylab = "EN",
       main = "Expected population size",
       cex.main = 1.7, # size of axis title
       xaxt="n",
       cex.axis = 1.2, # sze of axis label
       cex = 1.6, # size of symbols
       tck=-.015, # size of ticks in axis
       cex.lab = 1.6) # size of axis title
  axis(1, at = seq(1,max.seq, by = 1), cex.axis=1.2,
       labels = c(year.start:year.end), tck=-.015)
  ## Adding dotted line for lower and upper CI limits
  lines(EN.mf[n.start:n.end, 3], lty = "longdash", col = "orangered2")
  lines(EN.mf[n.start:n.end, 4], lty = "longdash", col = "orangered2")
  
  par(new = T)
  
  # plotting female
  plot(EN.f[n.start:n.end, 1], 
       col = "goldenrod1",
       ylim = c(0, max(EN.mf[, 4])),
       type = "b", 
       pch = 1,
       xlab = "",
       ylab = "",
       main = "",
       cex.main = 1.7,
       xaxt="n",
       yaxt="n",
       cex.axis = 1.2,
       cex = 1.6,
       tck=-.015,
       cex.lab = 1.6)
  ## Adding dotted line for lower and upper CI limits
  lines(EN.f[n.start:n.end, 3], lty = "dashed", col = "goldenrod1")
  lines(EN.f[n.start:n.end, 4], lty = "dashed", col = "goldenrod1")
  
  # create legend for Panel A
  legend("topleft", legend= c("Females + Males", "Females"), col = c('orangered2', 'goldenrod1'),
         cex=1.6, bg="transparent", pch = c(15,1),
         box.lty=0)
  legend("bottomright", text.font=2, legend= "A",
         cex=2.5, bg="transparent",
         box.lty=0)
  box()
  
  ## PANEL B = ENs for total population (female+male) and for females only ####
  # total population (female+male)
  plot(EN.mf[n.start:n.end, 1],
       col = "orangered2",
       ylim = c(0, max(EN.mf[, 4])),
       type = "b", 
       pch = 15,
       xlab = "Year",
       ylab = "EN",
       main = "Expected population size",
       cex.main = 1.7, # size of axis title
       xaxt="n",
       cex.axis = 1.2, # sze of axis label
       cex = 1.6,
       tck=-.015,
       cex.lab = 1.6)
  axis(1, at = seq(1,max.seq, by = 1), cex.axis=1.2,
       labels = c(year.start:year.end), tck=-.015)
  ## Adding dotted line for lower and upper CI limits
  lines(EN.mf[n.start:n.end, 3], lty = "longdash", col = "orangered2")
  lines(EN.mf[n.start:n.end, 4], lty = "longdash", col = "orangered2")
  
  par(new = T)
  
  # plotting male
  plot(EN.m[n.start:n.end, 1], 
       col = "turquoise3",
       ylim = c(0, max(EN.mf[, 4])),
       type = "b", 
       pch = 2,
       xlab = "",
       ylab = "",
       main = "",
       cex.main = 1.7,
       xaxt="n",
       yaxt="n",
       cex.axis = 1.2,
       cex = 1.6,
       tck=-.015,
       cex.lab = 1.6)
  ## Adding dotted line for lower and upper CI limits
  lines(EN.m[n.start:n.end, 3], lty = "dashed", col = "turquoise3")
  lines(EN.m[n.start:n.end, 4], lty = "dashed", col = "turquoise3")
  
  # create legend for Panel B
  legend("topleft", legend= c("Females + Males", "Males"), col = c('orangered2', 'turquoise3'),
         cex=1.6, bg="transparent", pch = c(15,2),
         box.lty=0)
  legend("bottomright", text.font=2, legend= "B",
         cex=2.5, bg="transparent",
         box.lty=0)
  box()
  
  ## PANEL C = sighting probability (p) for female ####
  plot(p.f[n.start:n.end, 1],
       col = "goldenrod1",
       ylim = c(0, max(p.f[, 4])),
       type = "b", 
       pch = 1,
       xlab = "Year",
       ylab = "p",
       main = "Sighting probabilities",
       cex.main = 1.7, # size of axis title
       xaxt="n",
       cex.axis = 1.2, # sze of axis label
       cex = 1.6,
       tck=-.015,
       cex.lab = 1.6)
  axis(1, at = seq(1,max.seq, by = 1), cex.axis=1.2,
       labels = c(year.start:year.end), tck=-.015)
  ## Adding dotted line for lower and upper CI limits
  lines(p.f[n.start:n.end, 3], lty = "longdash", col = "goldenrod1")
  lines(p.f[n.start:n.end, 4], lty = "longdash", col = "goldenrod1")
  
  # create legend for Panel C
  legend("topleft", legend= "Females", col = "goldenrod1",
         pch = 1, cex=1.6, bg="transparent",
         box.lty=0)
  legend("bottomright", text.font=2, legend= "C",
         cex=2.5, bg="transparent",
         box.lty=0)
  box()
  
  ## PANEL D = sighting probability (p) for male ####
  plot(p.m[n.start:n.end, 1],
       col = "turquoise3",
       ylim = c(0, max(p.f[, 4])),
       type = "b", 
       pch = 1,
       xlab = "Year",
       ylab = "p",
       main = "Sighting probabilities",
       cex.main = 1.7, # size of axis title
       xaxt="n",
       cex.axis = 1.2, # sze of axis label
       cex = 1.6,
       tck=-.015,
       cex.lab = 1.6)
  axis(1, at = seq(1,max.seq, by = 1), cex.axis=1.2,
       labels = c(year.start:year.end), tck=-.015)
  ## Adding dotted line for lower and upper CI limits
  lines(p.m[n.start:n.end, 3], lty = "longdash", col = "turquoise3")
  lines(p.m[n.start:n.end, 4], lty = "longdash", col = "turquoise3")
  
  # create legend for Panel C
  legend("topleft", legend= "Males", col = "turquoise3",
         pch = 1, cex=1.6, bg="transparent",
         box.lty=0)
  legend("bottomright", text.font=2, legend= "D",
         cex=2.5, bg="transparent",
         box.lty=0)
  box()
  
  ## PANEL E: survival probabilities (phis)
  # phis for females
  phi.f
  last.row <- rep(NA, ncol(phi.f))
  phi.f.full <- rbind(phi.f, last.row)
  ## Creating a line for point estimates. The y-axis goes from 0 to the
  ## highest upper CI limit.
  plot(phi.f.full[n.start:n.end, 1], 
       col = "goldenrod1",
       ylim = c(0, 1),
       type = "b", 
       pch = 1,
       xlab = "Year",
       ylab = "phi",
       main = "Survival probabilities",
       cex.main = 1.7,
       xaxt="n",
       cex.axis = 1.2,
       cex = 1.6,
       cex.lab = 1.6)
  axis(1, at = seq(1,max.seq, by = 1), cex.axis=1.2,
       labels = c(year.start:year.end), tck=-.015)
  
  ## Adding dotted line for lower and upper CI limits
  lines(phi.f.full[n.start:n.end, 3], lty = "dashed", col = "goldenrod1")
  lines(phi.f.full[n.start:n.end, 4], lty = "dashed", col = "goldenrod1")
  
  par(new = T)
  
  # phis for males
  phi.m.full <- rbind(phi.m, last.row)
  plot(phi.m.full[n.start:n.end, 1],
       col = "turquoise3",
       ylim = c(0, 1),
       type = "b", 
       pch = 2,
       xlab = "",
       ylab = "",
       xaxt="n",
       cex = 1.6,
       yaxt="n")
  lines(phi.m.full[n.start:n.end, 3], lty = "dotted", col = "turquoise3")
  lines(phi.m.full[n.start:n.end, 4], lty = "dotted", col = "turquoise3")
  
  # create legend for Panel E
  legend("bottomleft", legend= c("Females", "Males"), col = c("goldenrod1", "turquoise3"),
         pch = c(1,2), cex=1.6, bg="transparent",
         box.lty=0)
  legend("bottomright", text.font=2, legend= "E",
         cex=2.5, bg="transparent",
         box.lty=0)
  box()
  
  
  ## PANEL F: per capita recruitment rate (rho)
  # phis for females
  rho.f
  last.row <- rep(NA, ncol(rho.f))
  rho.f.full <- rbind(rho.f, last.row)
  ## Creating a line for point estimates. The y-axis goes from 0 to the
  ## highest upper CI limit.
  plot(rho.f.full[n.start:n.end, 1], 
       col = "goldenrod1",
       ylim = c(0, 1),
       type = "b", 
       pch = 1,
       xlab = "Year",
       ylab = "rho",
       main = "Per capita recruitment rates",
       cex.main = 1.7,
       xaxt="n",
       cex.axis = 1.2,
       cex = 1.6,
       cex.lab = 1.6)
  axis(1, at = seq(1,max.seq, by = 1), cex.axis=1.2,
       labels = c(year.start:year.end), tck=-.015)
  ## Adding dotted line for lower and upper CI limits
  lines(rho.f.full[n.start:n.end, 3], lty = "dashed", col = "goldenrod1")
  lines(rho.f.full[n.start:n.end, 4], lty = "dashed", col = "goldenrod1")
  
  par(new = T)
  
  # rhos for males
  rho.m.full <- rbind(rho.m, last.row)
  plot(rho.m.full[n.start:n.end, 1],
       col = "turquoise3",
       ylim = c(0, 1),
       type = "b", 
       pch = 2,
       xlab = "",
       ylab = "",
       xaxt="n",
       cex = 1.6,
       yaxt="n")
  lines(rho.m.full[n.start:n.end, 3], lty = "dotted", col = "turquoise3")
  lines(rho.m.full[n.start:n.end, 4], lty = "dotted", col = "turquoise3")
  
  # create legend for Panel F
  legend("topleft", legend= c("Females", "Males"), col = c("goldenrod1", "turquoise3"),
         pch = c(1,2), cex=1.6, bg="transparent",
         box.lty=0)
  legend("bottomright", text.font=2, legend= "F",
         cex=2.5, bg="transparent",
         box.lty=0)
  box()
  
}
