# loading libraries
library(lattice)
library(latticeExtra)
library(plyr)
library(tidyr)
library(directlabels)
setwd("~/Documents/SciLifeLab/Resources/Models/GEKKO/cyano/")
source("~/Documents/SciLifeLab/Resources/R_scripts/custom.theme.R")
source("~/Documents/SciLifeLab/Resources/R_scripts/custom.panel.functions.R")


# +++++++++++++ LOAD MODEL DATA ++++++++++++++++++++++++++++++++++++++++++++++++

# load modeling data according to following names
dat <- read.csv("result.csv") %>% subset(., mu > 0.01 & mu <1)


# +++++++++++++ PLOT FRACTIONS AND CONCENTRATIONS ++++++++++++++++++++++++++++++
#
# plot growth rates
plots <- list(xyplot(light+mu*100 ~ time,
  dat[1:9, ],
  par.settings=custom.lattice, type="l", lwd=2,
  xlab="time", ylab="% max light / µ",
  panel=function(x, y, ...) {
    panel.grid(h=-1, v=-1, col=grey(0.9))
    panel.xyplot(x, y, ...)
  }
))

# create plots of different 'process' variables and collect in list
plots[2:4] <- lapply(c("a", "c", "v"), function(i) {
  direct.label(
    method="visualcenter",
    xyplot(concentration ~ time, 
      subset(dat, variable == i), 
      groups=component, 
      xlab="time", ylab=paste0("[",i,"]"),
      type="l", lwd=2, 
      scales=list(alternating=FALSE),
      par.settings=custom.lattice, 
      panel=function(x, y, ...) {
        panel.grid(h=-1, v=-1, col=grey(0.9))
        panel.xyplot(x, y, ...)
      }
    )
  )
})
  

svg("plot_simulation.svg", width=6, height=5.8)
print(plots[[1]], split=c(1,1,2,2), more=TRUE)
print(plots[[2]], split=c(1,2,2,2), more=TRUE)
print(plots[[3]], split=c(2,1,2,2), more=TRUE)
print(plots[[4]], split=c(2,2,2,2))
#grid.text(c("A", "B", "C"), x=c(0.05, 0.5, 0.05), y=c(0.95, 0.95, 0.6), gp=gpar(fontsize=14))
dev.off()

