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

# load steady state model data
dat <- read.csv("result_steady_state_DN.csv")[-c(1, 23, 24)] %>% subset(., time >= 1)

dat <- gather(dat, key=component, value=concentration, -hv, -mu, -time)
dat <- separate(dat, component, c('variable', 'component'))

# ------------------------------------------------------------------------------
# same approach for loading the dynamic model data
dat2 <- read.csv("result_dynamic_DN.csv")[-1] %>% subset(., time >= 1)

dat2 <- gather(dat2, key=component, value=concentration, -hv, -mu, -time, -sigma)
dat2 <- separate(dat2, component, c('variable', 'component'))


# +++++++++++++ PLOT FRACTIONS AND CONCENTRATIONS ++++++++++++++++++++++++++++++
#
# plot growth rates
mu.plots <- lapply(c("dat", "dat2"), function(i) {
  doubleYScale(
    xyplot(hv ~ time, get(i)[unique(get(i)$time), ],
      par.settings=custom.lattice,
      type="l", lwd=2,
      ylim=c(0, 120),
      xlab="time", ylab="% light",
      panel=function(x, y, ...) {
        panel.grid(h=-1, v=-1, col=grey(0.9))
        panel.xyplot(x, y, ...)
      }
    ), 
    xyplot(mu ~ time, get(i)[unique(get(i)$time), ], type="l", lwd=2, ylim=c(0, 0.12)), 
    use.style=TRUE, add.axis=TRUE, add.ylab2=TRUE
  )
})


# generalized function to create plots of different 'process' variables
plot.var <- function(dat, var, comp=c("mai", "rib", "lhc", "pset", "cbm", "lpb"), 
  ylim=NULL, ylog=FALSE) {
  direct.label(
    method="visualcenter",
    xyplot(concentration ~ time, 
      data=subset(dat, variable==var & component %in% comp), 
      scales=list(y=list(limits=ylim, log=ylog)),
      groups=component,
      par.settings=custom.lattice,
      xlab="time", ylab=paste0("[", var, "]"),
      type="l", lwd=2, 
      panel=function(x, y, ...) {
        panel.grid(h=-1, v=-1, col=grey(0.9))
        panel.xyplot(x, y, ...)
      }
    )
  )
}


#svg("plot_simulation.svg", width=6, height=8)
print(mu.plots[[1]], split=c(1,1,2,3), more=TRUE)
print(mu.plots[[2]], split=c(2,1,2,3), more=TRUE)
print(plot.var(dat , var="c"), split=c(1,2,2,3), more=TRUE)
print(plot.var(dat2, var="c"), split=c(2,2,2,3), more=TRUE)
print(plot.var(dat , var="v", ylog=10), split=c(1,3,2,3), more=TRUE)
print(plot.var(dat2, var="v", ylog=10), split=c(2,3,2,3))
grid.text(c("steady state", "dynamic"), x=c(0.25, 0.75), y=c(0.98, 0.98), gp=gpar(fontsize=12))
#dev.off()

