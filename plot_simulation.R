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
dat <- read.csv("result.csv")[-c(1, 23, 24)] %>% subset(., time >= 1)

dat <- gather(dat, key=component, value=concentration, -hv, -mu, -time)
dat <- separate(dat, component, c('variable', 'component'))

# ------------------------------------------------------------------------------
# same approach for loading the dynamic model data
dat2 <- read.csv("result_dynamic.csv")[-1] %>% subset(., time >= 1)

dat2 <- gather(dat2, key=component, value=concentration, -hv, -time, -sigma)
dat2 <- separate(dat2, component, c('variable', 'component'))


# +++++++++++++ PLOT FRACTIONS AND CONCENTRATIONS ++++++++++++++++++++++++++++++
#
# plot growth rates
plots=list(
  doubleYScale(
    xyplot(hv ~ time, dat[1:100, ],
      type="l", lwd=2,
      ylim=c(0, 120),
      xlab="time", ylab="% max light",
      panel=function(x, y, ...) {
        panel.grid(h=-1, v=-1, col=grey(0.9))
        panel.xyplot(x, y, ...)
      }
    ), 
    xyplot(mu ~ time, dat[1:100, ], type="l", lwd=2, ylim=c(0, 0.12)), 
    use.style=TRUE, add.axis=TRUE, add.ylab2=TRUE
  )
)

# create plots of different 'process' variables and collect in list
plots[2:4] <- lapply(c("a", "v", "c"), function(i) {
  direct.label(
    method="visualcenter",
    xyplot(concentration ~ time, 
      subset(dat, variable == i), 
      groups=component, #ylim=c(0, 0.32),
      xlab="time", ylab=paste0("[",i,"]"),
      type="l", lwd=2, 
      scales=list(alternating=FALSE),
      panel=function(x, y, ...) {
        panel.grid(h=-1, v=-1, col=grey(0.9))
        panel.xyplot(x, y, ...)
      }
    )
  )
})

# additional plots for dynamic simulation: overlay of

#svg("plot_simulation_steadystate.svg", width=6, height=5.8)
print(plots[[1]], split=c(1,1,2,2), more=TRUE)
print(plots[[2]], split=c(1,2,2,2), more=TRUE)
print(plots[[3]], split=c(2,1,2,2), more=TRUE)
print(plots[[4]], split=c(2,2,2,2))
#grid.text(c("A", "B", "C"), x=c(0.05, 0.5, 0.05), y=c(0.95, 0.95, 0.6), gp=gpar(fontsize=14))
#dev.off()

