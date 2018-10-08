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

# list of raw model result files to load
files <- c(
  "result_steady_state.csv",
  "result_dynamic.csv",
  "result_dynamic_RIB_10.csv",
  "result_dynamic_RIB_20.csv",
  "result_dynamic_RIB_40.csv"
)

# Generalized function to load and combine data from multiple result tables
dat <- lapply(files, function(filename) {
  d <- read.csv(filename)[-1, ]
  d <- d[grep("X|slk", colnames(d), invert=TRUE)]
  d <- gather(d, key=component, value=concentration, -time)
  d <- separate(d, component, c('variable', 'component'))
  d$simulation <- gsub("result_|\\.csv", "", filename) %>% factor(., unique(.))
  d
}) %>% ldply


# +++++++++++++ PLOTTING +++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# change default colors to gradual ones
custom.lattice$superpose.line$col <- c(grey(0.5), heat.colors(7)[2:5])
custom.lattice$superpose.polygon$col <- c(grey(0.5), heat.colors(7)[2:5])


# generalized function to create plots of different 'process' variables
plot.var <- function(dat, var, groups, draw.panel.key=FALSE,
  comp=c("mai", "rib", "lhc", "pset", "cbm", "lpb"), 
  ylim=NULL, ylog=FALSE) {
  
  # first part of the plot is light intensity
  plot1 <- xyplot(c(0, concentration, 0) ~ c(min(time), time, max(time)), 
    subset(dat, variable=="hv" & simulation=="steady_state"),
    par.settings=custom.lattice,
    col="#A6A6A641", border=0,
    xlab="time", ylab="% light",
    ylim=c(0, 100), panel=panel.polygon)
  
  # second part of the plot are all concentrations
  plot2 <- xyplot(concentration ~ time, 
    data=subset(dat, variable==var & component %in% comp),
    groups=get(groups),
    par.settings=custom.lattice, 
    ylab=var,
    scales=list(y=list(limits=ylim, log=ylog)),
    type="l", lwd=2, #pch=19, cex=0.3,
    panel=function(x, y, ...) {
      panel.grid(h=-1, v=-1, col=grey(0.9))
      panel.xyplot(x, y, ...)
      if (is.na(comp)) panel.text(0.9, 0.9, var) else
        panel.text(0.9*max(x), 0.9*max(y), labels=toupper(comp))
      if (draw.panel.key) {
        panel.key(as.character(dat[[groups]]) %>% unique,
          points=FALSE, lines=TRUE, ...)
      }
    }
  )
  
  # combine in double Y scale plot
  doubleYScale(plot1, plot2, use.style=FALSE, add.axis=TRUE, add.ylab2=TRUE)
}

#svg("plot_simulation.svg", width=6, height=8)
print(plot.var(dat, var="mu", comp=NA, groups="simulation", ylim=c(0,0.12)), split=c(1,1,2,3), more=TRUE)
print(plot.var(dat, var="bm", comp=NA, groups="simulation", ylim=c(0,50)), split=c(2,1,2,3), more=TRUE)
print(plot.var(dat, var="c", comp="rib", groups="simulation", ylim=c(0,0.5)), split=c(1,2,2,3), more=TRUE)
print(plot.var(dat, var="c", comp="lhc", groups="simulation", ylim=c(0,0.5)), split=c(2,2,2,3), more=TRUE)
print(plot.var(dat, var="c", comp="cbm", groups="simulation", ylim=c(0,0.5)), split=c(1,3,2,3), more=TRUE)
print(plot.var(dat, var="c", comp="pset", groups="simulation", ylim=c(0,0.5), draw.panel.key=TRUE), split=c(2,3,2,3))
#grid.text(c("steady state", "dynamic"), x=c(0.25, 0.75), y=c(0.98, 0.98), gp=gpar(fontsize=12))
#dev.off()
