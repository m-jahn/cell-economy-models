# loading libraries
library(lattice)
library(latticeExtra)
library(tidyverse)
library(directlabels)
library(Rtools)
setwd("~/Documents/SciLifeLab/Resources/Models/GEKKO/cyano/")


# +++++++++++++ LOAD MODEL DATA ++++++++++++++++++++++++++++++++++++++++++++++++

# list of raw model result files to load
files <- list.files() %>% grep("result_", ., value = TRUE) %>% .[c(1:5)]

# Generalized function to load and combine data from multiple result tables
dat <- lapply(files, function(filename) {
  d <- read_csv(filename)[-1, ]
  d <- d[grep("X|slk", colnames(d), invert = TRUE)]
  d <- gather(d, key = component, value = concentration, -time)
  d <- separate(d, component, c('variable', 'component'), extra = "drop")
  d$simulation <- gsub("result_|\\.csv", "", filename) %>% factor(., unique(.))
  d$mu <- rep_len(subset(d, variable == "mu")$concentration, nrow(d))
  d
}) %>% dplyr::bind_rows() %>% mutate(simulation = factor(simulation))


# +++++++++++++ PLOTTING +++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# change default colors to gradual ones
heat.custom <- c(grey(0.5), heat.colors(length(files)+3)[(1:length(files))+1])
custom.lattice <- custom.lattice()
custom.lattice$superpose.line$col <- heat.custom
custom.lattice$superpose.polygon$col <- heat.custom
custom.lattice$superpose.symbol$col <- heat.custom


# generalized function to create plots of different 'process' variables
plot.var <- function(dat, yvar, xvar = "time", groups, draw.panel.key = FALSE,
  comp = c("mai", "rib", "lhc", "pset", "cbm", "lpb"), 
  ylim = NULL, ylog = FALSE) {
  
  # first part of the plot is light intensity
  plot1 <- xyplot(c(0, concentration, 0) ~ c(get(xvar)[1], get(xvar), tail(get(xvar),1)), 
    subset(dat, variable == "hv" & simulation == simulation[1]),
    par.settings = custom.lattice,
    col = "#A6A6A641", border = 0,
    ylab = "% light", xlab = xvar,
    ylim = c(0, 100), panel = panel.polygon)
  
  # second part of the plot are all concentrations
  plot2 <- xyplot(concentration ~ get(xvar), 
    data = subset(dat, variable == yvar & component %in% comp),
    groups = get(groups),
    par.settings = custom.lattice, 
    ylab = yvar,
    scales = list(y = list(limits = ylim, log = ylog)),
    type = "l", lwd = 2, #pch = 19, cex = 0.3,
    panel = function(x, y, ...) {
      panel.grid(h = -1, v = -1, col = grey(0.9))
      panel.xyplot(x, y, ...)
      if (!is.na(comp)) {
        panel.key(labels = toupper(comp), cex = 1, col = grey(0.5),
          points = FALSE, lines = FALSE)
      }
      if (draw.panel.key) {
        panel.key(as.character(dat[[groups]]) %>% unique, cex = 0.6,
          points = FALSE, lines = TRUE, ...)
      }
    }
  )
  
  # combine in double Y scale plot
  doubleYScale(plot1, plot2, use.style = FALSE, add.axis = TRUE, add.ylab2 = TRUE)
}

# plot change in final biomass
plot.bm <- doubleYScale(add.ylab2 = TRUE,
  xyplot(concentration ~ simulation, 
    subset(dat, variable == "bm" & time == tail(unique(time), 1)),
    groups = simulation,
    scales = list(x = list(rot = 25)),
    par.settings = custom.lattice, 
    panel = function(x, y, ...) {
      panel.grid(h = -1, v = -1, col = grey(0.9))
      panel.barplot(x, y, border = custom.lattice$superpose.symbol$col, 
        ewidth = 0.2, ...)
    }
  ),
  xyplot(concentration/concentration[1]*100 ~ simulation,
    subset(dat, variable == "bm" & time == tail(unique(time), 1)),
    type = NA, ylab = "% biomass")
)


#svg("plot_simulation.svg", width = 7, height = 10)
print(plot.var(dat, yvar = "c", xvar = "time", comp = "rib", groups = "simulation", ylim = c(0,0.5)), split = c(1,1,2,4), more = TRUE)
print(plot.var(dat, yvar = "c", xvar = "time", comp = "lhc", groups = "simulation", ylim = c(0,0.5)), split = c(2,1,2,4), more = TRUE)
print(plot.var(dat, yvar = "c", xvar = "time", comp = "cbm", groups = "simulation", ylim = c(0,0.5)), split = c(1,2,2,4), more = TRUE)
print(plot.var(dat, yvar = "c", xvar = "time", comp = "pset", groups = "simulation", ylim = c(0,0.5)), split = c(2,2,2,4), more = TRUE)
print(plot.var(dat, yvar = "u", xvar = "time", comp = "rib", groups = "simulation", ylim = c(0,0.5)), split = c(1,3,2,4), more = TRUE)
print(plot.var(dat, yvar = "mu", xvar = "time", comp = NA, groups = "simulation", ylim = c(0,0.15)), split = c(2,3,2,4), more = TRUE)
print(plot.var(dat, yvar = "bm", xvar = "time", comp = NA, groups = "simulation", ylim = c(1,10^4), draw.panel.key = TRUE, ylog = 10), split = c(1,4,2,4), more = TRUE)
print(plot.bm, split = c(2,4,2,4))
#dev.off()


# plot the growth rate 'advantage' of ribosomal reserve versus no reserve
#svg("plot_growth_advantage.svg", width = 6, height = 4)
dat %>% filter(variable == "mu") %>%
  
  group_by(time) %>% mutate(deltaMu = concentration/concentration[1]) %>%

  xyplot(deltaMu ~ time, .,
    par.settings = custom.lattice, groups = simulation,
    type = "l", lwd = 2, 
    panel = function(x, y, ...) {
      panel.grid(h = -1, v = -1, col = grey(0.9))
      panel.abline(h = 1, col = grey(0.6), lty = 2, lwd = 1.5)
      panel.xyplot(x, y, ...)
      panel.key(as.character(.[["simulation"]]) %>% unique, cex = 0.6,
          points = FALSE, lines = TRUE, ...)
    }
  )
#dev.off()
