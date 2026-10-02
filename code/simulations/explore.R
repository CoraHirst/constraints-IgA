### setup
#install.package("here")
library(here)

#install.package("ggplot2")
library(ggplot2)

#install.package("gridExtra")
library(gridExtra)

#install.packages("tidyverse")
library(tidyverse)

#install.packages("deSolve")
library(deSolve)

#load functions
MM_function = function(M_s, M_ns) {(Vmax*M_s)/(KM+M_ns+M_s)}


#choose plotting colors
my_colors = c("darkred", "maroon", "darkblue", "forestgreen", "darkorange", "goldenrod") #choose color range
num_colors = length(my_colors) #set number of colors needed for plotting
my_colors = colorRampPalette(my_colors)(num_colors) #increase number of colors staying in color range

####  heatmap way 1: lets look at the plots for a few different Ap/Thalf combos

# generation of new response
Aps = 10^c(1, 2, 3) #peak antibody titer will be the same as Mns for 90 days
#initial ab titer
A0 = 1 # initial specific ab concentration
tp = 30 #time of peak titer -- days

# steady state nonspecific antibodies
t_halfs = 3*10^c(1,2,3) #half lives of mucosal ab in days - 1 month, 3 months, 1 year, 10 years :) 


heatmap_list = lapply(Aps, function(Ap){ #element
  lapply(t_halfs, function(t_half){ #cols
    # setting t.half dependent parms for nonspecific population
    c = Ap*(log(2)/min(t_halfs))# calculate c such that the minimum M_ns = Ap = 100 
    d_mu = log(2)/t_half #set d_mu
    Mns = c/d_mu # steadys state nonspecific
    
    # setting t.half dependent parms for specific abs
    r0 = (2*log(Ap/A0)/tp) + d_mu #initial/greatest exponential expansion rate
    alpha = (r0 - d_mu)/tp # rate of decay of exponential growth rate
    
    # set time vector parms
    Tmax = Tmax # total days
    points_per_day = 24 # every hour
    time = seq(0, Tmax, length.out = Tmax*points_per_day) #generate time vector
    
    ### apply to calculate Ms titers over time
    Ms_soln = data.frame("time" = time, "t_half" = t_half, "Ms" = sapply(time, function(t){
      if(t<=tp){A0*exp((r0 - d_mu)*t - (alpha/2)*t^2)}
      else {Ap*exp(-d_mu*(t-tp))}
    }), "Mns" = Mns, "Ap" = Ap)
    
    ### apply to calculate Ls flux over time
    Ms_soln$phi_s = sapply(Ms_soln$Ms, function(Ms){MM_function(M_s = Ms, M_ns = Mns)})
    Ms_soln$phi_ns = sapply(Ms_soln$Ms, function(Ms){(Vmax*Mns)/(KM+Ms+Mns)})
    return(Ms_soln)
  }) %>% do.call(rbind, .)
}) 


plots_list = lapply(heatmap_list, function(MM_soln){

#ggplot
plot = ggplot() + geom_line(data = MM_soln, aes(x = time, y = (phi_s), group = factor(t_half), col = factor(t_half), linetype = "phi_s"), lwd = 1.5) +
  geom_line(data = MM_soln, aes(x = time, y = (phi_ns), group = factor(t_half), col = factor(t_half), linetype = "phi_ns"), lwd = 1.5) + 
  scale_color_manual(values = my_colors) +
  #coord_cartesian(ylim = c(0, 3), xlim = c(0, 15)) +
  scale_linetype_manual(values = c("phi_s" = "solid", "phi_ns" = "dashed"), labels = c("phi_s" = expression(phi[s]), "phi_ns" = expression( phi[ns]))) +
  labs(x = "time, days", y = bquote(phi[Ab] ~ "day"^-1), col = bquote("T"[1/2] ~ ", days"^-1), linetype = "Specificity", title = "flux of specific and non-specific antibody into the lumen over the course of infection") +
  theme_bw()

  return(plot)
})

do.call(grid.arrange, plots_list)


###

### Heatmap way 2

# generation of new response -fix the rate of generation c for the models (that doesn't seem fair right?)
Ap = 100 #peak antibody titer will be the same as Mns for 90 days
#initial ab titer
A0 = 1 # initial specific ab concentration
tp = 30 #time of peak titer -- days

# steady state nonspecific antibodies
t_halfs = c(30,90,400,4000) #half lives of mucosal ab in days - 1 month, 3 months, 1 year, 10 years :) 
d_mus = log(2)/t_halfs #half lives in days
names(d_mus) = t_halfs
c = Ap*(log(2)/min(t_halfs))# calculate c such that the minimum M_ns = Ap = 100
Mns_vec = c/d_mus #vector of approximate ss solutions
names(Mns_vec) = t_halfs

# transport
Vmax = 10^3 # order of magnitude of peak response # we might vary this later
k2 = 1/(2/24) #rate of transport - 1/T of transport which is 1/2hours or 1/2/24ths of a day
RT = Vmax/k2 #total number (concentration) of receptors, Vmax = k2RT
KM = Ap/3 #choosing KM such that at peak in the absence of Mns, Ms would achieve 75% of max flux at peak
k1 = k2/KM #k2/min(Mns_vec) #rate of antibody binding such that KM is on the order of 10

### parms for heatmap 
Aps = 10^c(1, 2, 3) #peak antibody titer will be the same as Mns for 90 days

t_halfs = 3*10^c(1,2,3) #half lives of mucosal ab in days - 1 month, 3 months, 1 year, 10 years :) 

heatmap_list = lapply(Aps, function(Ap){ #element
  lapply(t_halfs, function(t_half){ #cols
    # setting t.half dependent parms for nonspecific population
    d_mu = log(2)/t_half #set d_mu
    Mns = c/d_mu # steadys state nonspecific
    
    # setting t.half dependent parms for specific abs
    r0 = (2*log(Ap/A0)/tp) + d_mu #initial/greatest exponential expansion rate
    alpha = (r0 - d_mu)/tp # rate of decay of exponential growth rate
    
    # set time vector parms
    Tmax = Tmax # total days
    points_per_day = 24 # every hour
    time = seq(0, Tmax, length.out = Tmax*points_per_day) #generate time vector
    
    ### apply to calculate Ms titers over time
    Ms_soln = data.frame("time" = time, "t_half" = t_half, "Ms" = sapply(time, function(t){
      if(t<=tp){A0*exp((r0 - d_mu)*t - (alpha/2)*t^2)}
      else {Ap*exp(-d_mu*(t-tp))}
    }), "Mns" = Mns, "Ap" = Ap)
    
    ### apply to calculate Ls flux over time
    Ms_soln$phi_s = sapply(Ms_soln$Ms, function(Ms){MM_function(M_s = Ms, M_ns = Mns)})
    Ms_soln$phi_ns = sapply(Ms_soln$Ms, function(Ms){(Vmax*Mns)/(KM+Ms+Mns)})
    return(Ms_soln)
  }) %>% do.call(rbind, .)
}) 


plots_list = lapply(heatmap_list, function(MM_soln){
  
  #ggplot
  plot = ggplot() + geom_line(data = MM_soln, aes(x = time, y = (phi_s), group = factor(t_half), col = factor(t_half), linetype = "phi_s"), lwd = 1.5) +
    geom_line(data = MM_soln, aes(x = time, y = (phi_ns), group = factor(t_half), col = factor(t_half), linetype = "phi_ns"), lwd = 1.5) + 
    scale_color_manual(values = my_colors) +
    #coord_cartesian(ylim = c(0, 3), xlim = c(0, 15)) +
    scale_linetype_manual(values = c("phi_s" = "solid", "phi_ns" = "dashed"), labels = c("phi_s" = expression(phi[s]), "phi_ns" = expression( phi[ns]))) +
    labs(x = "time, days", y = bquote(phi[Ab] ~ "day"^-1), col = bquote("T"[1/2] ~ ", days"^-1), linetype = "Specificity", title = "flux of specific and non-specific antibody into the lumen over the course of infection") +
    theme_bw()
  
  return(plot)
})

do.call(grid.arrange, plots_list)
