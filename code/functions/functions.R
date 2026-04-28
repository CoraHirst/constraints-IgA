############ Competition Model Popualtions ###########

#####################################################
# Antibody concentrations in the respiratory mucosa #
#####################################################

#non-specific ab at steady state: 
M_ns = function(c,d_mu) {c/d_mu} #c is constant flux due to generation, d_mu is decay rate in the mucosa

#specific antibodies during infection
M_s = function(A0, t_p, r, t) {A0/(1+exp(-r*(t-t_p)))} #sigmoidal growth, A0 is max antibody, t_p is half the time to peak, r is initial growt

#####################################################
# Antibody concentrations in the lumen #
#####################################################

#michaelis menten kinetics - desolve ode model
IgA_competition.model <- function(t, y, parms)  # Single partial immune class (RPS)                       
{with(as.list(parms),	# allows the parameter file parms to be as a list
      {  
        y = pmax(y,0)         # avoid negative values (not ideal but ok)
        
        # map the state variables
        L_s = y[1]         # specific Ab in the lumen
       
        # make empty variables for the derivatives
        dL_s = 0 #change in specific Ab in the lumen
        
        #####################################
        #####################################
        #	 calculate the derivatives
        dL_s = M_s(A0, t_p, r, t)*(Vmax/(K + M_s(A0, t_p, r, t) + M_ns(c = c, d = d_mu))) - d_l*L_s
          
        # output the derivatives
        dy=c(dL_s) #vector of derivatives 
        
        return(list(dy)) #lists vector and returns
      }
  ) 
} 

#michaelis menten kinetics - desolve ode model
total_transported.model <- function(t, y, parms)  # Single partial immune class (RPS)                       
{with(as.list(parms),	# allows the parameter file parms to be as a list
      {  
        y = pmax(y,0)         # avoid negative values (not ideal but ok)
        
        # map the state variables
        tot.IgA = y[1]         # specific Ab in the lumen
        
        # make empty variables for the derivatives
        dtot.IgA = 0 #change in specific Ab in the lumen
        
        #####################################
        #####################################
        #	 calculate the derivatives
        dtot.IgA = (M_s(A0, t_p, r, t) + M_ns(c = c, d = d_mu))*(Vmax/(K + M_s(A0, t_p, r, t) + M_ns(c = c, d = d_mu))) - d_l*tot.IgA
        
        # output the derivatives
        dy=c(dtot.IgA) #vector of derivatives 
        
        return(list(dy)) #lists vector and returns
      }
) 
} 

############ Competition Model Popualtions ###########

#####################################################
# Antibody concentrations in the respiratory mucosa #
secondary_response.model <- function(t, y, parms)  # Single partial immune class (RPS)                       
{with(as.list(parms),	# allows the parameter file parms to be as a list
      {  
        y = pmax(y,0)         # avoid negative values (not ideal but ok)
        
        # map the state variables
        tot.IgA = y[1]         # specific Ab in the lumen
        
        # make empty variables for the derivatives
        dtot.IgA = 0 #change in specific Ab in the lumen
        
        #####################################
        #####################################
        #	 calculate the derivatives
        dtot.IgA = (M_s(A0, t_p, r, t) + M_ns(c = c, d = d_mu))*(Vmax/(K + M_s(A0, t_p, r, t) + M_ns(c = c, d = d_mu))) - d_l*tot.IgA
        
        # output the derivatives
        dy=c(dtot.IgA) #vector of derivatives 
        
        return(list(dy)) #lists vector and returns
      }
) 
} 