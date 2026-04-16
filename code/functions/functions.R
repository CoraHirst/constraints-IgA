# michaelis menten kinetics with competition
specific_flux = function(Vmax, K, M_s, M_ns) { 
  M_s*Vmax/(K + M_s + M_ns)} #define flux

# calculate steady state fraction of nonspecific IgA function
Mns_SteadyState = function(A0, tau, f, Vmax_0, d_mu, K) {
  gamma = A0/tau + f
  SS = (gamma - Vmax_0 - d_mu*K + sqrt((Vmax_0+d_mu*K-gamma)^2 + 4*d_mu*gamma*K))/(2*d_mu)
  return(SS)}

# ODE form
Mns.odes.model <- function(t, y, parms) 
{with(as.list(parms),	# allows the parameter file parms to be as a list
      {  
        # map the state variables
        Mns = y[1]
       
        ######################################
        gamma = A0/tau + f # calculate flux into mucosa
        
        ######################################
        #	 calculate the derivatives
        dMns = gamma - d_mu*Mns - (Vmax_0/(K+Mns))*Mns
        
        # output the derivatives
        dy=c(dMns)
        return(list(dy)) 
      }
    )
  } 


# Concentration of mucosal specific IgA during new response
M_specific = function(A0, r, rho, d_mu, t) {(A0*exp(-d_mu*t))/(1+exp(-r*(t-rho)))}