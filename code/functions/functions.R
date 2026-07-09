############ Competition Model Popualtions ###########

#####################################################
# Antibody concentrations in the lumen #
#####################################################
# Antibody concentrations in the respiratory mucosa #
primary_response.model <- function(t, y, parms)  # Single partial immune class (RPS)                       
{with(as.list(parms),	# allows the parameter file parms to be as a list
      {  
        y = pmax(y,0)         # avoid negative values (not ideal but ok)
        
        # map the state variables
        Ms = y[1]         # Antibodies in mucosa
        
        # make empty variables for the derivatives
        dMs = 0      # change in specific Ab in the mucosa
        
        #####################################
        # define time dependent functions   #
        #####################################
        r = function(t) {if(t <= tp) {r0 - alpha*t}
          else 0} #define time-varying r(t) - increase that slows and is then dominated by decay rate
        #####################################
        ## calculate the derivatives
        # cell populations
        dMs = r(t)*Ms - d_mu*Ms # specific Ab production and decay
        
        # output the derivatives
        dy=c(dMs) #vector of derivatives 
        
        return(list(dy)) #lists vector and returns
      }
) 
} 

# nonspecific antibodies dynamics
equilibrium_ns.model <- function(t, y, parms)  # Single partial immune class (RPS)                       
{with(as.list(parms),	# allows the parameter file parms to be as a list
      {  
        y = pmax(y,0)         # avoid negative values (not ideal but ok)
        
        # map the state variables
        Mns = y[1]         # Antibodies in mucosa
        Rns = y[2]
        Lns = y[3] 
        
        # make empty variables for the derivatives
        dMns = 0      # change in specific Ab in the mucosa
        dRns = 0 
        dLns = 0
        
        ################## Time varying var ############
        Re = RT - Rns
        #####################################
        ## calculate the derivatives
        # cell populations
        dMns = c - k1*Re*Mns - d_mu*Mns  # specific Ab production and decay
        dRns = k1*Re*Mns - k2*Rns
        dLns = k2*Rns - d_L*Lns
        # output the derivatives
        dy=c(dMns, dRns, dLns) #vector of derivatives 
        
        return(list(dy)) #lists vector and returns
      }
) 
} 

# michaelis menten kinetics - desolve ode model
IgA_competition.model <- function(t, y, parms)  # Single partial immune class (RPS)                       
{with(as.list(parms),	# allows the parameter file parms to be as a list
      {  
        y = pmax(y,0)         # avoid negative values (not ideal but ok)
        
        # map the state variables
        Ms = y[1]
        Mns = y[2]
        Rs = y[3]
        Rns = y[4]
        Ls = y[5]
        Lns = y[6]
        
        # make empty variables for the derivatives
        dMs = NaN
        dMns = NaN
        dRs = NaN
        dRns = NaN
        dLs = NaN
        dLns = NaN
        
        #####################################
        # define time dependent functions   #
        #####################################
        # exponential growth rate of Ms at time t
        r = function(t) {if(t < tp) {r0 - alpha*t}
          else 0} #define time-varying r(t) - increase that slows and is then dominated by decay rate
        
        # calculate Re at time t
        Re = RT - Rs - Rns # total number of receptors - those bound to specific ab - those bound to nonspecific ab
        
        #####################################
        #	 calculate the derivatives
        dMs = r(t)*Ms - k1*Ms*Re - d_mu*Ms  # specific Ab production and decay
        dMns = c - k1*Re*Mns - d_mu*Mns 
        dRs = k1*Ms*Re - k2*Rs
        dRns = k1*Re*Mns - k2*Rns
        dLs = k2*Rs - d_L*Ls
        dLns = k2*Rns - d_L*Lns
        
        # output the derivatives
        dy=c(dMs, dMns, dRs, dRns, dLs, dLns) #vector of derivatives 
        
        return(list(dy)) #lists vector and returns
      }
  ) 
} 



############ Recall Responses Model Popualtions ###########

#####################################################
# Antibody concentrations in the respiratory mucosa #
secondary_response.model <- function(t, y, parms)  # Single partial immune class (RPS)                       
{with(as.list(parms),	# allows the parameter file parms to be as a list
      {  
        y = pmax(y,0)         # avoid negative values (not ideal but ok)
        
        # map the state variables
        B = y[1]         # B cells 
        P_s = y[2]       # Antibody secreting cells (specific)
        P_q = y[3]       # quiescent antibody secreting cells (specific)
        M_s = y[4]       # specific Ab in mucosa
        L_s = y[5]       # specific Ab in lumen
        
        
        # make empty variables for the derivatives
        dB = 0      # change in specific Ab in the lumen
        dP_s = 0    # change in specific ASCs
        dP_q = 0    # change in quiescent ASCs
        dM_s = 0    # change in specific antibodies in mucosa
        dL_s = 0    # change in specific antibodies in the lumen
        
        #####################################
        # define V function 
        V.function = function(t_infection, t) {
          if(t>=min(t_infection) & t <= max(t_infection)){
            V = 1}
          else(V = 0)
        return(V)}
        # set V
        V = V.function(t_infection, t)
        #####################################
        ## calculate the derivatives
        # cell populations
        dB = r_b*(1-diff_frac)*B*V - diff_frac*B*V -d_B*B #B cells (memory and naive together)
        dP_s = diff_frac*B*V - quiessence_k*P_s - d_mu*P_s + induction_k*P_q*V # specific ASCs
        dP_q = quiessence_k*P_s - induction_k*P_q*V - d_q*P_q # quiescent sepecific ASCs
        # Specific antibodies        
        dM_s = rho*P_s - (Vmax/(K + M_s + M_ns(c = c, d = d_mu)))*M_s - d_Ab*M_s #mucosal
        dL_s = (Vmax/(K + M_s + M_ns(c = c, d = d_mu)))*M_s - d_l*L_s #lumenal
        
        # output the derivatives
        dy=c(dB, dP_s, dP_q, dM_s, dL_s) #vector of derivatives 
        
        return(list(dy)) #lists vector and returns
      }
  ) 
} 