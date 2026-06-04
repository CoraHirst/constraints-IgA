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
        M_s = y[1]
        M_ns = y[2]
        R_s = y[3]
        R_ns = y[4]
        L_s = y[5]
        L_ns = y[6]
        
        # make empty variables for the derivatives
        dM_s = NaN
        dM_ns = NaN
        dR_s = NaN
        dR_ns = NaN
        dL_s = NaN
        dL_ns = NaN
        
        #####################################
        #calculate Re 
        R_e = R_T - R_s - R_ns
        #####################################
        #	 calculate the derivatives
        dM_s = -k*M_s*R_e - d_mu*M_s
        dM_ns = -k*M_ns*R_e - d_mu*M_ns
        dR_s = k*M_s*R_e - r*R_s
        dR_ns = k*M_ns*R_e - r*R_ns
        dL_s = r*R_s - d_L*L_s
        dL_ns = r*R_ns - d_L*L_ns
        # output the derivatives
        dy=c(dM_s, dM_ns, dR_s, dR_ns, dL_s, dL_ns) #vector of derivatives 
        
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