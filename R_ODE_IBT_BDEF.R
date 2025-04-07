# Load the deSolve package
rm(list = ls());gc()
library(deSolve)

setwd('D:\\data\\Theoretical studies\\feedback\\')


minv <- 10^-10
# Define the ODE system
ode_system <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    if(N<=minv){N <- 0}
    if(C<=minv){C <- 0}
    if(Z<=minv){Z <- 0}
    if(R<=minv){R <- 0}
      
    F1 <- mu_n * (R^bn)  * N / (1.0 + b1 * N)
    F2 <- mu_z * (R^bz)  * C / (1.0 + b2 * C)
    
    dN <- q * (N0 - N) - F1 * C
    dC <- F1 * (1.0 - C) * C - F2 * Z
    dZ <- (tz * F2 - dz) * Z

    # Example: dXA/dt = XA - Gd * X
    dR <- Im - (Im / Rm) * R -(theta * N^2 + theta * eb * N + ec * (1.0 + (theta - 1.0) * eb^2 / (4.0 * ec))) * R/Rm
    
    # Return the list of derivatives
    return(list(c(dN, dC, dZ, dR)))
  })
}


##############################################################################################################
## Weak feedback scenario
# Define initial conditions and parameters
parameters <-     c(theta = 25, bn= 0.05, bz = 0.09, b1 = 3.3, b2 = 2.5, mu_n = 0.6, mu_z = 0.3, dz = 0.01, Rm = 30.0, Im = 1.0, 
                    tz = 0.4, eb = -0.9, ec = 1.0, N0 = 0.55, q = 0.2)
# # extinction parameter
# normal parameter: theta = 100, bn= 0.2, bz = 0.2, 
state <- c(N = 1.0, C = .15, Z = 1.0, R = 15)  # Example initial values

times <- seq(0, 100*365, by = 0.01)
out.burn <- ode(y = state, times = times, func = ode_system, parms = parameters, method = 'rk4')

# Solve the ODE system
state.new <- out.burn[NROW(out.burn),-1]

times <- seq(0, 10*365, by = 0.01)

out <- ode(y = state.new, times = times, func = ode_system, parms = parameters, method = 'rk4')

# Plot the results
windows()
plot(out, subset=time>1000)

out.final <- out

zimpact=parameters['mu_z']*(out.final[,'R']^parameters['bz'])*(1/(out.final[,'C']+parameters['b2']))
c_impact_n=parameters['mu_n']*(out.final[,'R']^parameters['bn'])*(1/(out.final[,'N']+parameters['b1']))

out.final <- data.frame(out.final, CiN=c_impact_n, ZiC=zimpact)


write.csv(out.final, 'out_weakfb.csv', row.names = F)


##############################################################################################################
## Strong feedback scenario
# Define initial conditions and parameters
parameters <-  c(theta = 100, bn= 0.2, bz = 0.2, b1 = 3.3, b2 = 2.5, mu_n = 0.6, mu_z = 0.3, dz = 0.01, Rm = 30.0, Im = 1.0, 
                    tz = 0.4, eb = -0.9, ec = 1.0, N0 = 0.55, q = 0.2)
# # extinction parameter
# normal parameter: theta = 100, bn= 0.2, bz = 0.2, 
state <- c(N = 1.0, C = .15, Z = 1.0, R = 15)  # Example initial values

times <- seq(0, 100*365, by = 0.01)
out.burn <- ode(y = state, times = times, func = ode_system, parms = parameters, method = 'rk4')

# Solve the ODE system
state.new <- out.burn[NROW(out.burn),-1]

times <- seq(0, 10*365, by = 0.01)

out <- ode(y = state.new, times = times, func = ode_system, parms = parameters, method = 'rk4')

# Plot the results
windows()
plot(out, subset=time>1000)


out.final <- out

zimpact=parameters['mu_z']*(out.final[,'R']^parameters['bz'])*(1/(out.final[,'C']+parameters['b2']))
c_impact_n=parameters['mu_n']*(out.final[,'R']^parameters['bn'])*(1/(out.final[,'N']+parameters['b1']))

out.final <- data.frame(out.final, CiN=c_impact_n, ZiC=zimpact)

write.csv(out.final, 'out_strongfb.csv', row.names = F)


# Define the time points at which to solve the ODE
times.burn <- seq(0, 36500, by = 0.1)

# Define initial conditions and parameters
parameters <-     c(theta = 25, bn= 0.05, bz = 0.2, mu_n = 0.6, mu_z = 0.3, b1 = 3.3, b2 = 2.5, dz = 0.01, Rm = 30.0, Im = 1.0, 
                    tz = 0.4, eb = -0.9, ec = 1.0, N0 = 0.55, q = 0.2)


# # Extinction scenario
# normal parameter: theta = 100, bn= 0.2, bz = 0.2, 
state <- c(N = 1.0, C = .15, Z = 1.0, R = 15)  # Example initial values

times <- seq(0, 10*365, by = 0.01)
out.burn <- ode(y = state, times = times, func = ode_system, parms = parameters, method = 'rk4')

# Solve the ODE system
state.new <- out.burn[NROW(out.burn),-1]

times <- seq(0, 1*365, by = 0.01)

out <- ode(y = state.new, times = times, func = ode_system, parms = parameters, method = 'rk4')

# Plot the results
windows()
plot(out, subset=time<365)

out.final <- out.burn

zimpact=parameters['mu_z']*(out.final[,'R']^parameters['bz'])*(1/(out.final[,'C']+parameters['b2']))
c_impact_n=parameters['mu_n']*(out.final[,'R']^parameters['bn'])*(1/(out.final[,'N']+parameters['b1']))

out.final <- data.frame(out.final, CiN=c_impact_n, ZiC=zimpact)

write.csv(out.final, 'outburn_extinct.csv', row.names = F)

###################################################################################
## Ploting ver 2
ind.t <- 1:(365*2)*100
weakfb <- read.csv('out_weakfb.csv', header = T)[ind.t,]
strongfb <- read.csv('out_strongfb.csv', header = T)[ind.t,]

weakfb[,-1] <- apply(weakfb[,-1],2,function(x){x-mean(x)})
strongfb[,-1] <- apply(strongfb[,-1],2,function(x){x-mean(x)})

yrg <- apply(rbind(weakfb,strongfb),2,range,na.rm=T)

exfb <- read.csv('outburn_extinct.csv', header = T)
exfb <- exfb[1:(10*1)*100,]
#exfb[,-1] <- apply(exfb[,-1],2,function(x){x-mean(x)})

win.graph(20,40);par(mfcol=c(5,1),mar=c(4,4,1,1))
plot(R~time,weakfb,type='l',ylim=yrg[,'R'],ylab='Phytoplankton species richness (R)',xlab='',lwd=2, lty=1, col='grey')
lines(R~time,strongfb, type='l',ylim=yrg[,'R'],lwd=2)

plot(C~time,weakfb,type='l',ylim=yrg[,'C'],ylab='Phytoplankton biomass (C)',xlab='',lwd=2, lty=1, col='grey')
lines(C~time,strongfb, type='l',ylim=yrg[,'C'],lwd=2)

plot(Z~time,weakfb,type='l',ylim=yrg[,'Z'],ylab='Zooplankton biomass (Z)',xlab='',lwd=2, lty=1, col='grey')
lines(Z~time,strongfb, type='l',ylim=yrg[,'Z'],lwd=2)

plot(CiN~time,weakfb,type='l',ylim=yrg[,'CiN'],ylab='Phytoplankton nutrient uptake rate ',xlab='',lwd=2, lty=1, col='grey')
lines(CiN~time,strongfb, type='l',ylim=yrg[,'CiN'],lwd=2)

plot(ZiC~time,weakfb,type='l',ylim=yrg[,'ZiC'],ylab='Zooplankton effects on phytoplankton',xlab='',lwd=2, lty=1, col='grey')
lines(ZiC~time,strongfb, type='l',ylim=yrg[,'ZiC'],lwd=2)
legend('topright',c('Scenario:\nWeak feedback','Strong feedback'),lty=c(1,1),col=c('grey',1),bty='n')


win.graph(20,40);par(mfcol=c(4,1),mar=c(4,4,1,1))

plot(C~time,exfb,type='l',ylab='Phytoplankton biomass (C)',xlab='',lwd=2, lty=1,ylim=c(0,0.12))

plot(R~time,exfb,type='l',ylab='Phytoplankton species richness (R)',xlab='',lwd=2, lty=1)

plot(Z~time,exfb,type='l',ylab='Zooplankton biomass (Z)',xlab='',lwd=2, lty=1)

plot(N~time,exfb,type='l',ylab='Nutrient (N)',xlab='',lwd=2, lty=1)

