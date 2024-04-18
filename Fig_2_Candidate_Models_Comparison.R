graphics.off() 
rm(list=ls(all=TRUE)) 

# Load the module 
source("DBDA2E-utilities.R")
library(rjags)
library(runjags)
library(VGAM)
library(bayesmeta)

#------------------------------------------------------------------------------
#THE DATA.
mydataframe = read.csv( file="BATHROOM_OTOOLE_2012.csv" )
C = as.numeric(mydataframe[, "Count"])  # MPN
V = as.numeric(mydataframe[, "Volume"]) # mL
N = length(C)

#PGA model------------------------------------------------------------------------------

# Package the data for shipping to JAGS:
dataList = list(
  C = C ,
  V = V ,
  N = N
)

# THE MODEL
modelstring ="
model {
for( i in 1 : N ) {

#Likelihood 
C[i] ~ dpois(lambda_1mL[i]*V[i])
lambda_1mL[i] = lambda_100mL[i]/100
lambda_100mL[i] ~ dgamma(shape,rate)

}

# Prior
shape ~ dgamma(0.01, 0.01)
rate ~ dgamma(0.01, 0.01)

}
" # close quote for modelstring
writeLines(modelstring,con="model.txt")
#------------------------------------------------------------------------------
# INTIALIZE THE CHAINS.
# Let JAGS do it
#------------------------------------------------------------------------------
# RUN THE CHAINS
require(rjags)
parameters = c("shape","rate") 
adaptSteps = 5000         # Number of steps to "tune" the samplers.
burnInSteps = 10^4         # Number of steps to "burn-in" the samplers.
nChains = 4               # Number of chains to run.
numSavedSteps=10^5       # Total number of steps in chains to save.
thinSteps=1               # Number of steps to "thin" (1=keep every step).
nPerChain = ceiling( ( numSavedSteps * thinSteps ) / nChains ) # Steps per chain.
# Create, initialize, and adapt the model:
jagsModel_PGA = jags.model( "model.txt" , data=dataList ,
                        n.chains=nChains , n.adapt=adaptSteps )
# Burn-in:
cat( "Burning in the MCMC chain...\n" )
update( jagsModel_PGA , n.iter=burnInSteps )
# The saved MCMC chain:
cat( "Sampling final MCMC chain...\n" )
mcmcCoda_PGA = coda.samples( jagsModel_PGA , variable.names=parameters ,
                         n.iter=nPerChain , thin=thinSteps )

#------------------------------------------------------------------------------
# EXAMINE THE RESULTS

# Convert coda-object codaSamples to matrix object for easier handling.
mcmcChain_PGA = as.matrix( mcmcCoda_PGA )
chainLength = NROW(mcmcChain_PGA)

# Display diagnostics of chain, for specified parameters:
parameterNames = varnames(mcmcCoda_PGA)
for ( parName in parameterNames ) {
  diagMCMC( codaObject=mcmcCoda_PGA , parName=parName)
}


#PLN model------------------------------------------------------------------------------

# Package the data for shipping to JAGS:
dataList = list(
  C = C ,
  V = V ,
  N = N
)

# THE MODEL.
modelstring ="
model {
for( i in 1 : N ) {

#Likelihood 
C[i] ~ dpois(lambda_1mL[i]*V[i])
lambda_1mL[i] = lambda_100mL[i]/100
lambda_100mL[i] ~ dlnorm(muOfLog,1/sigmaOfLog^2)

}

# Prior
muOfLog ~ dunif( -10 , 10)
sigmaOfLog ~ dexp(0.1)

}
" # close quote for modelstring
writeLines(modelstring,con="model.txt")
#------------------------------------------------------------------------------
# INTIALIZE THE CHAINS.
# Let JAGS do it
#------------------------------------------------------------------------------
# RUN THE CHAINS
require(rjags)
parameters = c("muOfLog","sigmaOfLog") 
adaptSteps = 5000         # Number of steps to "tune" the samplers.
burnInSteps = 10^4        # Number of steps to "burn-in" the samplers.
nChains = 4               # Number of chains to run.
numSavedSteps=10^5       # Total number of steps in chains to save.
thinSteps=1               # Number of steps to "thin" (1=keep every step).
nPerChain = ceiling( ( numSavedSteps * thinSteps ) / nChains ) # Steps per chain.
# Create, initialize, and adapt the model:
jagsModel_PLN = jags.model( "model.txt" , data=dataList ,
                        n.chains=nChains , n.adapt=adaptSteps )
# Burn-in:
cat( "Burning in the MCMC chain...\n" )
update( jagsModel_PLN , n.iter=burnInSteps )
# The saved MCMC chain:
cat( "Sampling final MCMC chain...\n" )
mcmcCoda_PLN = coda.samples( jagsModel_PLN , variable.names=parameters ,
                         n.iter=nPerChain , thin=thinSteps )

#------------------------------------------------------------------------------
# EXAMINE THE RESULTS

# Convert coda-object codaSamples to matrix object for easier handling.
mcmcChain_PLN = as.matrix(mcmcCoda_PLN)
chainLength = NROW(mcmcChain_PLN)

# Display diagnostics of chain, for specified parameters:
parameterNames = varnames(mcmcCoda_PLN)
for ( parName in parameterNames ) {
  diagMCMC( codaObject=mcmcCoda_PLN , parName=parName)
}


#PLO model------------------------------------------------------------------------------

# Package the data for shipping to JAGS:
dataList = list(
  C = C ,
  V = V ,
  N = N
)

# THE MODEL.
modelstring ="
model {
for( i in 1 : N ) {

#Likelihood PLN 
C[i] ~ dpois(lambda_1mL[i]*V[i])
lambda_1mL[i] = lambda_100mL[i]/100
lambda_100mL[i] ~ dlomax(alpha,sigma)

}

# Prior
alpha ~ dunif(0, 1000)
sigma ~ dunif(0, 1000)

}
" # close quote for modelstring
writeLines(modelstring,con="model.txt")
#------------------------------------------------------------------------------
# INTIALIZE THE CHAINS.
# Let JAGS do it
#------------------------------------------------------------------------------
# RUN THE CHAINS
require(rjags)
parameters = c("alpha","sigma") 
adaptSteps = 5000         # Number of steps to "tune" the samplers.
burnInSteps = 10^4        # Number of steps to "burn-in" the samplers.
nChains = 4               # Number of chains to run.
numSavedSteps=10^5        # Total number of steps in chains to save.
thinSteps=1               # Number of steps to "thin" (1=keep every step).
nPerChain = ceiling( ( numSavedSteps * thinSteps ) / nChains ) # Steps per chain.
# Create, initialize, and adapt the model:
jagsModel_PLO = jags.model( "model.txt" , data=dataList ,
                        n.chains=nChains , n.adapt=adaptSteps )
# Burn-in:
cat( "Burning in the MCMC chain...\n" )
update( jagsModel_PLO , n.iter=burnInSteps )
# The saved MCMC chain:
cat( "Sampling final MCMC chain...\n" )
mcmcCoda_PLO = coda.samples( jagsModel_PLO , variable.names=parameters ,
                         n.iter=nPerChain , thin=thinSteps )

#------------------------------------------------------------------------------
# EXAMINE THE RESULTS

# Convert coda-object codaSamples to matrix object for easier handling.
mcmcChain_PLO = as.matrix( mcmcCoda_PLO)
chainLength = NROW(mcmcChain_PLO)

# Display diagnostics of chain, for specified parameters:
parameterNames = varnames(mcmcCoda_PLO)
for ( parName in parameterNames ) {
  diagMCMC( codaObject=mcmcCoda_PLO, parName=parName)
}


# Plot CCDF ------------------------------------------------------------------------------

tiff(file = paste(".tiff", sep = ""), width = 2000, height = 2000, 
     units = "px",res = 1000, pointsize = 7)

options(scipen = 0)
point_x = 1
point_y = 1
myTicks_x<-c(1.0E+0, 1.0E+1, 1.0E+2, 1.0E+3, 1.0E+4, 1.0E+5, 1.0E+6, 1.0E+7)
myTicks_y<-c(0.01, 0.1, 1)
par(mar=c(4.1, 4.1, 1.1, 1.1))

plot(point_x, point_y, col="white", pch = 19,
     xlab = expression(paste(italic("E. coli"), " (MPN/100 mL)")),
     ylab = expression(italic("P") * group("(", list(italic("X") >= italic("x")), ")")),  # Mathematical notation for y-axis label 
     xaxt="n", 
     yaxt="n", 
     xlim = range(myTicks_x),ylim=range(myTicks_y),
     font.lab = 2, log='xy')

axis(side = 1, at = myTicks_x)
axis(side = 2, at = myTicks_y)


# Best-fit PGA ------------------------------------------------------------------------------
PGA_data=rgamma(n=500000, shape=(median(mcmcChain_PGA[,"shape"])), rate=median(mcmcChain_PGA[,"rate"]))
PGA_x=sort(PGA_data)
PGA_y = 1-ecdf(PGA_data)(sort(PGA_data) )
lines(PGA_x,PGA_y, col="#007BFF", lwd=1,lty=2)


# Best-fit PLN ------------------------------------------------------------------------------
PLN_data=rlnorm(n=500000, meanlog=(median(mcmcChain_PLN[,"muOfLog"])), sdlog=median(mcmcChain_PLN[,"sigmaOfLog"]))
PLN_x=sort(PLN_data)
PLN_y = 1-ecdf(PLN_data)(sort(PLN_data) )
lines(PLN_x,PLN_y, col="#28A745", lwd=1,lty=1)

# Best-fit PLO ------------------------------------------------------------------------------
PLO_data=rlomax(n=500000, shape=(median(mcmcChain_PLO[,"alpha"])), scale=median(mcmcChain_PLO[,"sigma"]))
PLO_x=sort(PLO_data)
PLO_y = 1-ecdf(PLO_data)(sort(PLO_data) )
lines(PLO_x,PLO_y, col="#FF5733", lwd=1,lty=4)

# Calculate ECDF values and shift them to get "greater than or equal to"
point_x = sort(C)
ecdf_values = ecdf(C)(point_x)
point_y = 1 - ecdf_values + 1/length(C)
points(point_x, point_y ,col = "black", bg = "black", pch = 23, lwd = 0.3, cex = 0.75)


legend("topright",                     # Position of legend within plot
       legend = c("PGA", 
                  "PLN", 
                  "PLO"),  # Labels in legend
       col = c("#007BFF", "#28A745", "#FF5733"), # Colors of lines/points in legend
       lty = c(2, 1, 4),               # Line types
       lwd = c(1, 1, 1),               # Line widths
       cex = 0.8,                      # Font size for legend
       bty = "n"                       # Remove box around legend
)

dev.off()
