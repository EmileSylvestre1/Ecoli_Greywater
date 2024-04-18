graphics.off() 
rm(list=ls(all=TRUE)) 

#------------------------------------------------------------------------------
#THE DATA.
mydataframe = read.csv( file="BATHROOM_OTOOLE_2012.csv" )
C = as.numeric(mydataframe[, "Count"]) #MPN
V = as.numeric(mydataframe[, "Volume"]) #mL
N = length(C)

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

# Data point ------------------------------------------------------------------------------

# Parameters Jahne et al. (2017)
#Wahsing machine - Wash: meanlog=0.69, sdlog=3.61
#Wahsing machine - Rinse:meanlog=0.00, sdlog=3.20
#Boathroom:meanlog=4.87, sdlog=1.88

lnorm_data=rlnorm(n=500000, meanlog=4.87, sdlog=1.88) #(to change manually)
lnorm_x=sort(lnorm_data)
lnorm_y = 1-ecdf(lnorm_data)(sort(lnorm_data) )
lines(lnorm_x,lnorm_y, col="#6A1B9A", lwd=1,lty=1)

# Calculate ECDF values and shift them to get "greater than or equal to"
point_x = sort(C)
ecdf_values = ecdf(C)(point_x)
point_y = 1 - ecdf_values + 1/length(C)
points(point_x, point_y ,col = "black", bg = "black", pch = 23, lwd = 0.5, cex = 0.75)

#Empirical CCDF
lines(point_x, point_y, col="#767676", lwd=0.5, lty=1)  # Adjust the color, line width, and line type as needed

legend("topright",                     # Position of legend within plot
       legend = c("LN (approximation)",
                  "Empirical"),  # Labels in legend
       col = c("#6A1B9A", "#767676"), # Colors of lines/points in legend
       lty = c(1, 1),               # Line types
       lwd = c(1, 1),               # Line widths
       cex = 0.8,                      # Font size for legend
       bty = "n"                       # Remove box around legend
)

dev.off()
