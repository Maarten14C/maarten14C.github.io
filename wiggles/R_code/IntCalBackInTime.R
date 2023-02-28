if(! 'rintcal' %in% installed.packages())
	install.packages('rintcal')
require(rintcal)
cc <- ccurve()

res <- 10e3
x <- seq(0, 3100, length=res)
y <- approx(cc[,1], cc[,2], x)$y
er <- .01*y; er[er<10] <- 10

width <- 9

for(i in 1:length(x)) {
  png(paste0("backintime/img_", 1e6+i, ".png"), height=15, width=15, res=300, units="cm")
  
  x.rng <- c(er[i]*width, -1*(er[i]*width))
  y.rng <- 1.2*c(max(er)*width, -1*max(er)*width)
  calibrate(y[i], er[i], cal.lim=x[i]+x.rng, C14.lim=y[i]+x.rng, legend1.loc=NA, legend2.loc=NA)
  segments(x[i], y[i], x[i], -99999, lty=2, col=2)
  segments(100e3, y[i], x[i], y[i], lty=2, col=2)
  dev.off()
}

# this one is going to be quite big, with so many frames and with many changes between frames
# system("ffmpeg -framerate 30 -pattern_type glob -y -i 'backintime/*.png' -c:v libx264 backintime6.mp4")

