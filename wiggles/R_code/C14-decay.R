C14n <- 1e4	# initial amount of 14C atoms
steps <- 50	# amount of years per calculation step
MaxAge <- 6e4	# maximum calculation age
thalf <- 5568	# half-life 14C

set.seed(90)
x <- runif(C14n, 0, 1)	# create random x-coors for atoms
y <- runif(C14n, 0, 1) # create random y-coors for atoms
survive <- rexp(C14n, log(2)/thalf) # exponential decay
o <- order(survive) 
survive <- survive[o] # list of created decay events
C14 <- cbind(survive, x, y) # location and timing of decays

ages <- seq(0, MaxAge, by=steps) # time windows

for(i in 1:length(ages))
 {
  png(paste0("decay_pngs/img_", 1e5+i, ".png"), units="cm", res=300, height=15, width=15)

  layout(matrix(c(1,2),ncol=1), heights=c(.6,.4))
  par(mar=c(0,1,.1,.1), mgp=c(2,1,0), xaxt="n", yaxt="n", bty="n")
  plot(0, type="n", xlab="", ylab="", main="", xlim=c(0,1), ylim=c(0,1))
  rect(0,0,1,1, col="black", border=NA)
  temp <- C14[min(which(C14[,1]>=ages[i])):nrow(C14),2:3]
  points(temp, pch=20, col="yellow", cex=.15)

  par(xaxt="s", yaxt="n", mar=c(3.2,1,0,.5))
  plot(ages,dexp(ages,log(2)/thalf), xlim=c(0,MaxAge),
   xlab="C-14 age", ylab="", main="", type="l", col="red", yaxs="i")
  lines(c(ages[i],ages[i]), c(0,dexp(ages[i],log(2)/thalf)), col="green", lwd=3)
  points(ages[i],0, pch=19, col="green")
  if(length(temp) > 0)
    legend("bottomleft", legend=nrow(temp), bty="n")
  if(round(i/10) == round(i)/10)
    cat(i, "")
  dev.off()
 }

# for an animated gif, using ImageMagick. Sometimes doesn't work, and the resulting animation doesn't have a very high colour resolution. The second command enhances the animated gif. Both ImageMagick and gifsicle will have to be installed.
#system("convert -loop 0 -delay 15 *.pdf temp.gif")
#system("gifsicle --optimize --colors 20 < temp.gif > C14-decay.gif")

# .mp4 plays well on chrome and MS-Powerpoint:
system("ffmpeg -framerate 20 -pattern_type glob -y -i 'decay_pngs/*.png' -c:v libx264 decay.mp4")

# webm plays well on firefox:
system("ffmpeg -framerate 20 -pattern_type glob -y -i 'decay_pngs/*.png' -c:v libvpx-vp9 -deadline best decay.webm")

