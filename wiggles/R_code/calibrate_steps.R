require(IntCal)
cc <- ccurve()

y <- 2400
er <- 50

xmin <- 2100
xmax <- 2800
ymin <- 2000
ymax <- 2700

ex <- 150

# make 3 movies, one with very slow steps (fps), one medium fast, one fast

cc.x <- function(x) {
  mu <- approx(cc[,1], cc[,2], x)$y
  sigma <- approx(cc[,1], cc[,3], x)$y
  return(cbind(mu, sigma))
}

prob <- function(y, er, mu, sigma)
  return(dnorm(mu, y, sqrt(er^2 + sigma^2)) / dnorm(y, y, sqrt(er^2 + sigma^2)))

yseq <- (y-(5*er)):(y+(5*er))
yprob <- dnorm(yseq, y, er)
draw.yprob <- xmax - (ex * yprob/max(yprob))
y.pol <- cbind(c(xmax, draw.yprob, xmax), c(min(yseq), yseq, max(yseq)))
cc.pol <- cbind(c(cc[,1], rev(cc[,1])), c(cc[,2]+cc[,3], rev(cc[,2]-cc[,3])))

xseq <- c(seq(2150, 2300, by=50), seq(2320, 2760, by=20))

draw.skeleton <- function() {
  plot(0, type="n", xlim=c(xmax, xmin), ylim=c(ymin, ymax), xlab="cal BP", ylab="C14 BP", bty="l", xaxs="i", yaxs="i")
  polygon(y.pol, col=rgb(0,0,0,.2), border=rgb(0,0,0,.5))
  polygon(cc.pol, col=rgb(0,.5,0,.5), border=rgb(0,.5,0,.5))
  segments(xmax, y+4*er, xmax-ex, y+4*er) # secondary x axis for C14 date
  segments(xmax-25, ymin, xmax-25, ymin+ex) # secondary x axis for calib date
  segments(xmax-10, y-er, xmax-10, y+er, col=2, lwd=2)
  points(xmax-10, y, col=2, pch=19)
}

done.points <- function(i) {
  if(i > 0) {
    xs <- xseq[1:i]
    probs <- prob(y, er, cc.x(xs)[,1], cc.x(xs)[,2])
    points(xs, ymin + ex * probs, col=4, pch=20, cex=1.5)
    lines(xs, ymin + ex * probs, col=4)
  }
}

if(!dir.exists("calib.pngs"))
  dir.create("calib.pngs")

k <- 0
for(i in 1:length(xseq)) {
  k <- k+1
  #Sys.sleep(1)

  png(paste0("calib.pngs/img_", 1e5+k, ".png"), units="cm", res=300, height=15, width=15)
  draw.skeleton()
  done.points(i-1)
  points(xseq[i], ymin, col=4, pch=3, cex=1.5)
  dev.off()

  k <- k+1
  #Sys.sleep(1)
  png(paste0("calib.pngs/img_", 1e5+k, ".png"), units="cm", res=300, height=15, width=15)
  draw.skeleton()
  done.points(i-1)
  points(xseq[i], ymin, col=4, pch=3, cex=1.5)
  segments(xseq[i], ymin, xseq[i], cc.x(xseq[i]), col=4, lty=2, lwd=1)
  dev.off()

  k <- k+1
  #Sys.sleep(1)
  png(paste0("calib.pngs/img_", 1e5+k, ".png"), units="cm", res=300, height=15, width=15)
  draw.skeleton()
  done.points(i-1)
  points(xseq[i], ymin, col=4, pch=3, cex=1.5)
  segments(xseq[i], ymin, xseq[i], cc.x(xseq[i]), col=4, lty=2, lwd=1)
  points(xseq[i], cc.x(xseq[i])[,1], col=4, pch=20, cex=1.5)
  dev.off()

  k <- k+1
  #Sys.sleep(1)
  png(paste0("calib.pngs/img_", 1e5+k, ".png"), units="cm", res=300, height=15, width=15)
  draw.skeleton()
  done.points(i-1)
  points(xseq[i], ymin, col=4, pch=3, cex=1.5)
  segments(xseq[i], ymin, xseq[i], cc.x(xseq[i]), col=4, lty=2, lwd=1)
  points(xseq[i], cc.x(xseq[i])[,1], col=4, pch=20, cex=1.5)
  segments(xseq[i], cc.x(xseq[i])[,1], xmax, cc.x(xseq[i])[,1], col=4, lty=2, lwd=1)
  dev.off()

  k <- k+1
  #Sys.sleep(1)
  png(paste0("calib.pngs/img_", 1e5+k, ".png"), units="cm", res=300, height=15, width=15)
  draw.skeleton()
  done.points(i-1)
  points(xseq[i], ymin, col=4, pch=3, cex=1.5)
  segments(xseq[i], ymin, xseq[i], cc.x(xseq[i]), col=4, lty=2, lwd=1)
  points(xseq[i], cc.x(xseq[i])[,1], col=4, pch=20, cex=1.5)
  segments(xseq[i], cc.x(xseq[i])[,1], xmax, cc.x(xseq[i])[,1], col=4, lty=2, lwd=1)
  points(xmax - ex * prob(y, er, cc.x(xseq[i])[,1], cc.x(xseq[i])[,2]), cc.x(xseq[i])[,1], col=4, pch=20, cex=1.5)
  dev.off()

  k <- k+1
  #Sys.sleep(1)
  png(paste0("calib.pngs/img_", 1e5+k, ".png"), units="cm", res=300, height=15, width=15)
  draw.skeleton()
  done.points(i-1)
  points(xseq[i], ymin, col=4, pch=3, cex=1.5)
  segments(xseq[i], ymin, xseq[i], cc.x(xseq[i]), col=4, lty=2, lwd=1)
  points(xseq[i], cc.x(xseq[i])[,1], col=4, pch=20, cex=1.5)
  segments(xseq[i], cc.x(xseq[i])[,1], xmax, cc.x(xseq[i])[,1], col=4, lty=2, lwd=1)
  points(xmax - ex * prob(y, er, cc.x(xseq[i])[,1], cc.x(xseq[i])[,2]), cc.x(xseq[i])[,1], col=4, pch=20, cex=1.5)
  points(xseq[i], ymin + ex * prob(y, er, cc.x(xseq[i])[,1], cc.x(xseq[i])[,2]), col=4, pch=20, cex=1.5)
  dev.off()
}

# .mp4 plays well on chrome and MS-Powerpoint:
system("ffmpeg -framerate .8 -pattern_type glob -y -i 'calib_pngs/*.png' -c:v libx264 calibrate_steps.mp4")

# webm plays well on firefox. -pix_fmt set to yuv420p to better translate colours between rgb pngs and yuv webm
system("ffmpeg -framerate .8 -pattern_type glob -y -i 'calib_pngs/*.png' -b:v 800k -pix_fmt yuv420p calibrate_steps.webm")

