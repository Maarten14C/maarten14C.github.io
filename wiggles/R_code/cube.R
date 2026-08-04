
if(!"plot3D" %in% installed.packages())
  install.packages("plot3D")
if(!"tuneR" %in% installed.packages())
  install.packages("tuneR")

library(plot3D)
library(tuneR)
source("~/Dropbox/ffmpeg.R") # code to produce video and audio files

graphite <- function(n=10000, m=30*60, start=1, lambda=8033, quick=FALSE, folder="block", width=2160, height=3024, res=300, theta=43, theta.move=5, phi=20, phi.move=2, cex.max=.4, cex.diff=.4, cex.min=0.1) {
  message("n: ", n, ", m: ", m)
  c14 <- rexp(n, 1/lambda) # C14 particles and their lifespans
  ord <- order(c14)
  x <- runif(n); y <- runif(n); z <- runif(n) # their locations within the block
  t <- seq(0, 60e3, length=m) # time
 
  # make the cube move around smoothly in a sleeping-8 shape
  u <- seq(0, 2*pi, length.out=m)
  theta <- theta + theta.move*sin(u)
  phi <- phi + phi.move*sin(2*u + pi/4)

  # work on clean pngs
  if(!dir.exists(folder))
    dir.create(folder)
  file.remove(list.files(folder, pattern="\\.png$", full.names=TRUE))

  yellows <- colorRampPalette(c("yellow", "greenyellow"))(100)
   
  for(i in start:m) {
    mm <- 10^ceiling(log10(m)) # move value up to next order of magnitude
    png(file.path(folder, paste0("img_", mm+i, ".png")), 
      height=height, width=width, res=res, , type="cairo", units="px")
    par(fig=c(0,1,.3,1), mar=c(0,1,1,1), new=FALSE)

    alive <- which(c14 > t[i])
    pop <- c()
    if(i > 1)
      pop <- which(c14 <= t[i] & c14 > t[i-1])

    plist <- scatter3D(x, y, z, col=NA, cex=0,
      theta=theta[i], phi=phi[i]) # gives the rotation

    # 'depths' of points within the cube, from viewer's point
    depth <- (1-x)^2 + y^2 # distance from viewer
    cex <- cex.max - pmax(cex.min, cex.diff * depth) # icon size smaller when further away, down to a limit

    points3D(x[alive], y[alive], z[alive],
      pch=19, cex=cex[alive],
      xlim=c(0,1), ylim=c(0,1), zlim=c(0,1),
      col=yellows, colvar=depth[alive],
      theta=theta[i], phi=phi[i], axes=FALSE, colkey=FALSE, bty="bl", edge.col=NA)

    # C14 atoms that decayed in this step
    if(length(pop) > 0)
      points3D(x[pop], y[pop], z[pop],
        pch=19, cex=4*cex[pop], add=TRUE,
        xlim=c(0,1), ylim=c(0,1), zlim=c(0,1), col="white",
        theta=theta[i], phi=phi[i], axes=FALSE, colkey=FALSE, bty="bl", edge.col=NA)

    # limbs
    segments3D(0,0,0, 1,0,0, col="black", add=TRUE) # bottom, front left
    segments3D(1,0,0, 1,1,0, col="black", add=TRUE) # bottom, front right
    segments3D(0,0,0, 0,1,0, col="grey15", add=TRUE) # bottom, back left
    segments3D(0,1,0, 1,1,0, col="grey15", add=TRUE) # bottom, back right
    segments3D(0,0,1, 0,1,1, col="black", add=TRUE) # top, back left
    segments3D(0,1,1, 1,1,1, col="black", add=TRUE) # top, back right
    segments3D(0,0,0, 0,0,1, col="black", add=TRUE) # vertical, left
    segments3D(1,1,0, 1,1,1, col="black", add=TRUE) # vertical, right
    segments3D(0,1,0, 0,1,1, col="grey15", add=TRUE) # vertical, back
    segments3D(0,0,1, 1,0,1, col="grey70", add=TRUE) # top, front left
    segments3D(1,0,1, 1,1,1, col="grey70", add=TRUE) # top, front right
    segments3D(1,0,0, 1,0,1, col="grey70", add=TRUE) # vertical, front

    # timeline
    par(fig=c(0,1,0,.3), mar=c(4,4,2,2), mgp=c(2, .7, 0), new=TRUE)
    survivors <- sapply(t[1:i], function(tt) sum(c14 > tt))
    plot(t[1:i], survivors, type="l", lwd=2, xlim=range(t), ylim=c(0, n), 
      xlab=expression(""^14*C~BP), ylab=expression(""^14*C~particles),
      bty="l", xaxt="n", xaxs="i", yaxt="n", yaxs="i")
    atx <- seq(0, 60e3, by=10e3)
    axis(1, at=atx, labels=format(atx, big.mark=",", scientific=FALSE))
    aty <- pretty(0:n, 3)
    axis(2, at=aty, labels=format(aty, big.mark=",", scientific=FALSE))
    points(t[i], survivors[i], pch=19)
    text(t[i], survivors[i]+500, survivors[i], pos=4, adj=c(0,0), cex=0.7)

    dev.off()
  }
  return(list(c14=c14, x=x, y=y, z=z, t=t))
}



pops <- function(m, fps=30, sr=44100, name="c14_1m.wav") {
  audio <- rep(0, ceiling(m / fps * sr))
  click <- sin(seq(0, 10*pi, length.out=200)) * exp(-seq(0, 5, length.out=200))

  for(i in 2:m) {
    pop <- which(c14$c14 <= c14$t[i] & c14$c14 > c14$t[i-1])
    if(length(pop) == 0)
      next

    frame.time <- (i-1) / fps
    start <- round(frame.time * sr)

    for(j in seq_along(pop)) {
      # spread clicks slightly within frame
      offset <- sample(0:(sr/fps - 1), 1)
      pos <- start + offset

      if(pos + length(click) - 1 <= length(audio))
        audio[pos:(pos + length(click) - 1)] <-
          audio[pos:(pos + length(click) - 1)] + click
    }
  }

  # normalise
  audio <- audio / max(abs(audio))

  wave <- Wave(left=round(audio * 32767), samp.rate=sr, bit=16)
  writeWave(wave, name)
}



set.seed(42)

m <- 30 * 60
n <- 1e4
c14 <- graphite(n, m)
ffmpeg("block", "block_1m.mp4", 30, 0, 0)

pops(m, name="block/c14_1m.wav")

# combine the audio and video:
combine.av("block", "block_1m.mp4", "c14_1m.wav", "graphite_audio_1m.mp4")
