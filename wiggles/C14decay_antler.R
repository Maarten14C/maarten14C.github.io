mypackages <- installed.packages()
if(!"png" %in% mypackages)
  install.packages("png")
if(!"tuneR" %in% mypackages)
  install.packages("tuneR")
if(!"terra" %in% mypackages)
  install.packages("terra")

source("~/Dropbox/ffmpeg.R")

#img <- readPNG("~/Dropbox/animations/antler_nobg.png") # b& w version to define a mask
img <- png::readPNG("~/Dropbox/animations/antler_thorsten-nobg.png") # low-res version to speed things up
mask <- img[, , 4] > 0.5 # define a mask - the shape of the antler
#mask <- mask[1:nrow(mask),]
mask_rast <- terra::rast(mask)

# make background NA
mask_dist <- mask_rast
values(mask_dist)[values(mask_dist)] <- NA # set all 'TRUE' cells to NA
d <- as.matrix(distance(mask_dist), wide=TRUE) # calculate distance from antler's edge
inside.pix <- which(mask, arr.ind=TRUE)
img_h <- nrow(img)
img_w <- ncol(img)
aspect <- img_h / img_w

# variables
set.seed(42)
n <- 1e4 # initial number of C14 atoms
border <- .6 # larger values cause more empty border regions
xmax <- 60e3 # maximum C14 age
right <- xmax # where to plot the righthand limit of the antler
folder <- "~/antler" # where to put the pngs and mp4
m <- 30*60 # number of frames
width <- 1800 # width of the pngs
height <- 2/3 * width
res <- 200 # png resolution. This affects e.g. font size
rgb.c14 <- c(.9, .05, .05) # colours of the C14 dots. Bright red

# antler dimensions
bottom <- n/2 # bottom part of the antler image
top <- n # top part
left <- xmax/7 # lefthand limit of the antler

# C14 particles and their lifespans
c14 <- rexp(n, 1/8033) 
t <- seq(0, 60e3, length=m) # time

# don't want select dots too close to the border/edge of the antler
w <- d[cbind(inside.pix[,1], inside.pix[,2])]
w <- w^border
sel <- inside.pix[sample(nrow(inside.pix), n, replace=FALSE, prob=w),]

x <- (sel[,2] - runif(n)) / ncol(mask) # add some jitter
y <- 1 - (sel[,1] - runif(n)) / nrow(mask)

# for plotting
dot.alpha <- runif(n, .4, 1) # transparency of the dots
#yellowgreens <- rgb(0, .45, .45, alpha=dot.alpha) # shades of yellowgreen
carboncolours <- rgb(rgb.c14[1], rgb.c14[2], rgb.c14[3], alpha=dot.alpha) # blood red
cexs <- 1.5*runif(n, .05, .12) # not all dots are equal (in size)

antler <- function(m=1e3, Folder=folder, as.png=TRUE) {
  # work on clean pngs
  if(!dir.exists(Folder))
    dir.create(Folder)
  file.remove(list.files(Folder, pattern="\\.png$", full.names=TRUE))
  
  mm <- 10^ceiling(log10(m)) # move starting value in file name up to next order of magnitude
  
  for(i in 1:m) {
    if(as.png)
      png(file.path(folder, paste0("img_", mm+i, ".png")), 
        height=height, width=width, res=res, type="cairo", units="px")
    par(mar=c(5,5,2,2))
    plot(0, type="n", xlim=c(0, xmax), ylim=c(0, n),
      bty="n", xaxt="n", yaxt="n", xaxs="i", yaxs="i",
      xlab=expression(""^14*C~BP), ylab=expression(""^14*C~particles))
    atx <- seq(0, xmax, by=10e3)
    axis(1, at=atx, labels=format(atx, big.mark=",", scientific=FALSE))
    aty <- pretty(0:n, 3)
    axis(2, at=aty, labels=format(aty, big.mark=",", scientific=FALSE))

    # calculate dimensions of the to-be-plotted antler png
    usr <- par("usr")
    pin <- par("pin")
    xscale <- pin[1] / diff(usr[1:2])
    yscale <- pin[2] / diff(usr[3:4])
    img_width_plot <- right - left
    img_height_plot <- img_width_plot * (img_h/img_w) * (xscale/yscale)
    bottom <- top - img_height_plot

    rasterImage(img, left, bottom, right, top) # plot the antler

    alive <- which(c14 > t[i]) # which C14 atoms have survived until at least this step
    pop <- c() # which ones decay within this step's time window
    if(i > 1)
      pop <- which(c14 <= t[i] & c14 > t[i-1])
    survivors <- sapply(t[1:i], function(tt) sum(c14 > tt))
  
    points(left+x[alive]*(right-left), bottom+y[alive]*(top-bottom), 
      cex=cexs[alive], pch=19, col=carboncolours[alive])

    if(length(pop) > 0) { # make the white spots look like explosions
      points(left+x[pop]*(right-left), bottom+y[pop]*(top-bottom),
        cex=5*cexs[pop], pch=19, col=rgb(1,1,1,.25))
      points(left+x[pop]*(right-left), bottom+y[pop]*(top-bottom),
        cex=4*cexs[pop], pch=19, col=rgb(1,1,1,.25))
      points(left+x[pop]*(right-left), bottom+y[pop]*(top-bottom),
        cex=3*cexs[pop], pch=19, col="white")
    }

    lines(t[1:i], survivors, lwd=2) # draw the exponential decay line
    points(t[i], survivors[i], pch=19) # where are we now
    text(t[i], survivors[i]+500, survivors[i], pos=4, adj=c(0,0), cex=0.7) # how many

    if(as.png)
      dev.off()
  }
    
  return(list(c14=c14, x=x, t=t))
}


pops <- function(m, fps=30, sr=44100, name="antler_c14_1m.wav") {
  n.samples <- ceiling(m / fps * sr)

  click <- sin(seq(0, 10*pi, length.out=200)) * exp(-seq(0, 5, length.out=200))

  left  <- rep(0, n.samples)
  right <- rep(0, n.samples)
  for(i in 2:m) {
    pop <- which(C14$c14 <= C14$t[i] & C14$c14 >  C14$t[i-1])

    if(length(pop) == 0)
      next

    frame.time <- (i - 1) / fps
    start <- round(frame.time * sr)

    for(j in pop) {
      offset <- sample(0:(sr / fps - 1), 1)
      pos <- start + offset
      if(pos + length(click) - 1 > n.samples)
        next

      # x controls stereo position
      xrank <- rank(x) / length(x)
      gain.left  <- sqrt(1 - xrank[j])
      gain.right <- sqrt(xrank[j])

      idx <- pos:(pos + length(click) - 1)

      left[idx]  <- left[idx] + gain.left * click
      right[idx] <- right[idx] + gain.right * click
    }
  }

  fade.sec <- 1 # fade in
  fade.n <- fade.sec * sr
  ramp <- seq(0, 1, length.out=fade.n)
  left[1:fade.n]  <- left[1:fade.n]  * ramp
  right[1:fade.n] <- right[1:fade.n] * ramp

  mx <- max(abs(c(left, right)))
  if(mx > 0) {
    left  <- left/mx
    right <- right/mx
  }

  wave <- tuneR::Wave(left=round(left * 32767), right=round(right * 32767), samp.rate=sr, bit=16)
  tuneR::writeWave(wave, name)
}




intro <- function(frames=30*3, Folder="antler_intro", as.png=TRUE) {
  # work on clean pngs
  if(!dir.exists(Folder))
    dir.create(Folder)
  file.remove(list.files(Folder, pattern="\\.png$", full.names=TRUE))

  mm <- 10^ceiling(log10(frames)) # move starting value in file name up to next order of magnitude

  axis.cols <- rgb(0, 0, 0, (1:frames)/frames)

  for(i in 1:frames) {
    if(as.png)
      png(file.path(Folder, paste0("img_", mm+i, ".png")),
        height=height, width=width, res=res, type="cairo", units="px")
    par(mar=c(5,5,2,2))
    plot(0, type="n", xlim=c(0, xmax), ylim=c(0, n), bty="n", xaxt="n", yaxt="n", xaxs="i", yaxs="i",
      xlab=expression(""^14*C~BP), ylab=expression(""^14*C~particles), col.lab=axis.cols[i])
    atx <- seq(0, xmax, by=10e3)
    axis(1, at=atx, labels=format(atx, big.mark=",", scientific=FALSE), col=axis.cols[i], col.axis=axis.cols[i])
    aty <- pretty(0:n, 3)
    axis(2, at=aty, labels=format(aty, big.mark=",", scientific=FALSE), col=axis.cols[i], col.axis=axis.cols[i])

    # calculate dimensions of the to-be-plotted antler png
    usr <- par("usr")
    pin <- par("pin")
    xscale <- pin[1] / diff(usr[1:2])
    yscale <- pin[2] / diff(usr[3:4])
    img_width_plot <- right - left
    img_height_plot <- img_width_plot * (img_h/img_w) * (xscale/yscale)
    bottom <- top - img_height_plot

    rasterImage(img, left, bottom, right, top)

    a <- i/frames
    cols <- rgb(rgb.c14[1], rgb.c14[2], rgb.c14[3], dot.alpha*a)
  #   cols <- rgb(.9, .05, .05, dot.alpha*a)
    points(left+x*(right-left), bottom+y*(top-bottom),
      cex=cexs*(i/frames), pch=19, col=cols)

    if(as.png)
      dev.off()
  }
}

C14 <- antler(m=m, as.png=TRUE)

ffmpeg("antler", "antler_1m.mp4", 30, 0, 0)

intro()
ffmpeg("antler_intro", "antler_intro.mp4", 30, 0, 0)
system("ffmpeg -y -f lavfi -i anullsrc=r=44100:cl=stereo -i antler_intro/antler_intro.mp4 -c:v copy -c:a aac -shortest   antler_intro/antler_intro_audio.mp4") # add silent audio to the intro

# create the audio for the main part
pops(m, name="antler/antler_1m_stereo.wav")

# combine the audio and video for the main part:
combine.av("antler", "antler_1m.mp4", "antler_1m_stereo.wav", "antler_audio_1m_stereo.mp4")

# now glue the intro to the rest: (needs a file listing the names of the mp4 files to join)
system("ffmpeg -y -f concat -safe 0 -i antler/files.txt -c copy antler/antler_complete.mp4")
