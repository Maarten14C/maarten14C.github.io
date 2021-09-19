# uniform distribution

unif <- F
xmin <- 0
xmax <- 55e3
x <- runif(1, 0, 55e3)
move <- 500
n <- 1500
if(!dir.exists("unif_pngs"))
  dir.create("unif_pngs")
energy <- dunif(x, xmin, xmax)
set.seed(10)

if(unif)
for(i in 2:n) {
  png(paste0("unif_pngs/img_", 1e4+i, ".png"), units="cm", res=300, height=10, width=18)
  plot(y <- seq(xmin, xmax, length=200), dunif(y, xmin, xmax), ylim=c(0, 2e-5), col="darkgreen", type="l", yaxs="i", xlab="cal BP", ylab="prob", bty="l")
  prop <- x[i-1] + rnorm(1, 0, move)
  energy[i] <- dunif(prop, xmin, xmax)
  if(energy[i] == 0)
    x[i] <- x[i-1] else
      x[i] <- prop
  rug(x, col=gray(.5))
  tmp <- density(x)
  lines(tmp$x, .05*tmp$y, col=gray(.5))
  segments(x[i], 0, x[i], energy[i], lwd=2, col=2)
  dev.off()
}
cat("\n")
if(unif)
  system("ffmpeg -framerate 30 -pattern_type glob -y -i 'unif_pngs/*.png' -c:v libx264 runif2.mp4")

####

# normal/Gaussian, n MCMC

mn <- 2450
er <- 50
x <- 2600
U <- dnorm(x, mn, er)
energy <- U

move <- 5
n <- 1500
if(!dir.exists("norm_pngs"))
  dir.create("norm_pngs")

set.seed(30)

norm <- F
if(norm)
for(i in 2:n) {
  png(paste0("norm_pngs/img_", 1e4+i, ".png"), units="cm", res=300, height=10, width=18)
  plot(y <- seq(mn-5*er, mn+5*er, length=200), dnorm(y, mn, er), ylim=c(0, .01), col="darkgreen", type="l", yaxs="i", xlab="cal BP", ylab="prob", bty="l")
  prop <- rnorm(1, x[i-1], move)
  energy[i] <- dnorm(prop, mn, er)
  if(energy[i]/energy[i-1] > runif(1))
    x[i] <- prop else {
      x[i] <- x[i-1]
      energy[i] <- energy[i-1]
      }
  rug(x, col=gray(.5))
  tmp <- density(x)
  lines(tmp$x, .005*tmp$y/max(tmp$y), col=gray(.5))
  segments(x[i], 0, x[i], energy[i], lwd=2, col=2)
  dev.off()
}
cat("\n")
if(norm)
  system("ffmpeg -framerate 30 -pattern_type glob -y -i 'norm_pngs/*.png' -c:v libx264   rnorm.mp4")

####

# gamma

mn <- 10
sh <- 1.5
x <- 20
energy <- dgamma(x, sh, sh/mn)

move <- 1
n <- 1500
if(!dir.exists("gamma_pngs"))
  dir.create("gamma_pngs")

set.seed(1)

gamma <- F
if(gamma)
for(i in 2:n) {
  png(paste0("gamma_pngs/img_", 1e4+i, ".png"), units="cm", res=300, height=10, width=18)
  plot(y <- seq(0, 3*mn, length=200), dgamma(y, sh, sh/mn), ylim=c(0, .08), col="darkgreen", type="l", yaxs="i", xlab="yr/cm", ylab="prob", bty="l")
  prop <- max(0, x[i-1]+rnorm(1,0,.25))
  U <- dgamma(prop, sh, sh/mn)
  if(U/energy[i-1] > runif(1)) {
    x[i] <- prop
    energy[i] <- U
  } else {
      x[i] <- x[i-1]
      energy[i] <- energy[i-1]
    }
  rug(x, col=gray(.5))
  tmp <- density(x)
  lines(tmp$x, .04*tmp$y/max(tmp$y), col=gray(.5))
  segments(x[i], 0, x[i], energy[i], lwd=2, col=2)
  dev.off()
}
cat("\n")
if(gamma)
  system("ffmpeg -framerate 30 -pattern_type glob -y -i 'gamma_pngs/*.png' -c:v libx264 rgamma.mp4")


## beta

mn <- .5
s <- 10
x <- .9
energy <- dbeta(x, s*mn, s * (1-mn))

move <- 1
n <- 1500
if(!dir.exists("beta_pngs"))
  dir.create("beta_pngs")

set.seed(1)

beta <- T
if(beta)
for(i in 2:n) {
  png(paste0("beta_pngs/img_", 1e4+i, ".png"), units="cm", res=300, height=10, width=18)
  plot(y <- seq(0, 1, length=200), dbeta(y, s*mn, s * (1-mn)), ylim=c(0, 2.7), col="darkgreen", type="l", yaxs="i", xlab="memory", ylab="prob", bty="l")
  prop <- x[i-1]+rnorm(1,0,.01)
  U <- dbeta(prop, s*mn, s * (1-mn))
  if(U/energy[i-1] > runif(1)) {
    x[i] <- prop
    energy[i] <- U
  } else {
      x[i] <- x[i-1]
      energy[i] <- energy[i-1]
    }
  rug(x, col=gray(.5))
  tmp <- density(x)
  lines(tmp$x, 1.2*tmp$y/max(tmp$y), col=gray(.5))
  segments(x[i], 0, x[i], energy[i], lwd=2, col=2)
  dev.off()
}
cat("\n")
if(beta)
  system("ffmpeg -framerate 30 -pattern_type glob -y -i 'beta_pngs/*.png' -c:v libx264 rbeta.mp4")
