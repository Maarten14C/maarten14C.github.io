require(coffee)

its <- 2e3

# a very short run, no thinning, with a clear burnin process
set.seed(20);sim.strat();strat(burnin=0, thinning=1, internal.thinning=1, its=its)

cc <- ccurve()

if(!dir.exists("strat_pngs"))
  dir.create("strat_pngs")	

pb <- txtProgressBar(1, round(its/100), style=3)

ex <- 150

for(j in 1:its) {
  if(j %% 100 == 0)
    setTxtProgressBar(pb, j)

  png(paste0("strat_pngs/img_", 1e4+j, ".png"), units="cm", res=300, height=10, width=15)
  layout(matrix(1:2, nrow=2), heights=c(3, 7))
  par(mar=c(3,3,1,1), mgp=c(1.7,.7,0), bty="l")
  plot(info$Us[1:j], type="l", xlim=c(0,info$Tr), ylim=rev(range(info$Us)), xlab="MCMC its", ylab="energy")

  plot(0, type="n", xlim=c(5600, 4e3), ylim=c(5,0), xlab="cal BP", ylab="position")
  #draw.dates(info$dets[,2], info$dets[,3], info$dets[,4], mirror=FALSE, hpd.col=NA, y.lab="position", col=rgb(0,0,1,.05), threshold=.02, y.lim=c(5,0), normalise=TRUE)
  for(i in 1:5) {
    tmp <- info$output[1:j,i] # the modelled years of interest
    segments(tmp, i, tmp, i-.08, col=gray(.5))
    if(length(unique(tmp)) > 5) {
      dns <- density(tmp)
      dns <- cbind(c(min(dns$x), dns$x, max(dns$x)), i - .5*c(0, dns$y/max(dns$y), 0))
      polygon(dns, border=gray(.5), col=NA)
      #lines(dns$x, i-.2*dns$y/max(dns$y), col=gray(.5))
      #segments(min(dns$x), i, max(dns$x), i, col=gray(.5))
    }

    # draw calib probs - slower than using draw.dates, but the heights coincide better with the subsequent drawing of the red segments. Interpolate to every 5 yr, not every 1 yr
    cc.prob <- dnorm(cc[,2], info$dets[i,2], sqrt(info$dets[i,3]^2 + cc[,3]^2))
    calib <- cbind(cc[,1], cc.prob)
    calib <- calib[which(calib[,2] > .0001),]
    calib <- approx(calib[,1], calib[,2], seq(min(calib[,1]), max(calib[,1]), by=5))
    calib <- cbind(calib$x, calib$y)
    pol <- cbind(c(min(calib[,1]), calib[,1], max(calib[,1])), i - ex*c(0,calib[,2],0))
    polygon(pol, col=rgb(0,0,1,.1), border=rgb(0,0,1,.5))

    # find calibrated probabilities of the modelled years
    cc.y <- approx(cc[,1], cc[,2], tmp[length(tmp)])$y
    cc.er <- approx(cc[,1], cc[,3], tmp[length(tmp)])$y
    prob <- dnorm(cc.y, info$dets[i,2], sqrt(info$dets[i,3]^2 + cc.er^2))
    segments(tmp[length(tmp)], i, tmp[length(tmp)], i-(ex*prob), col=2, lwd=2)
  }
  dev.off()
}

system("ffmpeg -framerate 30 -pattern_type glob -y -i 'strat_pngs/*.png' -c:v libx264 strat3.mp4")
