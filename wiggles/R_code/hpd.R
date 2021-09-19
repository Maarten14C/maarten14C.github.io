require(IntCal)
cc <- ccurve()

y <- 2400
er <- 50

xmin <- 2100
xmax <- 2800

if(!dir.exists("hpd_pngs"))
  dir.create("hpd_pngs")

prob <- function(y, er, mu, sigma)
  return(dnorm(mu, y, sqrt(er^2 + sigma^2)))

prob.x <- function(x) {
  mu <- approx(cc[,1], cc[,2], x)$y
  sigma <- approx(cc[,1], cc[,3], x)$y
  return(prob(y, er, mu, sigma))
}

steps <- .1
rounding <- 0

xseq <- seq(xmin, xmax, by=steps)
probs <- prob.x(xseq)
sumprobs <- sum(probs)
probs <- probs/sumprobs # sum to 1

o <- order(probs, decreasing=TRUE) # rank to find ages within the prob range
ranked <- cbind(xseq[o], cumsum(probs[o]))

probseq <- c(seq(0.01, .999, by=.003), rep(.68, 100), rep(.95, 100), rep(.99, 100))
probseq <- probseq[order(probseq)]
#probseq <- .5
for(j in 1:length(probseq)) {
  png(paste0("hpd_pngs/img_", 1e5+j, ".png"), units="cm", res=300, height=12, width=15)
  thisprob <- probseq[j]
  ranked <- cbind(ranked[,1:2], ranked[,2] < thisprob) # add a column which shows years within hpd ranges
  ranked <- ranked[order(ranked[,1]),] # order according to calendar age again

  from <- which(diff(ranked[,3]) > 0) # find borders
  to <- which(diff(ranked[,3]) < 0)

  plot(xseq, probs, type="l", xlim=c(xmax, xmin), ylim=c(0, 1.1*max(probs)), xlab="cal BP", ylab="", bty="l", xaxs="i", yaxs="i", col=4)
  min.hpd <- c()
  max.hpd <- c()
  for(i in 1:length(from)) {
    min.hpd[i] <- approx(ranked[(from[i]:(from[i]+1)),2], ranked[((from[i]-1):from[i]),1], thisprob, rule=0)$y
    max.hpd[i] <- approx(ranked[(to[i]:(to[i]+1)),2], ranked[((to[i]-1):to[i]),1], thisprob, rule=0)$y
    hpdseq <- seq(min.hpd[i], max.hpd[i], by=steps)
    hpdprobs <- prob.x(hpdseq)
    pol <- cbind(c(min(hpdseq), hpdseq, max(hpdseq)), c(0, hpdprobs/sumprobs, 0))
    polygon(pol, col=rgb(0,0,0,.2), border=NA)
  }
  line.height <- mean(prob.x(max.hpd[1]))/sumprobs
  abline(h=line.height, lwd=1)
  text(xmin+100, line.height+.00002, labels=paste0(round(100*thisprob, rounding), "%"))
  dev.off()
}

# .mp4 plays well on chrome and MS-Powerpoint:
system("ffmpeg -framerate 15 -pattern_type glob -y -i 'hpd_pngs/*.png' -c:v libx264 hpd.mp4")

# webm plays well on firefox. -pix_fmt set to yuv420p to better translate colours between rgb pngs and yuv webm
system("ffmpeg -framerate 15 -pattern_type glob -y -i 'hpd_pngs/*.png' -b:v 800k -pix_fmt yuv420p hpd.webm")
