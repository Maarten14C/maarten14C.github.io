### finish manual!!!

## hiatus and prior txt file full of problems. Perhaps best:
# 1) provide alternative settings as option in Bacon(), do not attempt to adapt _priors.txt.
# 2) get rid of difficult first line _priors.txt and instead introduce hiatus.depths as extra line.
# 3) file.remove("Cores/...") if prior file messes things up
# but now old non-default prior settings are kept in subsequent runs...
# or smthng like if(length(acc.shape)+length(acc.mean)+length(mem.strength)... > 0)... within Bacon.priors() to warn user about 1) non-existing prior file, 2) reading existing prior file with different settings than default, 3) reading existing prior file with no changed settings, or 4) you have provided prior settings different than those in the existing prior file. With details as to the differences?

# temporarily (?) added max.y option to deal with invalid ylim of some acc.mem plots

# produce an automatic _settings.txt file as in clam, with cal.curve(s) used etc.

# Bacon.hist should also give hpd ranges

# check if numbers in dets file are numeric (e.g. 1 could be | or l, 0 could be o)

## line 25 assumes at least 1 cc!=0...

## adding a section break and distinct acc.rate priors doesn't work, second values get clipped in priors file

## acc.mean should always be != acc.shape, otherwise error. Implement this, e.g. with warning

## put legend to greyscales in relevant graphs (as option)

## if hiatus, plot acc.posts of the individual sections

### allow for hiatuses in posterior Bacon model plots
### allow for BC/AD?

### All identifiers will start with Bacon
### some handy graph settings. bty draws (no) boxes outside plots, mar&mgp do spacing between axes&margins 
### par(bty="n", mar=c(4,3,1,1), mgp=c(1.5, .7, 0))

### for postscript avoid semi-transparency and do: greyscale=function(gr) grey(seq(1-gr,1, .2))

### load data, run bacon and analyse outputs
Bacon <- function(core="MSB2K", res=5, unit="cm", runname="", MinYr=-5, MaxYr=50000, d.min=c(), d.max=c(), ssize=1000, normal=FALSE, calcurve1="IntCal09", calcurve2="Marine09", calcurve3="SHCal04", calcurve4="ConstCal", deltaR=0, deltaSTD=0, t.a=3, t.b=4, cc=0, cutoff=.001, th0=c(), exx=15, default.acc=c(2, 10), acc.shape=c(), acc.mean=c(), default.mem=c(4, 0.7), mem.strength=c(), mem.mean=c(), default.hiatus=c(1, 100), hiatus.depths=c(), hiatus.shape=c(), hiatus.mean=c(), bins=50, dark=5, yrs=10, prob=.95, rounded=1, rev.axis=FALSE, store=TRUE, ask=TRUE, run=TRUE, postbomb=0, draw.res=200, yr.res=10, calc.every=1, burnin=min(200, ssize), plot.pdf=TRUE, greyscale=function(x) rgb(0,0,0,x), postscript=FALSE, C14.col=rgb(0,0,1,.35), C14.border=rgb(0,0,1,.5), cal.col=rgb(0,.5,.5,.35), cal.border=rgb(0,.5,.5,.5), normalise.dists=TRUE)
  {
    ### define file and calibration curve(s)
    det.file <- paste("Cores/", core, "/", core, ".dat", sep="")
    dets <- read.table(det.file, header=TRUE)
    if(length(th0)==0)
      th0 <- round(rnorm(2, dets[1,2], dets[1,3]))
    if(postbomb == 0 && ((ncol(dets)==4 && min(dets[,2]) < 0) || 
      ncol(dets) == 9 && min(dets[dets[,9] > 0,2]) < 0)) ## assumes at least 1 cc!=0...
	stop("\nWarning, you have negative C14 ages so should select a postbomb curve")
    if(postbomb!=0) if(postbomb=="NH1") postbomb <- 1 else if(postbomb=="NH2") postbomb <- 2 else 
      if(postbomb=="NH3") postbomb <- 3 else if(postbomb=="SH") postbomb <- 4

    ### read priors file, produce default one if it does not exist
    priors <- Bacon.priors(core, res, dets, d.min, d.max, default.acc, acc.shape, acc.mean, default.mem, mem.strength, mem.mean, default.hiatus, hiatus.depths, hiatus.shape, hiatus.mean)

    ### adapt depth range (and if needed remove or add dummy dets) if told to extrapolate
    if(length(d.min) > 0 || length(d.max) > 0) 
      if(d.min < priors$d.min || d.max > priors$d.max)
	cat("\n  Extrapolating beyond the dated levels\n")
    d <- seq(floor(priors$d.min), ceiling(priors$d.max), by=res) 
    if(max(d) < priors$d.max) d <- c(d, max(d)+res)
    if(min(d) < min(dets[,4]))
      {
        dets <- rbind(dets[1,], dets)
        dets[1,3:4] <- c(1e4, min(d))
      } 
    if(max(d) > max(dets[,4]))
      {
        dets <- rbind(dets, dets[nrow(dets),])
        dets[nrow(dets),3:4] <- c(1e4, max(d))
      }

    ### clean up old unfinished runs, produce files
    old.runs <- list.files(paste("Cores/", core, sep=""), pattern=".out.all")
    if(length(old.runs) > 0)
      for(i in old.runs) file.remove(paste("Cores/", core, "/", i, sep=""))
    prefix <- paste("Cores/", core, "/", core, runname, "_", length(d), sep="")
    bacon.file <- file(paste(prefix, ".bacon", sep=""), "w")
    if(!file.exists(outfile <- paste(prefix, ".out", sep="")))
      file.create(outfile)
      
    ### store values for future manipulations
    info <- list(core=core, det.file=det.file, prefix=prefix, runname=runname, dets=dets, d=d, K=length(d), Dc=res, bacon.file=bacon.file, ssize=ssize, normal=normal, calcurves=c(calcurve1, calcurve2, calcurve3, calcurve4), deltaR=deltaR, deltaSTD=deltaSTD, t.a=t.a, t.b=t.b, cc=1, acc.shape=priors$acc.shape, acc.mean=priors$acc.mean, mem.mean=priors$mem.mean, mem.strength=priors$mem.strength, hiatus.shape=priors$hiatus.shape, hiatus.mean=priors$hiatus.mean, hiatus.depths=priors$hiatus.depths, yrs=yrs, bins=bins, dark=dark, prob=prob, rounded=rounded, res=res, unit=unit, date=date(), th0=th0, MinYr=MinYr, MaxYr=MaxYr, postbomb=postbomb, calc.every=calc.every)
    if(store) info <<- info

    ### run bacon if initial graphs seem OK
    write.Bacon.file(info)
    question <- substr("",1,1)
    if(!run) ask <- FALSE
    if(ask)
      {
        ### plot initial data and priors
        layout(matrix(if(length(priors$hiatus.depths) > 0) c(1,2,3,4,4,4) else c(1,2,3,3),
          nrow=2, byrow=TRUE), heights=c(.3,.7))
        par(mar=c(3,3,1,1), mgp=c(1.5,.7,.0), bty="l", yaxs="i")
        PlotAccPrior(priors$acc.shape, priors$acc.mean)
        PlotMemPrior(priors$mem.strength, priors$mem.mean, res)
        if(length(priors$hiatus.depths) > 0)
          PlotHiatusPrior(priors$hiatus.shape, priors$hiatus.mean, priors$hiatus.depths)
        calib.plot(read.table(det.file, header=TRUE), info=info, yr.res=yr.res, rev.axis=rev.axis, d.lim=range(d), exx=exx, cutoff=cutoff, normal=normal, normalise.dists=normalise.dists, t.a=t.a, t.b=t.b, deltaR=deltaR, deltaSTD=deltaSTD, C14.col=C14.col, C14.border=C14.border, cal.col=cal.col, cal.border=cal.border)
        legend("topleft", core, bty="n", cex=1.5)
        question <- readline(cat("  Run", core, "with", length(d), "sections? (y/n) "))
      }
      if((ask==FALSE && run==TRUE) || substr(question,1,1)=="y")
        {
          system(paste("bin/bacon ", prefix, ".bacon ", outfile, " ", ssize, sep=""))
		  scissors(burnin, info)
          Bacon.PlotAgeModels(info, dets=dets, dark=dark, bins=bins, draw.res=draw.res, yr.res=yr.res, prob=prob, rounded=rounded, rev.axis=rev.axis, store=store, exx=exx, cutoff=cutoff, normalise.dists=normalise.dists, greyscale=greyscale, C14.col=C14.col, C14.border=C14.border, cal.col=cal.col, cal.border=cal.border)
		  if(plot.pdf)
	        dev.copy2pdf(file=paste(prefix, ".pdf", sep=""))
        } else
	      if(ask==TRUE) cat("  OK. Please adapt settings.\n\n")
    closeAllConnections()
  }



### set the parameters for the prior distributions
Bacon.priors <- function(core, res, dets, d.min, d.max, default.acc, acc.shape, acc.mean, default.mem, mem.strength, mem.mean, default.hiatus, hiatus.depths, hiatus.shape, hiatus.mean) 
  {
    Priors <- paste("Cores/", core, "/", core, "_priors.txt", sep="")
    if(file.exists(Priors) && file.info(Priors)$size > 0)
      {
	priors <- read.csv(Priors, sep="", header=FALSE)[,-1] # read priors file but not the row headers
	
	d <- unique(as.numeric(priors[1,]))
	pr <- list(d.min=min(d), d.max=max(d))
	if(length(d) > 2) # any defined elbows/sections?
	  pr$hiatus.depths <- sort(d[-c(min(d), max(d))], decreasing=TRUE)

	pr$acc.shape <- unique(as.numeric(priors[2,!is.na(priors[2,])]))
	pr$acc.mean <- unique(as.numeric(priors[3,!is.na(priors[3,])]))
	pr$mem.strength <- unique(as.numeric(priors[4,!is.na(priors[4,])]))
	pr$mem.mean <- unique(as.numeric(priors[5,!is.na(priors[5,])]))
	pr$hiatus.shape <- default.hiatus[1]
	pr$hiatus.mean <- default.hiatus[2]

	if(length(pr$hiatus.depths) > 0)
	  if(nrow(priors) > 5) #  && ncol(priors[(nrow(dets)-1):nrow(dets),]) > 1) ### does the file give hiatus information?
	    {
	      pr$hiatus.shape <- as.numeric(priors[6,])[-1]
	      pr$hiatus.mean <- as.numeric(priors[7,])[-1]
	    }
      } else
      {
	cat("\n  No priors file provided, I will write one.\n\n")
	pr <- list(d.min=min(dets[,4]), d.max=max(dets[,4]), acc.shape=default.acc[1], acc.mean=default.acc[2], mem.strength=default.mem[1], mem.mean=default.mem[2], hiatus.depths=hiatus.depths, hiatus.shape=default.hiatus[1], hiatus.mean=default.hiatus[2])	
      }

    ### update priors if any were provided as optional constants in Bacon(...) (i.e. changed from the defaults)
    if(length(d.min) > 0) pr$d.min <- d.min # will still use all dates if dmin > min(dets[,4])
    if(length(d.max) > 0) pr$d.max <- d.max # same for dmax < max(dets[,4])
    if(length(acc.shape) > 0) pr$acc.shape <- acc.shape
    if(length(acc.mean) > 0) pr$acc.mean <- acc.mean
    if(length(mem.strength) > 0) pr$mem.strength <- mem.strength
    if(length(mem.mean) > 0) pr$mem.mean <- mem.mean
    if(length(hiatus.shape) > 0) pr$hiatus.shape <- hiatus.shape
    if(length(hiatus.mean) > 0) pr$hiatus.mean <- hiatus.mean
    if(length(hiatus.depths) > 0) pr$hiatus.depths <- hiatus.depths

    ### update the priors file
    blah <- c("depths\t\t", pr$d.min, "\t", pr$hiatus.depths, "\t", pr$d.max, "\nacc.shape\t", unique(pr$acc.shape), "\nacc.mean\t", unique(pr$acc.mean), "\nmem.strength\t", unique(pr$mem.strength), "\nmem.mean\t", unique(pr$mem.mean))
    if(length(pr$hiatus.depths) > 0)
      blah <- c(blah, "\nhiatus.shape\t", unique(pr$hiatus.shape), "\nhiatus.mean\t", unique(pr$hiatus.mean))
    cat(blah, "\n", file=Priors)    
    closeAllConnections()
    pr
  }


### calcurve0=constcal, calcurve1="IntCal09", calcurve2="Marine09", calcurve3="SHCal04", calcurve4="ConstCal" 
write.Bacon.file <- function(info, bacon.file=info$bacon.file, dets=info$dets, postbomb=info$postbomb)
  {
    cat("## Ran on", info$date, "\n\n", file=bacon.file)
    cat("Cal 0 : ConstCal;\nCal 1 : ", info$calcurves[1], ", ", postbomb, 
      ";\nCal 2 : ", if(info$calcurves[2]=="Marine09") "IntCal09Marine" else info$calcurves[2], ## perhaps better rename to Marine09 inside cal.h and input.c
      ";\nCal 3 : ", info$calcurves[3], ";",
      if(info$calcurves[4]!="ConstCal")
	paste("\nCal 4 : GenericCal, Curves/", info$calcurves[4], ";", sep=""), 
      sep="", file=bacon.file)
    cat("\n\n##   id.   yr    std   depth  deltaR  deltaSTD     t.a   t.b   cc", file=bacon.file)
    if(ncol(dets) == 4)
      {
        cat("\nDet 0 : ", as.character(dets[1,1]), " ,  ", dets[1,2], ",  ",
          dets[1,3], ",  ", dets[1,4], ",  ", info$deltaR, ",  ", info$deltaSTD,
          ",  ", info$t.a, ",  ", info$t.b, ",  ", info$cc, ";", sep="", file=bacon.file)
        if(nrow(dets)>1)
          for(i in 2:nrow(dets))
            cat("\nDet ", i-1, " : ",  as.character(dets[i,1]),
            " , ", dets[i,2], ", ", dets[i,3], ", ", dets[i,4],
            ";", sep="", file=bacon.file)
      } else
    for(i in 1:nrow(dets))
        cat("\nDet ",i-1, " : ",  as.character(dets[i,1]), " ,  ",
          dets[i,2], ",  ", dets[i,3], ",  ", dets[i,4], ",  ",
          dets[i,5], ",  ", dets[i,6], ",  ", dets[i,7], ",  ",
          dets[i,8], ",  ", dets[i,9], ";", sep="", file=bacon.file)

    if(length(info$hiatus.depths) > 0) ### hiatus(es) inferred
      {
        cat("\n  Hiatus set at depth(s)", info$hiatus.depths, "\n")
        cat("\n\n### Depths and priors for fixed hiatuses, in descending order",
          "\n##### cm  alpha beta      ha     hb", file=bacon.file)
        for(i in 1:length(info$hiatus.depths))
          cat("\nHiatus ", i-1, ":  ", info$hiatus.depth[i], ",  ", info$acc.shape[i],
            ",  ", info$acc.shape[i]/info$acc.mean[i], ",  ", info$hiatus.shape[i],
            ",  ", info$hiatus.shape[i]/info$hiatus.mean[i], ";", sep="", file=bacon.file)
      }

    ### final parameters
    blah <- paste("\n\n##      K   MinYr   MaxYr   th0   th0p    w.a   w.b   alpha  beta",
      "\nBacon 0: ", ifelse(info$normal, "FixNor", "FixT"), ", ", length(info$d),
      ",  ", info$MinYr, ",  ", info$MaxYr, ",  ", info$th0[1], ",  ", info$th0[2], 
      ",  ", info$mem.strength*info$mem.mean, ",  ", info$mem.strength*(1-info$mem.mean), 
      ",  ", info$acc.shape[1], ",  ", info$acc.shape[1]/info$acc.mean[1], ", ", min(info$d), 
      ", ", max(info$d), ";\n", sep="")
    cat(blah, file=bacon.file)
  }



Bacon.PlotAgeModels <- function(info, dets=info$dets, d.lab="Depth", yr.lab="cal BP", bins=100, dark=5, prob=.95, rounded=0, draw.res=200, yr.res=10, rev.axis=FALSE, yr.lim=c(), store=TRUE, exx=15, cutoff=1e-4, plotMAP=FALSE, C14.col=rgb(0,0,1,.35), C14.border=rgb(0,0,1,.5), cal.col=rgb(0,.5,.5,.35), cal.border=rgb(0,.5,.5,.5), greyscale=function(x) rgb(0,0,0,x), main=c(), calc.every=info$calc.every, max.y=c(), normalise.dists=TRUE)
  {
    ### Load the twalk output and determinations
    cat("Reading twalk output:", paste(info$prefix, ".out", sep=""), "\n")
    info <- Bacon.AnaOut(paste(info$prefix, ".out", sep=""), "recacc.dat", info)

    ### leave for debugging purposes only
    Ana2(info)

    layout(matrix(if(length(info$hiatus.depths)>0) c(1:4,rep(5,4)) else c(1:3, rep(4,3)),
      nrow=2, byrow=TRUE), heights=c(.3,.7))
    par(mar=c(3,3,1,1), mgp=c(1.5,.7,.0), bty="l", yaxs="i")
    Bacon.PlotLogPost(info) # convergence information
    PlotAccPost(info, info$acc.shape, info$acc.mean)
    PlotMemPost(info, info$core, info$K, "", info$mem.strength, info$mem.mean, ds=1, Dc=info$Dc, max.y=max.y)
    if(length(info$hiatus.depths)>0)
      PlotHiatusPost(info, info$h.shape, info$h.mean)
    par(yaxs="r")
    
    calib.plot(read.table(info$det.file, header=TRUE), info=info, rev.axis=rev.axis, d.lim=range(info$d), d.lab=d.lab, yr.lab=yr.lab, yr.lim=yr.lim, plot.dists=FALSE)
    legend("topleft", ifelse(length(main)==0, info$core, main), bty="n", cex=1.5)
    d <- sort(c(seq(min(info$d), max(info$d), length=draw.res), info$hiatus.depths))
    hpds <- depth.ghost(d, info, rev.axis=rev.axis, dark=dark, yr.res=yr.res, prob=prob, greyscale=greyscale)
    calib.plot(read.table(info$det.file, header=TRUE), info=info, exx=exx, yr.res=yr.res, cutoff=cutoff, normal=info$normal, t.a=info$t.a, t.b=info$t.b, deltaR=info$deltaR, deltaSTD=info$deltaSTD, C14.col=C14.col, C14.border=C14.border, cal.col=cal.col, cal.border=cal.border, new.plot=FALSE, normalise.dists=normalise.dists)
    d <- sort(c(seq(min(info$d), max(info$d), by=calc.every), info$hiatus.depths))
    hpds <- depth.ghost(d, info, rev.axis=rev.axis, dark=dark, Plot=FALSE, prob=prob)
     
    info$hpds <- cbind(d, hpds)
    write.table(cbind(d[1:(length(d))], round(hpds, rounded)), paste(info$prefix, ".hpd", sep=""),
      col.names=c("depth", "min", "max", "MAP"), quote=FALSE, row.names=FALSE, sep="\t")

    ## calculate extremes of highest posterior density ranges for every depth
    rng <- hpds[,2]-hpds[,1]
    min.rng <- d[which(rng==min(rng))]
    max.rng <- d[which(rng==max(rng))]
    if(length(min.rng)==1) min.rng <- paste(" yr at", min.rng, info$unit) else
      min.rng <- paste(" yr between", min(min.rng), "and", max(min.rng), info$unit)
    if(length(max.rng)==1) max.rng <- paste(" yr at", max.rng, info$unit) else
      max.rng <- paste(" yr between", min(max.rng), "and", max(max.rng), info$unit)
    cat("Mean ", 100*prob, "% confidence intervals at ", round(mean(rng),0), " yr, min. ",
      min(rng), min.rng, ", max. ", max(rng), max.rng, "\n\n", sep="")

    if(store) info <<- info
  }



### for proxy.ghost
DepthsOfScore <- function(value, dat)
 {
   d <- c()
   for(i in 1:(nrow(dat)-1))
    {
      valueRange <- dat[i:(i+1),2]
      if(min(valueRange) <= value && max(valueRange) >= value)
        {
          slope <- (dat[i,2] - dat[i+1,2]) / (dat[i,1] - dat[i+1,1])
          intercept <- dat[i,2] - (slope*dat[i,1])
          if(slope==0) d[i-1] <- dat[i,1]
          d <- sort(c(d, (value - intercept) / slope ))
        }
      }
    unique(d)
  }


### to plot greyscale/ghost graphs of proxy values
proxy.ghost <- function(proxy=1, core=info$core, proxy.label=c(), res=100, yr.res=5, bins=info$bins, dark=1, d=info$d, rev.axis=FALSE, plotMAP=FALSE, yr.lim=c(), yr.lab="cal BP", proxy.lim=c(), xax="i", draw.box="l", sep="\t")
  {
    if(length(info$Tr)==0)
      stop("\nPlease first run Bacon.PlotAgeModels(info)\n\n")
    proxies <- read.csv(paste("Cores/", core, "/", core, "_proxies.csv", sep=""), header=TRUE)
    if(length(proxy.label)==0) proxy.label <- names(proxies)[proxy+1]
    proxy <- cbind(proxies[,1], proxies[,proxy+1])
    proxy <- proxy[!is.na(proxy[,2]),]
    proxy <- proxy[which(proxy[,1] <= max(d)),]
    proxy <- proxy[which(proxy[,1] >= min(d)),]
    if(length(unique(proxy[,2]))==1)
      stop("\nThese proxy values remain constant throughout the core, and cannot be plotted in grey-scale!\n\n")
    proxyseq <- seq(min(proxy[,2]), max(proxy[,2]), length=res)
    out <- list(yrseq=c(), binned=c(), maxs=c())
    ds <- c()
    d.length <- array(1, dim=c(res, 2))

    for(i in 1:res)
      {
        tmp  <- DepthsOfScore(proxyseq[i], proxy)
        ds <- c(ds, tmp)
        if(i > 1) d.length[i,1] <- d.length[(i-1),2]+1
        d.length[i,2] <- d.length[i,1]+length(tmp)-1
        if(length(tmp)==0) d.length[i,] <- d.length[i-1,]
      }
    cat("\nCalculating histograms... \n")
    Bacon.hist(ds, info, bins, Plot=FALSE)
    yr.min <- c()
    yr.max <- c()
    for(i in 1:length(hists))
      {
        yr.min <- min(yr.min, hists[[i]]$th0)
        yr.max <- max(yr.max, hists[[i]]$th1)
      }
    yr.seq <- seq(yr.min, yr.max, by=yr.res)

    all.counts <- array(0, dim=c(length(hists), length(yr.seq)))
    for(i in 1:length(hists))
      all.counts[i,] <- approx(seq(hists[[i]]$th0, hists[[i]]$th1, length=hists[[i]]$n), hists[[i]]$counts, yr.seq)$y
    all.counts[is.na(all.counts)] <- 0
    max.counts <- array(0, dim=c(res, length(yr.seq)))
    for(i in 1:res)
      for(j in 1:length(yr.seq))
        max.counts[i,j] <- max(all.counts[d.length[i,1]:d.length[i,2],j])
    if(length(yr.lim)==0)
    if(xax=="r")
      yr.lim <- range(pretty(c(1.04*max(yr.seq), .96*min(yr.seq)))) else # ugly but xaxs doesn't work with image!
      yr.lim <- range(yr.seq)[2:1]
    if(length(proxy.lim)==0) proxy.lim <- range(proxyseq)
    if(rev.axis)
    image(proxyseq, yr.seq, max.counts, xlim=proxy.lim, ylim=yr.lim[2:1], col=grey(seq(1,1-dark,-.02)), ylab=yr.lab, xlab=proxy.label, xaxs=xax) else
    image(yr.seq, proxyseq, t(max.counts), xlim=yr.lim, ylim=proxy.lim, col=grey(seq(1,1-dark,-.02)), xlab=yr.lab, ylab=proxy.label, xaxs=xax)
    box(bty=draw.box)
  }

############### adapt this to find "best" age-models using very small hpds (1%)?

### to plot greyscale/ghost graphs of age-depth model
depth.ghost <- function(d, info, bins=info$bins, prob=info$prob, dark=info$dark, rev.axis=FALSE, plotMAP=FALSE, yr.res=5, Plot=TRUE, range.col=1, best.col=2, greyscale=function(x) rgb(0,0,0,x))
  {
    hpds <- array(NA, dim=c(length(d), 2))
    Bacon.hist(d, info, bins)
    d.jumps <- diff(d)[1]

    peak <- 0
    for(i in 1:(length(d)))
      peak <- max(peak, hists[[i]]$counts)

    best <- c() ######

    for(i in 1:length(d))
      {
        yrs <- seq(hists[[i]]$th0, hists[[i]]$th1, length=bins)
        hst <- approx(yrs, hists[[i]]$counts, seq(min(yrs), max(yrs), by=yr.res))
	if(length(hst$y[which(hst$y > 0)]) < 2) 
	  {
	    best[i] <- hst$x[which(hst$y == max(hst$y))][1]
	    hpds[i,] <- range(hst$x)
	  } else
	    {
	      best[i] <- weighted.mean(hst$x, hst$y) ######
	      #best[i] <- hpd(hst, .05, info$yrs) ###### PREFERRED BUT DOES NOT WORK YET
	      tmp <<- hst
	      hpds[i,] <- hpd(hst, prob)
	      gr <- seq(0, min(1, dark*max(hst$y)/peak/diff(yrs[1:2])), length=50)
	      if(Plot)
		if(rev.axis)
		  image(hst$x, c(d[i]-(d.jumps/2), d[i]+(d.jumps/2)), t(t(hst$y)),
		    add=TRUE, col=greyscale(gr)) else
		  image(c(d[i]-(d.jumps/2), d[i]+(d.jumps/2)), hst$x, t(hst$y),
		    add=TRUE, col=greyscale(gr))
	    }
      }

best <<- best

    MAP <- which(info$Us==min(info$Us))[1] # the iteration with the best model
    MAP <- c(info$output[MAP,1], info$output[MAP,1]+
      cumsum(diff(info$d)*as.numeric(info$output[MAP,2:info$K])))
    MAP <- approx(info$d, MAP, d, rule=1)$y

    if(Plot)
      if(rev.axis)
	{
	  lines(hpds[,1], d, col=range.col, lty=3) # can't do hiatus yet
	  lines(hpds[,2], d, col=range.col, lty=3)
	  if(plotMAP) lines(MAP, d, col=best.col)
	} else
	{
	  lines(d, hpds[,1], col=range.col, lty=3)
	  lines(d, hpds[,2], col=range.col, lty=3)
	  if(plotMAP) lines(d, MAP, col=best.col)
	}
      cbind(hpds,  MAP)
  }



### calculate age distributions of depth(s)
Bacon.hist <- function(d, info, bins, tempfile="hists.dat", xlab="cal BP", ylab="Frequency", xlim=c(), Plot=TRUE)
  {
    infile <- noquote(paste(info$prefix, ".out", sep=""))
    tmp <- file(tempfile, "w")
    cat("bin/hist", infile, info$Tr, info$dets[1,4], info$Dc, info$K, bins, tempfile, file=tmp)
    for(i in d) cat(" ", i, file=tmp)
    cat("\n", file=tmp)
    system(readChar(tempfile, 1e5))
    source(tempfile)

    closeAllConnections()

    if(length(d)==1 && Plot==TRUE)
      {
	hst <- hists[[1]]
	if(length(xlim)==0) xlim <- c(hst$th1, hst$th0)
	pol <- cbind(c(hst$th0, seq(hst$th0, hst$th1, length=bins), hst$th1), c(0, hst$counts, 0))
	plot(0, type="n", xlim=c(hst$th1, hst$th0), ylim=c(0, max(pol[,2])), xlab=xlab, ylab=ylab, yaxs="i")
	polygon(pol, col=rgb(0,0,0,.5), border=rgb(0,0,0,.5))
      }
  }



### These help analyse the twalk output
### We will leave them for debugging purposes only
source("Ana.R")
# source("AutoCorr.R") # is now in Ana.R



### Time series of the log of the posterior
Bacon.PlotLogPost <- function(info, from=0, to=info$Tr)
  plot(from:(to-1), -info$Us[(from+1):to], type="l",
    ylab="Log of Objective", xlab="Iteration", main="")



### adapt to allow for hiatuses!
### Plots several models and MAP
Bacon.PlotPostModels <- function(info, from=1, to=info$Tr, length=100,
  MAP=sample(which(info$Us==min(info$Us)),1), rev.axis=FALSE)
  {
    for(i in seq(from, to, length=length))
      lines(info$d, c(info$output[i,1], info$output[i,1]+ cumsum(diff(info$d)*as.numeric(info$output[i,2:info$K]))), col="grey")
    if(rev.axis==FALSE)
      lines(info$d,
      c(info$output[MAP,1],info$output[MAP,1]+cumsum(diff(info$d)*as.numeric(info$output[MAP,2:info$K]))), col="red") else
      lines(      c(info$output[MAP,1],info$output[MAP,1]+cumsum(diff(info$d)*as.numeric(info$output[MAP,2:info$K]))), info$d, col="red")
  }



### cut away first bunch of iterations if MCMC burnin still visible
### NB will adapt the original .out file!
scissors <- function(burnin, set=info)
  {
    output <- read.table(paste(set$prefix, ".out", sep=""))
    if(burnin>= nrow(output)) 
	  stop("\nCannot remove that many iterations, there would be none left!\n\n", call.=FALSE)
    output <- output[burnin:nrow(output),]
    write.table(output, paste(set$prefix, ".out", sep=""), col.names=F, row.names=F)
    info$output <<- output
  }



### thin iterations by given proportion (if too much autocorrelation is visible in the MCMC series)
### NB will adapt the original .out file!
razor <- function(ratio=.1, set=info)
  {
    output <- read.table(paste(set$prefix, ".out", sep=""))
    if(ratio >= 1)
      stop("\nCannot remove that many iterations, there would be none left!\n\n", call.=FALSE)
    ratio <- sample(nrow(output), ratio*nrow(output))
    output <- output[-ratio,]
    write.table(output, paste(set$prefix, ".out", sep=""), col.names=F, row.names=F)
    info$output <<- output
  }



### plot chronological probability curves for events in the core
AgesOfEvents <- function(yrmin, yrmax, window, move, set=info, nicepars=TRUE, plot.steps=FALSE)
  {
    outfile <- paste(set$prefix, "_", window, "_probs.txt", sep="")
    file.create(outfile)
    MCMCname <- paste(set$prefix, ".out", sep="")    
    probfile <- paste("Cores/", set$core, "/", set$core, "_events.txt", sep="")
    if(!file.exists(probfile))
      stop("\nFile with probabilities for events per depth not found! Check the manual\n\n") 
    probs <- read.table(probfile)
    if(!is.numeric(probs[1,1]))
      stop("\nFirst line of the _events.txt file should NOT contain titles; please remove them\n\n")
    if(min(probs[,1]) < min(info$d) || max(probs[,1]) > max(info$d)) 
      {
	cat("\nSome depths in the _events.txt file go beyond the age-model; I will remove them\n\n")
	file.rename(probfile, paste(probfile, "_backup", sep=""))
	probs <- probs[which(probs[,1] >= min(info$d)),]
	probs <- probs[which(probs[,1] <= max(info$d)),]
	write.table(probs, probfile, col.names=FALSE, row.names=FALSE, quote=FALSE)
      }

    system(paste("bin/events", yrmin, yrmax, move, window, outfile, MCMCname, nrow(set$output), set$K, set$d[1], set$Dc, probfile, nrow(probs)))
    probs <- read.table(outfile)
    if(plot.steps)
      {
	d.sort <- sort(rep(1:nrow(probs),2))
	d.sort <- cbind(d.sort[-length(d.sort)], d.sort[-1])
	probs <- cbind(c(min(probs[,1]), probs[d.sort[,1],1], max(probs[,1])), c(0,probs[d.sort[,2],2],0))
      } else
	probs <- cbind(c(min(probs[,1]), probs[,1], max(probs[,1])), c(0,probs[,2],0))
    if(nicepars) par(yaxs="i", bty="l")
    plot(probs, type="n", xlab="cal. yr BP", ylab="prob (%)", ylim=c(0, 1.1*max(probs[,2])), xlim=c(yrmax, yrmin))
    polygon(probs, col="grey")
    if(move > window) cat("\nAre you sure you want the window widths to be smaller than the moves?\n\n")
  }


# for the moment only considers outer hpd ranges
# or simply plot the ranges as points, not lines
hpd <- function(hst, prob)
#	if(length(hst$y[hst$y > 0]) < 1) 
		#range(hst$x) else
		{
			o <- order(hst$y, decreasing=TRUE)
			hst <- cbind(hst$x[o], cumsum(hst$y[o])/sum(hst$y))
			if(min(hst[,2]) <= prob)
				range(hst[which(hst[,2]<=prob),1]) else
				range(hst[which(hst[,2]==max(hst[,2])),1])  # not pretty
		}



#calcurve0=constcal, calcurve1="IntCal09", calcurve2="IntCal09Marine", calcurve3="SHCal04", calcurve4=new 


### calibrate C14 dates and calculate distributions for any calendar dates
bacon.calib <- function(dat, info, calcurve1="IntCal09", calcurve2="Marine09", calcurve3="SHCal04", calcurve4="ConstCal", yr.res=1, cutoff=1e-7, normal=info$normal, t.a=3, t.b=4, d.R=0, d.STD=0)
  {
    if(calcurve1=="IntCal09") cc1 <- read.table("Curves/3Col_intcal09.14C") else
      cc1 <- read.csv(paste("Curves/", calcurve1, ".14c", sep=""), header=FALSE, skip=11)[,1:3]
    if(calcurve2=="Marine09") cc2 <- read.table("Curves/3Col_marine09.14C") else
      cc2 <- read.csv(paste("Curves/", calcurve2, ".14C", sep=""), header=FALSE, skip=11)[,1:3]
    if(calcurve3=="SHCal04") cc3 <- read.table("Curves/SHCal04.14C") else
      cc3 <- read.table(paste("Curves/", calcurve3, ".14C", sep=""))[,1:3]
    if(calcurve4 != "ConstCal") 
      cc4 <- read.table(paste("Curves/", calcurve4, sep=""))[,1:3]

    ## use Gaussian or t (Christen and Perez Radiocarbon 2009) calibration
    calib <- list()
    if(info$postbomb != 0)
      {
	if(info$postbomb==1) bomb <- read.table("Curves/postbomb_NH1.14C", header=FALSE) else
	if(info$postbomb==2) bomb <- read.table("Curves/postbomb_NH2.14C", header=FALSE) else
	if(info$postbomb==3) bomb <- read.table("Curves/postbomb_NH3.14C", header=FALSE) else
	bomb <- read.table("Curves/postbomb_SH.14C", header=FALSE)[,1:3]
	cc1 <- rbind(cc1, bomb[nrow(bomb):1,], deparse.level=0)
	cc3 <- rbind(cc3, bomb[nrow(bomb):1,], deparse.level=0)
      }

    d.cal <- function(cc, rcmean, w2, ta=t.a, tb=t.b, nor=normal, Cutoff=cutoff)
      {
        if(nor)
          cal <- cbind(cc[,1], dnorm(cc[,2], rcmean, sqrt(cc[,3]^2+w2))) else
          cal <- cbind(cc[,1], (tb + ((rcmean-cc[,2])^2) / (2*(cc[,3]^2 + w2))) ^ (-1*(ta+0.5)))
        cal[,2] <- cal[,2]/sum(cal[,2])
        cal <- cal[which(cal[,2]>0),]
        acc <- which(cal[,2] >= Cutoff)
        if(length(acc) > 1) cal <- cal[acc,]
        cal
      }

    if(ncol(dat)==4)
      for(i in 1:nrow(dat))
        calib[[i]] <- d.cal(cc1, dat[i,2]-d.R, dat[i,3]^2+d.STD^2, t.a, t.b, normal, cutoff) else
      {

        for(i in 1:nrow(dat))
          {
            det <- as.numeric(dat[i,])
            if(det[9]==0)
              {
                x <- seq(det[2]-max(100,4*det[3]), det[2]+max(100,4*det[3]), by=yr.res)
                cc <- cbind(x, x, rep(0,length(x))) # dummy 1:1 curve
                calib[[i]] <- d.cal(cc, det[2]-det[5], det[3]^2+det[6]^2,
                  det[7], det[8], normal, cutoff)
              } else
              {
                if(det[9]==1) cc <- cc1 else if(det[9]==2) cc <- cc2 else
                  if(det[9]==3) cc <- cc3 else cc <- cc4
                calib[[i]] <- d.cal(cc, det[2]-det[5], det[3]^2+det[6]^2,
                  det[7], det[8], normal, cutoff)
              }
          }
      }
    calib
  }


### produce plots of the calibrated distributions
calib.plot <- function(dat, info, rev.axis=FALSE, d.lim=range(dat[,4]), yr.lim=c(), yr.res=1, d.lab="Depth", yr.lab="cal yr BP", exx=15, cutoff=1e-3, C14.col=rgb(0,0,1,.5), C14.border=rgb(0,0,1,.75), cal.col=rgb(0,.5,.5,.5), cal.border=rgb(0,.5,.5,.75), normal=info$normal, t.a=3, t.b=4, deltaR=0, deltaSTD=0, new.plot=TRUE, plot.dists=TRUE, normalise.dists=TRUE)
  {
    calib <- bacon.calib(dat, info, info$calcurves[1], info$calcurves[2], info$calcurves[3], info$calcurves[4], yr.res, cutoff, normal, t.a, t.b, deltaR, deltaSTD)
    exx <- length(min(d.lim):max(d.lim)) * exx/50
    if(length(yr.lim)==0)
      {
		lims <- array(0, dim=c(length(calib),2))
		for(i in 1:length(calib)) lims[i,] <- range(calib[[i]][,1])
		yr.lim <- range(lims)
      }
    if(new.plot)  
      if(rev.axis)
        plot(0, type="n", xlim=yr.lim[2:1], ylim=d.lim[2:1], xlab=yr.lab, ylab=d.lab, main="") else
	      plot(0, type="n", xlim=d.lim, ylim=yr.lim, xlab=d.lab, ylab=yr.lab, main="")
      
    if(plot.dists)
      for(i in 1:nrow(dat))
        {
          cal <- calib[[i]]
		  o <- order(cal[,1])
		  if(normalise.dists)
			cal <- cbind(cal[o,1], exx*cal[o,2]/sum(cal[,2])) else
			cal <- cbind(cal[o,1], exx*cal[o,2]/max(cal[,2]))
          cal <- cal[cal[,2]>0,]
	      pol <- cbind(c(dat[i,4]-cal[,2], dat[i,4]+rev(cal[,2])), c(cal[,1], rev(cal[,1])))
	      if(rev.axis) pol <- cbind(pol[,2], pol[,1])
		  if(ncol(dat)==4 || (ncol(dat) > 4 && dat[i,9] > 0))
	       { col <- C14.col; border <- C14.border } else 
		   { col <- cal.col; border <- cal.border }
			polygon(pol, col=col, border=border) 
        }
  }



### FUNCTIONS TO PLOT THE PRIOR DISTRIBUTIONS

### plot the prior for the accumulation rate
PlotAccPrior <- function(s, mn, main="", xlim=c(0, 3*max(mn)), unit="cm", xlab=paste("Acc. rate (yr/", unit, ")", sep=""), ylab="Density", add=FALSE, legend=TRUE)
  {
    o <- order(s, decreasing=TRUE)
    priors <- unique(cbind(s[o],mn[o])[,1:2])
    if(length(priors) == 2)
      {
        curve(dgamma(x, s, s/mn), col=3, lwd=2, xlim=xlim, xlab=xlab, ylab=ylab, add=add)
        txt <- paste("acc.shape: ", priors[1], "\nacc.mean: ", priors[2])
      } else
      {
        curve(dgamma(x, priors[1,1], priors[1,1]/priors[1,2]), col=3, lwd=2, xlim=xlim, xlab=xlab, ylab=ylab, add=add)
        for(i in 2:nrow(priors))
          curve(dgamma(x, priors[i,1], priors[i,1]/priors[i,2]), col=3, lwd=2, xlim=xlim, xlab=xlab, ylab=ylab, add=if(i==1) add else TRUE)
        txt <- paste("acc.shape: ", toString(priors[,1]), "\nacc.mean: ", toString(priors[,2]))
      }
    if(legend)
      legend("topright", txt, bty="n", cex=.8, text.col=2, adj=c(0,0))
  }



### plot the prior for the memory (= accumulation rate varibility between neighbouring sections)
PlotMemPrior <- function(s, mn, Dc, ds=1, xlab="Memory (ratio)", ylab="Density", main="", add=FALSE, legend=TRUE, unit="cm")
  {
    curve(dbeta(x, s*mn, s*(1-mn)), 0, 1, col=3, lwd=2, xlab=xlab, ylab=ylab, add=add)
    txt <- paste("mem.strength:", s, "\nmem.mean:", mn, "\nsteps Dc:", Dc, unit)
    if(legend)
      legend("topright", txt, bty="n", cex=.8, text.col=2, adj=c(0,0))
    if(s*(1-mn) <= 1) cat("\nWarning! Chosen memory prior might cause problems.\nmem.strength * (1 - mem.mean) should be smaller than 1\n ")
  }



### plot the prior for the hiatus length
PlotHiatusPrior <- function(s, mn, hiatus=c(), xlab="Hiatus size (yr)", ylab="Density", main="", xlim=c(0, 3*max(mn)), add=FALSE, legend=TRUE)
  {
    o <- order(s, decreasing=TRUE)
    priors <- cbind(s[o],mn[o])[,1:2]
    if(length(priors) == 2)
      {
	curve(dgamma(x, priors[1], priors[1]/priors[2]), col=3, lwd=2, xlim=xlim, xlab=xlab, ylab=ylab, add=add) 
	txt <- paste("h.shape: ", toString(priors[1]), "\nh.mean: ", toString(priors[2]), "\nd: ", toString(hiatus))
      } else
      for(i in 2:nrow(priors))
	{
	  curve(dgamma(x, priors[i,1], priors[i,1]/priors[i,2]), col=3, lwd=2, xlim=xlim, xlab=xlab, ylab=ylab, add=if(i==1) add else TRUE)
	  txt <- paste("h.shape: ", toString(priors[,1]), "\nh.mean: ", toString(priors[,2]), "\nd: ", toString(hiatus))
	}
    if(legend)
      legend("topright", txt, bty="n", cex=.8, text.col=2, adj=c(0,0))
  }



### plot the posterior (and prior) of the accumulation rate
PlotAccPost <- function(info, s=info$acc.shape, mn=info$acc.mean, main="", xlab=paste("Acc. rate (yr/", info$unit, ")", sep=""), ylab="Frequency")
  {
    hi <- 2:(info$K-1)
    for(i in info$hiatus.depths) hi <- hi[-max(which(info$d < i))]
    post <- c()
    for(i in hi) post <- c(post, info$output[[i]])
    post <- density(post)
    lim.y <- c(0, max(post$y, 1.05*dgamma((s-1)/(s/mn), s, s/mn)))
    lim.x <- range(0, post$x, 2*mn)
    plot(0, type="n", xlim=lim.x, xlab=xlab, ylim=lim.y, ylab="")
    polygon(post, col=grey(.8), border=grey(.4))
    PlotAccPrior(s, mn, add=TRUE, xlim=range(post$x), xlab="", ylab=ylab, main=main, unit=info$unit)
  }



### plot the posterior (and prior) of the memory
PlotMemPost <- function(info, corenam, K, main="", s=info$mem.strength, mn=info$mem.mean, xlab=paste("Memory R (1", info$unit, "autocorrelation)"), ylab="Density", ds=1, Dc, unit=info$unit, max.y=c())
  {
    post <- density(info$output[,info$n]^(1/info$Dc), from=0, to=1)
    if(length(max.y)==0) max.y <- max(post$y, dbeta((0:100)/100, s*mn, s*(1-mn)))
    plot(0, type="n", xlab=xlab, xlim=c(0,1), ylim=c(0, max.y), ylab="", main="")
    polygon(post, col=grey(.8), border=grey(.4))
    PlotMemPrior(s, mn, Dc, add=TRUE, unit=unit, xlab="", ylab=ylab, main=main)
 }



### plot the posterior (and prior) of the hiatus
PlotHiatusPost <- function(info, shape=info$hiatus.shape, mean=info$hiatus.mean, main="", xlim=c(0, 3*max(info$acc.mean)), xlab=paste("Hiatus size (yr)", sep=""), ylab="Frequency", minbreaks=10)
  {
    hi <- c()
    for(i in info$hiatus.depths) hi <- c(hi, max(which(info$d <= i)))
    post <- c()
    for(i in hi) post <- c(post, info$output[[i+1]])
    plot(0, type="n", main="", xlab=xlab, xlim=xlim, ylab=ylab)
    polygon(density(post), col=grey(.8), border=grey(.4))
    PlotHiatusPrior(info$hiatus.shape, info$hiatus.mean, info$hiatus.depths, add=TRUE, xlim=xlim, xlab="", ylab=ylab, main=main)
  }



Cores <- function() list.files("Cores/") 


### Comments start with three #
### We will use global variables only for Core independent things, eg. Bacon.NumCols etc.



### dustbin functions (obsolete but could be recycled some day)

# AgeDiff <- function(dmin, dmax, info, bins, res, tempfile="hists.dat", xlab="Age difference (yr)")
#   {
#     Bacon.hist(c(dmin, dmax), info, bins, tempfile)
#     source(tempfile)
#     dmin <- hists[[1]]
#     dmax <- hists[[2]]
#     lims <- range(dmin$th0, dmin$th1, dmax$th0, dmax$th1)
#     lims <- seq(min(lims), max(lims), res)
#     dmin <- approx(seq(dmin$th0, dmin$th1,, dmin$n), dmin$counts, lims)$y
#     dmax <- approx(seq(dmax$th0, dmax$th1,, dmax$n), dmax$counts, lims)$y
#     dmin[is.na(dmin)] <- 0
#     dmax[is.na(dmax)] <- 0
#     dmin <- dmin/sum(dmin)
#     dmax <- dmax/sum(dmax)
#     counts <- c()
#     for(i in 0:(length(lims)-1)) # the 'delta' jumps
#       counts[i+1] <- sum(dmin[1:(length(lims)-i)]*dmax[(i+1):length(lims)])
#     plot(lims-min(lims), counts,  type="s", yaxs="i", ylim=c(0, 1.1*max(counts)), xlab=xlab)
#     segments(0, 0, 0, counts[1])
#   }

# Bacon.Age <- function(d, info, it, Ra)
#   if(j <- max(which(Ra <= d)) == 1)
#     info$output[it,1] + (d-Ra[j])*as.numeric(info$output[it,(j+1)]) else
#     info$output[it,1] + sum(diff(Ra[1:j])*as.numeric(info$output[it,2:j])) +
#     (d-Ra[j])*as.numeric(info$output[it,(j+1)])
