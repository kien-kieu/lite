# Line Tessellation (LiTe) library
# Release 1.1
# Authors: Katarzyna Adamczyk and Kiên Kiêu.
# Copyright INRA 2006-2015.
# Interdeposit Certification: IDDN.FR.001.030007.000.R.P.2015.000.31235
# License: GPL v3.

mySeed <- 5
tau <- 1.7
alpha <- 100
nupdates <- 15000
width <- 2.0

myTes <- new(TTessel)
setLiTeSeed(mySeed)
myTes$setDomain(width,width)
myMod <- ModelSquaredAreas(tau=tau,alpha=alpha*tau^4)
myMod$set_ttessel(myTes)
mySim <- new(SMFChain,myMod,0.33,0.33)

dev.new()
devTes <- dev.cur()
dev.new()
devEng <- dev.cur()
if(interactive()) {
    buf <- readline("make both new graphic devices visible and press enter")
}
dev.set(devTes)
plot(myTes,asp=1)
steps <- seq(from=0,to=nupdates,by=100)
energies <- rep(NA,length(steps))
dev.set(devEng)
plot(steps,energies,type="n",xlim=range(steps),ylim=c(-200,500),
     xlab="iteration",ylab="energy")
for(i in seq(along=steps)[-1]) {
    buf <- mySim$step(100)
    energies[i] <- myMod$getValue()
    dev.set(devEng)
    segments(steps[i-1],energies[i-1],steps[i],energies[i])
    dev.set(devTes)
    plot(myTes,asp=1)
}

