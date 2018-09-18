#library(flsa)
library(pixmap)
library(flsa)

### make the letters; a 32x32 immage
fl = matrix(numeric(32*32), ncol=32)
### make the S
fl[3:5,3:14]=1
fl[3:17,3:5]=1
fl[15:17, 3:14]=1
fl[15:29, 3:5]=1

### make the U
fl[3:29, 19:21]=1
fl[27:29, 19:30]=1

flimtrue = pixmapGrey(fl)

### add noise to the image
set.seed(1)
noise = runif(32*32) * (runif(32*32)<0.5) 
flnoise = fl + noise
flimnoise = pixmapGrey(flnoise)

fldenoise = flsa(flnoise, lambda2=0.5, verbose=T)
flimdenoise = pixmapGrey(fldenoise[1,,])

pdf("FLTrue.pdf")
plot(flimtrue)
dev.off()

pdf("FLNoise.pdf")
plot(flimnoise)
dev.off()

pdf("FLDenoised.pdf")
plot(flimdenoise)
dev.off()

