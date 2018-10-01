require(png)
require(EBImage)
require(RCurl)
require(PET)

im <- readPNG(getBinaryURL("https://cdn3.iconfinder.com/data/icons/halloween-29/64/ghost-512.png"))[,,1]
rad = radon(t(im))$rData
# Normalize intensity values from 0-1
rad = normalize(rad)
display(rad)
display(t(im))




P <- phantom()
R <- radon(P)
ir <- iradon(rad, XSamples=257, YSamples=257)
display(ir$irData)

rm(P,R)