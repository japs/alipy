import pyfits
import numpy as np


#img = pyfits.open("SMARTS.fits")

#print img[0]

data = pyfits.getdata("images4/V.fits")

print data.shape
data = data[200:700,200:700]
data= np.rot90(data)

pyfits.writeto("images4/Vcrop.fits", data)