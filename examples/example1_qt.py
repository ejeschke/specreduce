import sys
import numpy as np
from astropy.io import fits
from astropy import modeling as md

from PySpectrograph.Models import RSSModel

import specreduce

from specreduce.interidentify import InterIdentify
from specreduce import spectools as st
from specreduce import WavelengthSolution

#from the command line, run:
#
#python test_specred_rss.py [image name] [line list]


function='poly'
order=3
rstep=1
nrows=1
mdiff=20
thresh=3
niter=5
dc=3
ndstep=50
dsigma=5
method='Zeropoint'
res=2
dres=0.2
filename=None
smooth=3
inter=True
subback=0
textcolor='green'
log = None

linelist=sys.argv[2]
slines, sfluxes = st.readlinelist(linelist)

hdu = fits.open(sys.argv[1])

data = hdu[1].data
xarr = np.arange(data.shape[1])


grating = hdu[0].header['GRATING']
slitname = hdu[0].header['MASKID']
slit = st.getslitsize(slitname)
grang = hdu[0].header['GR-ANGLE']
arang = hdu[0].header['AR-ANGLE']

xbin, ybin = hdu[0].header['CCDSUM'].split()
xbin = int(xbin)
ybin = int(ybin)
xpos = -0.2666
ypos = 0.0117
objid = None


rss = RSSModel.RSSModel(grating_name=grating.strip(), gratang=grang,
                        camang=arang, slit=slit, xbin=xbin, ybin=ybin,
                        xpos=xpos, ypos=ypos)

rss.gamma = 0

res = 1e7 * rss.calc_resolelement(rss.alpha(), -rss.beta())
dres = res / 10.0
wcen = 1e7 * rss.calc_centralwavelength()
R = rss.calc_resolution(wcen / 1e7, rss.alpha(), -rss.beta())

ws = st.useRSSModel(xarr, rss, md.models.Polynomial1D(order), gamma=rss.gamma)
ws.fit()
print ws.x
print ws.wavelength
print ws.coef

istart = int(data.shape[0]/2)
InterIdentify(xarr, data, slines, sfluxes, ws, mdiff=mdiff, rstep=rstep,
              function=function, order=order, sigma=thresh, niter=niter,
              res=res, dres=dres, dc=dc, ndstep=ndstep, istart=istart,
              method=method, smooth=smooth, filename=filename,
              subback=subback, textcolor=textcolor, log=log, verbose=True)
