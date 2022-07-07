#! /usr/bin/env python

from netCDF4 import Dataset
import matplotlib
#matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import array
import matplotlib.cm as mcm
from mpl_toolkits.basemap import Basemap
from collections import defaultdict
import cmocean as cm
from scipy.spatial import cKDTree
import glob
import collections
import struct
import datetime
import time
import sys
sys.path.append('/home/bzhao/python_utils')
#import read_utils
#import plot_utils
#import math_utils
#import read_utils
import data_utils
#import get_info
#from pylab import *
import functools


def plot_pole_new(LON,LAT,VAR,levels,setcolor,setnorm,titlestr,POLE,ptype,MER):

    if  (POLE=='N'):
        m = Basemap(projection='npstere',lon_0=0,boundinglat=45)
    if  (POLE=='S'):
        m = Basemap(projection='spstere',lon_0=180,boundinglat=-45)
    m.drawcoastlines()
    m.fillcontinents()
    m.drawcountries()
    plt.title(titlestr)
    x, y =m(LON,LAT)
    if ptype=='cont':
        m.contourf(x,y,VAR, levels, origin='lower',cmap=setcolor, norm=setnorm, extend='both')
        #plt.colorbar(orientation='vertical',extend='both',shrink=0.4)
    if ptype=='plot':
        plt.plot(x,y,'.')#marker='.',color='k')
    if ptype=='scatter':
        #plt.plot(x,y,'.')#marker='.',color='k')
        VAR[abs(VAR)>999]='nan'
        print 'min='+str(np.nanmin(VAR))
        print 'max='+str(np.nanmax(VAR))
        valmin=min(levels)
        valmax=max(levels)

        plt.scatter(x,y,8*VAR/VAR,VAR,marker='o',vmin=valmin,vmax=valmax,cmap='jet',linewidths=0)

    m.drawparallels(np.arange(-90.,120.,15.),labels=[1,0,0,0]) # draw parallels
    m.drawmeridians(np.arange(0.,420.,30.),labels=MER) # draw meridians

def plot_pole_2(LON,LAT,VAR, LON1,LAT1,VAR1, levels, levels1, setcolor,setnorm,titlestr,POLE,ptype,MER):

    if  (POLE=='N'):
        m = Basemap(projection='npstere',lon_0=0,boundinglat=85)
    if  (POLE=='S'):
        m = Basemap(projection='spstere',lon_0=180,boundinglat=-45)
    m.drawcoastlines()
    m.fillcontinents()
    m.drawcountries()
    plt.title(titlestr)
    x, y =m(LON1,LAT1)
    plt.contour(x,y,VAR1, levels1, origin='lower',colors='k', linewidths=2)
    x, y =m(LON,LAT)
    plt.contourf(x,y,VAR, levels, origin='lower',cmap=setcolor, norm=setnorm, extend='both')

    m.drawparallels(np.arange(-90.,120.,15.),labels=[1,0,0,0]) # draw parallels
    m.drawmeridians(np.arange(0.,420.,30.),labels=MER) # draw meridians




class saltwatertile:

    def __init__(self, file): 
         
       header = np.genfromtxt(file, dtype='i4', usecols=(0), max_rows=8)
       #print header
       self.atm = 'x'.join([str(x) for x in header[3:5]])
       self.ocn = 'x'.join([str(x) for x in header[6:]])
       self.nx, self.ny = header[6], header[7]
       print self.atm, self.ocn 
       tile=np.genfromtxt(file, dtype=[('type','i1'), ('area','f8'), ('lon','f8'),('lat','f8'), ('gi1','i4'),
                           ('gj1','i4'), ('gw1','f8'),
                           ('idum','i4'), ('gi2','i4'), ('gj2','i4'), ('gw2','f8')], skip_header=8)
       n1=0
       n2=0
       for n in range(1, tile.shape[0]+1, 1):
           if tile[n-1][0] == 0:
               n1 = n
               break
       #print n1
       for n in range(n1, tile.shape[0]+1, 1):
           if tile[n-1][0] != 0:
               n2 = n
               break
       #print n2
       icetile=tile[n1-1:n2-1]
       #print icetile.shape
       #print 'hhh: ',icetile[0][2], icetile[-1][2]
       self.size = icetile.shape[0]
       self.gi = icetile['gi2'][:]
       self.gj = icetile['gj2'][:]
       print 'hhh: ',self.size,self.gi[-1],self.gj[-1]
       #return icetile


def get_nearest(lon, lat, LON, LAT, rad):
    lon[lon>80.0]=lon[lon>80.0]-360.0
    xs, ys, zs = lon_lat_to_cartesian(lon.flatten(), lat.flatten())
    xt, yt, zt = lon_lat_to_cartesian(LON.flatten(), LAT.flatten())
    points_in = zip(xs, ys, zs)
    print len(points_in)
    tree = cKDTree(points_in)
    #find indices of the nearest neighbors in the flattened array
    #d, inds = tree.query(zip(xt, yt, zt), k = 1)
    #get interpolated 2d field
    #zout = LON.copy().flatten()
    points = zip(xt, yt,zt)
    #print len(points)
    d, inds = tree.query(points, k = 1)
    return inds
    

def nearest_interp_new(z, LON, LAT, inds):
    zout = z.flatten()[inds].reshape(LON.shape)
    #zout.shape = LON.shape
    return zout

def lon_lat_to_cartesian(lon, lat, R = 1):
    """
    calculates lon, lat coordinates of a point on a sphere with
    radius R
    """
    lon_r = np.radians(lon)
    lat_r = np.radians(lat)

    x =  R * np.cos(lat_r) * np.cos(lon_r)
    y = R * np.cos(lat_r) * np.sin(lon_r)
    z = R * np.sin(lat_r)
    return x,y,z


def get_grid(mapl_tripole): #reads lat lon for tripolar ocean grid 
    #print atm, ocn 
    ncfile=Dataset(mapl_tripole,"r")

    LON     = ncfile.variables['lon_centers'][:]
    LAT     = ncfile.variables['lat_centers'][:]
    numlev   = ncfile.variables['mask'][:]
    ncfile.close()

    return LON, LAT, numlev



im,jm=360,120
km = 30
kitd = 12
lat_cutoff = 60.0

#CICE dimensions

ncat = 5
nilyr = 4
nslyr = 1


p5 = 0.5
c1 = 1.0
c2 = 2.0
Tice = 273.15
Tf = -1.8
saltmax = 3.2
nsal = 0.407
msal = 0.573
depressT = 0.054
rhoi      = 917.0
rhow      = 1026.0
rhos      = 330.0
cp_ice    = 2106.
cp_ocn    = 4218.
Lsub      = 2.835e6
Lvap      = 2.501e6
Lfresh    = Lsub-Lvap 

Tmlt = np.zeros(nilyr+1, dtype='float64')
salin = np.zeros(nilyr+1, dtype='float64')
for k in range(nilyr):
   zn = (float(k+1)-p5) / float(nilyr)
   salin[k]=(saltmax/c2)*(c1-np.cos(np.pi*np.power(zn,nsal/(msal+zn))))
   Tmlt[k] = -salin[k]*depressT
salin[nilyr] = saltmax
Tmlt[nilyr] = -salin[nilyr]*depressT

#PIO_DIR='/gpfsm/dnb04/projects/p94/verification/PIOMAS/RAW/MONTHLY'
PIO_DIR='/discover/nobackup/bzhao/ObservationData/PIOMAS'

grid_file = PIO_DIR+'/grid.dat'
mask_file = PIO_DIR+'/io.dat_360_120.output'

Year = sys.argv[1]
Mon = int(sys.argv[2])
Day = int(sys.argv[3])
icein = sys.argv[4]
p = icein.split('/')
iceout = '/'.join(p[:-1]+[p[-1]+'_piomas_inserted_'+'e'+Year+'-'+str(Mon)+'-'+str(Day)])
tilefile = sys.argv[5]
fld = 'giceday'
print icein
print iceout

#data_file = PIO_DIR+'/iceprod.H'+Year
#data_file = PIO_DIR+'/tice0.H'+Year
data_file = PIO_DIR+'/'+fld+'.H'+Year

index=[[0 for _ in range(ncat-1)] for _ in range(ncat)]
ii=range(ncat)
dist = defaultdict(list) 
for j in ii:
   for i in ii:
      if i != j:
          dist[j].append((abs(i-j),i+1)) 
for j in dist:
    dist[j].sort(key=lambda x: (x[0], x[1]*(-1)))
for j in dist:
    index[j]= [x[1]-1 for x in dist[j]]     
    #print j, dist[j]
print index
#POLE='N'

xx=np.loadtxt(grid_file)
#print xx.shape
yy = np.reshape(xx,(im*jm*2))
#print yy.shape

clon = np.reshape(yy[:im*jm],(jm,im))
clat = np.reshape(yy[im*jm:],(jm,im))

xx = np.genfromtxt(mask_file, dtype='i2', delimiter=2)
cmask = np.reshape(xx,(jm,im))


#print clon.shape
#print clat.shape

zz=None

with open(data_file, 'rb') as f: 
  for n in range(365):
      data = np.fromfile(f,dtype='float32', count=im*jm*kitd)
      if datetime.datetime(1999, Mon, Day).timetuple().tm_yday == n+1:
         zz = data.copy()

#print recl
#print data.shape
hice = np.reshape(zz, (kitd,jm,im))


#print hice[0,:20]


piomas_ic=np.array([0.0, 0.26, 0.71, 
              1.46, 2.61, 4.23, 6.39, 9.10, 12.39, 16.24, 20.62, 25.49], dtype='float64')

fac1=((0.6445-0.485)/(1.085-0.485))
fac2=((1.391-1.085)/(2.035-1.085))
fac3=((2.47-2.035)/(3.42-2.035))
fac4=((4.567-3.42)/(5.31-3.42))


missing=-9999.0


hicepm = hice.transpose()*piomas_ic
hicepm = hicepm.transpose()
hicepm = np.sum(hicepm, axis=0)


sw = saltwatertile(tilefile)
LON, LAT, numlevels = get_grid(sys.argv[6])

with Dataset(icein) as src, Dataset(iceout, "w") as dst:
    # copy global attributes all at once via dictionary
    dst.setncatts(src.__dict__)
    # copy dimensions
    for name, dimension in src.dimensions.items():
        dst.createDimension(
            name, (len(dimension) if not dimension.isunlimited() else None))
    # copy all file data except for the excluded
    for name, variable in src.variables.items():
        x = dst.createVariable(name, variable.datatype, variable.dimensions)
        dst[name][:] = src[name][:]
        # copy variable attributes all at once via dictionary
        dst[name].setncatts(src[name].__dict__)
    aicenout = dst['FR']
    vicenout = dst['VOLICE']
    vsnonout = dst['VOLSNO']
    tskinout = dst['TSKINI']
    eicenout = dst['ERGICE']
    esnonout = dst['ERGSNO']
    aicen = dst['FR'][:]
    vicen = dst['VOLICE'][:]
    vsnon = dst['VOLSNO'][:]
    tskin = dst['TSKINI'][:]
    eicen = dst['ERGICE'][:]
    esnon = dst['ERGSNO'][:]
    eicen = np.swapaxes(eicen,0,1)
    esnon = np.swapaxes(esnon,0,1)
    aicenpm5 = aicen.copy(order='C')
    vicenpm5 = vicen.copy(order='C')

    start = time.time()
    aicenpm = np.zeros((kitd, sw.ny, sw.nx), dtype='float64', order='C')
    hice_in = np.ma.array(hice[0])
    hice_in.mask = cmask == 0 
    lon_in_filtered = clon[np.logical_not(hice_in.mask)]
    lat_in_filtered = clat[np.logical_not(hice_in.mask)]
    nn = get_nearest(lon_in_filtered, lat_in_filtered, LON, LAT, 25)
    for n in range(kitd):
        hice_in = np.ma.array(hice[n])
        hice_in.mask = cmask == 0 
        hice_in_filtered = hice_in[np.logical_not(hice_in.mask)]
        hice_out = nearest_interp_new(hice_in_filtered, LON, LAT, nn)
        aicenpm[n,:,:] = hice_out
    end = time.time()
    print("Elapsed (interplating to tripolar) = %s" % (end - start))
    indi = sw.gi[:]-1 
    indj = sw.gj[:]-1   
    print indi.max(), indj.max()
    ind = np.array(range(sw.size))

    indis = indi[LAT[indj,indi]>lat_cutoff]
    indjs = indj[LAT[indj,indi]>lat_cutoff]
    inds = ind[LAT[indj,indi]>lat_cutoff]
    xx = piomas_ic[6:].reshape(6,1) 

    start = time.time()
    aicenpm5[0, inds] = aicenpm[1,indjs,indis]           + aicenpm[2,indjs,indis] * fac1     
    aicenpm5[1, inds] = aicenpm[2,indjs,indis]*(1.-fac1) + aicenpm[3,indjs,indis] * fac2
    aicenpm5[2, inds] = aicenpm[3,indjs,indis]*(1.-fac2) + aicenpm[4,indjs,indis] * fac3
    aicenpm5[3, inds] = aicenpm[4,indjs,indis]*(1.-fac3) + aicenpm[5,indjs,indis] * fac4
    aicenpm5[4, inds] = aicenpm[5,indjs,indis]*(1.-fac4) + np.sum(aicenpm[6:,indjs,indis],axis=0)
    vicenpm5[0, inds] = aicenpm[1,indjs,indis]*piomas_ic[1]           + aicenpm[2,indjs,indis]*piomas_ic[2]*fac1
    vicenpm5[1, inds] = aicenpm[2,indjs,indis]*piomas_ic[2]*(1.-fac1) + aicenpm[3,indjs,indis]*piomas_ic[3]*fac2
    vicenpm5[2, inds] = aicenpm[3,indjs,indis]*piomas_ic[3]*(1.-fac2) + aicenpm[4,indjs,indis]*piomas_ic[4]*fac3
    vicenpm5[3, inds] = aicenpm[4,indjs,indis]*piomas_ic[4]*(1.-fac3) + aicenpm[5,indjs,indis]*piomas_ic[5]*fac4
    vicenpm5[4, inds] = aicenpm[5,indjs,indis]*piomas_ic[5]*(1.-fac4) + np.sum(aicenpm[6:,indjs,indis]*xx, axis=0)
    exc = np.sum(aicenpm5,axis=0)-1.
    for n in range(ncat):
        aicenpm5[n,np.logical_and(exc > 5.e-7, aicenpm5[n,:]>exc)] -= exc[np.logical_and(exc > 5.e-7, aicenpm5[n,:]>exc)]  
        exc = np.sum(aicenpm5,axis=0)-1.

    hs = np.zeros((ncat,sw.size), dtype='float64') 
    hin = np.zeros((ncat,sw.size), dtype='float64') 
    qin = np.zeros((nilyr,ncat,sw.size), dtype='float64')
    qsn = np.zeros((nslyr,ncat,sw.size), dtype='float64')
    #hs[aicen>0.0] = vsnon[aicen>0.0]/aicen[aicen>0.0]
    qsn[:,vsnon>0.0] = esnon[:,vsnon>0.0]*nslyr/vsnon[vsnon>0.0] 
    qin[:,vicen>0.0] = eicen[:,vicen>0.0]*nilyr/vicen[vicen>0.0] 
    aicen[:,inds] = aicenpm5[:,inds] 

    yy = np.zeros((ncat, sw.size))
    yy[:, LAT[indj,indi]>lat_cutoff] = 1.0
    maska = np.logical_and(yy>0.0, np.logical_and(vicenpm5 > 0.0, vicen > 0.0))
    maskb = np.logical_and(yy>0.0, np.logical_and(vicenpm5 > 0.0, vicen == 0.0))
    maskc = np.logical_and(yy>0.0, np.logical_and(vicenpm5 == 0.0, vicen > 0.0))

    vicen[maska] = vicenpm5[maska] 
    eicen[:,maska] = qin[:,maska]*vicen[maska]/nilyr

    aicen[maskc] = 0.0
    vicen[maskc] = 0.0
    vsnon[maskc] = 0.0
    tskin[maskc] = Tice
    eicen[:,maskc] = 0.0
    esnon[:,maskc] = 0.0

       
    qi0 = -rhoi * (cp_ice*(Tmlt[:-1]-Tf) 
           + Lfresh*(1.-Tmlt[:-1]/Tf) - cp_ocn*Tmlt[:-1])
    qi0 = np.reshape(np.tile(qi0,ncat*sw.size),(sw.size,ncat,nilyr))  
    qi0 = qi0.transpose()
    qim = np.zeros((ncat,nilyr,sw.size), dtype='float64')
    tim = np.zeros((ncat,sw.size), dtype='float64')

    for n in range(ncat):
        for k in range(len(index[n])):
           qimq = qim[n,:,:] 
           timq = tim[n,:] 
           qinc = qin[:,index[n][k],:]
           timc = tskin[index[n][k],:]
           qimq[qimq==0.0] = qinc[qimq==0.0] 
           timq[timq==0.0] = timc[timq==0.0] 
           qim[n,:,:] = qimq
           tim[n,:] = timq
    for n in range(ncat):
        qimq = qim[n,:,:] 
        qinc = qi0[:,index[n][k],:]
        qimq[qimq==0.0] = qinc[qimq==0.0] 
        qim[n,:,:] = qimq
    tim[tim==0.0] = Tf + Tice
    qim = np.swapaxes(qim,0,1)
    vicen[maskb] = vicenpm5[maskb] 
    eicen[:,maskb] = qim[:,maskb]*vicen[maskb]/nilyr
    tskin[maskb] = tim[maskb]

    hs[aicen>0.0] = vsnon[aicen>0.0]/aicen[aicen>0.0]
    hin[aicen>0.0] = vicen[aicen>0.0]/aicen[aicen>0.0]
    hsn = (rhow - rhoi) / rhos * hin
    hs[hs>hsn] = hsn[hs>hsn]
    vsnon = hs*aicen 
    esnon[:,vsnon>0.0] = qsn[:,vsnon>0.0]*vsnon[vsnon>0.0]/nslyr 
    esnon[:,vsnon==0.0] = 0.0 

    inds = ind[np.logical_and(LAT[indj,indi]<60.0, LAT[indj,indi]>0.0)]
    aicen[:,inds] = 0.0
    vicen[:,inds] = 0.0
    eicen[:,:,inds] = 0.0
    esnon[:,:,inds] = 0.0
    vsnon[:,inds] = 0.0
    tskin[:,inds] = Tice

    aicenout[:] = aicen[:]
    tskinout[:] = tskin[:]
    vicenout[:] = vicen[:]
    vsnonout[:] = vsnon[:]
    eicenout[:] = np.swapaxes(eicen,0,1) 
    esnonout[:] = np.swapaxes(esnon,0,1)

    end = time.time()
    print("Elapsed (aggregating onto CICE ITD) = %s" % (end - start))
       
#lat_p = LAT[numlevels>0].flatten()




cmp = mcm.get_cmap('jet')
meridians=[1,0,1,1]

fig=plt.figure(figsize=(10,10), facecolor='w')
#fig.subplots_adjust(left=0.05, right=1.0, top=0.99, bottom=0.01,wspace=0.05,hspace=0.05)
ax=fig.add_axes([0.06, 0.0, 0.98, 0.98])

m = Basemap(projection='npstere',lon_0=0,boundinglat=45, resolution='l')
m.drawcoastlines()
m.fillcontinents()
m.drawcountries()
plt.title('PIOMAS: '+fld+' '+Year+'-'+str(Mon)+'-'+str(Day),y=1.05,size=25)

x, y =m(clon,clat)
#outside = (x <= m.xmin) | (x >= m.xmax) | (y <= m.ymin) | (y >= m.ymax)
#fbot = ma.masked_where(outside, fbot)
#m.pcolormesh(x,y,hice,cmap=cmp,vmin=0.0, vmax=4.0)
levl = 0 #m.pcolormesh(x,y,hice,cmap=cmp,vmin=-2, vmax=2)
m.pcolormesh(x,y,hicepm,cmap=cmp,vmin=0.0, vmax=5.0)
#m.pcolormesh(x,y,hice[levl],cmap=cmp,vmin=0.0, vmax=4.0)

#m.pcolormesh(x,y,hice,cmap=cmp,vmin=250, vmax=280)
#if POLE == 'N':
 #  m.plot(0.0,90.0,'ko',markersize=15, latlon=True)
m.drawparallels(np.arange(-90.,120.,15.),labels=[1,0,0,0]) # draw parallels
m.drawmeridians(np.arange(0.,420.,30.),labels=meridians) # draw meridians
plt.colorbar(orientation='vertical',extend='both',shrink=0.8)
#plt.tight_layout()
plt.show()
 

#print clon[0,:20]


