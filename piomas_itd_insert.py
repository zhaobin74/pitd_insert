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

def get_base_array(a):
    """Get the ndarray base which owns memory block."""
    if isinstance(a.base, np.ndarray):
        return get_base_array(a.base)
    return a

def share_base(a, b):
    return get_base_array(a) is get_base_array(b)


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
       print n1
       for n in range(n1, tile.shape[0]+1, 1):
           if tile[n-1][0] != 0:
               n2 = n
               break
       print n2
       icetile=tile[n1-1:]
       print icetile.shape
       print 'hhh: ',icetile[0][2], icetile[-1][2]
       self.size = icetile.shape[0]
       self.gi = icetile['gi2'][:]
       self.gj = icetile['gj2'][:]
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
    print len(points)
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


def get_grid(atm, ocn): #reads lat lon for tripolar ocean grid 
    ##ncfile=Dataset('/gpfsm/dnb42/projects/p17/gvernier/SAND_BOXES/PLOT_ODAS/DATA/grid_spec_720x410x40.nc', "r")
    #ncfile=Dataset('/discover/nobackup/yvikhlia/coupled/Forcings/Ganymed/a90x540_o720x410/INPUT/grid_spec.nc',"r")
    #ncfile=Dataset('/gpfsm/dnb02/projects/p23/bzhao/s2s3-duoc04/scratch/INPUT/grid_spec.nc',"r")
    if ocn=='1440x1080':
       ncfile = Dataset('/discover/nobackup/yvikhlia/coupled/Forcings/a'+atm+'_o'+ocn+'.newtile/INPUT/grid_spec.nc', "r")
    else:
       ncfile = Dataset('/discover/nobackup/yvikhlia/coupled/Forcings/a'+atm+'_o'+ocn+'/INPUT/grid_spec.nc', "r")
    LON     = ncfile.variables['x_T'][:]
    LAT     = ncfile.variables['y_T'][:]
    numlev     = ncfile.variables['num_levels'][:]
    ncfile.close()

    return LON, LAT, numlev



im,jm=360,120
km = 30
kitd = 12


#CICE dimensions

ncat = 5
nilyr = 4
nslyr = 1


p5 = 0.5
c1 = 1.0
c2 = 2.0
Tf = -1.8
Tice = 273.15
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
PIO_DIR='/gpfsm/dnb02/bzhao/ObservationData/PIOMAS'

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


sw = saltwatertile(tilefile)
LON, LAT, numlevels = get_grid(sw.atm, sw.ocn)

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
    aicen = dst['FR']
    vicen = dst['VOLICE']
    vsnon = dst['VOLSNO']
    eicen = dst['ERGICE']
    esnon = dst['ERGSNO']
    tskin = dst['TSKINI']
    aicenpm5 = aicen[:].copy()
    vicenpm5 = vicen[:].copy()
    #print('aicen, aicenpm5    ->', share_base(aicen, aicenpm5))  
    aicenpm = np.zeros((kitd, sw.ny, sw.nx), dtype='float64')
    hice_in = np.ma.array(hice[0])
    hice_in.mask = cmask == 0 
    lon_in_filtered = clon[np.logical_not(hice_in.mask)]
    lat_in_filtered = clat[np.logical_not(hice_in.mask)]
    nn = get_nearest(lon_in_filtered, lat_in_filtered, LON, LAT, 25)
    for n in range(kitd):
        hice_in = np.ma.array(hice[n])
        hice_in.mask = cmask == 0 
        hice_in_filtered = hice_in[np.logical_not(hice_in.mask)]
        #hice_out = nearest_interp(lon_in, lat_in, hice_in, LON, LAT)
        #hice_out = nearest_interp_new(lon_in, lat_in, hice_in, LON, LAT, 25)
        hice_out = nearest_interp_new(hice_in_filtered, LON, LAT, nn)

        #hice[n][numlevels>0] = hice_out[numlevels>0]
        #hice_out[numlevels==0] = 0
        aicenpm[n,:,:] = hice_out
    hs = np.zeros(ncat, dtype='float64')
    hi = np.zeros(ncat, dtype='float64')
    qin = np.zeros((ncat,nilyr), dtype='float64')
    qsn = np.zeros((ncat,nslyr), dtype='float64')
    start = time.time()
    for k in range(sw.size):
        if k % 10000 == 0:
           print ' ........ ', k 
           end = time.time()
           print("Elapsed (after compilation) = %s" % (end - start))
           start = end  
        i,j = sw.gi[k]-1,sw.gj[k]-1
        if LAT[j,i] > 60.0:
           aicenpm5[0, k] = aicenpm[1,j,i]           + aicenpm[2,j,i] * fac1
           aicenpm5[1, k] = aicenpm[2,j,i]*(1.-fac1) + aicenpm[3,j,i] * fac2
           aicenpm5[2, k] = aicenpm[3,j,i]*(1.-fac2) + aicenpm[4,j,i] * fac3
           aicenpm5[3, k] = aicenpm[4,j,i]*(1.-fac3) + aicenpm[5,j,i] * fac4
           aicenpm5[4, k] = aicenpm[5,j,i]*(1.-fac4) + np.sum(aicenpm[6:,j,i])
           vicenpm5[0, k] = aicenpm[1,j,i]*piomas_ic[1]           + aicenpm[2,j,i]*piomas_ic[2]*fac1
           vicenpm5[1, k] = aicenpm[2,j,i]*piomas_ic[2]*(1.-fac1) + aicenpm[3,j,i]*piomas_ic[3]*fac2
           vicenpm5[2, k] = aicenpm[3,j,i]*piomas_ic[3]*(1.-fac2) + aicenpm[4,j,i]*piomas_ic[4]*fac3
           vicenpm5[3, k] = aicenpm[4,j,i]*piomas_ic[4]*(1.-fac3) + aicenpm[5,j,i]*piomas_ic[5]*fac4
           vicenpm5[4, k] = aicenpm[5,j,i]*piomas_ic[5]*(1.-fac4) + np.sum(aicenpm[6:,j,i]*piomas_ic[6:])
           if np.sum(aicenpm5[:,k])-1. > 5.e-7:
              exc = np.sum(aicenpm5[:,k])-1.
              for n in range(ncat):
                  if aicenpm5[n,k] > exc:
                     aicenpm5[n,k] -= exc  
                     break
           for n in range(ncat):
               hs[n] = vsnon[n,k]/aicen[n,k] if aicen[n,k] > 0.0 else 0.0
               qsn[n,:] = esnon[n,:,k]*nslyr/vsnon[n,k] if vsnon[n,k] > 0.0 else 0.0
               qin[n,:] = eicen[n,:,k]*nilyr/vicen[n,k] if vicen[n,k] > 0.0 else 0.0
               aicen[n,k] = aicenpm5[n,k] 
               if vicenpm5[n,k] > 0.0 and vicen[n,k] > 0.0:
                   vicen[n,k] = vicenpm5[n,k] 
                   eicen[n,:,k] = qin[n,:]*vicen[n,k]/nilyr
               elif vicenpm5[n,k] > 0.0 and vicen[n,k] == 0.0:
                   vicen[n,k] = vicenpm5[n,k] 
                   upper_open = True
                   lower_open = True
                   for nu in range(n+1,ncat):
                       if qin[nu,0] != 0.0:
                          upper_open = False
                          break
                   for nl in range(n-1,-1,-1):
                       if qin[nl,0] != 0.0:
                          lower_open = False
                          break
                   if upper_open and lower_open:
                       qi0 = -rhoi * (cp_ice*(Tmlt[:-1]-Tf) 
                              + Lfresh*(1.-Tmlt[:-1]/Tf) - cp_ocn*Tmlt[:-1])
                       eicen[n,:,k] = qi0 * vicen[n,k] / nilyr
                       tskin[n,k] = Tf + Tice
                   elif upper_open:
                       eicen[n,:,k] = qin[nl,:]*vicen[n,k]/nilyr
                       tskin[n,k] = tskin[nl,k]
                   else:
                       eicen[n,:,k] = qin[nu,:]*vicen[n,k]/nilyr
                       tskin[n,k] = tskin[nu,k]
               elif vicenpm5[n,k] == 0.0 and vicen[n,k] > 0.0:
                   aicen[n,k] = 0.0
                   vicen[n,k] = 0.0
                   vsnon[n,k] = 0.0
                   tskin[n,k] = Tice
                   eicen[n,:,k] = 0.0
                   esnon[n,:,k] = 0.0
               if aicen[n,k] > 0.0:
                   hin = vicen[n,k] / aicen[n,k]   
                   hsn = (rhow - rhoi) / rhos * hin
                   if hs[n] > hsn:
                      hs[n] = hsn
                   vsnon[n,k] = hs[n]*aicen[n,k] 
                   if vsnon[n,k] > 0.0:
                      esnon[n,:,k] = qsn[n,:]*vsnon[n,k]/nslyr 
                   else:
                      esnon[n,:,k] = 0.0
               else:
                   vsnon[n,k] = 0.0
                   esnon[n,:,k] = 0.0  
                          

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
m.pcolormesh(x,y,hice[levl],cmap=cmp,vmin=0.0, vmax=1.0)

#m.pcolormesh(x,y,hice,cmap=cmp,vmin=250, vmax=280)
#if POLE == 'N':
 #  m.plot(0.0,90.0,'ko',markersize=15, latlon=True)
m.drawparallels(np.arange(-90.,120.,15.),labels=[1,0,0,0]) # draw parallels
m.drawmeridians(np.arange(0.,420.,30.),labels=meridians) # draw meridians
plt.colorbar(orientation='vertical',extend='both',shrink=0.8)
#plt.tight_layout()
plt.show()
 

#print clon[0,:20]


