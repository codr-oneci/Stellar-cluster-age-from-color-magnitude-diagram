import numpy as np
import sep
# additional setup for reading the test image and displaying plots
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib import rcParams
#matplotlib inline
import isochrones
from isochrones.priors import ChabrierPrior
from isochrones.mist import MIST_Isochrone
from isochrones import get_ichrone
#from isochrones.dartmouth import Dartmouth_Isochrone

#prima Padova 20Myr, Pas 20Myr, Ultima 140Myr
#Read the Padova isochrones i and r magnitudes

f_iso_20 = open("Padova_20myr.txt", "r")
data_iso_20=[]
#print(f_iso_20.readline())
for line in f_iso_20:
    data_iso_20.append(line.split("\t"))
iso_20_r=[]
iso_20_i=[]
iso_20_r_i=[]
for microlist in data_iso_20[1:-1]:
    str=microlist[0].split()
    iso_20_r.append(float(str[-3]))
    iso_20_i.append(float(str[-2]))
    iso_20_r_i.append(float(str[-3])-float(str[-2]))
#print(iso_20_i)
#print(iso_20_r)
#plt.plot(iso_20_r_i,iso_20_i)
#plt.show()


f_iso_40 = open("Padova_40myr.txt", "r")
data_iso_40=[]
#print(f_iso_40.readline())
for line in f_iso_40:
    data_iso_40.append(line.split("\t"))
iso_40_r=[]
iso_40_i=[]
iso_40_r_i=[]
#print(data_iso_80[1:-1][0][0].split())
for microlist in data_iso_40[1:-1]:
    str=microlist[0].split()
    iso_40_r.append(float(str[-3]))
    iso_40_i.append(float(str[-2]))
    iso_40_r_i.append(float(str[-3])-float(str[-2]))
#print(iso_40_i)
#print(iso_40_r)
#plt.plot(iso_40_r_i,iso_40_i)
#plt.show()

f_iso_60 = open("Padova_60myr.txt", "r")
data_iso_60=[]
#print(f_iso_60.readline())
for line in f_iso_60:
    data_iso_60.append(line.split("\t"))
iso_60_r=[]
iso_60_i=[]
iso_60_r_i=[]
#print(data_iso_80[1:-1][0][0].split())
for microlist in data_iso_60[1:-1]:
    str=microlist[0].split()
    iso_60_r.append(float(str[-3]))
    iso_60_i.append(float(str[-2]))
    iso_60_r_i.append(float(str[-3])-float(str[-2]))
#print(iso_60_i)
#print(iso_60_r)
#plt.plot(iso_60_r_i,iso_60_i)
#plt.show()



f_iso_80 = open("Padova_80myr.txt", "r")
data_iso_80=[]
#print(f_iso_80.readline())
for line in f_iso_80:
    data_iso_80.append(line.split("\t"))
iso_80_r=[]
iso_80_i=[]
iso_80_r_i=[]
#print(data_iso_80[1:-1][0][0].split())
for microlist in data_iso_80[1:-1]:
    str=microlist[0].split()
    iso_80_r.append(float(str[-3]))
    iso_80_i.append(float(str[-2]))
    iso_80_r_i.append(float(str[-3])-float(str[-2]))
#print(iso_80_i)
#print(iso_80_r)
#plt.plot(iso_80_r_i,iso_80_i)
#plt.show()



f_iso_100 = open("Padova_100myr.txt", "r")
data_iso_100=[]
#print(f_iso_100.readline())
for line in f_iso_100:
    data_iso_100.append(line.split("\t"))
iso_100_r=[]
iso_100_i=[]
iso_100_r_i=[]
#print(data_iso_80[1:-1][0][0].split())
for microlist in data_iso_100[1:-1]:
    str=microlist[0].split()
    iso_100_r.append(float(str[-3]))
    iso_100_i.append(float(str[-2]))
    iso_100_r_i.append(float(str[-3])-float(str[-2]))
#print(iso_100_i)
#print(iso_100_r)
#plt.plot(iso_100_r_i,iso_100_i)
#plt.show()

f_iso_120 = open("Padova_120myr.txt", "r")
data_iso_120=[]
#print(f_iso_120.readline())
for line in f_iso_120:
    data_iso_120.append(line.split("\t"))
iso_120_r=[]
iso_120_i=[]
iso_120_r_i=[]
#print(data_iso_80[1:-1][0][0].split())
for microlist in data_iso_120[1:-1]:
    str=microlist[0].split()
    iso_120_r.append(float(str[-3]))
    iso_120_i.append(float(str[-2]))
    iso_120_r_i.append(float(str[-3])-float(str[-2]))
#print(iso_120_i)
#print(iso_120_r)
#plt.plot(iso_120_r_i,iso_120_i)
#plt.show()

f_iso_140 = open("Padova_140myr.txt", "r")
data_iso_140=[]
#print(f_iso_140.readline())
for line in f_iso_140:
    data_iso_140.append(line.split("\t"))
iso_140_r=[]
iso_140_i=[]
iso_140_r_i=[]
#print(data_iso_80[1:-1][0][0].split())
for microlist in data_iso_140[1:-1]:
    str=microlist[0].split()
    iso_140_r.append(float(str[-3]))
    iso_140_i.append(float(str[-2]))
    iso_140_r_i.append(float(str[-3])-float(str[-2]))
#print(iso_140_i)
#print(iso_140_r)
#plt.plot(iso_140_r_i,iso_140_i)
#plt.show()

#The end of reading Padova Isochrones








def in_aperture(x,y):
    '''

    :param x: x coordinate of input point
    :param y: y coordinate of input point
    :return: value 1 if the point x,y is in the circular aperture defined by x0, y0, R
    '''
    x0=2048.0
    y0=2048.0
    R=2048*0.5
    if (x-x0)**2+(y-y0)**2>R**2:
        return 0
    else:
        return 1



rcParams['figure.figsize'] = [10., 8.]
#Process Bias and Dark frames data
#store dark and bias data in data_D and data_B
DarkFiles=["20200921.Dark.bin1x1.40.000secs00000310.fit","20200921.Dark.bin1x1.40.000secs00000311.fit","20200921.Dark.bin1x1.40.000secs00000312.fit",
           "20200921.Dark.bin1x1.40.000secs00000313.fit","20200921.Dark.bin1x1.40.000secs00000314.fit","20200921.Dark.bin1x1.40.000secs00000315.fit",
           "20200921.Dark.bin1x1.40.000secs00000316.fit","20200921.Dark.bin1x1.40.000secs00000317.fit","20200921.Dark.bin1x1.40.000secs00000318.fit",
           "20200921.Dark.bin1x1.40.000secs00000319.fit"]
ok=None
scaler=1.0/float(len(DarkFiles))
for D in DarkFiles:
    hdulist = fits.open(D)
    #hdulist = fits.open("data.001.fits")
    #data = fitsio.read("09-14-20.i'.bin1x1.90.000secs00000256.fits")
    if ok==None:
        data_D=scaler*hdulist[0].data
        ok=1
    else:
        data_D += scaler*hdulist[0].data

BiasFiles=["20200921.Bias.bin1x1.00000320.fit","20200921.Bias.bin1x1.00000321.fit","20200921.Bias.bin1x1.00000322.fit",
           "20200921.Bias.bin1x1.00000323.fit","20200921.Bias.bin1x1.00000324.fit","20200921.Bias.bin1x1.00000325.fit",
           "20200921.Bias.bin1x1.00000326.fit","20200921.Bias.bin1x1.00000327.fit","20200921.Bias.bin1x1.00000328.fit",
           "20200921.Bias.bin1x1.00000329.fit"]
ok=None
scaler=1.0/float(len(BiasFiles))
for B in DarkFiles:
    hdulist = fits.open(B)
    #hdulist = fits.open("data.001.fits")
    #data = fitsio.read("09-14-20.i'.bin1x1.90.000secs00000256.fits")
    if ok==None:
        data_B=scaler*hdulist[0].data
        ok=1
    else:
        data_B += scaler*hdulist[0].data

#process flat field data
RedFlatFiles=["20201001.FlatField.r.bin1x1.00000045.fit","20201001.FlatField.r.bin1x1.00000046.fit","20201001.FlatField.r.bin1x1.00000047.fit",
              "20201001.FlatField.r.bin1x1.00000048.fit","20201001.FlatField.r.bin1x1.00000049.fit","20201001.FlatField.r.bin1x1.00000050.fit",
              "20201001.FlatField.r.bin1x1.00000051.fit","20201001.FlatField.r.bin1x1.00000052.fit","20201001.FlatField.r.bin1x1.00000053.fit",
              "20201001.FlatField.r.bin1x1.00000054.fit","20201001.FlatField.r.bin1x1.00000075.fit","20201001.FlatField.r.bin1x1.00000076.fit",
              "20201001.FlatField.r.bin1x1.00000077.fit","20201001.FlatField.r.bin1x1.00000078.fit","20201001.FlatField.r.bin1x1.00000079.fit",
              "20201001.FlatField.r.bin1x1.00000080.fit","20201001.FlatField.r.bin1x1.00000081.fit","20201001.FlatField.r.bin1x1.00000082.fit",
              "20201001.FlatField.r.bin1x1.00000083.fit","20201001.FlatField.r.bin1x1.00000084.fit"]
ok=None
scaler=1.0/float(len(RedFlatFiles))
for RFF in RedFlatFiles:
    hdulist = fits.open(RFF)
    #hdulist = fits.open("data.001.fits")
    #data = fitsio.read("09-14-20.i'.bin1x1.90.000secs00000256.fits")
    if ok==None:
        data_RFF=scaler*hdulist[0].data
        ok=1
    else:
        data_RFF += scaler*hdulist[0].data
data_RFF=data_RFF/np.median(data_RFF) #normalization of flat field to median of 1.0 ADU
data_RFF=np.repeat(np.repeat(data_RFF,4,axis=0),4,axis=1)

InfraredFlatFiles=["20201001.FlatField.i.bin1x1.00000021.fit","20201001.FlatField.i.bin1x1.00000035.fit","20201001.FlatField.i.bin1x1.00000036.fit",
                   "20201001.FlatField.i.bin1x1.00000037.fit","20201001.FlatField.i.bin1x1.00000038.fit","20201001.FlatField.i.bin1x1.00000039.fit",
                   "20201001.FlatField.i.bin1x1.00000040.fit","20201001.FlatField.i.bin1x1.00000041.fit","20201001.FlatField.i.bin1x1.00000042.fit",
                   "20201001.FlatField.i.bin1x1.00000043.fit","20201001.FlatField.i.bin1x1.00000044.fit","20201001.FlatField.i.bin1x1.00000065.fit",
                   "20201001.FlatField.i.bin1x1.00000066.fit","20201001.FlatField.i.bin1x1.00000067.fit","20201001.FlatField.i.bin1x1.00000068.fit",
                   "20201001.FlatField.i.bin1x1.00000069.fit","20201001.FlatField.i.bin1x1.00000070.fit","20201001.FlatField.i.bin1x1.00000071.fit",
                   "20201001.FlatField.i.bin1x1.00000072.fit","20201001.FlatField.i.bin1x1.00000073.fit","20201001.FlatField.i.bin1x1.00000074.fit"]

ok=None
scaler=1.0/float(len(InfraredFlatFiles))
for IFF in InfraredFlatFiles:
    hdulist = fits.open(IFF)
    #hdulist = fits.open("data.001.fits")
    #data = fitsio.read("09-14-20.i'.bin1x1.90.000secs00000256.fits")
    if ok==None:
        data_IFF=scaler*hdulist[0].data
        ok=1
    else:
        data_IFF += scaler*hdulist[0].data
data_IFF=data_IFF/np.median(data_IFF) #normalization of flat field
data_IFF=np.repeat(np.repeat(data_IFF,4,axis=0),4,axis=1)
#Process Red data
RedFiles=["20200921.r'.bin1x1.40.000secs00000290.fit","20200921.r'.bin1x1.40.000secs00000291.fit",
          "20200921.r'.bin1x1.40.000secs00000292.fit","20200921.r'.bin1x1.40.000secs00000293.fit",
          "20200921.r'.bin1x1.40.000secs00000294.fit","20200921.r'.bin1x1.40.000secs00000295.fit",
          "20200921.r'.bin1x1.40.000secs00000296.fit","20200921.r'.bin1x1.40.000secs00000297.fit",
          "20200921.r'.bin1x1.40.000secs00000298.fit","20200921.r'.bin1x1.40.000secs00000299.fit"]
# read image into standard 2-d numpy array
ok=None
scaler=1.0/float(len(RedFiles))
for R in RedFiles:
    hdulist = fits.open(R)
    #hdulist = fits.open("data.001.fits")
    #data = fitsio.read("09-14-20.i'.bin1x1.90.000secs00000256.fits")
    if ok==None:
        data_R=scaler*hdulist[0].data
        ok=1
    else:
        data_R += scaler*hdulist[0].data

#substract dark and bias data
data_R=data_R-data_B-data_D
#data_R=data_R/data_RFF
#bkg = sep.Background(data, mask=mask, bw=64, bh=64, fw=3, fh=3)
bkg = sep.Background(np.float32(data_R))
# get a "global" mean and noise of the image background:
print(bkg.globalback)
print(bkg.globalrms)

# evaluate background as 2-d array, same size as original image
bkg_image = bkg.back()
# bkg_image = np.array(bkg) # equivalent to above

# show the background
#plt.imshow(bkg_image, interpolation='nearest', cmap='gray', origin='lower')
#plt.colorbar();
#plt.show()

# evaluate the background noise as 2-d array, same size as original image
bkg_rms = bkg.rms()


# show the background noise
#plt.imshow(bkg_rms, interpolation='nearest', cmap='gray', origin='lower')
#plt.colorbar();

# subtract the background
data_sub_R = data_R - bkg

objects_R = sep.extract(data_sub_R, 1.5, err=bkg.globalrms)


from matplotlib.patches import Ellipse

# plot background-subtracted image
fig, ax = plt.subplots()
m, s = np.mean(data_sub_R), np.std(data_sub_R)
#im = ax.imshow(data_sub_R, interpolation='nearest', cmap='gray',vmin=m-s, vmax=m+s, origin='lower')

# plot an ellipse for each object
cluster_ids_R=[]
for i in range(len(objects_R)):
    if in_aperture(objects_R['x'][i],objects_R['y'][i]):
        e = Ellipse(xy=(objects_R['x'][i], objects_R['y'][i]),
                width=6*objects_R['a'][i],
                height=6*objects_R['b'][i],
                angle=objects_R['theta'][i] * 180. / np.pi)
        e.set_facecolor('none')
        e.set_edgecolor('red')
        ax.add_artist(e)
        cluster_ids_R.append(i)

#plt.show()


flux_R, fluxerr_R, flag_R = sep.sum_circle(data_sub_R, objects_R['x'], objects_R['y'],
                                     3.0, err=bkg.globalrms, gain=1.0)

# show the first 10 objects results:
for i in cluster_ids_R:
    print("object {:d}: flux = {:f} +/- {:f}".format(i, flux_R[i], fluxerr_R[i]))

#End of processing Red data



#Process I data
InfraredFiles=["20200921.i'.bin1x1.40.000secs00000300.fit","20200921.i'.bin1x1.40.000secs00000301.fit","20200921.i'.bin1x1.40.000secs00000302.fit",
               "20200921.i'.bin1x1.40.000secs00000303.fit","20200921.i'.bin1x1.40.000secs00000304.fit","20200921.i'.bin1x1.40.000secs00000305.fit",
               "20200921.i'.bin1x1.40.000secs00000306.fit","20200921.i'.bin1x1.40.000secs00000307.fit","20200921.i'.bin1x1.40.000secs00000308.fit",
               "20200921.i'.bin1x1.40.000secs00000309.fit",]
# read image into standard 2-d numpy array
ok=None
scaler=1.0/float(len(InfraredFiles))
# read image into standard 2-d numpy array

for I in InfraredFiles:
    hdulist = fits.open(I)
    if ok==None:
        data_I=scaler*hdulist[0].data
        ok=1
    else:
        data_I += scaler*hdulist[0].data

data_I=data_I-data_B-data_D
#data_I=data_I/data_IFF
#bkg = sep.Background(data, mask=mask, bw=64, bh=64, fw=3, fh=3)
bkg = sep.Background(np.float32(data_I))
# get a "global" mean and noise of the image background:
print(bkg.globalback)
print(bkg.globalrms)

# evaluate background as 2-d array, same size as original image
bkg_image = bkg.back()
# bkg_image = np.array(bkg) # equivalent to above

# show the background
#plt.imshow(bkg_image, interpolation='nearest', cmap='gray', origin='lower')
#plt.colorbar();
#plt.show()

# evaluate the background noise as 2-d array, same size as original image
bkg_rms = bkg.rms()


# show the background noise
#plt.imshow(bkg_rms, interpolation='nearest', cmap='gray', origin='lower')
#plt.colorbar();

# subtract the background
data_sub_I = data_I - bkg

objects_I = sep.extract(data_sub_I, 1.5, err=bkg.globalrms)

# how many objects were detected
len(objects_I)

from matplotlib.patches import Ellipse
# plot background-subtracted image
fig, ax = plt.subplots()
m, s = np.mean(data_sub_I), np.std(data_sub_I)
im = ax.imshow(data_sub_I, interpolation='nearest', cmap='gray',
               vmin=m-s, vmax=m+s, origin='lower')

# plot an ellipse for each object
cluster_ids_I=[]
for i in range(len(objects_I)):
    if in_aperture(objects_I['x'][i],objects_I['y'][i]):
        e = Ellipse(xy=(objects_I['x'][i], objects_I['y'][i]),
                width=6*objects_I['a'][i],
                height=6*objects_I['b'][i],
                angle=objects_I['theta'][i] * 180. / np.pi)
        e.set_facecolor('none')
        e.set_edgecolor('red')
        ax.add_artist(e)
        cluster_ids_I.append(i)
print("arata elipse")
ax.set_xlabel("Pixel coordinate")
plt.show()



flux_I, fluxerr_I, flag_I = sep.sum_circle(data_sub_I, objects_I['x'], objects_I['y'],
                                     3.0, err=bkg.globalrms, gain=1.0)

# show the first 10 objects results:
for i in cluster_ids_I:
    print("object {:d}: flux = {:f} +/- {:f}".format(i, flux_I[i], fluxerr_I[i]))


#End of processing I data


#Source matching between frames

matches=[]

for id_R in cluster_ids_R:
    ideal_match=None
    smallest_distance=None
    for id_I in cluster_ids_I:
        new_dist=(objects_R['x'][id_R]-objects_I['x'][id_I])**2+(objects_R['y'][id_R]-objects_I['y'][id_I])**2
        if smallest_distance==None:
            smallest_distance=new_dist
            ideal_match=id_I
        if new_dist<smallest_distance:
            smallest_distance=new_dist
            ideal_match=id_I
    matches.append((id_R,ideal_match))

print(matches)

mag_I=[]
mag_R=[]
R_I=[]
for star in matches:
    R=-2.5*np.log10(flux_R[star[0]])
    I=-2.5*np.log10(flux_I[star[1]])
    mag_R.append(R)
    mag_I.append(I)
    R_I.append(R-I)

I_mag_correction=22.38 #used Stellarium to find the brightest star in the field:NGC 7654 766 with 8.38mag (V) corresponding to -14 in my plot

plt.plot(R_I, np.array(mag_I)+I_mag_correction, '.')
plt.title("Color Magnitude Diagram for M52",fontsize=13)
plt.xlabel("r'-i' [mag]",fontsize=13)
plt.ylabel("i' [mag]",fontsize=13)
plt.gca().invert_yaxis()
plt.show()

distance_M52=1400 #1.4 kiloparsecs
mag_I_abs=np.array(mag_I)+I_mag_correction+5-5*np.log10(distance_M52)
plt.plot(R_I, mag_I_abs, '.')
plt.title("Color Magnitude Diagram for M52",fontsize=13)
plt.xlabel("r'-i' [mag]",fontsize=13)
plt.ylabel("Absolute i' [mag]",fontsize=13)
plt.gca().invert_yaxis()
plt.show()

#fig, axis = plt.subplots(figsize=(8,6))
#axis.plot(R_I,mag_I_abs,'r.',alpha=0.1)
#axis.invert_yaxis()
#axis.set_xlabel("r'-i'")
#axis.set_ylabel("i'")
#plt.show()

#isozhrone analysis begins here




#iso=isochrones.isochrone.get_ichrone("mist", bands=['i','r'], tracks=False)
#N=1000 #stars in artificial cluster
#masses = ChabrierPrior().sample(N)
#age=9.72
#feh=-0.11
#tracks.generate(mass, age, feh, return_dict=True)
#df = iso.generate(masses, age, feh, distance=1700, AV=0.15) # 1700 pc
#df = df.dropna()
#print(len(df)) # about half of the original simulated stars are nans
#df.head()
#print(df.BP)
#df['BP-RP'] = df.BP - df.RP




# same plotting code as earlier
#Apply color index correction by matching curves:
R_I=np.array(R_I)+0.4956-0.203
fig, axis = plt.subplots(figsize=(8,6))
axis.plot(R_I,mag_I_abs,'b.',alpha=0.2)
axis.invert_yaxis()
axis.set_xlabel("r'-i'", fontsize=13)
axis.set_ylabel("i'",fontsize=13)
axis.set_title('Isochrone Fitting on M52 CMD',fontsize=13)

# now I plot the isochrone
axis.plot(iso_20_r_i,iso_20_i,label="20 Myr")
axis.plot(iso_40_r_i,iso_40_i,label="40 Myr")
axis.plot(iso_60_r_i,iso_60_i,label="60 Myr")
axis.plot(iso_80_r_i,iso_80_i,label="80 Myr")
axis.plot(iso_100_r_i,iso_100_i,label="100 Myr")
axis.plot(iso_120_r_i,iso_120_i,label="120 Myr")
axis.plot(iso_140_r_i,iso_140_i,label="140 Myr")
axis.legend()
plt.show()
