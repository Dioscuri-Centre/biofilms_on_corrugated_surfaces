# this program finds the trajectories of small, sub-difractive features such as small black or white dots in BR

# positions, sz, channels, and times can be selected as a range of values
# to select only one position use e.g. -p 1 1

#pixelsize=0.4172 # used for initial feature position calculation and for the final output trajectories
pixelsize=0.21 # 40x, nikon 1
dt=1/3600. # inter-frame time separation in [h]

import shutil
import numpy as np
from mmpyreader import mmpyreader
from pycromanager import Bridge
from argparse import ArgumentParser
from scipy.optimize import curve_fit
import os
import time
import json

def get_positions_of_features(name):
    global pixelsize
    data = np.genfromtxt(name, delimiter=',', skip_header = 1)
    return (1./pixelsize)*data[:,5:7]

def gaussian2d(xy,amp,cx,cy,s2,bkg):
    x,y=xy
    fun= bkg + amp*np.exp(-((x-cx)**2 + (y-cy)**2)/(2*s2) )
    return np.ravel(fun)

reset_gauss=True
p0=0

def gauss_xy(img):
    wy,wx=img.shape
    bkg=np.mean(img)
    global p0,reset_gauss
    if (reset_gauss):
        #p0 = np.array([np.max(img)-bkg, wx/2,wy/2,2,bkg])
        p0 = np.array([img[wy//2,wx//2]-bkg, wx/2,wy/2,2,bkg])
        reset_gauss=False
    x = np.linspace(0, wx, wx)
    y = np.linspace(0, wy, wy)
    x, y = np.meshgrid(x, y) # is the order of x and y okay?
    try:
        popt, pcov = curve_fit(gaussian2d,(x,y), np.ravel(img),p0)
    except RuntimeError:
        popt=p0
        reset_gauss=True
    if (popt[1]<0 or popt[1]>wx or popt[2]<0 or popt[2]>wy): # if cm outside the box
        popt[1],popt[2]=wx/2,wy/2

    if (reset_gauss==False):
        p0=popt
    return pixelsize*np.array([popt[1],popt[2]])  # x,y [in um] of the centre of mass

def xy_coordinates(movie,x1,y1):
    shift=np.array([x1,y1])*pixelsize
    global reset_gauss
    reset_gauss=True
    return np.array([gauss_xy(frame)+shift for frame in movie])
    
def process(dir_name,poss_dir_name):
    global store,mm,builder,output_name

    
    stack = reader.load_folder(dir_name,True)


    axes = stack._metadata.getOrderedAxes()
    print(axes)
    axnb = list(stack._mmstore.getNextIndex(ax) for ax in axes)
    print(axnb)

    no_channels=axnb[axes.index('channel')]
    no_times=axnb[axes.index('time')]
    no_zs=axnb[axes.index('z')]
    no_positions=axnb[axes.index('position')]

    if (no_channels!=1 or no_zs!=1):
        error
    #print(no_channels, no_times,no_zs,no_positions)

    for pi in range(0,no_positions):        
        poss=get_positions_of_features(poss_dir_name+"/Results"+str(pi)+".csv")
        print(poss)

        print("loading ",no_times," frames for position: ",pi)

        #movie,meta=stack.get_substack(P=0,C=0,Z=0)
        movie=[]
        for t in range(0,no_times):
            imgMM_src = stack.get_img(C=0,P=pi,T=t,Z=0,unpack=False)
            h,w = imgMM_src.getHeight(),imgMM_src.getWidth()
            img = imgMM_src.getRawPixels().reshape(h,w) #returns a numpy array
            movie.append(img)
            imgMM_src._close()
        movie=np.array(movie)
        print(movie.shape)


        for ni in range(0,len(poss)):
            x1,y1=int(poss[ni,0])-10,int(poss[ni,1])-10
            x2,y2=int(poss[ni,0])+10,int(poss[ni,1])+10
            print("processing region ",ni,": x1=",x1,", y1=",y1,"\tx2=",x2,",y2=",y2)
            out=xy_coordinates(movie[:,y1:y2,x1:x2],x1,y1)
            print(out)
            
            with open(output_name, 'ba') as fout:
                np.save(fout,np.array([pi,ni,x1,y1,x2,y2]))
                np.save(fout,out) # velocity field
                np.save(fout,movie[0,y1:y2,x1:x2]) # first frame from the movie

    

parser = ArgumentParser()
parser.add_argument("iname",help="directory with MM data (TIFF files)")
parser.add_argument("coords",help="directory with Results*.csv files for all positions")
parser.add_argument("oname",help="output file, typically with an *.npb extension")
#parser.add_argument("-t","--time",nargs='+', type=int,help="time range [frame_first,frame_last]")
#parser.add_argument("-p","--position",nargs='+', type=int,help="position range [position_first,position_last]")
#parser.add_argument("-z","--zposition",nargs='+', type=int,help="z position range [z_position_first,z_position_last]")
#parser.add_argument("-c","--channel",nargs='+', type=int,help="channel range [first,last]")
args = parser.parse_args()

input_dir_name=os.path.abspath(args.iname)
input_coords_dir_name=os.path.abspath(args.coords)
output_name=os.path.abspath(args.oname)

print("Input directory = "+input_dir_name)
print("Input directory for coordinates (Results*.csv) = "+input_coords_dir_name)
print("Output file = "+output_name)
fout=open(output_name, 'bw')
fout.close()

reader = mmpyreader.MMpyreader(verbose=True)
mm=reader.mm
builder = mm.data().getCoordsBuilder()
print("bridge opened correctly")


tt1=time.time()

if (os.path.isdir(input_dir_name) is False):
    print("input image is not a directory")
    exit()
if (os.path.isdir(input_coords_dir_name) is False):
    print("input coordinates is not a directory")
    exit()

process(input_dir_name,input_coords_dir_name)

#summary_file.close()
print("finished...")
tt2=time.time()
print("time elapsed: ",str(tt2-tt1))
#reader.bridge.close()


