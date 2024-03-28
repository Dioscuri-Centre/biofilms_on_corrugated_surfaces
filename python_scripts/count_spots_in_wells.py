# this program finds bright spots (bacteria in FL mode)
# in regions determined based on Klaudia's manually determined positions of wells
# returned coordinates are in pixels

# positions, sz, channels, and times can be selected as a range of values
# to select only one position use e.g. -p 1 1

pixelsize=0.4172 # for Nikon1, 20x
#pixelsize=0.332 # for Nikon2, 20x


import shutil
import numpy as np
from mmpyreader import mmpyreader
#from scipy.ndimage import zoom
from pycromanager import Bridge
from argparse import ArgumentParser
import os
import time
import json
import cv2 as cv
import scipy.linalg
from scipy.ndimage import gaussian_filter
from skimage.measure import label,regionprops,regionprops_table
import skimage

def get_positions_of_wells(name):
    global pixelsize
    data = np.genfromtxt(name, delimiter=',', skip_header = 1)
    #print(data)
    return np.reshape((1./pixelsize)*data[:,1:3],(len(data)//4,4,-1))


def process_tile(src,threshold,di=0,dj=0):
    h,w=src.shape
    img=src.astype("float32")
#    x=np.tile(np.linspace(0,w-1,w),h)
#    y=np.repeat(np.linspace(0,w-1,w),h)
#    A = np.c_[x, y, np.ones(w*h)]
#    C,_,_,_ = scipy.linalg.lstsq(A, np.reshape(img,-1))    # coefficients
    #print(C)
#    Z = C[0]*x + C[1]*y + C[2]
#    if (verbose>2):
#        plt.rcParams["figure.figsize"] = (8,8)
#        plt.imshow(img)
#        plt.show()
#    img=img-np.reshape(Z,(h,w))

    img=gaussian_filter(img,1)
    img=255*img/np.max(img)

    #imgmin=np.min(img)
    #imgmax=np.max(img)
    #img=(img-imgmin)/(imgmax-imgmin)
    #img=255-img/16
    #img=blur(img)
    #img[img<0.6]=0
    
    thr,img_peaks = cv.threshold(img,threshold,255,cv.THRESH_BINARY)
    label_im=label(img_peaks)
    peaks=regionprops(label_im)
    
    xys=[]
    for p in peaks:
        if (p.area>1 and p.area<1000):
            y,x=p.centroid
            y=int(y)+di
            x=int(x)+dj
            xys.append([x,y,p.area])
    return xys


def detect_bright_spots(src, pos):
    x1=int(np.min([pos[0,0],pos[3,0]]))
    y1=int(np.min([pos[0,1],pos[1,1]]))
    x2=int(np.max([pos[1,0],pos[2,0]]))
    y2=int(np.max([pos[2,1],pos[3,1]]))

    img=src[y1:y2,x1:x2].astype("float32")
    b=20
    img_larger=src[y1-b:y2+b,x1-b:x2+b].astype("float32")
#    h,w=img.shape
    img=gaussian_filter(img,2)   # blur a bit to reduce noise
    img0=gaussian_filter(img_larger,b)[b:-b,b:-b]  
    img=1.*img/(img0)  # remove bkg by averaging with a gaussian profile

    xys=skimage.feature.peak_local_max(img,min_distance=1,threshold_abs =1.05)
    for i in range(0,len(xys)):
        xys[i]=xys[i]+[y1,x1]

    return xys, src[y1:y2,x1:x2], x1, y1
    #return xys, img, x1, y1

def process_2d_slice(pi, ci, ti, zi):
    global store,mm,builder, stack, output_dir_name, summary_file, no_zs
    #if (no_zs>1):
    coords = builder.z(zi).time(ti).channel(ci).stagePosition(pi).build()
    #else:
    #    coords = builder.time(ti).channel(ci).stagePosition(pi).build()
    path0=os.path.join(input_coords_name,"Results"+str(pi)+".csv")
    if (stack._mmstore.hasImage(coords) and os.path.exists(path0)):
            
        print("processing: ",pi,ti,ci,zi)

        imgMM_src = stack.get_img(C=ci,P=pi,T=ti,Z=zi,unpack=False)
        h,w = imgMM_src.getHeight(),imgMM_src.getWidth()
        img = imgMM_src.getRawPixels().reshape(h,w) #returns a numpy array
        imgMM_src._close()

        poss=get_positions_of_wells(path0)
        for ni in range(0,len(poss)):
            print("analysing well",ni)
            spotsxy, image, j0,i0 = detect_bright_spots(img,poss[ni])
            print(len(spotsxy))
            path=os.path.join(output_dir_name,"spots"+str(pi)+".npy")
            with open(path, 'ab') as fout:
                np.save(fout,np.array([pi,ti,ci,zi,ni,i0,j0]))
                np.save(fout,spotsxy)
                np.save(fout, image)
            summary={"posfile":path0,
                     "spots":path, "P":pi,"C":ci,"T":ti,"Z":zi,"N":ni}
            json.dump(summary,summary_file)
            summary_file.write("\n")
            summary_file.flush()
    else:
        if os.path.exists(path0):
            print("not found: ",pi,ti,ci,zi)
        else:
            print("coordinates file not found: ",path0)
    

parser = ArgumentParser()
parser.add_argument("iname",help="input directory with TIFF file")
parser.add_argument("coords",help="input directory with coordinates (each in sep. folder named 0...119")
parser.add_argument("oname",help="output directory")
parser.add_argument("-t","--time",nargs='+', type=int,help="time range [frame_first,frame_last]")
parser.add_argument("-p","--position",nargs='+', type=int,help="position range [position_first,position_last]")
parser.add_argument("-z","--zposition",nargs='+', type=int,help="z position range [z_position_first,z_position_last]")
parser.add_argument("-c","--channel",nargs='+', type=int,help="channel range [first,last]")
parser.add_argument("-ptcz","--coordinates",nargs="+",type=int,help="dimesions for each coordinate: P T C Z")
args = parser.parse_args()

input_dir_name=os.path.abspath(args.iname)
output_dir_name=os.path.abspath(args.oname)
input_coords_name=os.path.abspath(args.coords)
print("Input image directory = "+input_dir_name)
print("Input coordinates directory = "+input_coords_name)
print("Output directory = "+output_dir_name)

# create directory and open summary file for writing
if (not os.path.isdir(output_dir_name)):
    os.mkdir(output_dir_name)

summary_file=open(os.path.join(output_dir_name, "summary.json"),"w")

reader = mmpyreader.MMpyreader(verbose=True)
mm=reader.mm
builder = mm.data().getCoordsBuilder()
print("bridge opened correctly")

stack = reader.load_folder(input_dir_name,True)

axes = stack._metadata.getOrderedAxes()
print(axes)
axnb = list(stack._mmstore.getNextIndex(ax) for ax in axes)
print(axnb)

if ('channel' in axes):
    no_channels=axnb[axes.index('channel')]
else:
    no_channels=1
if ('time' in axes):
    no_times=axnb[axes.index('time')]
else:
    no_times=1
if ('z' in axes):
    no_zs=axnb[axes.index('z')]
else:
    no_zs=1
if ('position' in axes):
    no_positions=axnb[axes.index('position')]
else:
    no_positions=1


if (args.coordinates):
    no_positions=args.coordinates[0]
    no_times=args.coordinates[1]
    no_channels=args.coordinates[2]
    no_zs=args.coordinates[3]


print(no_channels, no_times,no_zs,no_positions)

t1,t2=0,no_times-1
if (args.time):
    t1=max(args.time[0],0)
    t2=min(args.time[1],no_times-1)
    print("only time frames from ",t1," to ",t2)

p1,p2=0,no_positions-1
if (args.position):
    p1=max(args.position[0],0)
    p2=min(args.position[1],no_positions-1)
    print("only positions from ",p1," to ",p2)

z1,z2=0,no_zs-1
if (args.zposition):
    z1=max(args.zposition[0],0)
    z2=min(args.zposition[1],no_zs-1)
    print("only z positions from ",z1," to ",z2)

c1,c2=0,no_channels-1
if (args.channel):
    c1=max(args.channel[0],0)
    c2=min(args.channel[1],no_channels-1)
    print("only channels from ",c1," to ",c2)


tt1=time.time()



# I decided to do it like this so that the algorithm works even if some tiff files are missing, and only say one 
# file for a particular position from a number of positions is in the directory.
for pi in range(p1,p2+1):
    f=open(os.path.join(output_dir_name,"spots"+str(pi)+".npy"), 'wb')
    f.close()
    for ti in range(t1,t2+1):
        for zi in range(z1,z2+1):
            for ci in range(c1,c2+1):
                process_2d_slice(pi,ci,ti,zi)

summary_file.close()
print("finished...")
tt2=time.time()
print("time elapsed: ",str(tt2-tt1))
#reader.bridge.close()


