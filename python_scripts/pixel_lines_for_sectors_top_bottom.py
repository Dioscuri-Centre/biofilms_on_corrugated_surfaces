# THIS VERSION IS FOR two strains G and R
# and for creating two lines: top and bottom line, for each well

# this program calculates average pixel intensities and std devs along a line
# the line is determined based on Klaudia's manually determined positions of wells

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

def generate_correction_field_from_coeffs():
    
    # for Nikon #1 experiment 2024_01_12 fig 4
    C=[ 1.  ,       -1.38800465 , 0.99781413 ,-0.40403354 ,-0.1889882 , -1.9386761 ,  6.38929358 ,-0.16358318  ,1.29683599 ,-3.62559667]
    w,h=1920,1200 

    x=np.tile(np.linspace(0,1,w),h)-0.5
    y=np.repeat(np.linspace(0,1,h),w)-0.5
    y2,xy,x2 = y*y,x*y,x*x
    y3,y2x,yx2,x3 = y*y*y,x*y*y,x*x*y,x*x*x
    Z = C[9]*x3 + C[8]*yx2 + C[7]*y2x + C[6]*y3 + C[5]*x2 + C[4]*y2 + C[3]*xy + C[2]*x + C[1]*y + C[0]
    return np.reshape(1./Z,(h,w))    

def flat_field_correction_apply(img,corr):
    img2=img.astype('float32')*corr
    img=img2.astype('int16')
    return img

def get_positions_of_wells(name):
    global pixelsize
    data = np.genfromtxt(name, delimiter=',', skip_header = 1)
    #print(data)
    return np.reshape((1./pixelsize)*data[:,1:3],(len(data)//4,4,-1))

def get_pixels_along_line(img, x1, y1, x2, y2):
    dx=x2-x1
    dy=y2-y1
    l=np.sqrt(dx**2+dy**2)
    dx/=l
    dy/=l
    pix=[img[int(y1+dy*i),int(x1+dx*i)] for i in range(0,int(l))]
    return np.array(pix)

def sectors_to_pixelline(img, poss, ni): # create pixel line data
    pixellines=[]    
    for i in range(0,6):
        c=0.0+0.1*i/5
        x1=(poss[ni,0,0]*c+(1-c)*poss[ni,3,0])
        y1=(poss[ni,0,1]*c+(1-c)*poss[ni,3,1])
        x2=x1+ (poss[ni,1,0]-poss[ni,0,0])
        y2=y1+ (poss[ni,1,1]-poss[ni,0,1])
        pixelline=get_pixels_along_line(img,x1,y1,x2,y2)
        pixellines.append(pixelline)

    pixellines=np.array(pixellines)
    pixelav_top=np.mean(pixellines,axis=0)

    pixellines=[]    
    for i in range(0,6):
        c=0.9+0.1*i/5
        x1=(poss[ni,0,0]*c+(1-c)*poss[ni,3,0])
        y1=(poss[ni,0,1]*c+(1-c)*poss[ni,3,1])
        x2=x1+ (poss[ni,1,0]-poss[ni,0,0])
        y2=y1+ (poss[ni,1,1]-poss[ni,0,1])
        pixelline=get_pixels_along_line(img,x1,y1,x2,y2)
        pixellines.append(pixelline)

    pixellines=np.array(pixellines)
    pixelav_bottom=np.mean(pixellines,axis=0)
    
    return pixelav_top, pixelav_bottom

def process_2d_slice(pi, ci, ti, zi):
    global store,mm,builder, stack, output_dir_name, summary_file, no_zs, args, flat_field_correction
    #if (no_zs>1):
    coords = builder.z(zi).time(ti).channel(ci).stagePosition(pi).build()
    #else:
    #    coords = builder.time(ti).channel(ci).stagePosition(pi).build()
    path0=os.path.join(input_coords_name,"Results"+str(pi)+".csv")
    if (stack._mmstore.hasImage(coords) and os.path.exists(path0)):
            
        print("processing: ",pi,ti,ci,zi,end='\r')

        imgMM_src = stack.get_img(C=ci,P=pi,T=ti,Z=zi,unpack=False)
        h,w = imgMM_src.getHeight(),imgMM_src.getWidth()
        img = imgMM_src.getRawPixels().reshape(h,w) #returns a numpy array
        imgMM_src._close()

        if (args.flatfield):
            img = flat_field_correction_apply(img,flat_field_correction)

        poss=get_positions_of_wells(path0)
        for ni in range(0,len(poss)):
            print("processing: ",pi,ti,ci,zi," analysing well ",ni,end='\r')
            if (pi<30):
                pixelav_top, pixelav_bottom = sectors_to_pixelline(img,poss,ni)
            else:
                pixelav_bottom, pixelav_top = sectors_to_pixelline(img,poss,ni) # swap for the wells on the other side of the channel
            path=os.path.join(output_dir_name,"pixelline"+str(pi)+".npy")
            with open(path, 'ab') as fout:
                np.save(fout,np.array([pi,ti,ci,zi,ni]))
                np.save(fout,pixelav_top)
                np.save(fout,pixelav_bottom)
            summary={"posfile":path0,
                     "pixels":path, "P":pi,"C":ci,"T":ti,"Z":zi,"N":ni}
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
parser.add_argument('-f', '--flatfield', action='store_true', help="use the coefficients hard-coded in the script to do flat-field correction")  # on/off flag for flatfield correction
args = parser.parse_args()

input_dir_name=os.path.abspath(args.iname)
output_dir_name=os.path.abspath(args.oname)
input_coords_name=os.path.abspath(args.coords)
print("Input image directory = "+input_dir_name)
print("Input coordinates directory = "+input_coords_name)
print("Output directory = "+output_dir_name)

if (args.flatfield):
    print("flat field corrections with hard-coded coefficients")
    flat_field_correction=generate_correction_field_from_coeffs()

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
    f=open(os.path.join(output_dir_name,"pixelline"+str(pi)+".npy"), 'wb')
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


