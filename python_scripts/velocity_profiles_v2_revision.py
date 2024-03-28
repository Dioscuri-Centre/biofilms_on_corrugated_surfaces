# this program calculates the velocity profile in given regions of 1s videos
# the regions are determined based on Witek's manually determined positions of wells
# velocities are returned in pixels/h for the original, unbinned image

# positions, sz, channels, and times can be selected as a range of values
# to select only one position use e.g. -p 1 1

pixelsize=0.4172 # used only for well position calculation
dt=1/3600. # inter-frame time separation in [h]

import shutil
import numpy as np
from mmpyreader import mmpyreader
from pycromanager import Bridge
from argparse import ArgumentParser
import os
import time
import json

def get_positions_of_wells(name):
    global pixelsize
    if (os.path.exists(name)):
        data = np.genfromtxt(name, delimiter=',', skip_header = 1)
    #print(data)
        return np.reshape((1./pixelsize)*data[:,1:3],(len(data)//4,4,-1))
    else:
        print("file with coordinates not found: ",name)
        exit()

def blur(im):
    kernel = np.array([1.0,2.0,2.0,2.0,1.0]) # Here you would insert your actual kernel of any size
    blurred = np.apply_along_axis(lambda x: np.convolve(x, kernel, mode='same'), 0, im)
    blurred= np.apply_along_axis(lambda x: np.convolve(x, kernel, mode='same'), 1, blurred)
    return blurred

def velocity_field_from_pair(im0, im1):
    i0=blur(im0)
    i1=blur(im1)
    if (len(i0)!=len(i1) or len(i0[0])!=len(i1[0])):
        raise ValueError("dimensions of i0, i1 unequal")

    dx=16

    w=len(i0[0])
    h=len(i0)

    di=i1-i0
    didy=np.apply_along_axis(lambda x: np.convolve(x, [1/2.,0,-1/2.], mode='same'), 0, i0)
    didx=np.apply_along_axis(lambda x: np.convolve(x, [1/2.,0,-1/2.], mode='same'), 1, i0)
    ny,nx=int(h/dx), int(w/dx)
    vx=np.zeros((ny,nx))
    vy=np.zeros((ny,nx))
    for ii in range(0,ny):
        for jj in range(0,nx):
            i,j=ii*dx,jj*dx
            aa=np.array([(didx[i:i+dx,j:j+dx]).reshape(-1),(didy[i:i+dx,j:j+dx]).reshape(-1)])
            sol=np.linalg.lstsq(aa.transpose(),di[i:i+dx,j:j+dx].reshape(-1), rcond=1e-9)
            vx[int(i/dx),int(j/dx)]=-sol[0][0]        
            vy[int(i/dx),int(j/dx)]=sol[0][1]        

    x,y = np.meshgrid(np.linspace(dx/2.,w-dx/2.,nx),np.linspace(h-dx/2.,dx/2.,ny))
    return (w,h,x,y,vx,vy)



##def get_pixels_along_line(img, x1, y1, x2, y2):
##    dx=x2-x1
##    dy=y2-y1
##    l=np.sqrt(dx**2+dy**2)
##    dx/=l
##    dy/=l
##    pix=[img[int(y1+dy*i),int(x1+dx*i)] for i in range(0,int(l))]
##    return np.array(pix)
##
##def sectors_to_pixelline(img, poss, ni): # create pixel line data
##    pixellines=[]    
##    for i in range(0,11):
##        c=0.4+0.2*i/10
##        x1=(poss[ni,0,0]*c+(1-c)*poss[ni,3,0])
##        y1=(poss[ni,0,1]*c+(1-c)*poss[ni,3,1])
##        x2=x1+ (poss[ni,1,0]-poss[ni,0,0])
##        y2=y1+ (poss[ni,1,1]-poss[ni,0,1])
##        pixelline=get_pixels_along_line(img,x1,y1,x2,y2)
##        pixellines.append(pixelline)
##
##    pixellines=np.array(pixellines)
##    pixelav=np.mean(pixellines,axis=0)
##    pixelstd=np.std(pixellines,axis=0)
##    return pixelav, pixelstd

def velocity_field(movie):
    global dt
    vxs=[]
    vys=[]
    varvxs=[]
    varvys=[]
    coords=None
    dis=[1,2,4,8,16,32]
    for di in dis:
        print("di=",di)
        vecfields=[]
        for i in np.arange(0,len(movie)-di,2):
            w,h,x,y,vx,vy = velocity_field_from_pair(movie[i],movie[i+di])
            coords=(x,y)
    #        vx=vx/(di)
    #        vy=vy/(di)
            vecfields.append((vx,vy))

        vecfields=np.array(vecfields)
        vxs.append(np.mean(vecfields[:,0],axis=0))
        vys.append(np.mean(vecfields[:,1],axis=0))
        varvxs.append(np.var(vecfields[:,0],axis=0))
        varvys.append(np.var(vecfields[:,1],axis=0))
    vxs=np.array(vxs)
    vys=np.array(vys)
    varvxs=np.array(varvxs)
    varvys=np.array(varvys)



    vv=np.transpose(np.sqrt(vxs**2+vys**2),(1,2,0)) # absolute velocities

    vx=np.empty((len(vxs[0]),len(vxs[0,0])))
    vy=np.empty((len(vxs[0]),len(vxs[0,0])))
    for i in range(0,len(vxs[0])):
        for j in range(0,len(vxs[0,0])):
            try:
                k=np.argmax(vv[i,j][vv[i,j]<0.5]) # select the largest velocity that is smaller than a given value (for linearity)
            except ValueError:
                k=5
            vx[i,j]=vxs[k,i,j]/(dt*dis[k])
            vy[i,j]=vys[k,i,j]/(dt*dis[k])

    x,y=coords
    print(x.shape,y.shape,vx.shape,vy.shape)
    return np.array([x,y,vx,vy])


    
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
        poss=get_positions_of_wells(poss_dir_name+"/Results"+str(pi)+".csv")
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
            x1,y1=int(poss[ni,0,0])-25,int(poss[ni,0,1])-40
            x2,y2=int(poss[ni,2,0])+25,int(poss[ni,2,1])+40
    #                else: # for reversed wells
    #                    x1,y1=poss[ni,0,0],poss[ni,0,1]
    #                    x2,y2=poss[ni,2,0],poss[ni,2,1] 
            print("processing region ",ni,": x1=",x1,", y1=",y1,"\tx2=",x2,",y2=",y2)
            out=velocity_field(movie[:,y1:y2,x1:x2])
    #        print(out)
            
    ##        path=os.path.join(output_name,"vel_profile"+str(pi)+"_"+str(ni)+".npy")
            with open(output_name, 'ba') as fout:
                np.save(fout,np.array([pi,ni,x1,y1,x2,y2]))
                np.save(fout,out) # velocity field
                np.save(fout,movie[0,y1:y2,x1:x2]) # first frame from the movie

#            summary={"posfile":path0,
#                     "pixels":path, "P":pi,"C":ci,"T":ti,"Z":zi,"N":ni}
#           json.dump(summary,summary_file)
#            summary_file.write("\n")
#            summary_file.flush()

    

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


