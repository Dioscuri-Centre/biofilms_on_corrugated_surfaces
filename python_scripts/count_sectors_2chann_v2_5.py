# this is the version that seems to work best for 2-channel videos

import numpy as np
from argparse import ArgumentParser
import os
import time
import json

#_threshold=80 

def normalize(a,mina,maxa):
    #y=(a-mina)/(maxa-mina) old version
    y=np.log(a-mina+1)/np.log(maxa-mina)
    return np.convolve(y,np.ones(10)/10.,"valid") #  default 10; 5 inceases the number of sectors slightly and 20 decreases it but does not increase relative differences

def analyse_sectors(pix1, pix2):
    # nonzero=(pix1+pix2>0.2) old version
    nonzero=(pix1+pix2>0.5) #  0.25 does not change much
    #psel1=pix1[nonzero]
    #psel2=pix2[nonzero]
    sel=(pix2>pix1)*nonzero
    binary=np.full(len(pix1),0)
    binary[nonzero]=1
    binary[sel]=2
    n1=np.sum(sel)
    n2=np.sum(nonzero)-n1
    #plt.plot(binary)
    #plt.plot((pix1+pix2>0.1))
    #plt.show()
    d=np.diff(binary)
    secs=np.concatenate(([0],np.argwhere(d!=0)[:,0],[len(d)]))

    # -------the following code is adapted from Alistair Miles, alimanfoo/find_runs.py-------
    n=len(binary)
    loc_run_start = np.empty(n, dtype=bool)
    loc_run_start[0] = True
    np.not_equal(binary[:-1], binary[1:], out=loc_run_start[1:])
    run_starts = np.nonzero(loc_run_start)[0]
    run_values = binary[loc_run_start] # find run values
    run_lengths = np.diff(np.append(run_starts, n)) # find run lengths
    # ----------------
    
    ns1=len(run_lengths[run_values==1])  # number of sectors of type 1
    ns2=len(run_lengths[run_values==2])  # number of sectors of type 2

    return len(secs)-1, n1, n2, ns1,ns2, secs


parser = ArgumentParser()
parser.add_argument("iname") # input must be a two-colour OME TIFF
parser.add_argument("oname")
#parser.add_argument("-t","--threshold",nargs='+', type=float,help="threshold for pixel intensity")
parser.add_argument('-b', '--bottom', action='store_true', help="use the bottom pixelline")  # on/off flag for flatfield correction
args = parser.parse_args()

input_dir_name=args.iname
output_file_name=args.oname

#if (args.threshold):
#    _threshold=args.threshold[0]
#    print("_threshold=",_threshold)


summary=[]
n=0
with open(os.path.join(input_dir_name,'summary.json'),"r") as f:
    for line in f:
        summary.append(json.loads(line))
        n=n+1
#print(summary[0]["pixels"])

# re-order files so that two colours are even and odd elements of summary
#tmp=(summary[0::2],summary[1::2])
#summary=np.concatenate(tmp)

pixnames=list( dict.fromkeys([summary[i]["pixels"] for i in range(0,n)]) )
diffNs=list( dict.fromkeys([summary[i]["N"] for i in range(0,n)]) )
diffCs=list( dict.fromkeys([summary[i]["C"] for i in range(0,n)]) )
diffTs=list( dict.fromkeys([summary[i]["T"] for i in range(0,n)]) )


if (len(diffCs) != 2):
    raise NameError("Number of channels is not two!")
no_N=len(diffNs) # number of wells per field of view
print("Number of wells per field of view: ",no_N)


# first scan all wells and times to determine min and max signal in each channel:
chmin=[[1e6,1e6,1e6,1e6] for i in range(0,len(diffTs))]
chmax=[[0,0,0,0] for i in range(0,len(diffTs))]

for pi in range(0,len(pixnames)):
    print(pixnames[pi])
    npos=sum([s['pixels']==pixnames[0] for s in summary])
    with open(os.path.join(input_dir_name,'pixelline'+str(pi)+".npy"), 'rb') as f:
        for i in range(0,npos):
            coords=np.load(f) # pi,ti,ci,zi,ni
            if coords[0] != pi:
                raise NameError('P does not agree')
            ti,ci,zi,ni=coords[1],coords[2],coords[3],coords[4]
            pixelav=np.load(f)
            pixelstd=np.load(f)
            if (args.bottom): pixelav=pixelstd # use the bottom pixelline, assuming it is written where normally std dev line would be
            min,max=np.min(pixelav),np.max(pixelav)
            if (min<chmin[ti][ci]):
                chmin[ti][ci]=min
            if (max>chmax[ti][ci]):
                chmax[ti][ci]=max

print("chmins=",chmin)
print("chmaxs=",chmax)
#exit()


fout=open(output_file_name,"w")
fout.write("pi,ti,zi,ni,n1,n2,ns1,ns2,Nsecs,xpos\n")
for pi in range(0,len(pixnames)):
    print(pixnames[pi])
    npos=sum([s['pixels']==pixnames[0] for s in summary])
    with open(os.path.join(input_dir_name,'pixelline'+str(pi)+".npy"), 'rb') as f:
        for i in range(0,npos//(2*no_N)):
            data=[[],[]]
            for j in range(0,2*no_N):
                coords=np.load(f) # pi,ti,ci,zi,ni
                if coords[0] != pi:
                    print(coords[0]," ",pi,end='\r')
                    raise NameError('P does not agree')
                ti,ci,zi,ni=coords[1],coords[2],coords[3],coords[4]
                #print(ni,ti,zi,ci)
                pixelav=np.load(f)
                pixelstd=np.load(f)
                if (args.bottom): pixelav=pixelstd # use the bottom pixelline, assuming it is written where normally std dev line would be
                data[j//no_N].append([coords,pixelav])  # add data to each channel (j//no_N), assuming the last index in the summary file is ni
            #print(data[0][0][0]," ",data[0][1][0])
            #print(data[1][0][0]," ",data[1][1][0])
            # data[c][n][coords,pixelav]: c = channel, n = position within FOV

            for j in range(0,no_N):
                pi,ti,ci1,zi,ni=data[0][j][0]
                _,_,ci2,_,_=data[1][j][0]
                pix1,pix2=normalize(data[0][j][1],chmin[ti][ci1],chmax[ti][ci1]), normalize(data[1][j][1],chmin[ti][ci2],chmax[ti][ci2])
                #print(ti," ",zi," ",ni)

                nosecs, n1,n2, ns1,ns2, xpos = analyse_sectors(pix1,pix2)
                print(ti," ",zi," ",ni,": ",nosecs,end='\r')
                fout.write(str(pi)+","+str(ti)+","+str(zi)+","+str(ni)+","+str(n1)+","+str(n2)+","+str(ns1)+","+str(ns2)+","+str(nosecs))
                for k in range(0,len(xpos)):
                    fout.write(","+str(xpos[k]))            
                fout.write("\n")
                fout.flush()

            #print(pixelstd)
    #pixel_name=summary[i]["pixels"]
    #if pixel_name != old_pixel_name:
    #    old_pixel_name=pixel_name
fout.close()      

