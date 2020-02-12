#!/usr/bin/env python

import numpy as n
import glob
import h5py
import matplotlib.pyplot as plt
import cfi_config as c

#alpha_norm               Dataset {56452}
#braggs                   Dataset {56452, 3}
#dcos                     Dataset {56452, 2}
#dh                       Dataset {SCALAR}
#dop_errs                 Dataset {56452}
#dops                     Dataset {56452}
#dt                       Dataset {SCALAR}
#heights                  Dataset {56452}
#lats                     Dataset {56452}
#link                     Dataset {56452}
#lons                     Dataset {56452}
#rgs                      Dataset {30}
#t                        Dataset {56452}
#times                    Dataset {48}
#v                        Dataset {2, 48, 30}
#v_resid                  Dataset {56452}
#ve                       Dataset {2, 48, 30}

keys=["alpha_norm","braggs","dcos","dh","dop_errs","dops",
      "dt","heights","lats","link","lons","rgs","t","rgs",
      "t","times","v","v_resid","ve"]


class mmaria_data:
    def __init__(self,dname,debug=False):
        self.fl = glob.glob("%s/*.h5"%(dname))
        self.fl.sort()
        self.mint=[]
        self.maxt=[]
        self.debug=debug
        
        for f in self.fl:
            if self.debug:
                print("reading %s"%(f))
            h=h5py.File(f,"r")
            t=h["t"].value
            self.mint.append(n.min(t))
            self.maxt.append(n.max(t))
            h.close()
        self.mint=n.array(self.mint)
        self.maxt=n.array(self.maxt)        

    def get_bounds(self):
        return([n.min(self.mint),n.max(self.maxt)])

    def read_data(self,t0,t1):
        """
        Read all meteor radar network data between these times (unix)
        """
        # first file
        first_idx=n.where( (self.mint < t0) & (self.maxt > t0)) [0]
        print(first_idx)
        if len(first_idx) == 0:
            print("no data")
            return(None)
        
        # last file
        last_idx=n.where( (self.mint < t1) & (self.maxt > t1)) [0]
        print(last_idx)
        if len(last_idx) == 0:
            print("no data")
            return(None)

        alpha_norm=n.zeros([0],dtype=n.float32)
        braggs=n.zeros([0,3],dtype=n.float32)
        dcos=n.zeros([0,2],dtype=n.float32)
        dh=1.5
        dop_errs=n.zeros([0],dtype=n.float32)
        dops=n.zeros([0],dtype=n.float32)
        dt=0
        heights=n.zeros([0],dtype=n.float32)
        lats=n.zeros([0],dtype=n.float32)
        lons=n.zeros([0],dtype=n.float32)        
        link=n.zeros([0],dtype="<U60")
        rgs=n.zeros([30],dtype=n.float32)
        t=n.zeros([0],dtype=n.float32)
        times=n.zeros([0],dtype=n.float32)
        v=n.zeros([2,0,30],dtype=n.float32)
        v_resid=n.zeros([0],dtype=n.float32)
        ve=n.zeros([2,0,30],dtype=n.float32)

        for i in range(first_idx,last_idx+1):
            h=h5py.File(self.fl[i],"r")
            didx=n.where( ((h["t"].value) > t0) & ((h["t"].value) < t1))[0]
            t = n.concatenate((t,h["t"].value[didx]))
            alpha_norm = n.concatenate((alpha_norm,h["alpha_norm"].value[didx]))
            braggs = n.concatenate((braggs,h["braggs"].value[didx,:]))

            dcos = n.concatenate((dcos,h["dcos"].value[didx,:]))
            dh=h["dh"].value
            dop_errs = n.concatenate((dop_errs,h["dop_errs"].value[didx]))
            dops = n.concatenate((dops,h["dops"].value[didx]))
            dt=h["dt"].value
            heights=n.concatenate((heights,h["heights"].value[didx]))
            lats=n.concatenate((lats,h["lats"].value[didx]))
            lons=n.concatenate((lons,h["lons"].value[didx]))
            link=n.concatenate((link,h["link"].value[didx]))
            v_resid=n.concatenate((v_resid,h["v_resid"].value[didx]))
            
            # mean horizontal wind model
            rgs=h["rgs"].value
            times=n.concatenate((times,h["times"].value))
            v=n.concatenate((v,h["v"].value),axis=1)
            ve=n.concatenate((ve,h["ve"].value),axis=1)
            
            if self.debug:
                print("file idx %d"%(i))
                
        return({"t":t,
                "alpha_norm":alpha_norm,
                "braggs":braggs,
                "dcos":dcos,
                "dh":dh,
                "dop_errs":dop_errs,
                "dops":dops,
                "dt":dt,
                "heights":heights,
                "lats":lats,
                "lons":lons,
                "link":link,
                "rgs":rgs,
                "v_resid":v_resid,
                "times":times,
                "v":v,
                "ve":ve})



# directory with all mmaria network data
md=mmaria_data(c.data_directory)
# what is the data bounds (first and last time stamp)
print(md.get_bounds())

# read all meteor radar data between these two timestamps
d=md.read_data(1514774804,1514974804)
plt.pcolormesh(d["times"],d["rgs"]/1e3,n.transpose(d["v"][0,:,:]),vmin=-100,vmax=100)
plt.xlabel("Time (unix)")
plt.ylabel("Altitude (km)")
plt.colorbar()
plt.show()
