import time
import datetime
import numpy as n
import h5py
import matplotlib.pyplot as plt

# our internal modules
import mmaria_read as mr
import cfi_dac as cfi
import cfi_config as c
import mean_wind_est as mw

md=mr.mmaria_data(c.data_directory)#for many files in a directory
b=md.get_bounds()

def avg_hor_acfs(dcos_thresh=0.8,
                 mean_wind_time_avg=4*3600.0,
                 h0=90.0,
                 dh=5,
                 ds_h=25.0,
                 dtau=300.0,
                 s_h=n.arange(25.0,500.0,25.0),
                 ds_z=1,
                 years=[2018,2019,2020],                 
                 months=[5,6,7],
                 name="summer_hacf",
                 remove_mean=False,
                 n_days=31):
    n_lags=len(s_h)
    all_acfs=[]
    all_errs=[]
    
    n_avg=0.0
    
    for year in years:
        for month in months:
            d0=datetime.date(year,month,1)
            t0=time.mktime(d0.timetuple())
    
            for day in range(n_days):
                d=md.read_data(t0=t0+day*24*3600.0,t1=t0+(day+1)*24*3600.0)
                n_meas=len(d["t"])
                print("n_meteors %d"%(n_meas))
                if n_meas > 100:
                    if remove_mean:
                        times,times_h,v,ve,rgs,lat0,lon0,dt,dh=mw.mean_wind(meas=d, 
                                                                            dt=mean_wind_time_avg,
                                                                            dh=1.0,
                                                                            max_alt=105,
                                                                            min_alt=78,
                                                                            dcos_thresh=dcos_thresh,
                                                                            ofname="res/tmp.h5",
                                                                            data='dict')
                
                    meas=cfi.get_meas(meas_file=d,
                                      mean_rem=remove_mean,
                                      plot_dops=False,
                                      dcos_thresh=dcos_thresh,
                                      mean_wind_file="res/tmp.h5",
                                      data='mmaria')
            
                    #title='Vertical AFC with data from date {}.{}.{} to {}.{}.{}'.format(y,m0,day,y1,m1,day)
#                        return(h0,dtau,ds_h,acfs,errs,shs,s_h,names)

                    ih0,idtau,dis_h,acfs,errs,ishs,si_h,names=cfi.hor_acfs(meas,
                                                              h0=h0,
                                                              dh=dh,
                                                              ds_z=ds_z,
                                                              ds_h=ds_h,
                                                              s_h=s_h,
                                                              dtau=dtau,
                                                              title=name)
                    all_acfs.append(acfs)#+=acfs
                    all_errs.append(errs)
                    n_avg+=1.0
    all_acfs=n.array(all_acfs)
    all_errs=n.array(all_errs)
    err_vars=n.zeros([len(s_h),6])
    acfs=n.zeros([len(s_h),6])
    for i in range(len(s_h)):
        for ci in range(6):
            err_vars[i,ci]=n.var(all_acfs[:,i,ci])
            ws=0.0
            for mi in range(int(n_avg)):
                w=1.0/all_errs[mi,i,ci]
                ws+=w
                acfs[i,ci]+=w*all_acfs[mi,i,ci]
            acfs[i,ci]=(1.0/ws)*acfs[i,ci]
                    
    #print(all_acfs.shape)
    colors=["C0","C1","C2","C3","C4","C5"]
    
    cfi.plot_hor_acfs(shs=s_h,
                  names=names,
                  acfs=acfs,
                  ds_z=ds_z,
                  dtau=dtau,
                  ds_h=ds_h,
                  err_vars=err_vars,
                  colors=colors,
                  n_avg=n_avg)
    
    
    plt.savefig("C:/Users/OleK/Master_thesis/figs/fig_%s.png"%(name))
    ho=h5py.File("res/%s.h5"%(name),"w")
    ho["acf"]=acfs
    ho["s_h"]=s_h
    ho["err_var"]=err_vars/n.sqrt(n_avg)
    ho["h0"]=h0
    ho["dtau"]=dtau
    ho["ds_h"]=ds_h
    ho.close()
    
    
    
    
    # 2018 to 2020
# summer: 5,6,7
# winter: 11,12,1
# fall: 8,9,10
# spring: 2,3,4
#y=2019
#y1=2019
#m0=5
#m1=8
#day=1

avg_hor_acfs(dcos_thresh=0.8,
             mean_wind_time_avg=4*3600.0,
             h0=90.0,
             dh=5,
             ds_h=25.0,
             dtau=300.0,
             s_h=n.arange(25.0,500.0,25.0),
             ds_z=1,
             years=[2018,2019,2020],
             months=[5,6,7],
             name="summer_hacf",
             remove_mean=False,
             n_days=31)

#avg_ver_acfs(dcos_thresh=0.8,
#             mean_wind_time_avg=4*3600.0,
#             h0=80.0,
#             ds_h=100.0,
#             dtau=300.0,
#             s_z=n.arange(0.0,20.0,1.0),
#             years=[2018,2019,2020],
#             months=[11,12,1],
#             name="winter_vacf",
#             n_days=31)
