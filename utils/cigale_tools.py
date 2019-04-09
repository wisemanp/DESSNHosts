# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import subprocess
import glob
import matplotlib.pyplot as plt
import seaborn as sns
import logging
import os
import argparse
import multiprocessing
import astropy.io.fits as fits
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy import units as u
from astroquery.irsa_dust import IrsaDust
from shutil import copyfile

bands = ['g','r','i','z']

def plot_host_cmd(df, miika_class,msize=1,alpha=0.1,plot=True):
    from astropy.cosmology import FlatLambdaCDM
    from astropy.cosmology import Planck15 as cosmo
    sn_names = []
    sn_app_g,sn_app_r = [],[]
    sn_abs_g,sn_abs_r = [],[]
    sn_zs = []
    agn_names = []
    agn_app_g,agn_app_r= [],[]
    agn_abs_g,agn_abs_r = [],[]
    agn_zs=[]
    for i in df.index:
        if df['SPECZ'].loc[i]>0:
            mclass = miika_class[miika_class['TRANSIENT_NAME']==df['TRANSIENT_NAME'].loc[i]]

            if mclass['Est'].values ==1:
                sn_app_g.append(df['MAG_AUTO_G'].loc[i])
                sn_app_r.append(df['MAG_AUTO_R'].loc[i])
                sn_abs_g.append(df['MAG_AUTO_G'].loc[i]-5*(np.log10(cosmo.luminosity_distance(z=df['SPECZ'].loc[i]).value*10**6)-1))
                sn_abs_r.append(df['MAG_AUTO_R'].loc[i]-5*(np.log10(cosmo.luminosity_distance(z=df['SPECZ'].loc[i]).value*10**6)-1))
                sn_names.append(df['TRANSIENT_NAME'].loc[i])
                sn_zs.append(df['SPECZ'].loc[i])
            else:
                agn_app_g.append(df['MAG_AUTO_G'].loc[i])
                agn_app_r.append(df['MAG_AUTO_R'].loc[i])
                agn_abs_g.append(df['MAG_AUTO_G'].loc[i]-5*(np.log10(cosmo.luminosity_distance(z=df['SPECZ'].loc[i]).value*10**6)-1))
                agn_abs_r.append(df['MAG_AUTO_R'].loc[i]-5*(np.log10(cosmo.luminosity_distance(z=df['SPECZ'].loc[i]).value*10**6)-1))
                agn_names.append(df['TRANSIENT_NAME'].loc[i])
                agn_zs.append(df['SPECZ'].loc[i])
    print ('%s objects have redshift, out of %s total objects'%(len(sn_zs)+len(agn_zs),len(df)))
    agn_zs = np.array(agn_zs)
    sn_zs = np.array(sn_zs)
    if plot:
        f,ax = plt.subplots()
        ax.scatter(np.array(sn_abs_r),np.array(sn_app_g)-np.array(sn_app_r),marker='.',s=(msize*6*sn_zs)**2,color='g',alpha=alpha,label='SN Hosts')
        ax.scatter(np.array(agn_abs_r),np.array(agn_app_g)-np.array(agn_app_r),marker='.',s=(msize*6*agn_zs)**2,color='r',alpha=alpha,label='AGN Hosts')
        leg=ax.legend()
        for lh in leg.legendHandles:
            lh.set_alpha(1)
        xs = np.linspace(-26,-18,100)
        ys = -0.375*xs -7.7
        ax.plot(xs,ys,'-k')
        ax.invert_xaxis()
        ax.set_xlabel('$M_r$',size=18)
        ax.set_ylabel('$g - r$',size=18)
        ax.set_xlim(-18,-27)
        ax.set_ylim(-0.5,4)
        f.subplots_adjust(bottom=0.13,top=0.98)


        f2,ax2 = plt.subplots()
        ax2.scatter(np.array(sn_abs_r),np.array(sn_app_g)-np.array(sn_app_r),marker='.',s=(msize*6*sn_zs)**2,color='g',alpha=alpha,label='SN Hosts')
        leg=ax2.legend()
        for lh in leg.legendHandles:
            lh.set_alpha(1)
        xs = np.linspace(-26,-18,100)
        ys = -0.375*xs -7.7
        ax2.plot(xs,ys,'-k')
        ax2.invert_xaxis()
        ax2.set_xlabel('$M_r$',size=18)
        ax2.set_ylabel('$g - r$',size=18)
        ax2.set_xlim(-18,-27)
        ax2.set_ylim(-0.5,4)
        f2.subplots_adjust(bottom=0.13,top=0.98)


        f3,ax3 = plt.subplots()
        ax3.scatter(np.array(agn_abs_r),np.array(agn_app_g)-np.array(agn_app_r),marker='.',s=(msize*6*agn_zs)**2,color='r',alpha=alpha,label='AGN Hosts')
        leg=ax3.legend()
        for lh in leg.legendHandles:
            lh.set_alpha(1)
        xs = np.linspace(-26,-18,100)
        ys = -0.375*xs -7.7
        ax3.plot(xs,ys,'-k')
        ax3.invert_xaxis()
        ax3.set_xlabel('$M_r$',size=18)
        ax3.set_ylabel('$g - r$',size=18)
        ax3.set_xlim(-18,-27)
        ax3.set_ylim(-0.5,4)
        f3.subplots_adjust(bottom=0.13,top=0.98)
    galaxies = []
    for counter,mag in enumerate(sn_abs_r):
        col = sn_app_g[counter]-sn_app_r[counter]


        if col<-0.375*mag -7.7:
            pass
        else:
            galaxies.append(sn_names[counter])
    return galaxies


def prep_cigale_data(name_df = None,sn_name_fn='/home/wiseman/code/des_stacks/source_lists/all_transients.txt',dered=True,ml=False,fz = None):
    print ('sn_name_fn',sn_name_fn)
    if os.path.isfile(sn_name_fn):
        sn_names = np.genfromtxt(sn_name_fn,dtype=str,delimiter='\n')
        name_df = pd.DataFrame(sn_names,columns=['TRANSIENT_NAME'])
    print (name_df)
    latest = pd.read_csv('/media/data3/wiseman/des/coadding/results/sngals_deep_v4.csv')
    latest.reset_index(inplace=True,drop=True)
    print(latest.head(10))
    latest = latest.merge(name_df,on='TRANSIENT_NAME',how='inner')
    dlr1s = latest[(latest['DLR_RANK']==1)|(latest['DLR_RANK']==-1)]
    #add redshifts from SNSPECT
    snspect = pd.read_csv('/media/data3/wiseman/des/coadding/catalogs/snspect.csv',index_col=0)
    for i in range(len(dlr1s)):
        if dlr1s['SPECZ'].iloc[i]>0:
            pass
        else:

            sn = dlr1s.iloc[i]['TRANSIENT_NAME']
            spec_obs = snspect[snspect['TRANSIENT_NAME']==sn]
            done=False
            for z in spec_obs['Z_GAL'].values:
                if z>0:
                    dlr1s['SPECZ'].iloc[i] = z
                    print (dlr1s['TRANSIENT_NAME'].iloc[i])
                    print ('Added z of %s from GAL to %s'%(z,sn))
                    done=True
                    break
                done = False
            if done == False:
                for z in spec_obs['Z_SN'].values:
                    if z >0:
                        dlr1s['SPECZ'].iloc[i] = z
                        print (dlr1s['TRANSIENT_NAME'].iloc[i])
                        print ('Added z of %s from SN to %s'%(z,sn))
                        done=True
                        break
    if fz:
        force_redshifts = np.loadtxt(fz,dtype='str')
        for counter,sn in enumerate(force_redshifts[:,0]):
            snloc = dlr1s[dlr1s['TRANSIENT_NAME']==sn]
            dlr1s['SPECZ'].loc[snloc] = float(force_redshifts[counter,1])
            print ('Forced z of %s for %s as requested'%(force_redshifts[counter,1],sn))
    if ml ==True:
        miika_class = pd.read_csv('f/media/data1/pursiainen/agn_sn_classifier_all_des/31_agn_transient_plot_final.dat',sep='\t',skiprows=5,names=['TRANSIENT_NAME','Est','Per 0','Per 1'])
        miika_testset = pd.read_csv('/media/data1/pursiainen/agn_sn_classifier_all_des/31_agn_transient_train-and-test.dat',sep='\t',skiprows=5,names=['TRANSIENT_NAME','Typ','Est','Per 0','Per 1'])
        miika_class = miika_class.append(miika_testset)
        galaxies = plot_host_cmd(dlr1s,miika_class,plot=False)
        galdf = pd.DataFrame(galaxies,columns=['TRANSIENT_NAME'])


        allgals = dlr1s[dlr1s['SPECZ']>0].merge(galdf,on='TRANSIENT_NAME',how='inner')
    else:
        allgals = dlr1s[dlr1s['SPECZ']>0]
    allgals = allgals[['TRANSIENT_NAME','SPECZ','Y_IMAGE',
                       'RA','DEC',
                       'MAG_AUTO_G','MAGERR_AUTO_G',
                       'MAG_AUTO_R','MAGERR_AUTO_R',
                       'MAG_AUTO_I','MAGERR_AUTO_I',
                       'MAG_AUTO_Z','MAGERR_AUTO_Z',
                       'FLUX_AUTO_G','FLUXERR_AUTO_G',
                       'FLUX_AUTO_R','FLUXERR_AUTO_R',
                       'FLUX_AUTO_I','FLUXERR_AUTO_I',
                       'FLUX_AUTO_Z','FLUXERR_AUTO_Z',
                       'MAG_ZEROPOINT_G','MAG_ZEROPOINT_ERR_G',
                       'MAG_ZEROPOINT_R','MAG_ZEROPOINT_ERR_R',
                       'MAG_ZEROPOINT_I','MAG_ZEROPOINT_ERR_I',
                       'MAG_ZEROPOINT_Z','MAG_ZEROPOINT_ERR_Z']]
    print ('Going to work on %s galaxies!'%len(allgals))
    dered_suffix=''
    if dered==True:
        dered_suffix='_dered'
        coords = SkyCoord(allgals['RA'].values*u.deg,allgals['DEC'].values*u.deg)

        vals_dict  = {'g':[],'r':[],'i':[],'z':[]}
        for i in range(len(allgals)):
            ext = IrsaDust.get_extinction_table(coords[i])
            e = ext.to_pandas()
            e.set_index('Filter_name',drop=True,inplace=True)
            for b in bands:
                vals_dict[b].append(e['A_SandF'].loc['SDSS %s'%b])
        for b in bands:
            allgals['MAG_AUTO_%s'%b.capitalize()] -= vals_dict[b]

    for b in bands:
        allgals['FLUX_AUTO_mJy_%s'%b] = ''
        allgals['FLUX_AUTO_mJy_%s'%b] = ''
        allgals['FLUX_AUTO_mJy_%s'%b] =(10**3)*10**(3.56-(allgals['MAG_ZEROPOINT_%s'%b.capitalize()]/2.5))*allgals['FLUX_AUTO_%s'%b.capitalize()]
        allgals['FLUXERR_AUTO_mJy_%s'%b] = allgals['FLUX_AUTO_mJy_%s'%b].values*((2.303*allgals['MAG_ZEROPOINT_ERR_%s'%b.capitalize()]/2.5)**2 +\
    (allgals['FLUXERR_AUTO_%s'%b.capitalize()]/allgals['FLUX_AUTO_%s'%b.capitalize()])**2)**0.5
    for_cigale = allgals.rename(columns={'TRANSIENT_NAME':'id',
                                         'SPECZ':'redshift',
                                         'FLUX_AUTO_mJy_g':'decam_g',
                                         'FLUX_AUTO_mJy_r':'decam_r',
                                         'FLUX_AUTO_mJy_i':'decam_i',
                                         'FLUX_AUTO_mJy_z':'decam_z',
                                         'FLUXERR_AUTO_mJy_g':'decam_g_err',
                                         'FLUXERR_AUTO_mJy_r':'decam_r_err',
                                         'FLUXERR_AUTO_mJy_i':'decam_i_err',
                                         'FLUXERR_AUTO_mJy_z':'decam_z_err'
                                       })

    for_cigale = for_cigale[['id','redshift','decam_g','decam_r','decam_i','decam_z',
                             'decam_g_err','decam_r_err','decam_i_err','decam_z_err']]

    for_cigale.to_csv('/media/data3/wiseman/des/cigale/cigale-v2018.0/pcigale/%s.dat'%(os.path.split(sn_name_fn)[-1].split('.')[0]+dered_suffix),sep=' ',index=False)

def cigale_config(data_fn):
    from configobj import ConfigObj
    config = ConfigObj('/media/data3/wiseman/des/cigale/cigale-v2018.0/pcigale/pcigale.ini',encoding='utf8')
    config['data_file'] = '%s'%data_fn
    config.write()
    return

def run_cigale(outdir):
    os.chdir('/media/data3/wiseman/des/cigale/cigale-v2018.0/pcigale/')
    proc = subprocess.Popen(['pcigale','check'],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    outs,errs = proc.communicate()
    print(outs)
    proc = subprocess.Popen(['pcigale','run'],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    outs,errs = proc.communicate()
    if not os.path.isdir('/media/data3/wiseman/des/cigale/cigale-v2018.0/pcigale/%s'%outdir):
        os.mkdir('/media/data3/wiseman/des/cigale/cigale-v2018.0/pcigale/%s'%outdir)
    copyfile('/media/data3/wiseman/des/cigale/cigale-v2018.0/pcigale/out/results.fits',
              '/media/data3/wiseman/des/cigale/cigale-v2018.0/pcigale/%s/results.fits'%outdir)
    copyfile('/media/data3/wiseman/des/cigale/cigale-v2018.0/pcigale/%s.dat'%outdir,'/media/data3/wiseman/des/cigale/cigale-v2018.0/pcigale/%s/cig_in.dat'%outdir)
