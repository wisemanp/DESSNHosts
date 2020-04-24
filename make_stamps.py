import astropy
import numpy as np
import seaborn as sns
import pandas as pd
import aplpy
import os
import matplotlib.pyplot as plt
import astropy.io.fits as fits
from astropy.table import Table
from matplotlib.colors import LogNorm
import glob
from astropy.coordinates import SkyCoord
from astropy import units as u
import subprocess
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve
from astropy.wcs import WCS
import configparser
import _pickle as cpickle
from shutil import copyfile

from des_stacks import des_stack as stack
from des_stacks.utils import stack_tools,source_tools
bands = ['g','r','i','z']
def get_sn_radec(ra,dec,field):
    
    f=open('/media/data3/wiseman/des/coadding/config/chiplims.pkl','rb')
    chiplims = cpickle.load(f)

    #################
    obj_field = field
    the_field = chiplims[obj_field]
    for ccd in the_field.keys():
        if the_field[ccd][0][0] > ra > the_field[ccd][2][0]:
            if the_field[ccd][0][1] < dec < the_field[ccd][1][1]:
                return ccd

    
class sn():
    def __init__(self,sn_name=None,ra=None,dec=None,field=None):
        self.sn_name = sn_name
        self.snra=ra
        self.sndec=dec
        self.f=field
        self.peak_fns = {}
        self._get_sn_dat()
        
        
    def _get_sn_dat(self):
        if self.sn_name:
            self.snra,self.sndec,self.f,self.y,self.chip = source_tools.get_sn_dat(self.sn_name)
            print (self.f)
        elif self.snra:
            self.chip = get_sn_radec(self.snra,self.sndec,self.f)
            self.chip=15
            
            self.y = 'none'
            print ('Found chip! %s'%self.chip)
        cuts = {}
        
        for b in bands :
            cuts[b] = stack_tools.get_cuts(self.f,b)
        sg,sr,si,sz = [stack.Stack(self.f, b, self.y, [str(self.chip)] ,
                                   'coadding',cuts[b]) for b in bands]
        for s in [sg,sr,si,sz]:
            s.cuts = cuts[s.band]
        self.desdir = '/media/data3/wiseman/des/coadding/5yr_stacks/CAP/%s/'%self.sn_name
        self.host_imgs= {}
        if not os.path.isdir(self.desdir):
            os.mkdir(self.desdir)
        
        for b in bands:
            cutstring = '%s_%s'%(cuts[b]['teff'],cuts[b]['psf'])
            img_fn = os.path.join(self.desdir,
                    'ccd_%s_%s_%s_clipweighted_sci.resamp.fits'%(self.chip,b,cutstring))
            print (b)
            print('Should have a file called: %s'%img_fn)
            if not os.path.isfile(img_fn):
                print("Can't find %s"%img_fn)
                if not os.path.isfile(os.path.join('/media/data3/wiseman/des/coadding/5yr_stacks',
                                     'MY%s'%self.y,self.f,'CAP',str(self.chip),
                                      'ccd_%s_%s_%s_clipweighted_sci.resamp.fits'%(self.chip,b,cutstring))):
                    print("Resample doesn't exist in CAP folder, trying to make a new one")   
    
                    # if there is no white image, make one
                    det_name = os.path.join(sg.out_dir,'MY%s'%self.y,self.f,'CAP',str(self.chip),
                                            '%s_%s_%s_riz.fits'%(self.y,self.f,self.chip))
                    
                    noff1,noff2 = 0,0
                    while True:
                        det_name,noff1,noff2 = stack_tools.resample_chip_for_cap(sg,sr,si,sz,self.chip,npix_off1=noff1,npix_off2 = noff2)
                        if noff1 == 0 and noff2 == 0:
                            break
                copyfile(os.path.join('/media/data3/wiseman/des/coadding/5yr_stacks',
                                         'MY%s'%self.y,self.f,'CAP',str(self.chip),
                                          'ccd_%s_%s_%s_clipweighted_sci.resamp.fits'%(self.chip,b,cutstring)),
                            img_fn)
                    
            self.host_imgs[b] = os.path.join(self.desdir,
                    'ccd_%s_%s_%s_clipweighted_sci.resamp.fits'%(self.chip,b,cutstring))

        
    
    
        #self.lightcurve = pd.read_csv('/media/data1/wiseman/des/DES13C1feu/DES13C1feu.dat',index_col= 0, sep= '\t')
        #self.lightcurve.dropna(inplace=True)
        
        #cp=configparser.ConfigParser()
        ## read the .ini file
        #cp.read('/media/data3/wiseman/des/coadding/config/snobs_params.ini')
        ## Make a list of years
        #year_mjd_lims =[int(lim.strip()) for lim in cp.get('year_mjd_lims','Y'+str(self.y)).split(',')]
        
        #self.lightcurve = self.lightcurve.loc[(self.lightcurve.index>int(year_mjd_lims[0]))&(self.lightcurve.index<int(year_mjd_lims[1]))]   
        
    def plt_lc(self,graph,n): 
        
        xg,yg= self.lightcurve.index[n*1],self.lightcurve['f_g'].iloc[n*1]
        graph.scatter(xg,yg,color='b',marker='o',edgecolor='b',s=100,zorder=2)
        xr,yr= self.lightcurve.index[n*1],self.lightcurve['f_z'].iloc[n*1]
        graph.scatter(xr,yr,color='r',marker='o',edgecolor='r',s=100,zorder=2)
        xi,yi= self.lightcurve.index[n*1],self.lightcurve['f_i'].iloc[n*1]
        graph.scatter(xi,yi,color='g',marker='o',edgecolor='g',s=100,zorder=2)
        xz,yz= self.lightcurve.index[n*1],self.lightcurve['f_r'].iloc[n*1]
        graph.scatter(xz,yz,color='purple',marker='o',edgecolor='purple',s=100,zorder=2)

        self.lightcurve.plot(y=['f_g','f_r','f_i','f_z'],color=['b','purple','g','r'],ax=graph,linewidth= 2.5)
        graph.axis('off')
        graph.legend().set_visible(False)
        #graph.vlines(56690,x,45)
        #plt.savefig('graph_%s.png'%x,transparent = True)
        #plt.clf()
        
                   
    def insert_star(self,img,star,rad,b,n):
        host = fits.open(img)
        host_head = host[0].header       
        w = WCS(host_head)
        pixcrd = w.wcs_world2pix(np.array([[self.snra,self.sndec]]),1)[0]
        pixra,pixdec = int(pixcrd[0]),int(pixcrd[1])
        host[0].data[pixdec-rad:pixdec+rad+1,pixra-rad:pixra+rad+1]=host[0].data[pixdec-rad:pixdec+rad+1,pixra-rad:pixra+rad+1]+star
        host.writeto(self.host_imgs[b][:-5]+'%s.fits'%n,overwrite=True)
        return self.host_imgs[b][:-5]+'%s.fits'%n   
    
    def make_image_host(self):
        imgs = []
        for b in ['i','r','g']:#['z','i','r']:
            imgs.append(self.host_imgs[b])
            #imgs.append('/media/data3/wiseman/des/coadding/random/DES16C3cje_%s.fits'%b)
        aplpy.make_rgb_cube([imgs[0],
                             imgs[1],
                             imgs[2]],
                            os.path.join(self.desdir,'cube_host.fits'))
        aplpy.make_rgb_image(os.path.join(self.desdir,'cube_host.fits'),
                             os.path.join(self.desdir,'rgb_host.png'),
                             stretch_g='log',stretch_r='log',stretch_b='log',
                             vmin_r=1.,vmax_r=5000.,
                             vmin_b=1.,vmax_b=3000.,
                             vmin_g=1.,vmax_g=5000.)
        
        fig = plt.figure(figsize=(7,7))
        f = aplpy.FITSFigure(os.path.join(self.desdir,'cube_host_2d.fits'),
                             figure=fig,subplot =[0.02,0.02,0.95,0.95])
        h = fits.getheader(os.path.join(self.desdir,'cube_host_2d.fits'))
        f.show_rgb(os.path.join(self.desdir,'rgb_host.png'),vertical_flip=True)
        hor_line = np.array([[self.snra-0.00027,self.snra+0.00027],[self.sndec,self.sndec]])
        ver_line = np.array([[self.snra,self.snra],[self.sndec-0.00027,self.sndec+0.00027]])
        f.show_lines([ver_line,hor_line],color='c',linewidth=3)

        f.recenter(self.snra,self.sndec,radius=0.0027/2)  
        f.tick_labels.hide()
        f.ticks.hide()
        f.add_scalebar(0.0027/2)
        #f.show_circles(34.97413016, -4.95281899,1/3600,color='c')
        f.scalebar.set_color('c')
        f.scalebar.set_label('5"')
        f.scalebar.set_font_size(20)
        f.scaelbar.set_linewidth(4)
        plt.annotate(self.sn_name,xy=(0.1,0.85),xycoords='axes fraction',color='c',size=25)
        #plt.annotate(self.sn_name,xy=(0.1,0.7),xycoords='axes fraction',c='c',size=13)
        plt.savefig(os.path.join(self.desdir,'%s_rgb_host.png'%self.sn_name))
        print ('Saving image at %s'%os.path.join(self.desdir,'%s_rgb_host.png'%self.sn_name))
        
        
    def make_image(self,n):
        imgs = []
        for b in ['z','i','r']:
            imgs.append(self.host_imgs[b][:-5]+'%s.fits'%n )
            
        aplpy.make_rgb_cube([imgs[0],
                             imgs[1],
                             imgs[2]],
                            os.path.join(self.desdir,'cube_%s.fits'%n))
        aplpy.make_rgb_image(os.path.join(self.desdir,'cube_%s.fits'%n),os.path.join(self.desdir,'rgb_%s.png'%n),
                             vmin_r=-1.,vmax_r=15850.,vmin_b=-1.,vmax_b=15850.,vmin_g=-1.,vmax_g=15850.)
        
        fig = plt.figure(figsize=(7,7))
        f = aplpy.FITSFigure(os.path.join(self.desdir,'cube_%s_2d.fits'%n),
                             figure=fig,subplot =[0.02,0.02,0.95,0.95])
        f.show_rgb(os.path.join(self.desdir,'rgb_%s.png'%n))
        f.recenter(self.snra,self.sndec,radius=0.015)  
        f.tick_labels.hide()
        f.ticks.hide()
        g = fig.add_axes([0.03,0.7,0.9,0.25])
        self.plt_lc(g,n)
        g.axis('off')
        plt.savefig(os.path.join(self.desdir,'rgb_plus_lc_%s.png'%n))
        plt.clf()
        return os.path.join(self.desdir,'rgb_plus_lc_%s.png'%n)


if __name__ =="__main__":
    miika = np.loadtxt('/home/wiseman/code/des_stacks/source_lists/miika_5yr.txt',dtype='str')
    ccs = np.loadtxt('/home/wiseman/code/des_stacks/source_lists/spec_CC.txt',dtype='str')
    hashes = '#'*100
    #miika_hosts = dlr1s.merge(pd.DataFrame(miika,columns=['TRANSIENT_NAME']),how='inner',on='TRANSIENT_NAME')
    for ret in miika:
        print(hashes)
        print(ret)
        SN = sn(ret)
    
        SN.make_image_host()
