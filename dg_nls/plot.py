import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import os
import pickle
os.environ["PROJ_LIB"] = "C:\\ProgramData\\Anaconda3\\Library\\share"
from mpl_toolkits.basemap import Basemap
plt.close('all')

matplotlib.rcParams['text.usetex']=True
matplotlib.rcParams['mathtext.fontset'] = 'cm'
plt.rc('font', **{'family': 'serif', 'serif': ['cmr10']})

titlefont = {'fontsize':12}
labelfont = {'fontsize':10}
tickfont = {'fontsize':8}
marker_size = 4
marker_color = 'c'

xdim = 101
ydim = 101
x = np.linspace(-0.5,0.5,xdim)
y = np.linspace(-0.5,0.5,ydim)
x,y = np.meshgrid(x,y)

lev = 301
plt.figure(figsize=[5+3/8,5+3/8])    
plt.subplot(221)
t=0
f = pickle.load(open( "ftle_0.pickle", "rb" ))
plt.pcolormesh(x,y,-f,shading='gouraud')
plt.annotate('A', xy=(0.02, 0.91), xycoords='axes fraction',color='white')
plt.subplot(222)
t=1
f = pickle.load(open( "ftle_1.pickle", "rb" ))
plt.pcolormesh(x,y,f,shading='gouraud')
plt.annotate('B', xy=(0.02, 0.91), xycoords='axes fraction',color='white')
plt.subplot(223)
t=2
f = pickle.load(open( "ftle_2.pickle", "rb" ))
plt.pcolormesh(x,y,f,shading='gouraud')
plt.annotate('C', xy=(0.02, 0.91), xycoords='axes fraction',color='white')
plt.subplot(224)
t=4
f = pickle.load(open( "ftle_4.pickle", "rb" ))
plt.pcolormesh(x,y,f,shading='gouraud')
plt.annotate('D', xy=(0.02, 0.91), xycoords='axes fraction',color='white')
plt.savefig('nls-Time_FTLE_Comparison.png', transparent=False, bbox_inches='tight',dpi=300)
plt.savefig('nls-Time_FTLE_Comparison.eps', transparent=False, bbox_inches='tight')
#plt.close('all')
    
    
    
    
    
    
    
    
    
    
    
    
    