import matplotlib
matplotlib.use('Agg')

import sys
import matplotlib.pyplot as plt
#from matplotlib.pyplot import figure
import numpy as np

"""
SMALL_SIZE = 8
MEDIUM_SIZE = 10
BIGGER_SIZE = 12

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
"""

#plt.rc('axes', titlesize=18)      # fontsize of the axes title
#plt.rc('axes', labelsize=14)      # fontsize of the x and y labels
#plt.rc('legend', fontsize=18)    # legend fontsize

plt.rc('axes', titlesize=24)      # fontsize of the axes title
plt.rc('axes', labelsize=20)      # fontsize of the x and y labels
plt.rc('legend', fontsize=20)     # legend fontsize
plt.rc('xtick', labelsize=16)    # fontsize of the tick labels
plt.rc('ytick', labelsize=16)    # fontsize of the tick labels

plt.rc('figure.subplot', left=0.2)
plt.rc('figure.subplot', right=0.95)
plt.rc('figure.subplot', bottom=0.13)
plt.rc('figure.subplot', top=0.88)




#plt.rc('font', size=24)          # controls default text sizes
#plt.rc('figure', titlesize=24)    # fontsize of the figure title


print(plt.style.available)
#plt.style.use('ggplot')
#plt.style.use('fivethirtyeight')
#plt.style.use('dark_background')

#sys.exit(0)

colors = {}

colors.update({"11": "#8b1a0e"})  # red       
colors.update({"12": "#FF4500"})  # orangered 
colors.update({"13": "#B22222"})  # firebrick 
colors.update({"14": "#DC143C"})  # crimson   

colors.update({"21": "#5e9c36"})  # green
colors.update({"22": "#006400"})  # darkgreen  
colors.update({"23": "#228B22"})  # forestgreen
colors.update({"24": "#808000"})  # olive      

colors.update({"31": "#8A2BE2"})  # blueviolet
colors.update({"32": "#00008B"})  # darkblue  

colors.update({"41": "#2F4F4F"})  # darkslategray




def plot(filename, ind):


    plt.title('Notitle')
    plt.xlabel('Time, fs')
    plt.ylabel('ACF/Integral, a.u.^2')

    clrs = ["11", "12", "21", "31", "41"]

    for n in [0.2, 0.4, 0.6, 0.8, 1.0]:
        t, dat = np.loadtxt("_acf-%2.1f.txt" % (n), usecols=(0, ind), unpack=True)
        indx = [0.2, 0.4, 0.6, 0.8, 1.0].index(n)
        cl = clrs[indx]
        plt.plot(t, dat, label='%2.1f' % (n), linewidth=4, color = colors[cl]) 
        plt.legend()
        plt.show()

    plt.savefig(filename, dpi=300)
    plt.close()


plot("_C_re.png", 5)
plot("_C_im.png", 6)
plot("_I_re.png", 7)


