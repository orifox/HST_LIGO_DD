#!/usr/bin/env python
# 2019
# L. Strolger, STScI

#
"""
prob_trigger.py [--options]
returns the likelihood that we should execute a trigger
options are:
   --verbose (shows all output, by default)
   --quiet (shows only the likelihood)
   --help
   --plotit (generates plots)
   --area = [float] (in deg2)
   --distance = [float] (in Mpc)
   --bprob = [float] (BNS probablity)
   --velocity = [float] (spectral expansion velocity, if known)
   --hours = [float] (time since burst, in hours)
Options can be combined. There are some default values for testing
"""

import os,sys,pdb,scipy,glob
from pylab import *

d0 = 150 #Mpc
q0 = 1000 #deg2
b0 = 1 #BNS prob
verbose=True
hours=12 ##from burst
velocity = 0
plotit=False

def sigmoid_sn(x,*p):
    x0, T, k, b = p
    b = sqrt(b**2)
    T = sqrt(T**2)
    y = (T / (1 + exp(-k*(x-x0))))-b
    return y



def prob_confusion(d0, q0, verbose=True, plotit=True):
    p0 = (1460, 1, -0.003, 0.)
    p1 = (1460, 1, 0.003, 0.)
    Ngal = q0 *pi**2/180**2*d0**2*0.14
    Prob = sigmoid_sn(Ngal, *p0)
    Nsn=Ngal/100*5/265*5
    Nlo=Nsn-Nsn/5*2
    Nhi=Nsn/5*10.-Nsn
    if verbose:print('Ngal=%.1f, Nsn=%.1f+%.1f-%.1f, Nprob=%.2f' %(Ngal, Nsn, Nhi, Nlo, Prob))

    if plotit:
        clf()
        xx = arange(0,5000,10) ## number of potential hosts
        ax1 = subplot(221)
        ax2 = subplot(222)
        ax3 = subplot(223)
        ax4 = subplot(224)
        ## limit of confusion with other SNe
        ax1.plot(xx, sigmoid_sn(xx, *p0), 'b-', label=r'P$_{\rm good}$')
        ax5 = ax1.twinx()
        ax5.plot(xx, sigmoid_sn(xx, *p1), 'b--', label=r'P$_{\rm bad}$')
        ## by volume
        volume = xx / 0.14 ## local grp has 2/14 Mpc^-3
        ax2.plot(volume/1e6, sigmoid_sn(xx, *p0), 'b-', label=r'P$_{\rm good}$')
        ax6 = ax2.twinx()
        ax6.plot(volume/1e6, sigmoid_sn(xx, *p1), 'b--', label=r'P$_{\rm bad}$')
        ## by distance
        omegas = linspace(2000, 20000, 5) #deg^2
        my_map = plt.get_cmap('RdBu_r')
        colors = r_[linspace(0.1,1,len(omegas)),linspace(0.1,1,len(omegas))]    
        my_colors = my_map(colors)
        for i,omega in enumerate(omegas):
            distance = sqrt(volume/omega*(180/pi)**2)
            if i==0:
                ax3.plot(distance, sigmoid_sn(xx, *p0), '-', color=my_colors[i], label=r'%2d $\times10^3$deg$^2$'%(omega/1e3))
            else:
                ax3.plot(distance, sigmoid_sn(xx, *p0), '-', color=my_colors[i], label=r'%2d'%(omega/1e3))
        ## by angular size
        distances = [30, 40, 80, 100, 150]#, 200]
        colors = r_[linspace(0.1,1,len(distances)),linspace(0.1,1,len(distances))]    
        my_colors = my_map(colors)
        for i,distance in enumerate(distances):
            omega = (180/pi)**2*volume/(distance**2)
            ax4.plot(omega/1e3, sigmoid_sn(xx, *p0),'-', color=my_colors[i], label='%2d Mpc'%distance)
        ax3.set_xlim(0,200)
        ax4.set_xlim(0,20)
        ax1.set_xlabel('Number of Normal Galaxies')
        ax2.set_xlabel(r'Volume (Gpc$^3$)')
        ax3.set_xlabel('Distance (Mpc)')
        ax4.set_xlabel(r'Sky localization (10$^3$ deg$^2$)')
        tight_layout()
        ax1.legend(frameon=False, loc=4)
        ax2.legend(frameon=False, loc=4)
        ax5.legend(frameon=False, loc=1)
        ax6.legend(frameon=False, loc=1)
        ax3.legend(frameon=False)
        ax4.legend(frameon=False)
        ax5.set_yticks([])
        ax6.set_yticks([])
        savefig('prob_ccsne.png')
    return(Prob)
    
    
def prob_spec(velocity, verbose=True, plotit=True):

    p0=(10000, 1., 1./2000, 0.)
    Prob = sigmoid_sn(velocity, *p0)
    if verbose: print('Prob_velocity = %.2f' %Prob)
    if plotit:
        clf()
        velocities = arange(0, 30000, 500)
        ax=subplot(111)
        ax.plot(velocities, sigmoid_sn(velocities, *p0), 'r-')
        ax.set_xlabel('Expansion Velocity')
        savefig('prob_spec.png')
    return(Prob)

def prob_time(hours_since_burst, verbose=True, plotit=True):

    Prob=exp(-hours/24/7.)
    if verbose: print('Prob_hours=%.2f' %Prob)

    if plotit:
        xx=arange(0,14,1/24.)
        ax=subplot(111)
        ax.plot(xx,exp(-xx/7.), 'r-')
        ax.set_xlabel('Trigger Submitted (Days)')
        savefig('prob_time.png')

    return(Prob)
        

if __name__=='__main__':
    import getopt
    try:
        opt, arg = getopt.getopt(
            sys.argv[1:],"v,h",
            longopts=["verbose", "quiet","help","plotit",
                      "distance=", "area=", "bprob=",
                      "velocity=","hours="])
    except getopt.GetoptError:
        print ("Error : incorrect option or missing argument.")
        print (__doc__)
        sys.exit(1)
    for o, a in opt:
        if o in ["-h", "--help"]:
            print (__doc__)
            sys.exit(0)
            
        elif o == "-v":
            verbose=True
        elif o == "--verbose":
            verbose=True
        elif o == "--quiet":
            verbose=False

        elif o =="--plotit":
            plotit=True

        elif o =="--distance":
            d0 = float(a)
        elif o =="--area":
            q0 = float(a)
        elif o =="--bprob":
            b0 = float(a)
        elif o =="--velocity":
            velocity = float(a)
        elif o =="--hours":
            hours = float(a)
            
    if verbose:
        print('Assuming distance = %.1f Mpc, area = %.1f deg2, BNS_prob=%.1f, OT_prob = 1.0' %(d0,q0,b0),)
        print('hours since burst= %.1f hours, spec_velocity = %.1f km/s (if zero, no spectral confirmation)' %(hours, velocity))
            

    
    P_ccsne=prob_confusion(d0, q0, verbose=verbose, plotit=plotit)
    if verbose: print('BNS_prob=%.1f' %b0)
    P_time=prob_time(hours, verbose=verbose, plotit=plotit)
    if velocity != 0.:
        P_spec=prob_spec(velocity, verbose=verbose,plotit=plotit)
    else:
        P_spec=1.0

    print (" Likelihood that we should trigger = %.1f" %(P_ccsne*b0*P_time*P_spec))
    
    
    
