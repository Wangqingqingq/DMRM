import numpy as np
import h5py
import matplotlib.pyplot as plt
from astropy.cosmology import Planck15

Volume = 100.0**3   # (Mpc)^3

def himassfunc(mhi,c,miu):
    Nbin = 10
    Nboot= 500
    mlow = np.min(mhi)
    mhig = np.max(mhi)
    step = (mhig-mlow)/float(Nbin)
    #print(mlow,mhig)    
    msbin= np.zeros(Nbin)
    ndens= np.zeros(Nbin)
    error= np.zeros(Nbin)
    for i in range(Nbin):
        ix       = (mhi>=i*step+mlow) & (mhi<(i+1)*step+mlow)
        print('Number:',len(mhi[ix]))
        msbin[i] = np.mean(mhi[ix])
        tmp      = np.zeros(Nboot)
        for j in range(Nboot):
          irand    = np.random.randint(low=0,high=(len(mhi)-1),size=len(mhi))
          resample = mhi[irand]
          iy       = (resample>=i*step+mlow) & (resample<(i+1)*step+mlow)
          tmp[j]   = float(len(resample[iy]))/(Volume*miu)
          #print(tmp[j])
        ndens[i] = float(len(mhi[ix]))/(Volume*miu)
        error[i] = np.std(tmp)
        print(msbin[i],ndens[i],error[i])
    ms    = np.zeros(len(msbin))
    for k in range(len(msbin)):
        ms[k] =10**msbin[k]
    lmbin = mass_to_l(ms,0.5)
    print(lmbin)
    print(ndens)
    plt.errorbar(lmbin,ndens,yerr=error,fmt=c,ecolor = c,markerfacecolor=c,elinewidth=2.5)
    #plt.yscale('log')
    #plt.xlim(mlow,mhig)
    #plt.ylim(1e-5,1.0)
    plt.xlabel(r'$\rm L_{HI}/L_{sun}$')
    plt.ylabel(r'$\rm n(L_{HI})$')
    plt.loglog()
    plt.savefig('HI_LmFunc.pdf')
    return {'Msbin':msbin,'Ndens':ndens,'Error':error}

def mass_to_l(m,z):
    lum  = np.zeros(len(m))
    dl   = Planck15.luminosity_distance(z)
    for i in range(len(m)):
        flux   = m[i]/(2.356*10**5*dl**2)
        lum[i] = flux*4*np.pi*dl**2
    return lum

def main():

    prefix = '../hih2_galaxy/'
    fname  = ['hih2_galaxy_033.hdf5','hih2_galaxy_040.hdf5','hih2_galaxy_050.hdf5','hih2_galaxy_067.hdf5','hih2_galaxy_099.hdf5']
    color  = ['tomato','darkorange','gold','lightgreen','lightskyblue']
    colorl = ['mistyrose','peachpuff','lightgoldenrodyellow','honeydew','lightblue']
#struc  = h5py.File(prefix+fname,'r')
    for i in range(2):
        struc  = h5py.File(prefix+fname[i],'r')
        data   = struc
        group_keys     = list(struc.keys())[15]
        #print("Keys: %s" % group_keys)
        m_hi_gd14_map  = (struc[group_keys])
        m_hi_gd14_mapl = np.array(m_hi_gd14_map)*1.2
        #for i in range(len(m_hi_gd14_map)):
        #    print(np.log10(m_hi_gd14_map[i]))
        ix             = np.array(m_hi_gd14_map)>1.0e+6  # This threshold need to be updated 
        ix2            = np.array(m_hi_gd14_mapl)>1.0e+6  # later according to the SKA2 detection limit
        mhi            = np.log10(np.array(m_hi_gd14_map[ix]))
        mhi2           = np.log10(np.array(m_hi_gd14_mapl[ix2]))
    himassfunc(mhi,color[1],1)
    himassfunc(mhi2,colorl[1],1.2)


if __name__=='__main__':
    main()
