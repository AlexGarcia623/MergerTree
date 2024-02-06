import sys
import os
import time
import h5py
import numpy as np
import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt
import illustris_python as il
import scipy.integrate as integrate

# set RC params
plt.rcParams['figure.facecolor']    = 'white'
mpl.rcParams['axes.facecolor']      = 'white'
mpl.rcParams['font.size']           = 20
mpl.rcParams['axes.linewidth']      = 1.5
mpl.rcParams['xtick.direction']     = 'in'
mpl.rcParams['ytick.direction']     = 'in'
mpl.rcParams['xtick.minor.visible'] = 'true'
mpl.rcParams['ytick.minor.visible'] = 'true'
mpl.rcParams['xtick.major.width']   = 1.5
mpl.rcParams['ytick.major.width']   = 1.5
mpl.rcParams['xtick.minor.width']   = 1.0
mpl.rcParams['ytick.minor.width']   = 1.0
mpl.rcParams['xtick.major.size']    = 7.5
mpl.rcParams['ytick.major.size']    = 7.5
mpl.rcParams['xtick.minor.size']    = 3.5
mpl.rcParams['ytick.minor.size']    = 3.5
mpl.rcParams['xtick.top']           = True
mpl.rcParams['ytick.right']         = True

os.environ['MANPATH']="/home/paul.torrey/local/texlive/2018/texmf-dist/doc/man:$MANPATH"
os.environ['INFOPATH']="/home/paul.torrey/local/texlive/2018/texmf-dist/doc/info:$INFOPATH"
os.environ['PATH']="/home/paul.torrey/local/texlive/2018/bin/x86_64-linux:/home/paul.torrey/local/texlive/2018/texmf-dist:$PATH"

mpl.rcParams['text.usetex']        = True
mpl.rcParams['font.family']        = 'serif'
mpl.rc('font',**{'family':'sans-serif','serif':['Times New Roman'],'size':15})
mpl.rc('text', usetex=True)

path = '/orange/paul.torrey/IllustrisTNG/Runs/'

run  = 'L35n2160TNG'

base_path = path + '/' + run
post_dir  = path + '/' + run + '/postprocessing'
tree_dir  = post_dir + '/trees/SubLink/'

# out_dir   = path + '/' + run + '/output'

snap = 99

SUBHALO_ID   = 526478 # TNG0053
file_num     = 0 
file_index   = None
file_located = False

# small = 1/500
# big   = 1/100

MAJOR_MERGER_FRAC = 1.0/200.0
MERGER_MASS_TYPE  = 4 # 0 -> Gas, 1 -> DM, 4 -> Stars


h  = 6.774E-01
xh = 7.600E-01
zo = 3.500E-01
mh = 1.6725219E-24
kb = 1.3086485E-16
mc = 1.270E-02

def ageFromScaleFactor(a,aEnd=1):
    h  = 6.774E-01 
    H0 = 100 * h# km/s/Mpc
    H0 /= 3.086e+19 # km/s/km = 1/s
    OmegaM0 = 0.30890000000000001
    OmegaLambda0 = 0.69110000000000005
    
    t = []
    
    if (len(a) > 0):
        expr = lambda a: 1 / np.sqrt(OmegaM0/a + OmegaLambda0*a**2)
        # Calculate t in s
        t = np.array([integrate.quad(expr, i, aEnd)[0]/H0 for i in a])
        # Convert t to Gyr
        t /= 3.15576e+16

    return t

# Convert to raw SubhaloID
SUBHALO_ID += int(snap * 1.00E+12)

# Locate the SubLink tree file that the Subhalo exists in
while not file_located:    
    tree_file = h5py.File( tree_dir + 'tree_extended.%s.hdf5' %file_num, 'r' )
    
    subID = np.array(tree_file.get('SubhaloIDRaw'))
    
    overlap = SUBHALO_ID in subID
    
    if (overlap):
        file_located = True
        file_index   = np.where( subID == SUBHALO_ID )[0][0]
    else:
        file_num += 1

print('Found Subhalo ID %s in file %s' %(SUBHALO_ID, file_num))
# With the file number in hand, let's learn about our galaxy!
tree_file = h5py.File( tree_dir + 'tree_extended.%s.hdf5' %file_num, 'r' )

redshifts = None

with h5py.File('/orange/paul.torrey/IllustrisTNG/Runs/L35n2160TNG/postprocessing/trees/LHaloTree/trees_sf1_099.0.hdf5','r') as f:
    for i in f:
        tree1 = f.get(i)
        redshifts = np.array(tree1.get('Redshifts'))
        break

snap_to_z = { snap: redshifts[i] for i, snap in enumerate(range(0,100)) }
z_to_snap = { round(redshifts[i],2): snap for i, snap in enumerate(range(0,100)) }

snap_to_scf = { snap: 1/(1+redshifts[i]) for i, snap in enumerate(range(0,100)) }

scf = 1/(redshifts + 1)
age = ageFromScaleFactor(scf)
        
rootID  = tree_file["SubhaloID"][file_index]
fp      = tree_file["FirstProgenitorID"][file_index]
snapNum = tree_file["SnapNum"][file_index]
mass    = np.log10( tree_file["SubhaloMassType"][file_index,MERGER_MASS_TYPE] * 1.00E+10 / h )
size    = tree_file["SubhaloHalfmassRad"][file_index] * (snap_to_scf[snapNum]) / h
subfind = tree_file["SubfindID"][file_index]
subhalo = tree_file["SubhaloID"][file_index]
Zgas    = tree_file['SubhaloGasMetallicity'][file_index]
Zgas_sfr= tree_file['SubhaloGasMetallicitySfr'][file_index]

OH       = Zgas * (zo/xh) * (1.00/16.00)
Zgas     = np.log10(OH) + 12
OH_SFR   = Zgas_sfr * (zo/xh) * (1.00/16.00)
Zgas_sfr = np.log10(OH_SFR) + 12

print( 'Subfind ID: %s' %subfind )
print( 'Subhalo ID: %s' %subhalo )

sizeList    = [size]
massList    = [mass]
snapList    = [snap]
ZgasList    = [Zgas]
ZgasSFRList = [Zgas_sfr]

all_mergers_mass    = []
all_mergers_snap    = []
all_mergers_mf      = []
all_mergers_Zgas    = []
all_mergers_ZgasSFR = []

wantedSnaps = [99, 67, 50, 33, 25, 21, 17]

print('z=0 Mass: %s' %mass)
print('z=0 Z: %s' %round(Zgas,2))

while fp != -1:
    fpIndex = file_index + (fp - rootID)
    
    Sfid    = tree_file["SubfindID"][fpIndex]
    mass    = np.log10( tree_file["SubhaloMassType"][fpIndex,MERGER_MASS_TYPE] * 1.00E+10 / h )
    fpSnap  = tree_file["SnapNum"][fpIndex]
    size    = tree_file["SubhaloHalfmassRad"][fpIndex] * (snap_to_scf[fpSnap]) / h
    Zgas    = tree_file['SubhaloGasMetallicity'][fpIndex]
    Zgas_sfr= tree_file['SubhaloGasMetallicitySfr'][fpIndex]
    
    OH       = Zgas * (zo/xh) * (1.00/16.00)
    Zgas     = np.log10(OH) + 12
    OH_SFR   = Zgas_sfr * (zo/xh) * (1.00/16.00)
    Zgas_sfr = np.log10(OH_SFR) + 12
    
    sizeList.append( size )
    massList.append( mass )
    snapList.append( fpSnap )
    ZgasList.append( Zgas )
    ZgasSFRList.append( Zgas_sfr )
    
#     if (fpSnap in wantedSnaps):
#         print( 'z=%s Mass: %s' %(round(snap_to_z[fpSnap],1),mass) )

    nextProgenitor = tree_file["NextProgenitorID"][fpIndex]
    while nextProgenitor != -1:
        npIndex = file_index + (nextProgenitor - rootID)

        npMass  = np.log10( tree_file["SubhaloMassType"][npIndex,MERGER_MASS_TYPE] * 1.00E+10 / h )
        npSnap  = tree_file["SnapNum"][npIndex]
        
        npZgas    = tree_file['SubhaloGasMetallicity'][npIndex]
        npZgasSFR = tree_file['SubhaloGasMetallicitySfr'][npIndex]
        
        npOH       = npZgas * (zo/xh) * (1.00/16.00)
        npZgas     = np.log10(npOH) + 12
        npOH_SFR   = npZgasSFR * (zo/xh) * (1.00/16.00)
        npZgasSFR  = np.log10(npOH_SFR) + 12
        
        if ((10**npMass / 10**mass) > MAJOR_MERGER_FRAC) and (npMass > 0):
            # plt.scatter( npSnap, npMass ,color='red' )
        
            all_mergers_mass.append( npMass )
            all_mergers_snap.append( npSnap )
            all_mergers_mf  .append( 10**npMass / (10**mass) )
                        
            all_mergers_Zgas.append( npZgas )
            all_mergers_ZgasSFR.append( npZgasSFR )
            
            if (npZgasSFR < 0):
                print('Zgas %s, ZgasSFR %s' %(npZgas, npZgasSFR))
            
            print('Progenitor at snap %s (z=%s) has mass %s' %(npSnap, round(snap_to_z[npSnap],4), npMass) )
            # print('\tfraction: %s' %(10**npMass / 10**mass))
            print('Progenitor at snap %s (z=%s) has Z %s' %(npSnap,round(snap_to_z[npSnap],4),round(npZgasSFR,2)))
        
        nextProgenitor = tree_file["NextProgenitorID"][npIndex]
    
    fp = tree_file["FirstProgenitorID"][fpIndex]


all_mergers_mass    = np.array(all_mergers_mass    )
all_mergers_snap    = np.array(all_mergers_snap    )
all_mergers_mf      = np.array(all_mergers_mf      )
all_mergers_Zgas    = np.array(all_mergers_Zgas    )
all_mergers_ZgasSFR = np.array(all_mergers_ZgasSFR )
    

fig, axs = plt.subplots( 2, 1, figsize=(8,7), sharex=True )    

axs[0].scatter( snapList, massList, label=r'${\rm Main~ Progenitor~ - TNG0053}$' )
axs[0].plot( snapList, np.log10(10**np.array(massList)/200), color='k', linestyle='--', alpha=0.5 )
axs[0].plot( snapList, np.log10(10**np.array(massList)/100), color='k', linestyle=':', alpha=0.5 )
axs[0].scatter( all_mergers_snap, all_mergers_mass, color='red', label=r'${\rm Merger~ Partner}$' )

axs[0].axvspan( 72, 84, color='orange', alpha=0.25 )

# leg = axs[0].legend(frameon=True,labelspacing=0.05,loc='upper left',
#                  handletextpad=0, handlelength=0, markerscale=-1,framealpha=1,edgecolor=(1, 1, 1, 0))

# for i in range(len(leg.get_texts())): leg.legendHandles[i].set_visible(False)

# colors = ['C0','red']
# for index, text in enumerate(leg.get_texts()):
#     text.set_color(colors[index])

if (MERGER_MASS_TYPE == 0):
    axs[0].set_ylabel(r'$\log\left( M_{\rm gas}~[M_\odot] \right)$')
elif (MERGER_MASS_TYPE == 1):
    axs[0].set_ylabel(r'$\log\left( M_{\rm DM}~[M_\odot] \right)$')
elif (MERGER_MASS_TYPE == 4):
    axs[0].set_ylabel(r'$\log\left( M_*~[M_\odot] \right)$')

# plt.axvline( all_mergers_snap[0], color='k', linestyle='--')
# plt.text( 0.3,0.2 ,r'${\rm Last~ major~ merger:~} z=%s$' %round(snap_to_z[all_mergers_snap[0]],4), transform=plt.gca().transAxes )
# plt.text( 0.3,0.15,r'${\rm Merger~ fraction:~} %s$' %(round(all_mergers_mf[0],3)), transform=plt.gca().transAxes )

for z in [0.0,0.5,1.0,2.0,3.01,4.01,5.0]:
    axs[0].axvline( z_to_snap[z], color='gray', linestyle=':', alpha=0.5 )
    axs[1].axvline( z_to_snap[z], color='gray', linestyle=':', alpha=0.5 )
    axs[0].text( z_to_snap[z] + 0.5, 6.0, r'$z = %s$' %round(z,1), rotation=90, alpha=0.5 )
    axs[1].text( z_to_snap[z] + 0.5, 5.5, r'$z = %s$' %round(z,1), rotation=90, alpha=0.5 )

axs[1].scatter( snapList, ZgasSFRList, label=r'${\rm Main~ Progenitor~ - TNG0053}$' )
axs[1].scatter( all_mergers_snap, all_mergers_ZgasSFR, color='red', label=r'${\rm Merger~ Partner}$' )

leg = axs[1].legend(frameon=True,labelspacing=0.05,loc='lower right',
                    handletextpad=0, handlelength=0, markerscale=-1,framealpha=1,edgecolor=(1, 1, 1, 0))

for i in range(len(leg.get_texts())): leg.legendHandles[i].set_visible(False)

colors = ['C0','red']
for index, text in enumerate(leg.get_texts()):
    text.set_color(colors[index])

axs[1].set_ylabel(r'$\log({\rm O/H}) + 12 ~{\rm (dex)}$')
axs[1].set_xlabel(r'${\rm Snapshot}$')

# for z in [0.5,1.0,2.0,3.01,4.01,5.0]:
#     axs[1].axvline( z_to_snap[z], color='gray', linestyle=':', alpha=0.5 )
#     axs[1].text( z_to_snap[z] + 0.5, 5.5, r'$z = %s$' %round(z,1), rotation=90, alpha=0.5 )
        

plt.tight_layout()

plt.subplots_adjust( hspace=0.01 )

plt.savefig('MassHistory.pdf',bbox_inches='tight')

plt.clf()

plt.scatter( snapList, sizeList, label=r'${\rm Main~ Progenitor~ - TNG0053}$' )

leg = plt.legend(frameon=True,labelspacing=0.05,loc='upper left',
                    handletextpad=0, handlelength=0, markerscale=-1,framealpha=1,edgecolor=(1, 1, 1, 0))

for i in range(len(leg.get_texts())): leg.legendHandles[i].set_visible(False)

colors = ['C0','red']
for index, text in enumerate(leg.get_texts()):
    text.set_color(colors[index])

plt.ylabel(r'${\rm Half-mass~Radius~(kpc)}$')
plt.xlabel(r'${\rm Snapshot}$')

for z in [0.5,1.0,2.0,3.01,4.01,5.0]:
    plt.axvline( z_to_snap[z], color='gray', linestyle=':', alpha=0.5 )
    plt.text( z_to_snap[z] + 0.5, 7.0, r'$z = %s$' %round(z,1), rotation=90, alpha=0.5 )

for snap in all_mergers_snap:
    plt.axvline( snap, color='red', linestyle='--' )
    
plt.tight_layout()

plt.savefig('SizeHistory.pdf')

plt.clf()

fig, axs = plt.subplots( 2, 1, figsize=(8,6), sharex=True,
                        gridspec_kw={'height_ratios':[1,0.5]})

axs[0].scatter( snapList, ZgasSFRList, label=r'${\rm Main~ Progenitor~ - TNG0053}$' )
axs[0].scatter( all_mergers_snap, all_mergers_ZgasSFR, color='red', label=r'${\rm Merger~ Partner}$' )

leg = axs[0].legend(frameon=True,labelspacing=0.05,loc='best',
                    handletextpad=0, handlelength=0, markerscale=-1,framealpha=1,edgecolor=(1, 1, 1, 0))

for i in range(len(leg.get_texts())): leg.legendHandles[i].set_visible(False)

colors = ['C0','red']
for index, text in enumerate(leg.get_texts()):
    text.set_color(colors[index])

axs[0].set_ylabel(r'$\log({\rm O/H}) + 12 ~{\rm (dex)}$')
# axs[0].set_xlabel(r'${\rm Snapshot}$')

for z in [0.5,1.0,2.0,3.01,4.01,5.0]:
    axs[0].axvline( z_to_snap[z], color='gray', linestyle=':', alpha=0.5 )
    axs[0].text( z_to_snap[z] + 0.5, 5.5, r'$z = %s$' %round(z,1), rotation=90, alpha=0.5 )
    
# axs[1].scatter( snapList, massList, label=r'${\rm Main~ Progenitor~ - TNG0053}$' )

nans_mask = ( all_mergers_ZgasSFR > 0  )

axs[1].scatter( all_mergers_snap, all_mergers_mf, color='red', label=r'${\rm Merger~ Partner}$' )

# leg = plt.legend(frameon=True,labelspacing=0.05,loc='upper left',
#                     handletextpad=0, handlelength=0, markerscale=-1,framealpha=1,edgecolor=(1, 1, 1, 0))

# for i in range(len(leg.get_texts())): leg.legendHandles[i].set_visible(False)

colors = ['C0','red']
for index, text in enumerate(leg.get_texts()):
    text.set_color(colors[index])

axs[1].set_ylabel(r'${M_{\rm partner} / M_{\rm main}}$')
axs[1].set_xlabel(r'${\rm Snapshot}$')

# plt.axvline( all_mergers_snap[0], color='k', linestyle='--')
# plt.text( 0.3,0.2 ,r'${\rm Last~ major~ merger:~} z=%s$' %round(snap_to_z[all_mergers_snap[0]],4), transform=plt.gca().transAxes )
# plt.text( 0.3,0.15,r'${\rm Merger~ fraction:~} %s$' %(round(all_mergers_mf[0],3)), transform=plt.gca().transAxes )

for z in [0.5,1.0,2.0,3.01,4.01,5.0]:
    axs[1].axvline( z_to_snap[z], color='gray', linestyle=':', alpha=0.5 )
    
for snaps in all_mergers_snap[~nans_mask]:
    axs[0].axvline( snaps, color='red', linestyle='--', alpha=0.3 )
    
axs[1].axhline( MAJOR_MERGER_FRAC, color='k')
axs[1].text( 40, MAJOR_MERGER_FRAC*1.25, r'${\rm Minimum~ Threshold}$' )

xmin, xmax = axs[1].get_xlim()
ymin, ymax = axs[1].get_ylim()

axs[1].fill_between( np.arange(-100,200), MAJOR_MERGER_FRAC, -5, color="gray", alpha=0.5)
axs[1].fill_between( np.arange(-100,200), MAJOR_MERGER_FRAC, -5, facecolor="none", edgecolor='k', hatch='X', alpha=0.99)

axs[1].set_xlim(xmin, xmax)
axs[1].set_ylim(ymin, ymax)

plt.tight_layout()

plt.subplots_adjust( hspace=0.05 )

plt.savefig('ZHistory.pdf',bbox_inches='tight')

tree_file.close()