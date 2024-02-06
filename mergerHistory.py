import h5py
import numpy as np


#############################################################################
########### CHANGE THE FOLLOWING BASED ON THE DIRECTORY STRUCTURE ###########
path = '/orange/paul.torrey/IllustrisTNG/Runs/'
run  = 'L35n2160TNG'

base_path = path + '/' + run
post_dir  = path + '/' + run + '/postprocessing'
tree_dir  = post_dir + '/trees/SubLink/'
#############################################################################


## Variables
MAJOR_MERGER_FRAC = 1.0/200.0 # Sets major merger fraction, here i have 1:200
MERGER_MASS_TYPE  = 4 # 0 -> Gas, 1 -> DM, 4 -> Stars

h  = 6.774E-01
xh = 7.600E-01
zo = 3.500E-01
mh = 1.6725219E-24
kb = 1.3086485E-16
mc = 1.270E-02


#############################################################################
# Convert to raw SubhaloID (this is the form the tree files are in)
snap        = 99 # This galaxy exists at redshift of 0
SUBHALO_ID  = 526478 # I am looking for this specific Subhalo
SUBHALO_ID += int(snap * 1.00E+12)

# Locate the SubLink tree file that the Subhalo exists in
file_num     = 0 
file_index   = None
file_located = False
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
#############################################################################


#############################################################################
### The following is not necessary... these files might not be on Rivanna ###
# Get information from the headers on the redshift, scalefactors, etc
redshifts = None
# this tree file can be any LHaloTree file
different_tree_file = '/orange/paul.torrey/IllustrisTNG/Runs/L35n2160TNG/postprocessing/trees/LHaloTree/trees_sf1_099.0.hdf5'
with h5py.File(different_tree_file,'r') as f:
    for i in f:
        tree1 = f.get(i)
        redshifts = np.array(tree1.get('Redshifts'))
        break
# This will help me convert back and forth betwen redshift, snapshot and scalefactor
snap_to_z   = { snap: redshifts[i] for i, snap in enumerate(range(0,100)) }
z_to_snap   = { round(redshifts[i],2): snap for i, snap in enumerate(range(0,100)) }
snap_to_scf = { snap: 1/(1+redshifts[i]) for i, snap in enumerate(range(0,100)) }

scf = 1/(redshifts + 1)
age = ageFromScaleFactor(scf)
#############################################################################


#############################################################################
# With the file number in hand, let's learn about our galaxy!
with h5py.File( tree_dir + 'tree_extended.%s.hdf5' %file_num, 'r' ) as tree_file:
    
    rootID  = tree_file["SubhaloID"][file_index] ## ID in the tree
    snapNum = tree_file["SnapNum"][file_index] ## Snapshot number of this subhalo
    fp      = tree_file["FirstProgenitorID"][file_index] ## ID of the first progenitor (most massive galaxy at next snaphot)
    
    ## Properties of subhalo
    mass    = np.log10( tree_file["SubhaloMassType"][file_index,MERGER_MASS_TYPE] * 1.00E+10 / h )
    size    = tree_file["SubhaloHalfmassRad"][file_index] * (snap_to_scf[snapNum]) / h
    Zgas    = tree_file['SubhaloGasMetallicity'][file_index]
    Zgas_sfr= tree_file['SubhaloGasMetallicitySfr'][file_index]

    ## Metallicity conversions
    OH       = Zgas * (zo/xh) * (1.00/16.00)
    Zgas     = np.log10(OH) + 12
    OH_SFR   = Zgas_sfr * (zo/xh) * (1.00/16.00)
    Zgas_sfr = np.log10(OH_SFR) + 12

    # Containers
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

    ## Interested in z=0, 0.5, 1, 2, 3, 4, 5
    wantedSnaps = [99, 67, 50, 33, 25, 21, 17]

    print('z=0 Mass: %s' %mass)
    print('z=0 Z: %s' %round(Zgas,2))

    while fp != -1: # While there is still a first progenitor
        fpIndex = file_index + (fp - rootID) # Get the index of that progenitor

        # Get the properties of the first progenitor
        mass    = np.log10( tree_file["SubhaloMassType"][fpIndex,MERGER_MASS_TYPE] * 1.00E+10 / h )
        fpSnap  = tree_file["SnapNum"][fpIndex]
        size    = tree_file["SubhaloHalfmassRad"][fpIndex] * (snap_to_scf[fpSnap]) / h
        Zgas    = tree_file['SubhaloGasMetallicity'][fpIndex]
        Zgas_sfr= tree_file['SubhaloGasMetallicitySfr'][fpIndex]

        OH       = Zgas * (zo/xh) * (1.00/16.00)
        Zgas     = np.log10(OH) + 12
        OH_SFR   = Zgas_sfr * (zo/xh) * (1.00/16.00)
        Zgas_sfr = np.log10(OH_SFR) + 12

        # Save the properties of the first progenitor
        sizeList.append( size )
        massList.append( mass )
        snapList.append( fpSnap )
        ZgasList.append( Zgas )
        ZgasSFRList.append( Zgas_sfr )
        
        #### If you want to only learn about the first progenitor(s) stop here ####

        nextProgenitor = tree_file["NextProgenitorID"][fpIndex]
        while nextProgenitor != -1: ## If there is a next progenitor at this redshift
            npIndex = file_index + (nextProgenitor - rootID) ## get next progenitor ID

            ## Get properties of next progenitor
            npMass    = np.log10( tree_file["SubhaloMassType"][npIndex,MERGER_MASS_TYPE] * 1.00E+10 / h )
            npSnap    = tree_file["SnapNum"][npIndex]
            npZgas    = tree_file['SubhaloGasMetallicity'][npIndex]
            npZgasSFR = tree_file['SubhaloGasMetallicitySfr'][npIndex]

            npOH       = npZgas * (zo/xh) * (1.00/16.00)
            npZgas     = np.log10(npOH) + 12
            npOH_SFR   = npZgasSFR * (zo/xh) * (1.00/16.00)
            npZgasSFR  = np.log10(npOH_SFR) + 12

            ## If this is a major merger (and real)
            if (10**npMass / 10**mass) > MAJOR_MERGER_FRAC and npMass > 0:
                
                ## save the data
                all_mergers_mass.append( npMass )
                all_mergers_snap.append( npSnap )
                all_mergers_mf  .append( 10**npMass / (10**mass) )

                all_mergers_Zgas.append( npZgas )
                all_mergers_ZgasSFR.append( npZgasSFR )

                print('Progenitor at snap %s (z=%s) has mass %s' %(npSnap, round(snap_to_z[npSnap],4), npMass) )
                print('Progenitor at snap %s (z=%s) has Z %s' %(npSnap,round(snap_to_z[npSnap],4),round(npZgasSFR,2)))

            nextProgenitor = tree_file["NextProgenitorID"][npIndex] ## Breaks infinite recursion

        fp = tree_file["FirstProgenitorID"][fpIndex] ## Breaks infinite recursion
#############################################################################