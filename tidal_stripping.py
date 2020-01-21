##########################################################################################
#
# MOVING RADIUS FINDER
# Info: 'moving_radius_finder' will track the mean radius,velocity, 
# and center of mass of a moving satellite. 
# 'Data' will be your hdf5 file
# 'Mass_Percentage' is the fraction of mass you want to count as
# the halo.
# 'Shell_Amounts' simply defines the number of sherical bins to be used
# for calculating the minimum radius containing XXXX fraction of the halo mass
#
# HALF MASS RADIUS INITIAL
# 'half_mass_radiusInitial' simply inputs the hdf5 file as 'Data', and eliminates
# the particles outside the initial half mass radius (input as HMR). To calculate
# the initial half mass radius simply enter a mass fraction of 0.5 into the initial
# hdf5 snapshot of 'moving_radius_finder' and return the mean radius.
#
# HALF MASS RADIUS ITER
# 'half_mass_radius_iter' acts in the same way as 'half_mass_radiusInitial' except
# takes the extra input parameter IDs, which is the particle IDs. The part IDs will 
# be in the output from the previus function, 
# ex: HMR_Initial = half_mass_radiusInitial(Data,HMR)
#     HMR_Next_Snapshot = half_mass_radiusIter(Data,HMR,HMR_Initial[:,3])
#     Since the part IDs are stored in the 4th column of the output from
#     half_mass_radiusInitial
#
##########################################################################################

import h5py
import numpy as np
import math


# Radius finder function
def moving_radius_finder(Data,Mass_Percentage,Shell_Amounts):
    """ A function for finding the approximated
    radius of a moving satellite object withing a dark matter halo.
    Works with hdf5 outputs from gadget snapshot files"""

    # Load in Coordinates
    Coords = np.array(Data.get('Coordinates'))
    PartID = np.array(Data.get('ParticleIDs'))
    Vel = np.array(Data.get('Velocities'))
    Num_Parts = len(Coords[:,0])
    x = Coords[:,0]
    y = Coords[:,1]
    z = Coords[:,2]

    # Find Most Dense Region
    x_ave = np.sum(x)/Num_Parts
    y_ave = np.sum(y)/Num_Parts
    z_ave = np.sum(z)/Num_Parts

    # Centralize Coordinates
    x = x - x_ave
    y = y- y_ave
    z = z- z_ave
    Data = np.column_stack((Coords,PartID))

    # Get Radius
    Radius = np.sqrt(x**(2) + y**(2) + z**(2))
    Data = np.column_stack((Data,Radius,Vel))
    max_radius = int(np.amax(abs(Radius)))
    Step = np.linspace(0,max_radius,Shell_Amounts)

    # Finding Desired Mass Fraction
    for i in range(len(Step)):
        Fraction = Radius <= Step[i]
        Current_Mass_Percentage = (np.sum(Fraction))/Num_Parts
        if Current_Mass_Percentage >= Mass_Percentage:
            Shell_Radius_Tol = Step[i]
            # Eliminate particles outside shell Radius
            Data = Data[Data[:,4] <= Shell_Radius_Tol]
            Parts_In_Shell = len(Data[:,4])
            Mean_Radius = np.sum(Data[:,4])/Parts_In_Shell
            Stdev = np.std(Data[:,4])
            x_n = np.sum(Data[:,0])/len(Data[:,1])
            y_n = np.sum(Data[:,1])/len(Data[:,1])
            z_n = np.sum(Data[:,2])/len(Data[:,1])
            CoM = [x_n,y_n,z_n]
            Vel = [np.sum(Data[:,5])/len(Data[:,1]),
            np.sum(Data[:,6])/len(Data[:,1]),
            np.sum(Data[:,7])/len(Data[:,1])]
            Stdev_Vel = [np.std(Data[:,5]),np.std(Data[:,6]),np.std(Data[:,7])]
            Max_Rad = Shell_Radius_Tol
            break

    return [Mean_Radius,Stdev,CoM,Vel,Stdev_Vel,Max_Rad];


# Radius finder function
def half_mass_radiusInitial(Data,HMR):
    """ A function for finding initial half-mass radius
    of a spherical dark matter halo.
    Works with hdf5 outputs from gadget snapshot files"""

    # Load in Coordinates
    Coords = np.array(Data.get('Coordinates'))
    PartID = np.array(Data.get('ParticleIDs'))
    Num_Parts = len(Coords[:,0])
    x = Coords[:,0]
    y = Coords[:,1]
    z = Coords[:,2]

    # Find Most Dense Region
    x_ave = np.sum(x)/Num_Parts
    y_ave = np.sum(y)/Num_Parts
    z_ave = np.sum(z)/Num_Parts

    # Centralize Coordinates
    x = x - x_ave
    y = y- y_ave
    z = z- z_ave
    Data = np.column_stack((Coords,PartID))

    # Get Radius
    Radius = np.sqrt(x**(2) + y**(2) + z**(2))
    Data = np.column_stack((Data,Radius))

    # Eliminate particles outside Half Mass Radius
    Data = Data[Data[:,4] <= HMR]
    #N_Parts = len(Data[:,4])

    return (Data);

# Radius finder function
def half_mass_radiusIter(Data,HMR,IDs):
    """ A function for iteritively finding amount of particles
    in an initial half-mass radius.
    Works with hdf5 outputs from gadget snapshot files"""

    # Load in Coordinates
    Coords = np.array(Data.get('Coordinates'))
    PartID = np.array(Data.get('ParticleIDs'))
    Data = [0,0,0,0] # To initialise for now, delete/fix later

    # Eliminate Particles not in ID list
    for i in range(len(PartID)):
        if PartID[i] in IDs:
            New_Row = [Coords[i,0],Coords[i,1],Coords[i,2],PartID[i]]
            Data = np.vstack([Data,New_Row])

    # Eliminate first row
    Data = np.delete(Data, (0), axis=0)
    Num_Parts = len(Data[:,0])
    x = Data[:,0]
    y = Data[:,1]
    z = Data[:,2]

    # Find Most Dense Region
    x_ave = np.sum(x)/Num_Parts
    y_ave = np.sum(y)/Num_Parts
    z_ave = np.sum(z)/Num_Parts

    # Centralize Coordinates
    x = x - x_ave
    y = y- y_ave
    z = z- z_ave

    # Get Radius
    Radius = np.sqrt(x**(2) + y**(2) + z**(2))
    Data = np.column_stack((Data,Radius))

    # Eliminate particles outside Half Mass Radius
    Data = Data[Data[:,4] <= HMR]


    return (Data);
