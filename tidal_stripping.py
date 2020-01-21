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
    """ A function for finding amount of particles
    in an initial half-mass radius.
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
    """ A function for finding amount of particles
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
