import math
import numpy as np
import pandas as pd
import MDAnalysis as mda
from numpy.linalg import norm
from MDAnalysis.analysis import rms, align
from scipy.spatial.distance import cdist
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA
from numpy import savetxt
import os



def universe(topology,trajectory):
    u = mda.Universe(topology,trajectory)
    return u

def water_check(universe,water_group,ligand_group,cutoff):
    water_lt = []
    cutoff_distance = cutoff  # Angstroms
    for ts in universe.trajectory:
        water_around = universe.select_atoms(f"group water and around {cutoff_distance} group ligand", water=water_group, ligand=ligand_group)
        water_around_num = len(np.unique(water_around.resids))
        water_lt.append(water_around_num)
    return np.array(water_lt) 
    
def hbonds_analysis(arg1,arg2,universe):
    hbonds = HBA(
      universe=universe,
      between=[arg1, arg2]
      )

    arg1_hydrogens_sel = hbonds.guess_hydrogens(arg1)
    arg1_acceptors_sel = hbonds.guess_acceptors(arg1)

    arg2_hydrogens_sel = hbonds.guess_hydrogens(arg2)
    arg2_acceptors_sel = hbonds.guess_acceptors(arg2)

    hbonds.hydrogens_sel = f"({arg1_hydrogens_sel}) or ({arg2_hydrogens_sel})"
    hbonds.acceptors_sel = f"({arg1_acceptors_sel}) or ({arg2_acceptors_sel})"

    hbonds.run()
    return hbonds    
    
def hbond_process(universe,arg1,arg2):
    hbonds = hbonds_analysis(arg1,arg2,universe) 
    with open(f'hbonds{arg1}{arg2}.csv', 'w') as hbfile:
        for hbond in hbonds.hbonds:
            time_frame = hbond[0]
            donor = int(hbond[1])
            hydrogen = int(hbond[2])
            acceptor = int(hbond[3])
            donor_id = u.atoms[[donor]].resids[0]
            donor_resname = u.atoms[[donor]].resnames[0]
            donor_name = u.atoms[[donor]].names[0]
            hydrogen_name = u.atoms[[hydrogen]].names[0]
            acceptor_id = u.atoms[[acceptor]].resids[0]
            acceptor_name = u.atoms[[acceptor]].names[0]
            acceptor_resname = u.atoms[[acceptor]].resnames[0]    
            hbfile.write(str(time_frame)+' '+str(donor_id)+str(donor_resname)+'@'+str(donor_name)+'_'+str(hydrogen_name)+'_'+str(acceptor_id)+str(acceptor_resname)+'@'+str(acceptor_name)+'\n')
    hbfile.close()
    
def distance_cal(universe,selection1,selection2,dis_cutoff):
    count = 0
    for ts in universe.trajectory:
        dis = norm(selection1.positions - selection2.positions)
        if dis <= dis_cutoff:
            count += 1
    return count    
    
def protein_distance_contact_map(arg1,arg2,universe,dis_cutoff):
    # Create a new 2D dataframe
    # arg1 and arg2 should be a tuple (first_resid, last_resid)
    data = pd.DataFrame(0, index=range(arg1[0],arg1[1]+1), columns=range(arg2[0],arg2[1]+1))
    data.index = [i for i in range(arg1[0],arg1[1]+1)]
    data.columns = [i for i in range(arg2[0],arg2[1]+1)]
    u = universe
    for res1 in range(arg1[0],arg1[1]+1): 
        count = 0
        for res2 in range(arg2[0],arg2[1]+1):
            res1_c = u.select_atoms(f"resid {res1} and name CA")           
            res2_c = u.select_atoms(f"resid {res2} and name CA")
            dis_count = distance_cal(u,res1_c,res2_c,dis_cutoff)
            data.at[int(res1),int(res2)] = data.at[int(res1),int(res2)] + dis_count
    data.to_csv('distance_contact.csv',index=True)        

def compute_aromatic_normal(ring_atoms):
    # Calculate the centroid of the ring
    centroid = np.mean(ring_atoms.positions, axis=0)

    # Perform PCA to find the normal vector of the ring plane
    coords_centered = ring_atoms.positions - centroid
    _, _, vh = np.linalg.svd(coords_centered)
    normal = vh[2, :]

    return normal

def calculate_angle(v1, v2):
    """Calculate the angle between two vectors in degrees."""
    cos_angle = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
    angle = np.arccos(np.clip(cos_angle, -1, 1))
    return np.degrees(angle)

def cation_pi_check(cation_pos,centroid,normal,distance,dis_cutoff=6.0,angle_min=45,angle_max=135):
    """ The default cutoff of distance is 6 angstrom, the angle range is > 135 or < 45 degree """
    if distance < dis_cutoff:
        # Compute angle
        vector = cation_pos - centroid
        angle = calculate_angle(normal, vector[0])

        if angle > angle_max or angle < angle_min:
            return True
        else:
            return False
            
def cation_pi_analysis(universe,cation,aromatics):
    cation = u.select_atoms(f'{cation} and name N1')  # Modify as needed
    aromatics = u.select_atoms(f'{aromatics} and name CG CD1 CD2 CE1 CE2 CZ')
    frame = 0
    with open(f'cation_pi.csv', 'w') as cpfile:
        for ts in u.trajectory:
            cation_pos = cation.positions
            for aromatic in aromatics.residues:
                ring_atoms = aromatic.atoms
                normal = compute_aromatic_normal(ring_atoms)
                # Calculate centroid of the ring for distance measurement
                centroid = np.mean(ring_atoms.positions, axis=0)
                distance = np.linalg.norm(cation_pos - centroid)
                if cation_pi_check(cation_pos,centroid,normal,distance):
                    cpfile.write(str(frame)+' '+str(aromatic.resid)+'-'+str(aromatic.resname)+'\n')
            frame += 1        
    cpfile.close()        
    
def theta_calc(A,O,B): # O is the center point of the angle
    OA = A - O
    OB = B - O
    theta = np.arccos(np.dot(OA,OB)/(norm(OA)*norm(OB)))
    return np.rad2deg(theta)
    
def theta_calc_trajectory(point_A,point_O,point_B,universe,theta_initial,abs_angle=True):
    theta_lt = []
    for ts in universe.trajectory:
        A = universe.select_atoms(point_A).center_of_geometry()
        O = universe.select_atoms(point_O).center_of_geometry()
        B = universe.select_atoms(point_B).center_of_geometry()
        if abs_angle == True:
            theta = theta_calc(A,O,B) - theta_initial
        elif abs_angle == False:    
            theta = theta_calc(A,O,B)
        theta_lt.append(theta)
    return np.array(theta_lt)   
    
def loop_angle(topology, crdfile, trajectory, point_A, point_B, point_O, traj_analysis=True, abs_angle=True):
    topfile = topology
    inpfile = crdfile
    u = mda.Universe(topfile, inpfile)

    pointO = u.select_atoms(point_O).center_of_geometry()
    pointA = u.select_atoms(point_A).center_of_geometry()
    pointB = u.select_atoms(point_B).center_of_geometry()

    theta_AOB = theta_calc(pointA, pointO, pointB)   
    if traj_analysis == True:
        u = universe(topology, trajectory)
        theta = theta_calc_trajectory(point_A, point_O, point_B, u, theta_AOB, abs_angle)
        savetxt('/loop_angle.csv', theta, delimiter=' ')
    elif traj_analysis == False:
        print(theta_AOB)


def aligner(universe,reference):
    """
    reference is the selections of relatively stable residues.
    
    """
    average = align.AverageStructure(universe,universe,select='protein and name CA',ref_frame=0).run()
    ref = average.results.universe
    aligned_traj = align.AlignTraj(universe
    ,ref
    ,select=f'resid {reference} and name CA'
    ,filename='aligned_traj.dcd'
    ,in_memory=False).run()

def rmsf_analysis(topfile,selection):
    """selection is the resids like 1 to 218  """
    print("Starting align the trajectory...")
    aligner() 
    u = mda.Universe(topfile,'aligned_traj.dcd') 
    c_alphas = u.select_atoms(f'resid {selection} and name CA')
    R = rms.RMSF(c_alphas).run()
    savetxt('rmsf.csv',R.results.rmsf,delimiter=' ')
      
    














































   
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
              
