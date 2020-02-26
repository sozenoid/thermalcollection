#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 14:55:29 2020

@author: macenrola
"""
def parse_xyz_block(xyz_block):
	"""
	PRE: Takes in an xyz block as a list of lines
	POST: Will return an array of array as np.array(np.array(), np.array(),...) each sub array contains [x,y,z] coord of a specific atom
	"""
	# let go of the atom names, the order is supposed to keep that information
	coords = np.array([np.array([float(y) for y in x.strip().split()[1:]]) for x in xyz_block[2:]], dtype=np.float64) # skip the two first lines that hold atom number and energy info
	return coords

def return_list_of_snapshots_and_atom_list(xyzfile):
	"""
	PRE: Takes in an xyz file
	POST: Returns a list of atoms and a list of snapshots of the trajectory as arrays
	"""
	natoms = -1
	snapshot_list=[]
	with open(xyzfile, 'rb') as r:
		tempxyz=[]
		for i, line in enumerate(r):
# =============================================================================
# 			if i==1000000: break
# =============================================================================
			if i==0:
				natoms=line.strip()
			if line.strip()==natoms and tempxyz!=[]:
				if i<1000:
					atm_list=[x.strip().split()[0] for x in tempxyz[2:]]
				xyz_array=parse_xyz_block(tempxyz)
				snapshot_list.append(xyz_array)
				tempxyz=[]
				nsnaps=len(snapshot_list)
				if nsnaps%10000==0: print "{}th snapshot processed".format(nsnaps)
			tempxyz.append(line)

	return atm_list, np.array(snapshot_list)[:]

def split_velocity_file(xyz_vel, natoms_large_frag=126):
	"""
	=========
	GO FOR IT; there was a bug in the file splitting part i'd say, probably better to keep trajectories short 
	=========
	PRE  : Takes in a velocity file formatted as an open babel xyz file, the atom breakdown is given in natoms, the smaller fragment is on top of the stack
	POST : returns two xyz_velocity files and adapts the headers for the two molecule components
	"""
	small_frag = xyz_vel+"-small_frag.xyz"
	large_frag = xyz_vel+"-large_frag.xyz"
	marker=None
	line_list=[]
	with open(small_frag, "wb") as ws:
		with open(large_frag, "wb") as wb:
			with open(xyz_vel, "rb") as r:
				for i, line in enumerate(r):
					if i==0:
						natoms = int(line)
					if i%10000==0: print i
					line_list.append(line)
					if i==0:
						marker=line
					if marker == line and i>0:
						ws.writelines("     {}\n".format(natoms-natoms_large_frag))
						ws.writelines(line_list[1])
						ws.writelines(line_list[2:natoms-natoms_large_frag+2])
						wb.writelines("     {}\n".format(natoms_large_frag))
						wb.writelines(line_list[1])
						wb.writelines(line_list[natoms-natoms_large_frag+2:-1])
						line_list=[line]
						
def convert_position_snapshots_to_velocity_snapshots(xyztraj, dt=1e-15):
	"""
	PRE  : Takes in an xyz trajectory where the snaps are separated by a small time step (for the velocity to be "instant")
	POST : will compute the velocity based on the time step size and also provide an estimate of the kinetic energy at each step
	"""
	print xyztraj
	atm_list, snapshot_pos= return_list_of_snapshots_and_atom_list(xyztraj) # in Angstroms
	print snapshot_pos.shape
	snapshots_vels = np.diff(snapshot_pos, axis=0) # in Angstrom/fs
	print snapshots_vels.shape
	dic_weights={'C':12,'H':1,'N':14, 'O':16, 'F':19.}
	conv=(1e-10/1e-15)**2*1.66e-27*6.022e23/4.184/1e3 # (Ang/fs)**2*amu*mol/kcal
	ekin = [float(0.5*conv*sum([y[0]**2*y[1] for y in zip(np.reshape(x, x.shape[0]*x.shape[1]), [dic_weights[i] for k in atm_list for i in [k]*3])])) for x in snapshots_vels]
	ekin = [x for x in ekin if np.abs(x)<1e3 ]
	with open(xyztraj+"-EKIN", "w") as w:
		for els in ekin:
			w.write(str(els)+"\n")

def plot_summary_EKIN(fcomplex, gcomplex, smallfrag, bigfrag):
	"""
	PRE : Takes in EKIN files that contain a list of kinetic energies in kcal/mol
	POST: Plot them along with their variance
	"""
	for f in zip([fcomplex, gcomplex, smallfrag, bigfrag], ["E(G2@CB[7])","E(G2) vacuum","E(G2) @CB[7]","G2@E(CB[7])"], [.2,.2,0.2,0.2]):
		with open(f[0], "rb") as r:
			width =1000
			ekin = [float(x.strip()) for x in r.readlines()]
			if f[1] == "E(G2) vacuum":
				ekin = ekin[10000:]
			ekin = [x  for x in ekin[:40000] if x<200]
			moving_av = np.convolve(ekin, np.ones(width)/width)
			plt.plot(moving_av[1002:-width], label="{0} ({1:3.1f}/{2:3.1f})".format(f[1], np.mean(ekin[2:]), np.var(ekin[2:])), linewidth=1, alpha=1)
			plt.plot(ekin[2:], linewidth=f[2], alpha=f[2])
	plt.title(gcomplex.split("/")[-1])
	plt.legend()
# =============================================================================
# 	plt.show()
# =============================================================================
	plt.savefig(gcomplex+"-curves.png")
	
if __name__ == "__main__":
	from matplotlib import pyplot as plt
	import numpy as np
	import sys
	if len(sys.argv) <3:
		print """Provide more argments (3) please. Use on GUEST FILES ONLY
		Use make_avg_coords.py slice fname or 
		make_avg_coords.py kin fname"""

	elif len(sys.argv) >= 3 and sys.argv[1] == "slice":
		print "will read the guest and compute the kinetic energy for the guest, the fragments and the complex"
		name = sys.argv[2]
		split_velocity_file(name)
	elif len(sys.argv) >= 3 and sys.argv[1] == "kin":
		print "will read the guest and compute the kinetic energy for the guest, the fragments and the complex"
		name = sys.argv[2]
		convert_position_snapshots_to_velocity_snapshots(name)
	elif len(sys.argv) == 3 and sys.argv[1] == "plot":
		print "will read the guest and compute the kinetic energy for the guest, the fragments and the complex"
		number = sys.argv[2].split('-')[0]
		fcomplex = "/home/macenrola/Documents/Thermal_energy_collector/known-guest-pdbs/OUT/COMPLEXES/{}-orig.pdb.xyz_a_0.0_s_0.0_0.0_docked.xyz-pos-1.xyz-EKIN".format(number)
		gcomplex = "/home/macenrola/Documents/Thermal_energy_collector/known-guest-pdbs/OUT/SINGLE_GUESTS/{}-orig.pdb.xyz-pos-1.xyz-EKIN".format(number)
		smallfrag = "/home/macenrola/Documents/Thermal_energy_collector/known-guest-pdbs/OUT/COMPLEXES/{}-orig.pdb.xyz_a_0.0_s_0.0_0.0_docked.xyz-pos-1.xyz-small_frag.xyz-EKIN".format(number)
		bigfrag = "/home/macenrola/Documents/Thermal_energy_collector/known-guest-pdbs/OUT/COMPLEXES/{}-orig.pdb.xyz_a_0.0_s_0.0_0.0_docked.xyz-pos-1.xyz-large_frag.xyz-EKIN".format(number)
			
		plot_summary_EKIN(fcomplex, gcomplex, smallfrag, bigfrag)
		