#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 17 18:17:06 2020

@author: macenrola
"""

def get_block_between_tags(fname, begin_tag, end_tag):
	"""
	PRE: Takes in a file and two tags, 
	POST: Will return the text between the first occurence of these two tags, between the first and the second occurence if the begin and end tag is the same
	"""
	with open(fname, "rt") as r: 
		lines=r.readlines()
	if begin_tag not in lines or end_tag not in lines:
		print("No cubic data, returning None")
		return None
	begin = lines.index(begin_tag)
	end = lines[begin+1:].index(end_tag)
	return lines[begin: begin+end+1]

def get_data_from_out(outfile):
	"""
	PRE  : The output comes from an input file with the route `#n PM6D3 Opt=(Restart) Freq=(Anharmonic, PrintDerivatives, InternalModes, cubic, HPModes)`
		The starting geometry should probably be minimised to at least Tight or very tight in Gaussian
	POST : Will return the cubic coefficients, the cartesian modes 
		Double check that the modes coefficients are given to 5 decimals
	"""
	tag_modes = " Harmonic frequencies (cm**-1), IR intensities (KM/Mole), Raman scattering\n"
	mode_block = get_block_between_tags(outfile, tag_modes, tag_modes)
# =============================================================================
# 	print "".join(mode_block)
# =============================================================================
	# Then gets the mode numbers and positions in the block based on the "A" sequencies in the block
	AA_indices = [i-1 for i, j in enumerate(mode_block) if "Frequencies ---" in j]
	spacing = AA_indices[1]-AA_indices[0]
	mode_dic={}
	for i in AA_indices:
# =============================================================================
# 		print "".join(mode_block[i-1:i-1+spacing])
# =============================================================================
		mode_num = mode_block[i-1].strip().split()
		meta = [x.strip().split() for x in mode_block[i+1:i+5]] # gets the frequency, reduced masses, force constants and IR intensities
		meta_vals = list(map(list, zip(*[x[x.index("---")+1:] for x in meta]))) # gets the values only, by mode
		mode_vals = list(map(list, zip(*[x.strip().split() for x in mode_block[i+6:i+spacing-1]]))) #gets the mode data along with the kind of atom
		for k, i in enumerate(mode_num):
			mode_dic[i] = (meta_vals[k], mode_vals[3+int(k)]) # adds the meta values and the mode values to the dictionary
		mode_dic["order"] = mode_vals[:3]
		
	cubic_dic = {}
	begin_cubic = " :        CUBIC FORCE CONSTANTS IN NORMAL MODES         :\n"
	end_cubic = " :       QUARTIC FORCE CONSTANTS IN NORMAL MODES        :\n"
	cubic_block = get_block_between_tags(outfile, begin_cubic, end_cubic)
	if cubic_block is not None:
	# =============================================================================
	# 	print "".join(cubic_block)
	# =============================================================================
		for i, line in enumerate(cubic_block[9:-5]):
			parts = line.strip().split()
			cubic_dic["{}-{}-{}".format(*parts[:3])] = parts[3:6]
		#print "".join(cubic_block)
	
		return mode_dic, cubic_dic
	else:
		return mode_dic




def get_mode_localisation_on_xylene(mode, mode_xyl_atoms):
	"""
	PRE: Takes in modes formatted as by get_data_from_out
	POST: Will print the sum of amplitudes corresponding to the xylene moiety, they are assumed to be located first in the list of atoms
	Observation: The xylene guest tends to have more movement than expected from its isolated form. It gets a larger share of kinetic energy from the NM distribution than when in the vacuum
	Conversely the CB gets less of that energy. The sum of amplitude of normed normal modes is 70 in the mxylene-h-single-CB6 system while 48 in vacuum
	Here are the results for energy localisation, fascinating
	mxylene-h-single-CB6.com_OUT.out.xyz.com_OUT.out.xyz.com_OUT.out.xyz.com_OUT.out
	69.9628961769
	mxylene-h-single-CB7.com_OUT.out.xyz.com_OUT.out.xyz.com_OUT.out.xyz.com_OUT.out
	66.8494180851
	mxylene-protonated-CB6.com_OUT.out.xyz.com_OUT.out.xyz.com_OUT.out.xyz.com_OUT.out
	71.2384169403
	mxylene-protonated-CB7.com_OUT.out.xyz.com_OUT.out.xyz.com_OUT.out.xyz.com_OUT.out
	68.9233001733
	mxylene-CB6.com_OUT.out.xyz.com_OUT.out.xyz.com_OUT.out.xyz.com_OUT.out
	64.3905017355
	mxylene-CB7.com_OUT.out.xyz.com_OUT.out.xyz.com_OUT.out.xyz.com_OUT.out
	62.4629103182
	"""
	xyl_localised_modes=[]
	for k in sorted(mode)[:-1]: # to avoid the order 
		LOC_ON_XYL_COEFF=sum([float(x)**2 for x in mod[k][1][:mode_xyl_atoms*3-1]])
		xyl_localised_modes.append((LOC_ON_XYL_COEFF, int(k), mod[k][0][0]))
		
	return sorted(xyl_localised_modes, key=lambda x: x[1])

def get_mode_locatlization_with_COMP_GUEST(comp_out_file,  guest_out_file):
	"""
	PRE  : 	Takes the outfiles from gaussian using at least: #n PM6D3 opt=(maxstep=999, maxcyc=999) maxdisk=100GB freq=(PrintDerivatives, InternalModes, HPModes)
			It is assumed that the guest atoms are presented first in the xyz coordinates of the complex
			The atom number of the guest will be taken from the 'order' field in the mode_dic
			
	POST : 	Will use the get_data_from_out method and show the degree of mode localization of the complex on the xylene guest
			
	"""
	mode_comp  =  get_data_from_out(comp_out_file) # get the complex modes
	mode_guest =  get_data_from_out(guest_out_file) # get the guest modes
	n_atoms_guest = int(mode_guest['order'][1][-1]) # get the guest atom number 
	guest_localised_modes=[]
	for k in sorted(mode_comp)[:-1]: # to avoid the order  
		LOC_ON_GUEST_COEFF=sum([float(x)**2 for x in mode_comp[k][1][:n_atoms_guest*3-1]]) # compute the localization number
		guest_localised_modes.append((LOC_ON_GUEST_COEFF, int(k), mode_comp[k][0][0])) #stores the localisation number, the mode number and the frequency
	print("{} atoms and {} modes".format(n_atoms_guest, len(mode_guest.keys())-1))
	res = sorted(guest_localised_modes, key=lambda x: x[0])
	sorted_coeffs = [x[0] for x in res] 
	total_loc = sum(sorted_coeffs) # total number of localization for the guest
	normal_dist = [0]*(len(sorted_coeffs)-n_atoms_guest*3-6)+[1]*(n_atoms_guest*3-6) # if only the 3N-6 modes of the guest were located on it 
	plt.plot(sorted_coeffs, label="guest int mode #={0:4.2f}".format(total_loc))
	plt.plot(normal_dist, label="guest 3N-6={0:4.2f}".format(sum(normal_dist)))
	plt.xlabel("mode #, sorted by guest loc coeff")
	plt.ylabel("mode loc coeff")
	plt.title(guest_out_file.split("/")[-1].split("-")[0])
	plt.legend()
	plt.savefig(guest_out_file+".png")
	plt.close()
	
	
	res = sorted(guest_localised_modes, key=lambda x: x[0])
	sorted_coeffs = [x[0] for x in res] 
	normal_dist = [0]*(len(sorted_coeffs)-n_atoms_guest*3-6)+[1]*(n_atoms_guest*3-6) # if only the 3N-6 modes of the guest were located on it 
	total_loc = get_approx_cv(mode_comp, sorted_coeffs)
	guest_loc = get_approx_cv(mode_guest, [1]* (n_atoms_guest*3-6))
	print("{0:2.2f} tot vs {1:2.2f} just guest === Ratio {2:1.2f}".format(total_loc, guest_loc, guest_loc/total_loc))
# =============================================================================
# 	plt.plot(sorted_coeffs, label="guest int mode #={0:4.2f}".format(total_loc))
# 	plt.plot(normal_dist, label="guest 3N-6={0:4.2f}".format(guest_loc))
# 	plt.xlabel("mode #, sorted by guest loc coeff")
# 	plt.ylabel("mode loc coeff")
# 	plt.title(guest_out_file.split("/")[-1].split("-")[0])
# 	plt.legend()
# 	plt.show()
# =============================================================================
	
def get_approx_cv(modic, coeffs):
	"""
	PRE  : Takes a mode dictionary and the coefficients of localisation, for a lone molecule, all coefficients should be one
	POST : Takes the frequencies in cm-1 and computes an approximative CV using the formula page 137 of mcquarrie 1976 statistical mechanics
	"""
	c=2.99e10
	h = 6.62607004e-34
	kb = 1.38064852e-23 
	T = 300
	approx_cv = 0
	for k, loc in zip(sorted(modic)[:-1], coeffs): # to avoid the order  
		freq = float(modic[k][0][0])
		stemp =freq*c*h/kb
		approx_cv += loc*(stemp/T)**2 * math.exp(stemp/T)/(math.exp(stemp/T)-1)**2
	return approx_cv

if __name__ == "__main__":
	import matplotlib.pyplot as plt
	import glob
	import math
	#print(get_data_from_out("/home/macenrola/Documents/Thesis/thermal_collector/VIBS/GUESTS/156-orig.pdb.xyz.com_OUT.out")['order'][1][-1])
	for f in sorted(glob.glob("/home/macenrola/Documents/Thesis/thermal_collector/VIBS/GUESTS/*orig.pdb.xyz.com_OUT.out")):
		print(f)
		num=f.split("/")[-1].split("-")[0]
		print(num)
# =============================================================================
# 		get_mode_locatlization_with_COMP_GUEST("/home/macenrola/Documents/Thesis/thermal_collector/VIBS/COMPLEXES/{}-orig.pdb.xyz_a_0.0_s_0.0_0.0_docked.xyz.com_OUT.out".format(num), "/home/macenrola/Documents/Thesis/thermal_collector/VIBS/GUESTS/{}-orig.pdb.xyz.com_OUT.out".format(num))
# =============================================================================

		try:
			get_mode_locatlization_with_COMP_GUEST("/home/macenrola/Documents/Thesis/thermal_collector/VIBS/COMPLEXES/{}-orig.pdb.xyz_a_0.0_s_0.0_0.0_docked.xyz.com_OUT.out".format(num), "/home/macenrola/Documents/Thesis/thermal_collector/VIBS/GUESTS/{}-orig.pdb.xyz.com_OUT.out".format(num))
		except:
			print("passing {}".format(f))
