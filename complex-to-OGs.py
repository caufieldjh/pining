#!/usr/bin/python
#pining_for_new_data.py
'''
A small script for converting proteins complexes as provided by the 
IntAct Complex Portal into sets of eggNOG OGs, while retaining
complex names. 
If provided with a raw edgelist (one line per edge, two tab-delimited
interactors, preferably eggNOG OGs), also identifies matching complexes.

Will be integrated into the main pining project files.
'''

import argparse
import itertools
import glob, gzip, os, sys, urllib2
from datetime import date
from tqdm import *

#Methods
def find_complexes_in_edgelist(edgelist_filename, og_cplxs):
	'''Given the filename of an edgelist (one line per edge, two 
	tab-delimited interactors, preferably eggNOG OGs), identifies all
	complexes with at least one potentially matching interaction 
	(assuming all-vs-all or "star" interaction model).
	This is just based on set membership; no weight is provided for
	edge adjacency or shortest paths involving complex members.'''
	
	output_filename = "complex_matches.txt"
	
	#Parse the edgelist input first, assigning each edge a numerical ID
	edges = {}
	i = 0
	with open(edgelist_filename) as edgelist_file:
		for line in edgelist_file:
			splitline = (line.rstrip()).split("\t")
			int_a = splitline[0]
			int_b = splitline[1]
			edges[i] = (int_a, int_b)
			i = i+1
			
	#Set up the score dict.
	#Scores are raw matches first, fraction of total matches second
	cplx_match_scores = {}
	for cplx in og_cplxs:
		cplx_match_scores[cplx] = (0, 0)
	
	#Iterate through complexes, parse them as sets of possible interactions,
	#then for each interaction, try to find it or its reciprocal within
	#the provided edgelist.
	for cplx in og_cplxs:
		raw_score = 0 #Total matches, even if it's the same interaction
		unique_score = 0 #Total unique matches
		matched_already = []
		#print(cplx)
		these_poss_interactions = set(itertools.combinations(og_cplxs[cplx], 2))
		#print(these_poss_interactions)
		total_poss_int = float(len(these_poss_interactions))
		#print("Total possible interactions: %s" % total_poss_int)
		for poss_interaction in these_poss_interactions:
			for edge in edges:
				if edges[edge] in (poss_interaction, reversed(poss_interaction)):
					#print("Found match: %s" % str(poss_interaction))
					raw_score = raw_score +1
					if poss_interaction not in matched_already:
						unique_score = unique_score +1
						matched_already.append(poss_interaction)
		if total_poss_int > 0:
			fract_score = unique_score / total_poss_int
		else:
			fract_score = 0
		cplx_match_scores[cplx] = (raw_score, fract_score)
	
	#write the output file
	with open(output_filename, 'w+b') as outfile:
		outfile.write("Cplx_ID\tRaw\tFract\n")
		for cplx in cplx_match_scores:
			outfile.write("%s\t%s\t%s\n" % (cplx, cplx_match_scores[cplx][0],
											cplx_match_scores[cplx][1]))
	
	return output_filename
	
def get_eggnog_maps(): 
	'''Download and unzip the eggNOG 4.5 ID conversion file
	and the appropriate NOG files.
	Filters file to just Uniprot IDs; the resulting file is the map file.
	This function does not map the target proteins in the user's 
	input file. That happens elsewhere.'''
	
	baseURL = "http://eggnogdb.embl.de/download/eggnog_4.5/"
	convfilename = "eggnog4.protein_id_conversion.tsv.gz"	
	#File contains ALL database identifiers and corresponding proteins
	
	convfilepath = baseURL + convfilename
	outfilepath = convfilename[0:-3]
	dl_convfile = 1	#If 1, we need to download
	if os.path.isfile(convfilename): 
		#Already have the compressed file, don't download
		print("Found compressed ID conversion file on disk: %s" % convfilename)
		decompress_convfile = 1
		dl_convfile = 0
	if os.path.isfile(outfilepath): 
		#Already have the decompressed file, don't download
		print("Found ID conversion file on disk: %s" % outfilepath)
		decompress_convfile = 0
		dl_convfile = 0
	
	if dl_convfile == 1:
		print("Downloading ID mapping file - "
			"this file is ~400 Mb compressed so this may take some time.")
		print("Downloading from %s" % convfilepath)
		response = urllib2.urlopen(convfilepath)
		compressed_file = open(os.path.basename(convfilename), "w+b") 
		#Start local compressed file
		chunk = 1048576
		while 1:
			data = (response.read(chunk)) #Read one Mb at a time
			compressed_file.write(data)
			if not data:
				print("\n%s file download complete." % convfilename)
				compressed_file.close()
				break
			sys.stdout.flush()
			sys.stdout.write(".")
		decompress_convfile = 1
		
	if decompress_convfile == 1:
		print("Decompressing map file. Lines written, in millions:")
		#Done in chunks since it's a large file
		with gzip.open(convfilename) as infile: 
			#Open that compressed file, read and write to uncompressed file
			outfile = open(outfilepath, "w+b")
			linecount = 0
			for line in infile:
				outfile.write(line)
				linecount = linecount +1
				if linecount % 100000 == 0:
						sys.stdout.write(".")
				if linecount % 1000000 == 0:
						sys.stdout.flush()
						sys.stdout.write(str(linecount/1000000))
			infile.close()
		newconvfilename = outfilepath
		outfile.close()
	
	#Download and decompress member NOG files (2 of them)
	#Using base NOGs and mammalian NOGs (maNOG) for now
	nogURL = baseURL + "data/NOG/"
	nogfilename = "NOG.members.tsv.gz"
	eunogURL = baseURL + "data/euNOG/"
	eunogfilename = "euNOG.members.tsv.gz" 
	all_nog_locations = [[nogURL, nogfilename], [eunogURL, eunogfilename]]
	
	for location in all_nog_locations:
		baseURL = location[0]
		memberfilename = location[1]
		memberfilepath = baseURL + memberfilename
		outfilepath = memberfilename[0:-3]
		if os.path.isfile(memberfilename): 
			print("\nFound compressed NOG membership file on disk: %s" % memberfilename)
			decompress_memberfile = 1
		if os.path.isfile(outfilepath): 
			print("\nFound NOG membership file on disk: %s" % outfilepath)
			decompress_memberfile = 0
		else:
			print("\nDownloading NOG membership file - this may take some time.")
			print("Downloading from %s" % memberfilepath)
			response = urllib2.urlopen(memberfilepath)
			compressed_file = open(os.path.basename(memberfilename), "w+b") 
			#Start local compressed file
			chunk = 1048576
			while 1:
				data = (response.read(chunk)) #Read one Mb at a time
				compressed_file.write(data)
				if not data:
					print("\n%s file download complete." % memberfilename)
					compressed_file.close()
					break
				sys.stdout.flush()
				sys.stdout.write(".")
			decompress_memberfile = 1
			
		if decompress_memberfile == 1:
			print("Decompressing NOG membership file %s" % memberfilename)
			#Done in chunks since it's a large file
			with gzip.open(memberfilename) as infile: #Open that compressed file, read and write to uncompressed file
				outfile = open(outfilepath, "w+b")
				linecount = 0
				for line in infile:
					outfile.write(line)
					linecount = linecount +1
					if linecount % 100000 == 0:
						sys.stdout.write(".")
					if linecount % 1000000 == 0:
						sys.stdout.flush()
						sys.stdout.write(str(linecount/1000000))
				infile.close()
			outfile.close()
			
	#Clean up by removing compressed files
	print("\nRemoving compressed files.")
	all_compressed_files = [convfilename, nogfilename, eunogfilename]
	for filename in all_compressed_files:
		if os.path.isfile(filename):
			os.remove(filename)
	
	#Load and filter the ID conversion file as dictionary
	print("Parsing ID conversion file. Lines read, in millions:")
	with open(convfilename[0:-3]) as infile:
		id_dict = {}	#Dictionary of eggNOG protein IDs with database IDs as keys
		#Gets filtered down to relevant database IDs (i.e., Uniprot IDs)
		linecount = 0
		for line in infile:
			linecount = linecount +1
			line_raw = ((line.rstrip()).split("\t"))	#Protein IDs are split for some reason; merge them
			one_id_set = [line_raw[0] + "." + line_raw[1], line_raw[2], line_raw[3]]
			if "UniProt_AC" in one_id_set[2]:
				id_dict[one_id_set[1]] = one_id_set[0]
			if linecount % 100000 == 0:
				sys.stdout.write(".")
			if linecount % 1000000 == 0:
				sys.stdout.flush()
				sys.stdout.write(str(linecount/1000000))
		infile.close()

	#Use filtered ID conversion input to map to NOG members
	print("\nReading NOG membership files.")
	all_nog_filenames = [nogfilename[0:-3], eunogfilename[0:-3]]
	nog_members = {}	
	#Dictionary of NOG ids with protein IDs as keys (need to split entries for each)
	nog_count = 0
	for filename in all_nog_filenames:
		temp_nog_members = {}	
		#We will have duplicates within each set but don't want to lose the information.
		print("Reading from %s" % filename)
		with open(filename) as infile:
			for line in infile:
				nog_count = nog_count +1
				line_raw = ((line.rstrip()).split("\t"))
				nog_id = line_raw[1]
				line_members = line_raw[5].split(",")
				for protein_id in line_members:			
					#The same protein could be in more than one OG at the same level
					if protein_id in temp_nog_members:
						temp_nog_members[protein_id] = temp_nog_members[protein_id] + "," + nog_id
					else:
						temp_nog_members[protein_id] = nog_id
			infile.close()
		nog_members.update(temp_nog_members)
	
	upids_length = str(len(id_dict))
	nogs_length = str(nog_count)
	proteins_length = str(len(nog_members))
	
	print("Mapping %s Uniprot IDs to %s NOGs through %s eggNOG protein IDs:" % (upids_length, nogs_length, proteins_length))
	upid_to_NOG = {}	#Conversion dictionary. Values are OGs, keys are UPIDs.
	mapped_count = 0	#upids mapped to nogs.
	for upid in id_dict:
		if id_dict[upid] in nog_members:
			upid_to_NOG[upid] = nog_members[id_dict[upid]]
			mapped_count = mapped_count +1
			if mapped_count % 100000 == 0:
				sys.stdout.write(".")
			if mapped_count % 1000000 == 0:
				sys.stdout.flush()
				sys.stdout.write(str(mapped_count/1000000))
		
	#Use this mapping to build map file, named "uniprot_og_maps_*.txt"
	print("Writing map file.")
	nowstring = (date.today()).isoformat()
	mapfilename = "uniprot_og_maps_" + nowstring + ".txt"

	mapfile = open(mapfilename, "w+b")
	for mapping in upid_to_NOG:
		mapfile.write(mapping + "\t" + upid_to_NOG[mapping] + "\n")	
		#Each line is a uniprot ID and an OG id
	mapfile.close()
	print("Map file complete.")

def load_cplxs(cplx_filename_list):
	'''Loads three items from IntAct Complex Portal files for each complex:
	The complex accession ID, the recommended name, and the list of
	all proteins in the complex.
	Returns a dict with accessions as keys; values are complex names
	and lists of UPIDs (without stoichiometry).'''
	
	cplxs = {}
	
	for filename in cplx_filename_list:
		with open(filename) as cplx_file:
			cplx_file.readline() #Skip header
			for line in cplx_file:
				splitline = (line.rstrip()).split("\t")
				ac = splitline[0]
				name = splitline[1]
				upids = []
				for value in splitline[4].split("|"):
					upids.append((value.split("("))[0]) #Remove stoich.
				cplxs[ac] = [name,upids]
	
	return cplxs

def make_output(cplxs, og_map, unmapped, outfilename):
	'''Produces tsv file of complex IDs, names, and OG members.
		Also provides dict of complexes and their members.'''
	
	outdict = {}
	
	with open(outfilename, 'w+b') as outfile:
		for cplx in cplxs:
			these_ogs = []
			for upid in cplxs[cplx][1]:
				if upid in unmapped:
					these_ogs.append(upid)
				else:
					these_ogs.append(og_map[upid])
			these_ogs = list(set(these_ogs)) #Remove duplicates
			lineout = "%s\t%s\t%s\n" % (cplx, cplxs[cplx][0], ",".join(these_ogs))
			outfile.write(lineout)
			
			outdict[cplx] = these_ogs
			
	return outdict
	
def map_prots_to_ogs(ids):
	'''Uses eggNOG files to create a dict of protein IDs and 
	corresponding OGs.
	Takes list of Uniprot IDs as input.'''
	
	prot_OG_maps = {} #Dictionary of all Uniprot ID to OG maps
	target_OG_maps = {} 
	#Dictionary to save protein-OG mapping specific for this input set
	
	prots_without_OG = [] #UPIDs of proteins w/o corresponding OGs
	
	mapfile_list = glob.glob("uniprot_og_maps_*")
	if len(mapfile_list) == 0:
		sys.exit("Can't find OG mapping file. Exiting...")
	else:
		mapfilename = mapfile_list[0]
	
	#Load proteins to OG map file
	print("Loading entries from map file...")
	with open(mapfilename) as mapfile:
		filelen = sum(1 for line in mapfile) -1
		mapfile.seek(0)
		pbar = tqdm(total=filelen)
		for line in mapfile:
			splitline = (line.rstrip()).split()
			upid = splitline[0]
			og = splitline[1]
			prot_OG_maps[upid] = og
			pbar.update(1)
			
	pbar.close()
	
	#Now look up input UPIDs
	for upid in ids:
		if upid in prot_OG_maps:
			matching_og = prot_OG_maps[upid]
		else:
			matching_og = upid	
			#If the protein doesn't map to an OG it retains its upid
			prots_without_OG.append(upid)
				
		target_OG_maps[upid] = matching_og
				
	return target_OG_maps, prot_OG_maps, prots_without_OG

#Main
def main():
	
	print(sys.argv[0])
	
	parser = argparse.ArgumentParser()
	parser.add_argument('--edgelist', help="a filename to use for the "
						"edgelist, as a list of two columns where "
						"each line is an edge and each item is an "
						"eggNOG OG ID (for each interactor)")
	args = parser.parse_args()
	
	if args.edgelist:
		edgelist_filename = args.edgelist
		print("Will use %s as edgelist input file." % edgelist_filename)
		
	intact_cplx_list = glob.glob('intact*.tsv')
	if len(intact_cplx_list) == 0:
		sys.exit("Found no complex input files.\n"
						"Please name all files so they begin with \'intact\'.")
	else:
		print("Found the following complex files:")
		for filename in intact_cplx_list:
			print(filename)
	
	#Load contents of complex files
	print("Loading complexes and their proteins.")
	prot_cplxs = load_cplxs(intact_cplx_list)
	cplx_count = len(prot_cplxs)
	print("Input contains %s complexes." % (cplx_count))
	
	#Get flat UPID list from complexes
	prot_ids = []
	for cplx in prot_cplxs:
		for upid in prot_cplxs[cplx][1]:
			prot_ids.append(upid)
	prot_ids = list(set(prot_ids)) #Remove duplicates
	prot_id_count = len(prot_ids)
	print("Input contains %s unique complex member IDs.\n" % (prot_id_count))
			
	#Locate Uniprot to OG map files
	#Download the corresponding files if needed
	print("Preparing orthologous group membership files.")
	
	mapping_file_list = glob.glob('uniprot_og_maps*.txt')
	if len(mapping_file_list) > 1:
		sys.exit("Found more than one Uniprot to eggNOG OG mapping file. Check for duplicates.")
	if len(mapping_file_list) == 0:
		print("No Uniprot to eggNOG mapping files found or they're incomplete. Rebuilding them.")
		get_eggnog_maps()
	if len(mapping_file_list) == 1:
		print("Found the Uniprot ID mapping file:")
		for filename in mapping_file_list:
			print(filename)
	
	print("Mapping target proteins to orthologous groups.")
	og_map, all_og_map, unmapped = map_prots_to_ogs(prot_ids)
	
	unmapped_count = len(unmapped)
	mapped_count = len(og_map) - unmapped_count
	all_count = len(all_og_map)
	all_og_count = len(set(all_og_map.values()))
	og_count = len(set(og_map.values()))
	
	print("Full OG map covers %s Uniprot protein IDs "
			"and %s OGs." % (all_count, all_og_count))
	print("Input proteins mapped to %s total OGs." % og_count)
	print("%s proteins from the input mapped to OGs.\n"
			"%s ids did not map to OGs."
			% (mapped_count, unmapped_count))
	#if len(unmapped) > 0:
	#	print("These ids did not map:")
	#for upid in unmapped:
	#	print(upid)
	
	print("Writing output table of complexes and OGs.")
	if len(intact_cplx_list) == 1:
		outfilename = intact_cplx_list[0][:-4] + "_OGs.tsv"
	else:
		outfilename = "complexes_as_OGs.tsv"
	og_cplxs = make_output(prot_cplxs, og_map, unmapped, outfilename)
	print("Wrote results to %s." % outfilename)
	
	if args.edgelist:
		print("Now finding complex matches in edgelist.")
		matched_complex_filename = find_complexes_in_edgelist(edgelist_filename, og_cplxs)
		print("Wrote results to %s." % matched_complex_filename)

if __name__ == "__main__":
	sys.exit(main())
