#!/usr/bin/python
#pining_for_new_data.py
'''
A tool for retrieving protein-protein interaction data sets, using
a list of proteins to search for.
As the number of returned interactions may be quite large,
IntAct and BioGRID are searched as local databases.

Produces combined, annotated versions of the network of the 
retrieved interactions.

'''
__author__= "Harry Caufield"
__email__ = "j.harry.caufield@gmail.com"

import argparse
import glob, gzip, operator, os, random, re, sys, urllib2, zipfile
from datetime import date

import pandas as pd
from pandas import pivot_table

#Constants and Options
directories = ["output","databases", "og_info"]

upid_match = re.compile('[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}')

nowstring = (date.today()).isoformat()

#Classes

#Functions

def load_prot(filename):
	'''
	Obtains the protein IDs to search for from a list of proteins
	along with other data fields, if present.
	Assumes one unique protein per Uniprot ID.
	
	Uses headings if provided, otherwise assigns names to them
	based on whether they are categorical or not.
	
	Returns a dict of protein IDs and a list of lists of associated
	observations. Also returns Pandas data frame.
	Also returns heading names of the observations if available,
	or just assigns generic ones if not provided, and
	determines if they are qualitative or not.
	This is crucial to how the observations will be treated.
	
	'''
	
	def is_name(string):
		#Tests if something looks like a group name or a value
		#Returns True if it's a name, False if not
		this_is_a_name = False
		
		try:
			float(string)
		except ValueError:
			this_is_a_name = True
		
		return this_is_a_name
	
	prot_dict = {} #Protein IDs are keys, observations are values
					#If no observations are present all values are "NA"
	
	groups = [] #Names of experimental groups and observation values
				#List of tuples of group names and types
				#Types are "qual" or "quant"
	
	data = [] #Raw data for setting up data frame
		
	with open(filename) as protfile:
		
		#Check the first line first to see if it looks like a heading
		#or a protein ID, in which case, store it
		
		firstline = protfile.readline()
		splitline = (firstline.rstrip()).split("\t")
		
		if re.match(upid_match, splitline[0]) != None:
			#Looks like a Uniprot ID.
			i = 1
			j = 1
			for heading in splitline[1:]:
				if is_name(heading):
					groups.append(("Group%s" % i, "qual"))
					i = i +1
				else:
					groups.append(("Value%s" % j, "quant"))
					j = j +1
			
			upid = splitline[0]
			prot_dict[upid] = splitline[1:]	
			
		else: #It's not a Uniprot ID, at least, so assume a header line
			for heading in splitline[1:]:
				groups.append((heading, "NA"))
		
		#Load data in file
		#Check on observation format first
		secondline = protfile.readline()
		splitline = (secondline.rstrip()).split("\t")
		i = 1
		new_groups = []
		for heading in groups:
			heading_name = heading[0]
			if is_name(splitline[i]): 
				new_groups.append((heading_name, "qual"))
			else:
				new_groups.append((heading_name, "quant"))
			i = i+1
			
		groups = new_groups
		
		upid = splitline[0]
		if upid not in prot_dict.keys():
			prot_dict[upid] = [splitline[1:]]
		else:
			prot_dict[upid].append(splitline[1:])
		
		#Now read rest of the file		
		for line in protfile:
			splitline = (line.rstrip()).split("\t")
			data.append(splitline)
			upid = splitline[0]
			if upid not in prot_dict.keys():
				prot_dict[upid] = [splitline[1:]]
			else:
				prot_dict[upid].append(splitline[1:])
	
	#Set up data frame
	group_names = ["id"]
	for group in groups:
		group_names.append(group[0])		
	prot_frame = pd.DataFrame(data, columns=group_names, dtype=float)

	return prot_dict, prot_frame, groups

def get_ppi_db(name):
	'''
	Download either IntAct or BioGRID local copies.
	Extracts to the 'databases' folder.
	'''
	
	if name == "IntAct":
		baseURL = "ftp://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/"
		dbfilename = "intact.zip"
		
	if name == "BioGRID":
		baseURL = "https://thebiogrid.org/downloads/archives/Latest%20Release/"
		dbfilename = "BIOGRID-ALL-LATEST.mitab.zip"
		
	dbfilepath = baseURL + dbfilename
	
	#Database names may change
	if name == "IntAct":
		outfilename = "intact.txt"
	elif name == "BioGRID":
		biogrid_db_list = glob.glob('BIOGRID-ALL-*.mitab.txt')
		if len(biogrid_db_list) > 0:
			outfilename = biogrid_db_list[0]
		
	dl_dbfile = True	#If true, we need to download, and this is the default
	if os.path.isfile(dbfilename): 
		#Already have the compressed file, don't download
		print("Found compressed database file on disk: %s" % dbfilename)
		decompress_dbfile = True
		dl_dbfile = False

	if os.path.isfile(outfilename): #A bit redundant for BioGRID
		#Already have the decompressed file, don't download
		print("Found database file on disk: %s" % outfilename)
		decompress_dbfile = False
		dl_dbfile = False
		
	if dl_dbfile:
		print("Downloading %s database file." % name)
		print("Downloading from %s" % dbfilepath)
		response = urllib2.urlopen(dbfilepath)
		#print(response.headers.items())
		
		compressed_file = open(os.path.basename(dbfilepath), "w+b")
		#Start local compressed file
		chunk = 1048576
		while 1:
			data = (response.read(chunk)) #Read one Mb at a time
			compressed_file.write(data)
			if not data:
				print("\n%s file download complete." % dbfilename)
				compressed_file.close()
				break
			sys.stdout.flush()
			sys.stdout.write(".")
		decompress_dbfile = True
		
	if decompress_dbfile:
		print("Decompressing %s file." % name)
		try:
			with zipfile.ZipFile(dbfilename, "r") as infile:
				for filename in infile.namelist():
					if filename not in ['intact_negative.txt']:
						#Ignore some files
						print(filename)
						outfile = open(filename, 'w+b')
						outfile.write(infile.read(filename))
						outfile.close()
						outfilename = filename
						break #Just want one file.
		except zipfile.BadZipfile as e:
			sys.exit("Something is wrong with this database file.\n"
						"Please remove it and re-download.\n"
						"Error: %s" % e)
					
	return outfilename
		
def get_eggnog_maps(): 
	'''
	Download and unzip the eggNOG 4.5 ID conversion file
	and the appropriate NOG files.
	Filters file to just Uniprot IDs; the resulting file is the map file.
	This function does not map the target proteins in the user's 
	input file. That happens elsewhere.
	'''
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
	
def get_eggnog_annotations():
	#Downloads and extracts the eggNOG NOG annotations. 
	
	baseURLs = ["http://eggnogdb.embl.de/download/latest/data/euNOG/",
				"http://eggnogdb.embl.de/download/latest/data/NOG/"]
	euannfilename = "euNOG.annotations.tsv.gz"	
	#The annotations for eukaryote-specific NOGs
	lucaannfilename = "NOG.annotations.tsv.gz"	
	#The annotations for other NOGs, but not bacteria-specific NOGs
	annfilenames = [euannfilename, lucaannfilename]
	
	this_url = 0
	for annfilename in annfilenames:
		annfilepath = baseURLs[this_url] + annfilename
		this_url = this_url +1
		outfilepath = annfilename[0:-3]
		if os.path.isfile(annfilename): 
			print("Found compressed annotation file on disk: " + annfilename)
		else:
			response = urllib2.urlopen(annfilepath)
			print("Downloading from " + annfilepath)
			compressed_file = open(os.path.basename(annfilename), "w+b") 
			#Start local compressed file
			chunk = 1048576
			while 1:
				data = (response.read(chunk)) 
				#Read one Mb at a time
				compressed_file.write(data)
				if not data:
					print("\n" + annfilename + " file download complete.")
					compressed_file.close()
					break
				sys.stdout.flush()
				sys.stdout.write(".")
		
		print("Decompressing annotation file.")
		with gzip.open(annfilename) as infile: 
			#Open that compressed file, read and write to uncompressed file
			file_content = infile.read()
			outfile = open(outfilepath, "w+b")
			outfile.write(file_content)
			infile.close()
		outfile.close()
		
	print("\nRemoving compressed files.")
	all_compressed_files = [euannfilename, lucaannfilename]
	for filename in all_compressed_files:
		os.remove(filename)

def map_ogs_to_annotations(ogs):
	#Uses eggNOG to create a dict of OG IDs and corresponding
	#annotations, as tuples of functional category (FuncCat) and
	#text summary description.
	#Only maps OGs passed to it.
	#Handles combined OGs (e.g., "ENOG1234,ENOG1235") and other IDs
	
	og_note_map = {} #Just annotations for OGs mapping to input proteins
	all_og_notes = {} #All annotations
	
	os.chdir(directories[2])
	
	annfile_list = glob.glob("*annotations.tsv")
	if len(annfile_list) == 0:
		sys.exit("Can't find OG annotation files. Exiting...")
	
	for annfilename in annfile_list:
		with open(annfilename) as annfile:
			for line in annfile:
				splitline = (line.rstrip()).split("\t")
				og = splitline[1]
				funccat = splitline[4]
				desc = splitline[5]
				all_og_notes[og] = (funccat, desc)
	
	for og in ogs:
		if og in all_og_notes:
			og_note_map[og] = all_og_notes[og]
		else:
			split_og = og.split(",")
			if len(split_og) > 1:
				funccats = []
				descs = []
				for this_og in split_og:
					this_funccat, this_desc = all_og_notes[this_og]
					funccats.append(this_funccat)
					descs.append(this_desc)
				funccats = list(set(funccats)) #Unique cats only
				funccats = "".join(funccats)
				descs = list(set(descs)) #Unique descs only
				descs = "|".join(descs)
				og_note_map[og] = (funccats,descs)
			else:
				og_note_map[og] = ("NA","NA")
					
	os.chdir("..")	
	
	return og_note_map, all_og_notes

def map_prots_to_ogs(ids):
	#Uses eggNOG files to create a dict of protein IDs and 
	#corresponding OGs.
	#Takes list of Uniprot IDs as input.
	
	prot_OG_maps = {} #Dictionary of all Uniprot ID to OG maps
	target_OG_maps = {} 
	#Dictionary to save protein-OG mapping specific for this input set
	
	prots_without_OG = [] #UPIDs of proteins w/o corresponding OGs
	
	os.chdir(directories[2])
	
	mapfile_list = glob.glob("uniprot_og_maps_*")
	if len(mapfile_list) == 0:
		sys.exit("Can't find OG mapping file. Exiting...")
	else:
		mapfilename = mapfile_list[0]
	
	#Load proteins to OG map file
	with open(mapfilename) as mapfile:
		for line in mapfile:
			splitline = (line.rstrip()).split()
			upid = splitline[0]
			og = splitline[1]
			prot_OG_maps[upid] = og
	
	#Now look up input UPIDs
	for upid in ids:
		if upid in prot_OG_maps:
			matching_og = prot_OG_maps[upid]
		else:
			matching_og = upid	
			#If the protein doesn't map to an OG it retains its upid
			prots_without_OG.append(upid)
				
		target_OG_maps[upid] = matching_og
			
	os.chdir("..")		
	
	return target_OG_maps, prot_OG_maps, prots_without_OG
	
def search_int_file(ids, filename, db, target_ogs, all_og_map):
	'''
	Searches a PSI-MI TAB format set of interactions (filename) 
	for the Uniprot IDs provided (ids) and for and IDs with the same
	OGs as the provided IDs, using target_ogs and the og map.
	Indexes file first, searches IDs in the index, then returns
	interactions at the specified indices.
	Also returns interactions in short form, as proteins, corresponding
	OGs, and counts/IDs of publications each interaction was observed
	in.
	'''
	
	ids = set(ids) #ids should be unique, plus sets are more efficient
	
	all_int = {} #Just the interactors from each interaction.
					#Keys are line numbers
					#Values are tuples of interactors
	all_og_int = {} #Same as all_int but values are tuples of OGs
					#the proteins belong to.
	
	complete_interactions = [] #List of all matching interactions
								#As split lists
								#Includes all fields provided in database
	prot_match_count = 0 #Count of interactions with at least one
							#match in the target set
	og_match_count = 0	#Count of interactions with at least one match
						#to an OG with members in the target set
	
	ppi_pubs = {} #Dictionary of tuples of protein IDs as keys
						#Values are lists of all pub_ids the interaction
						#or any interactions with OG membership shared
						#for both proteins has been observed.
						
	short_interactions = [] #Short form of interaction list, where each
							#interaction has the form:
							#[protA, protB, OG_A, OG_B, pub_count, pub_ids]
							
	if db == "previous":
		os.chdir(directories[0])
	else:
		os.chdir(directories[1])
		
	if db == "BioGRID":
		#BioGRID has its own ID format.
		#Will need to convert - work in progress.
		pass

	#Just load interactors from each interaction first
	with open(filename) as intactfile:
		#Count the lines in the file first to determine how many PPI
		filelen = sum(1 for line in intactfile) -1
		print("Searching %s interactions." % filelen)
		
		#Progbar
		prog_width = filelen / 10000
		sys.stdout.write("[%s]" % (" " * prog_width))
		sys.stdout.flush()
		sys.stdout.write("\b" * (prog_width+1))
		
		intactfile.seek(0)
		intactfile.readline() #Skip header
		
		i = 1
		#Load all interactor pairs as original IDs and as OGs
		for line in intactfile:
			splitline = (line.rstrip()).split("\t")
			#Only store if two Uniprot IDs
			if (splitline[0].split(":"))[0] == "uniprotkb" and \
				(splitline[1].split(":"))[0] == "uniprotkb":
				interactorA = (splitline[0].split(":"))[1]
				interactorB = (splitline[1].split(":"))[1]
				interaction = (interactorA, interactorB)
				all_int[i] = interaction
				
				if interactorA in all_og_map:
					og_A = all_og_map[interactorA]
				else:
					og_A = interactorA
				
				if interactorB in all_og_map:
					og_B = all_og_map[interactorB]
				else:
					og_B = interactorB
					
				og_interaction = (og_A, og_B)
				all_og_int[i] = og_interaction
			
			i = i +1
			
			if i % 10000 == 0:
				sys.stdout.flush()
				sys.stdout.write("#")
		
		#Now search interactor pairs for target interactors
		get_these_lines = set()
		
		for line_num in all_int:
			for interactor in all_int[line_num]:
				if interactor in ids:
					get_these_lines.add(line_num)
					prot_match_count = prot_match_count +1
					break
					
		for line_num in all_og_int:
			for interactor in all_og_int[line_num]:
				if interactor in target_ogs:
					get_these_lines.add(line_num)
					og_match_count = og_match_count +1
					break
		
		intactfile.seek(0)
		intactfile.readline() #Skip header 	
		
		j = 1
		
		#Add all interactions from the file
		for line in intactfile:
			if j in get_these_lines:
				complete_interactions.append((line.rstrip()).split("\t"))
			j = j+1
		
	#Now construct the compressed interaction dict
	for ppi in complete_interactions:
		interactorA = ppi[0]
		interactorB = ppi[1]
		pub_ids = ppi[8]
		if interactorA != interactorB: #Ignore self interactions
			interactorA = (ppi[0].split(":"))[1]
			interactorB = (ppi[1].split(":"))[1]
			interaction = (interactorA, interactorB)
			rec_interaction = (interactorB, interactorA)
			
			#Only want unique pub_ids
			added = False
			if interaction in ppi_pubs:
				if pub_ids not in ppi_pubs[interaction]:
					ppi_pubs[interaction].append(pub_ids)
					added = True
			if rec_interaction in ppi_pubs and not added:
				if pub_ids not in ppi_pubs[rec_interaction]:
					ppi_pubs[rec_interaction].append(pub_ids)
					added = True
			
			if not added:
				ppi_pubs[interaction] = [pub_ids]
	
	#Add PPI and pub IDs to short interaction list, with OGs
	for ppi in ppi_pubs:
		prot_A = ppi[0]
		prot_B = ppi[1]
		pub_ids = ppi_pubs[ppi]
		pub_count = str(len(pub_ids)) #Convert to string
		pub_ids = ",".join(ppi_pubs[ppi]) #Convert to string
		
		if prot_A in all_og_map:
			og_A = all_og_map[prot_A]
		else:
			og_A = prot_A
		
		if prot_B in all_og_map:
			og_B = all_og_map[prot_B]
		else:
			og_B = prot_B
		
		
		values = [prot_A, prot_B, og_A, og_B, pub_count, pub_ids]
		
		short_interactions.append(values)
											 
	os.chdir("..")
	
	return complete_interactions, prot_match_count, og_match_count, \
			short_interactions
	
def save_prot_list(filename, prot_dict, og_map, og_note_map):
	#Saves the input proteins with OG assignments and annotations
	
	os.chdir(directories[0])
	
	with open(filename, 'w') as outfile:
		#Write the header
		outfile.write("UPID\teggNOG_OG\tFuncCat\tDescription\n")
		
		for upid in prot_dict.keys():
			og = og_map[upid]
			funccat = og_note_map[og][0]
			desc = og_note_map[og][1]
			outfile.write("%s\t%s\t%s\t%s\n" % (upid, og, funccat, desc))
	
	os.chdir("..")
	
def save_interactions(filename, ppi, mode):
	#Saves a set of interactions to a file
	#This function takes a mode argument which may be one of:
	#"tab" - for full PSI-MI TAB format entries
	#"short" - for lists of interactors and publications only
	#These only determine the header
	
	if mode == "tab":
		header = load_psi_tab_header()
	if mode == "short":
		header = "prot_A\tprot_B\tOG_A\tOG_B\tPubCount\tPublications"
	
	os.chdir(directories[0])
	
	with open(filename, 'wb') as outfile:
		#Write the header
		outfile.write("%s\n" % header)
		for interaction in ppi:
			flatline = "\t".join(interaction)
			outfile.write("%s\n" % flatline)
		
	os.chdir("..")
			
def graph_interactions(ids, unique_ppi, unique_og_ppi, prot_dict, 
						og_map, pub_int_counts, og_pub_int_counts,
						con_limit, og_con_limit):
	'''
	Builds the network graph
	This just works with interactions between
	target interactors for now.
	ids is the set of target interactors.
	Produces high confidence graphs as well using counts of PPI
	observations (i.e., if a PPI is seen in multiple publications).
	Uses the OG map to also produce OG network variants.

	con_limit is the number of publications a PPI must be observed in
	to be included in the high confidence graph.
	
	Exports figures of the graph visualizations
	and file of each graph as an edge list, along with the count
	of publications observing each interaction
	(this is a bit muddled with OGs, for now)
	and the med.k values from all_protein_k.txt
	for each mouse strain and group.
	
	The full prot_dict is passed here for future use.
	'''
	
	def draw_graph(nxgraph, name, node_color):
		#Subfunction for drawing graph with pygraphviz.
		A = to_agraph(nxgraph)
		A.node_attr['color'] = node_color
		A.draw(name, prog="sfdp")
		
	ppi = [] #All target-only protein interactions as tuples
	hc_ppi = [] #High confidence interactions, target-only, based
				#on number of publication observed in
	og_int = [] #All target-only, OG-mapped interactions as tuples
				#These are not sets as we are interested in 
				#interactions seen multiple times.
	hc_og_int = [] #High confidence interactions, target-only,
					#OG-mapped, based on number of publication 
					#observed in
	exp_OG = [] #All target-only, OG-mapped interactions as tuples,
				#expanded to multiple interacting species
				#(so some interactions are predictions)
	exp_hc_og_int = [] #All high-condfidence, target-only, OG-mapped 
						#interactions as tuples, expanded to multiple 
						#interacting species 
						#(so some interactions are predictions)
	
	#Set up the value fields we will use based on prot_dict
	
	#Ensuring that the graph only includes PPI among target interactors
	for interaction in unique_ppi:
		if interaction[0] in ids and interaction[1] in ids:
			ppi.append((interaction[0], interaction[1]))
			og_int.append((og_map[interaction[0]], og_map[interaction[1]],
							pub_int_counts[interaction]))
			if pub_int_counts[interaction] >= con_limit:
				hc_ppi.append((interaction[0], interaction[1]))
				hc_og_int.append((og_map[interaction[0]], og_map[interaction[1]],
							pub_int_counts[interaction]))
	
	#The OG-expanded graphs may include non-target interactors
	#(both those shared in OGs and those interacting with those
	#shared in OGs)
	#for now, just because I haven't filtered any further
	for interaction in unique_og_ppi:
		og_interaction = []
		for interactor in interaction:
			try:
				og_interaction.append(og_map[interactor])
			except KeyError:
				og_interaction.append(interactor)
		exp_OG.append((og_interaction[0], og_interaction[1],
							og_pub_int_counts[interaction]))
		if og_pub_int_counts[interaction] >= og_con_limit:
			hc_og_int.append((og_interaction[0], og_interaction[1],
							og_pub_int_counts[interaction]))
		
	for input_type in ["protein", "high_confidence_protein", "OG",
						"high_confidence_OG","exp_OG",
						"exp_high_confidence_OG"]:
		G=nx.Graph()
		pubs = ""
		if input_type == "protein":
			for interaction in ppi:
				G.add_edge(interaction[0], interaction[1], 
							pubs = pub_int_counts[interaction])
			node_color = "red"
		if input_type == "high_confidence_protein":
			for interaction in hc_ppi:
				G.add_edge(interaction[0], interaction[1], 
							pubs = pub_int_counts[interaction])
			node_color = "red"
		if input_type in ["OG", "high_confidence_OG",
								"exp_OG", "exp_high_confidence_OG"]:
			if input_type == "OG":
				for interaction in og_int:
					G.add_edge(interaction[0], interaction[1], 
								pubs = interaction[2])
			elif input_type == "high_confidence_OG":
				for interaction in hc_og_int:
					G.add_edge(interaction[0], interaction[1], 
								pubs = interaction[2])
			elif input_type == "exp_OG":
				for interaction in exp_OG:
					G.add_edge(interaction[0], interaction[1], 
								pubs = interaction[2])
			elif input_type == "exp_high_confidence_OG":
				for interaction in exp_hc_og_int:
					G.add_edge(interaction[0], interaction[1], 
								pubs = interaction[2])
			node_color = "green"
		print("The %s graph has %d nodes, %d edges, and %d connected components." \
			% (input_type, nx.number_of_nodes(G), nx.number_of_edges(G),
				nx.number_connected_components(G)))
		if input_type in ["high_confidence_protein",
							"high_confidence_OG",
							"exp_high_confidence_OG"]:
			print("\tHigh confidence indicates PPI observed in at least "
					"%s different publications (%s for OG-expanded graphs)." 
					% (con_limit, og_con_limit))
		
		print("Building layout and drawing %s graph..." % input_type)
		data_fields = ["pubs"]
		
		os.chdir(outfiledir)
		nx.write_edgelist(G, "heart_%s_edgelist.csv" % input_type, 
							delimiter=",", data=data_fields)
							
		draw_graph(G, "heart_%s_graph.svg" % input_type, node_color)
		os.chdir("..")

def load_psi_tab_header():
	#Loads the PSI-MI TAB header format from a file.
	#Can also get this from IntAct database file 
	#but this way we don't need to read it.
	
	with open("psimitabheader.txt") as format_file:
		header = (format_file.readline()).rstrip()
	
	return header
	
def save_wide_data(filename, prot_frame, groups):
	#Saves a set of observations to a file in wide format
	#i.e. one set of observations per line rather than
	#one observation per line.
	#Uses Pandas data frame
	#Groups will tell us what kind of observations we have.
	
	os.chdir(directories[0])
	
	cols = []
	vals = []
	#Determine indices beyond the id
	for group in groups:
		if group[1] == "qual":
			cols.append(group[0])
		else:
			vals.append(group[0])
	
	prot_frame_wide = pivot_table(prot_frame, values=vals, 
									index='id', columns=cols)
	
	prot_frame_wide.to_csv(filename, sep='\t', na_rep='NA')
		
	os.chdir("..")
			
#Main
def main():
	prot_dict = {} #Keys are unique protein IDs.
					#Values are lists of lists of additional data,
					#if present.
	
	outputs = {} #Names of output files
	
	print(sys.argv[0])
	
	parser = argparse.ArgumentParser()
	parser.add_argument('-inputfile', help="name of a text file containing "
						"UniprotAC IDs of proteins to search for, "
						"one protein ID per line.")
	args = parser.parse_args()
	
	#Set up the output and database storage directories
	for directory in directories:
		if not os.path.isdir(directory):
			print("Setting up directory for %s." % directory)
			os.mkdir(directory)
	
	#Load protein ID input file
	if args.inputfile:
		try:
			have_obs = False
			
			print("Loading %s as input file." % args.inputfile)
			protfilename = args.inputfile
			prot_dict, prot_frame, groups = load_prot(protfilename)
			prot_ids = prot_dict.keys()
			group_count = len(groups)
			
			print("Loaded %s unique protein IDs." % len(prot_ids))
			
			if group_count > 0: 
				qual_groups = []
				quant_groups = []
				for group in groups:
					if group[1] == "qual":
						qual_groups.append(group[0])
					if group[1] == "quant":
						quant_groups.append(group[0])
				
				print("Input contains %s types of observations:" % group_count)
				print("%s (%s) are categorical and" % (len(qual_groups), qual_groups))
				print("%s (%s) are numerical.\n" % (len(quant_groups), quant_groups))
				
				if len(qual_groups) > 0 and len(quant_groups) > 0:
					have_obs = True
				
		except IOError as e:
			sys.exit("Can't find or open that file...\n%s" % e)
	else:
		sys.exit("No input file provided.")
	
	if have_obs:
		#Write the output protein data file
		#This is the input data in wide format, essentially
		print("Converting input data observations to wide form.")
		dataoutfilename = "%s_wide.txt" % protfilename[0:-4]
		save_wide_data(dataoutfilename, prot_frame, groups)
		outputs["Input data in wide form"] = dataoutfilename
	
	#Locate PPI database files
	#Download if necessary.
	os.chdir(directories[1])
	intact_db_list = glob.glob('intact*.txt')
	if len(intact_db_list) > 1:
		sys.exit("Found more than one copy of the IntAct database.\n"
						"Check for duplicates.")
						
	if len(intact_db_list) == 0:
		print("No IntAct database found. Downloading new copy.")
		intactfilename = get_ppi_db("IntAct")
		
	if len(intact_db_list) == 1:
		print("Found IntAct database file.")
		intactfilename = intact_db_list[0]
	
	intactfile_loc = os.path.join(directories[1], intactfilename)
		
	biogrid_db_list = glob.glob('BIOGRID-ALL-*.mitab.txt')
	if len(biogrid_db_list) > 1:
		sys.exit("Found more than one copy of the BioGRID database.\n"
					"Check for duplicates.")
					
	if len(biogrid_db_list) == 0:
		print("No BioGRID database found. Downloading new copy.")
		biogridfilename = get_ppi_db("BioGRID")
		
	if len(biogrid_db_list) == 1:
		print("Found BioGRID database file.")
		biogridfilename = biogrid_db_list[0]
	
	biogridfile_loc = os.path.join(directories[1], biogridfilename)
		
	os.chdir("..")
		
	#Locate Uniprot to OG map files
	#Download the corresponding files if needed
	print("Preparing orthologous group membership and annotation files.")
	
	os.chdir(directories[2])
	mapping_file_list = glob.glob('uniprot_og_maps*.txt')
	if len(mapping_file_list) > 1:
		sys.exit("Found more than one Uniprot to eggNOG OG mapping file. Check for duplicates.")
	if len(mapping_file_list) == 0:
		print("No Uniprot to eggNOG mapping files found or they're incomplete. Rebuilding them.")
		get_eggnog_maps()
		
	annotation_file_list = glob.glob('*annotations.tsv')
	if len(annotation_file_list) == 0:
		print("No eggNOG annotation files found or they're incomplete. Retrieving them.")
		get_eggnog_annotations()
			
	os.chdir("..")
	
	#Map proteins to OGs
	#og_map includes only target protein IDs as keys
	#all_og_map includes ALL protein IDs.
	#We use this second map to get OGs for interactors outside
	#the target set.
	
	print("Mapping target proteins to orthologous groups.")
	og_map, all_og_map, unmapped = map_prots_to_ogs(prot_ids)
	
	#Map OGs to their annotations
	these_ogs = set(og_map.values())
	og_note_map, all_og_notes = map_ogs_to_annotations(these_ogs)
	
	unmapped_count = len(unmapped)
	mapped_count = len(og_map) - unmapped_count
	all_count = len(all_og_map)
	all_og_count = len(set(all_og_map.values()))
	og_count = len(set(og_map.values()))
	
	print("Full OG map covers %s Uniprot protein IDs "
			"and %s OGs." % (all_count, all_og_count))
	print("Input proteins mapped to %s total OGs." % og_count)
	print("%s proteins from the input mapped to OGs.\n"
			"%s proteins did not map to OGs."
			% (mapped_count, unmapped_count))
	if len(unmapped) > 0:
		print("These proteins did not map:")
	for upid in unmapped:
		print(upid)
	
	#Now write the output protein list file
	#with one uniprot ID and its corresponding OG per line
	simpleoutfilename = "%s_prots_and_ogs.txt" % protfilename[0:-4]
	save_prot_list(simpleoutfilename, prot_dict, og_map, og_note_map)
	outputs["Proteins in the input and corresponding OGs"] \
			= simpleoutfilename
	
	#Look for filtered interaction output if we already have it
	#If we don't have an output file, use the IntAct set.
	ppioutfilename = "%s_matching_ppi.txt" % protfilename[0:-4]
	shortoutfilename = "%s_matching_short.txt" % protfilename[0:-4]
	ppioutfilepath = os.path.join(directories[0], ppioutfilename)
	
	if not os.path.isfile(ppioutfilepath):
		have_output = False
		print("Did not find previous output file.")
		print("Searching IntAct interactions for provided protein IDs.")
		db = "IntAct"
		complete_interactions, prot_match_count, og_match_count, \
			short_interactions \
			= search_int_file(prot_ids, intactfilename, db, these_ogs, all_og_map)
		#print("Searching BioGRID interactions for provided protein IDs.")
		#db = "BioGRID"
		#complete_interactions, prot_match_count, og_match_count, \
			#short_interactions \
			#= search_int_file(prot_ids, biogridfilename, db, these_ogs, all_og_map)
			
		##This should append interactions rather than replacing them
		##And then will need to check for duplicates.
	else:
		have_output = True
		print("Found existing interaction output file: %s" % ppioutfilename)
		print("Analyzing previous interaction output.")
		db = "previous"
		complete_interactions, prot_match_count, og_match_count, \
			short_interactions \
			= search_int_file(prot_ids, ppioutfilename, db, these_ogs, all_og_map)
	
	ppi_count = len(complete_interactions)
	short_count = len(short_interactions)
	
	print("\nFound %s total interactions involving target protein IDs "
			"or shared OGs." % ppi_count)
	print("%s interactions involve at least one target protein." % prot_match_count)
	print("%s interactions involve at least one matching OG." % og_match_count)
	print("%s unique protein interactions after compression." % short_count)
	
	#Save the matching set if we don't have one yet
	if not have_output:
		print("\nSaving matching interactions...")
		save_interactions(ppioutfilename, complete_interactions, "tab")
		outputs["Matching complete interactions in PSI-MI TAB format"] \
			= ppioutfilename
		print("\nSaving short interaction file with publication counts...")
		save_interactions(shortoutfilename, short_interactions, "short")
		outputs["Matching short interactions with OGs and publication counts"] \
			= shortoutfilename
	
	print("Complete. See output files:\n")
	for output_type in outputs:
		print("%s: %s\n" % (output_type, outputs[output_type]))
			
	sys.exit()
	
	#Testing cutoff
	
	#Get top lists of them
	#Get counts first
	taxo_counts = {key: len(value) for key, value in taxo_dict.items()}
	pub_counts = {key: len(value) for key, value in pub_dict.items()}
	taxo_int_counts = {key: len(value) for key, value in taxo_ints.items()}
	pub_int_counts = {key: len(value) for key, value in pub_ints.items()}

	high_taxo_dict = sorted(taxo_counts.items(), key=operator.itemgetter(1),
								reverse=True)[0:15]
	high_pub_dict = sorted(pub_counts.items(), key=operator.itemgetter(1),
								reverse=True)[0:15]
	high_taxo_int_dict = sorted(taxo_int_counts.items(), key=operator.itemgetter(1),
								reverse=True)[0:15]
	high_pub_int_dict = sorted(pub_int_counts.items(), key=operator.itemgetter(1),
								reverse=True)[0:15]
								
	print("\nTaxonomic breakdown:\n")
	print("Taxons\t\t\t\t\tTotal Interactions")
	for taxon in high_taxo_dict:
		split_taxon = taxon[0].split("\t")
		taxidA = (split_taxon[0].split("|"))[0]
		taxidB = (split_taxon[1].split("|"))[0]
		print("%s\tvs.\t%s\t\t\t%s" % (taxidA, taxidB, taxon[1]))
	
	'''
	This count isn't useful unless it's done with OG vs. OG interactions.
	'''	
	#print("\nInteractions seen in the most taxons:")
	#print("Interaction\t\t\tCount of Different Taxons")
	#for count in high_taxo_int_dict:
	#	print("%s\t\t%s" % (count[0], count[1]))
	
	print("\nPublication breakdown:\n")
	print("Publication\t\tPMID\t\tTotal Interactions")
	for pub in high_pub_dict:
		print("%s\t\t%s" % (pub[0], pub[1]))
	
	print("\nInteractions (nonredundant) seen in the most publications:")
	print("Interaction\t\t\tCount of Different Publications")
	for count in high_pub_int_dict:
		print("%s\t\t%s" % (count[0], count[1]))
		
	#Save the matching set if we don't have one yet
	if have_output == 0:
		outfilename = "matching_ppi.txt"
		print("\nSaving filtered interactions to %s." % outfilename)
		save_interactions(clean_interactions, outfilename)
	
	'''
	OG-expanded interaction networks below.
	'''
	
	#Get top lists of them
	#Get counts first
	og_taxo_counts = {key: len(value) for key, value in og_taxo_dict.items()}
	og_pub_counts = {key: len(value) for key, value in og_pub_dict.items()}
	og_taxo_int_counts = {key: len(value) for key, value in og_taxo_ints.items()}
	og_pub_int_counts = {key: len(value) for key, value in og_pub_ints.items()}

	og_high_taxo_dict = sorted(og_taxo_counts.items(), key=operator.itemgetter(1),
								reverse=True)[0:15]
	og_high_pub_dict = sorted(og_pub_counts.items(), key=operator.itemgetter(1),
								reverse=True)[0:15]
	og_high_taxo_int_dict = sorted(og_taxo_int_counts.items(), key=operator.itemgetter(1),
								reverse=True)[0:15]
	og_high_pub_int_dict = sorted(og_pub_int_counts.items(), key=operator.itemgetter(1),
								reverse=True)[0:15]
								
	print("\nTaxonomic breakdown:\n")
	print("Taxons\t\t\t\t\tTotal Interactions")
	for taxon in og_high_taxo_dict:
		split_taxon = taxon[0].split("\t")
		taxidA = (split_taxon[0].split("|"))[0]
		taxidB = (split_taxon[1].split("|"))[0]
		print("%s\tvs.\t%s\t\t\t%s" % (taxidA, taxidB, taxon[1]))
	
	print("\nInteractions seen in the most taxons:")
	print("Interaction\t\t\tCount of Different Taxons")
	for count in og_high_taxo_int_dict:
		print("%s\t\t%s" % (count[0], count[1]))
	
	print("\nPublication breakdown:\n")
	print("Publication\t\tPMID\t\tTotal Interactions")
	for pub in og_high_pub_dict:
		print("%s\t\t%s" % (pub[0], pub[1]))
	
	print("\nInteractions (nonredundant) seen in the most publications:")
	print("Interaction\t\t\tCount of Different Publications")
	for count in og_high_pub_int_dict:
		print("%s\t\t%s" % (count[0], count[1]))
		
	#Build and visualize ALL graphs now
	print("\n\nBuilding and visualizing graphs.")
	graph_interactions(prot_ids, unique, og_unique, prot_dict, all_og_map, 
						pub_int_counts, og_pub_int_counts, 
						con_limit, og_con_limit)
	
	print("All done.")

if __name__ == "__main__":
	sys.exit(main())
