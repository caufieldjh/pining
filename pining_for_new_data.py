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
	observations.
	Also returns heading names of the observations if available,
	or just assigns generic ones if not provided, and
	determines if they are qualitative or not.
	This is crucial to how the observations will be treated.
	
	Saves a file containing the list of interactors alone,
	as "[originalfilename]_prot_only.txt".
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

	outfilename = "%s_prot_only.txt" % filename[0:-4]
			
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
			upid = splitline[0]
			if upid not in prot_dict.keys():
				prot_dict[upid] = [splitline[1:]]
			else:
				prot_dict[upid].append(splitline[1:])
				
	#Set up output file
	
	os.chdir(directories[0])
	
	#Now write the output protein file
	#With one uniprot ID per line
	with open(outfilename, 'w') as outfile:
		#Write the header

		for upid in prot_dict.keys():
			outfile.write("%s\n" % upid)
	
	os.chdir("..")
	
	return prot_dict, groups

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
	
	dl_dbfile = 1	#If 1, we need to download
	if os.path.isfile(dbfilename): 
		#Already have the compressed file, don't download
		print("Found compressed database file on disk: %s" % dbfilename)
		decompress_dbfile = 1
		dl_dbfile = 0
	if os.path.isfile(outfilepath): 
		#Already have the decompressed file, don't download
		print("Found database file on disk: %s" % outfilepath)
		decompress_dbfile = 0
		dl_dbfile = 0
		
	if dl_dbfile == 1:
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
		decompress_dbfile = 1
		
	if decompress_dbfile == 1:
		print("Decompressing %s file." % name)
		with zipfile.ZipFile(dbfilename, "r") as infile:
			for filename in infile.namelist():
				print(filename)
				outfile = open(filename, 'w+b')
				outfile.write(infile.read(filename))
				outfile.close()
		
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
		
def map_prots_to_ogs(ids):
	#Uses EggNOG files to create a dict of protein IDs and 
	#corresponding OGs.
	#Takes list of Uniprot IDs as input.
	
	prot_OG_maps = {} #Dictionary of all Uniprot ID to OG maps
	target_OG_maps = {} 
	#Dictionary to save protein-OG mapping specific for this interaction set
	
	mapped_count = 0
	prots_without_OG = 0
	
	mapfile_list = glob.glob("uniprot_og_maps_*")
	if len(mapfile_list) == 0:
		sys.exit("Can't find OG mapping file. May need to start over.")
	else:
		mapfilename = mapfile_list[0]
	
	with open(mapfilename) as mapfile:
		for line in mapfile:
			splitline = (line.rstrip()).split()
			upid = splitline[0]
			og = splitline[1]
			prot_OG_maps[upid] = og
			
	for protein in ids:
		if protein in prot_OG_maps:
			matching_og = prot_OG_maps[protein]
			mapped_count = mapped_count +1
		else:
			matching_og = protein	
			#If the protein doesn't map to an OG it retains its upid
			prots_without_OG = prots_without_OG +1
				
		target_OG_maps[protein] = matching_og
	
	og_count = len(set(target_OG_maps.values()))
	
	print("Full OG map covers %s Uniprot protein IDs." % len(prot_OG_maps))		
	print("Mapped %s proteins to OGs. %s proteins did not map to OGs."
			% (mapped_count, prots_without_OG))
	print("Map includes %s OGs in total." % og_count)
			
	return target_OG_maps, prot_OG_maps
	
def search_int_file(ids, filename, all_og_map):
	#Searches a PSI-MI TAB format set of interactions (filename) 
	#for the Uniprot IDs provided (ids) and for and IDs with the same
	#OGs as the provided IDs, using the OG map.
	
	all_int = [] #All interactions in the IntAct set
				 #Treated as a set once all interactions added
	
	found_ppi = [] #Interactions containing searched interactors as UPIDs
	found_og_ppi = [] #Interactions containing searched interactors
						#as UPIDs, as well as interactions sharing
						#the same OGs. Still provided as UPIDs here.
	
	ids = set(ids) #ids should be unique, plus sets are more efficient
	
	ogs = []
	upids_not_mapped = [] 
	for upid in ids:
		try:
			if all_og_map[upid] not in ogs:
				ogs.append(all_og_map[upid])
		except KeyError:
			upids_not_mapped.append(upid)
			all_og_map[upid] = upid
	print("Protein IDs without corresponding OGs: %s"
			% ", ".join(upids_not_mapped))
	
	ogs = set(ogs)
	
	single_matches = 0 #Count of interactions involving 1 target protein
	double_matches = 0 #Count of interactions involving 2 target proteins
	single_og_matches = 0 #Interactions from shared OGs involving 1 OG member
	double_og_matches = 0 #Interactions from shared OGs involving 2 OG members
	
	with open(filename) as intactfile:
		#Count the lines in the file first to determine how many PPI
		#Then load into memory
		filelen = sum(1 for line in intactfile) -1
		print("Searching %s interactions." % filelen)
		prog_width = filelen / 10000
		
		#Progbar 1
		sys.stdout.write("[%s]" % (" " * prog_width))
		sys.stdout.flush()
		sys.stdout.write("\b" * (prog_width+1))
		
		intactfile.seek(0)
		intactfile.readline()
		
		for line in intactfile:
			splitline = (line.rstrip()).split("\t")
			all_int.append(tuple(splitline))
			
		all_int = set(all_int)
	
	j = 0
	
	for interaction in all_int:
		#Count all interactions containing target interactors
		j = j +1
		interactors = []
		for interactor in interaction[0:2]:
			try:
				interactors.append((interactor.split(":"))[1])
			except IndexError:
				pass
		interactors = tuple(interactors)
		
		#Lookup each interactor to find corresponding OG
		og_interactors = []
		for interactor in interactors:
			try:
				og_interactors.append(all_og_map[interactor])
			except KeyError:
				og_interactors.append(interactor)
				
		found_count = 0
		og_found_count = 0
		int_added = 0
		og_int_added = 0
		for prot_id in ids:
			#Check if one or both interactors is a target interactor
			if prot_id in interactors:
				found_count = found_count +1
			if found_count > 0 and int_added == 0:
				found_ppi.append(interaction)
				int_added = 1
				if found_count == 2:
					break
		for og in ogs:
			#Check if one or both interactors is a target interactor
			#but this time, with OGs
			if og in og_interactors:
				og_found_count = og_found_count +1
			if og_found_count > 0 and og_int_added == 0:
				found_og_ppi.append(interaction)
				og_int_added = 1
				if og_found_count == 2:
					break
		if og_found_count == 1:
			single_og_matches = single_og_matches +1
		elif found_count == 2:
			double_og_matches = double_og_matches +1
				
		if j % 10000 == 0:
			sys.stdout.flush()
			sys.stdout.write("#")
			
	return found_ppi, single_matches, double_matches, found_og_ppi, \
			single_og_matches, double_og_matches

def clean_ppi(interactions):
	#Removes interactions from the list if they are self interactions
	#only or involve non-protein interactors.
	#Also produces the set of unique interactions.
	
	clean_interactions = []
	self_interactions = []
	non_ppi = []
	
	unique = set() #All sets of interacting protein pairs
					#from those in the ppi
				#This is a set of tuples, with no repeats
				
	for interaction in interactions:
		interactors = []
		for interactor in interaction[0:2]:
			remove = 0
			if interactor == "-":
				self_interactions.append(interaction)
				remove = 1
				break
			elif interactor.split(":")[0] != "uniprotkb":
				non_ppi.append(interaction)
				remove = 1
				break
			interactors.append((interactor.split(":"))[1])
		interactors = tuple(interactors)
		if remove == 0:
			clean_interactions.append(interaction)
			if interactors not in unique:
				unique.add(interactors)
				
	return clean_interactions, self_interactions, non_ppi, unique
	
def describe_ppi(interactions):
	'''
	Gets dicts, of taxons and publications
	associated with interactions and interactors.
	Also gets the inverse: taxons and publications each 
	interaction is seen in - though in this case,
	interactions are just pairs of interactors in any arrangement
	(so A vs. B is the same as B vs. A).
	'''
	
	taxo_dict = {} #Taxonomic breakdown of interactions.
					#Keys are names and taxids.
					#Values are lists of interactions 
					#involving the source or pair of sources.
					#Cross-taxon interactions are counted separately.
	
	pub_dict = {} #Sources of interactions.
					#Keys are a tab-delimited string of author and PMID.
					#Values are lists of interactions
					#involving the source.
	
	taxo_ints = {} #Keys are tuples of interactors in an interaction
					#Values are names and taxids as used in taxo_dict.
					#Values are stored in a set as we only want uniques.
	
	pub_ints = {}  #Keys are tuples of interactors in an interaction
					#Values are author and PMID strings as in pub_dict.
					#Values are stored in a set as we only want uniques.
	
	for interaction in interactions:
		interactorA = (interaction[0].split(":"))[1]
		interactorB = (interaction[1].split(":"))[1]
		
		#Get taxonomic breakdown and taxon lists
		taxonA = (interaction[9].split(":"))[1]
		taxonB = (interaction[10].split(":"))[1]
		taxonAB = taxonA + "\t" + taxonB
		if taxonAB not in taxo_dict:
			taxo_dict[taxonAB] = [interaction]
		else:
			taxo_dict[taxonAB].append(interaction)
		
		#Add to taxon list for these interactors
		if (interactorA, interactorB) in taxo_ints or \
			(interactorB, interactorA) in taxo_ints:
			if (interactorA, interactorB) not in taxo_ints:
				taxo_ints[(interactorA, interactorB)] = set([taxonAB])
			else:
				taxo_ints[(interactorA, interactorB)].add(taxonAB)
		else:
			taxo_ints[(interactorA, interactorB)] = set([taxonAB])
		
		#Get publication breakdown
		author = interaction[7]
		pub_ids = interaction[8].split("|")
		for value in pub_ids:
			splitvalue = value.split(":")
			if splitvalue[0] == "pubmed":
				pmid = splitvalue[1]
		pub_key = author + "\t" + pmid
		if pub_key not in pub_dict:
			pub_dict[pub_key] = [interaction]
		else:
			pub_dict[pub_key].append(interaction)
	
		#Add to pub list for these interactors
		if (interactorA, interactorB) in pub_ints or \
			(interactorB, interactorA) in pub_ints:
			if (interactorA, interactorB) not in pub_ints:
				pub_ints[(interactorA, interactorB)] = set([pub_key])
			else:
				pub_ints[(interactorA, interactorB)].add(pub_key)
		else:
			pub_ints[(interactorA, interactorB)] = set([pub_key])
		
	return taxo_dict, pub_dict, taxo_ints, pub_ints
			
def save_interactions(ppi, outfilename):
	#Saves a set of interactions to a file
	tfile = outfilename
	
	os.chdir(outfiledir)
	with open(tfile, 'wb') as outfile:
		#Write the header
		outfile.write("%s\n" % psimitabheader)
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
			
#Main
def main():
	prot_dict = {} #Keys are unique protein IDs.
					#Values are lists of lists of additional data,
					#if present.
	
	print(sys.argv[0])
	
	parser = argparse.ArgumentParser()
	parser.add_argument('-inputfile', help="name of a text file containing "
						"UniprotAC IDs of proteins to search for, "
						"one protein ID per line.")
	args = parser.parse_args()
	
	psi_tab_header = load_psi_tab_header()
	
	#Set up the output and database storage directories
	for directory in directories:
		if not os.path.isdir(directory):
			print("Setting up directory for %s." % directory)
			os.mkdir(directory)
	
	#Load protein ID input file
	if args.inputfile:
		try:
			print("Loading %s as input file." % args.inputfile)
			protfilename = args.inputfile
			prot_dict, groups = load_prot(protfilename)
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
				
		except IOError as e:
			sys.exit("Can't find or open that file...\n%s" % e)
	else:
		sys.exit("No input file provided.")
	
	#Locate PPI database files
	#Download if necessary.
	os.chdir(directories[1])
	intact_db_list = glob.glob('intact*.txt')
	if len(intact_db_list) > 1:
		#Intact provides a list of negative interactions in the same file
		if 'intact_negative.txt' in intact_db_list:
			print("Found IntAct database file.")
		else:
			sys.exit("Found more than one copy of the IntAct database.\n"
						"Check for duplicates.")
	if len(intact_db_list) == 0:
		print("No IntAct database found. Downloading new copy.")
		get_ppi_db("IntAct")
	if len(intact_db_list) == 1:
		print("Found IntAct database file.")
		
	biogrid_db_list = glob.glob('BIOGRID-ALL-*.mitab.txt')
	if len(biogrid_db_list) > 1:
		sys.exit("Found more than one copy of the BioGRID database.\n"
					"Check for duplicates.")
	if len(biogrid_db_list) == 0:
		print("No BioGRID database found. Downloading new copy.")
		get_ppi_db("BioGRID")
	if len(biogrid_db_list) == 1:
		print("Found BioGRID database file.")
		
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
	
	sys.exit("Complete for now. Databases and OG map file are ready to use.")
	
	
	#Testing cutoff
	
	
	#Map proteins to OGs
	#og_map includes only target protein IDs as keys
	#all_og_map includes ALL protein IDs
	print("Mapping target proteins to orthologous groups.")
	og_map, all_og_map = map_prots_to_ogs(prot_ids)
	
	#Look for filtered interaction output if we already have it
	#If we don't have an output file, use the IntAct set.
	outfilename = "matching_ppi.txt"
	ppi_ofile_list = glob.glob(os.path.join(outfiledir, outfilename))
	if len(ppi_ofile_list) >1:
		sys.exit("Found multiple potential interaction input files. Exiting...")
	elif len(ppi_ofile_list) == 0 :
		have_output = 0
		print("Did not find previous output file.")
		print("Searching IntAct interactions for provided protein IDs.")
		match_ppi, single, double, og_match_ppi, og_single, og_double \
			= search_int_file(prot_ids, intactfile_loc, all_og_map)
	elif len(ppi_ofile_list) == 1:
		have_output = 1
		print("Found existing interaction output file: %s" % ppi_ofile_list[0])
		print("Analyzing previous interaction output.")
		match_ppi, single, double, og_match_ppi, og_single, og_double \
			= search_int_file(prot_ids, ppi_ofile_list[0], all_og_map)
	
	print("\nFound %s total interactions involving target protein IDs." % len(match_ppi))
	print("%s interactions involve one target protein." % single)
	print("%s interactions involve two target proteins." % double)
	
	#Filter self interactions and non-protein interactions
	#Also determine how many interactions are unique 
	#(remove reciprocals and compress across publications)
	print("\nFiltering out self and non-protein interactions...")
	clean_interactions, self_interactions, non_ppi, unique = clean_ppi(match_ppi)
	print("Found %s PPI." % len(clean_interactions))
	print("Found %s self interactions." % len(self_interactions))
	print("Found %s non-protein interactions." % len(non_ppi))
	print("%s interactions are unique (by protein ID) across all " 
			"publications in this list." % len(unique))
	
	#Describe the interactions
	taxo_dict, pub_dict, taxo_ints, pub_ints = describe_ppi(clean_interactions)
	
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
	
	print("\n\nWith all interactors sharing OG membership...")
	print("\nFound %s total interactions involving target protein IDs." % len(og_match_ppi))
	print("%s interactions involve one target protein." % og_single)
	print("%s interactions involve two target proteins." % og_double)
	
	#Filter self interactions and non-protein interactions
	#Also determine how many interactions are unique 
	#(remove reciprocals and compress across publications)
	print("\nFiltering out self and non-protein interactions...")
	og_clean_interactions, og_self_interactions, og_non_ppi, og_unique \
			= clean_ppi(og_match_ppi)
	print("Found %s PPI." % len(og_clean_interactions))
	print("Found %s self interactions." % len(og_self_interactions))
	print("Found %s non-protein interactions." % len(og_non_ppi))
	print("%s interactions are unique (by protein ID) across all " 
			"publications in this list." % len(og_unique))
	
	#Describe the interactions
	og_taxo_dict, og_pub_dict, og_taxo_ints, og_pub_ints = describe_ppi(og_clean_interactions)
	
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
		
	#Save the matching set if we don't have one yet
	if have_output == 0:
		outfilename = "matching_og_ppi.txt"
		print("\nSaving filtered interactions to %s." % outfilename)
		save_interactions(og_clean_interactions, outfilename)
	
	
	#Build and visualize ALL graphs now
	print("\n\nBuilding and visualizing graphs.")
	graph_interactions(prot_ids, unique, og_unique, prot_dict, all_og_map, 
						pub_int_counts, og_pub_int_counts, 
						con_limit, og_con_limit)
	
	print("All done.")

if __name__ == "__main__":
	sys.exit(main())
