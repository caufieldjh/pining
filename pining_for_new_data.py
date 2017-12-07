#!/usr/bin/python
#pining_for_new_data.py
'''
A tool for retrieving protein-protein interaction data sets, using
a list of proteins to search for.
As the number of returned interactions may be quite large,
IntAct and BioGRID are searched as local files.

Produces combined, annotated versions of the network of the 
retrieved interactions in several different formats.

Saves interactions and corresponding orthologous groups to a graph
database (Neo4j). Tries to start the Neo4j server if it isn't already
running. Deletes all existing records in the database at each run.

Also visualizes graph in Cytoscape, though only if it is running.


'''
__author__= "Harry Caufield"
__email__ = "j.harry.caufield@gmail.com"

import argparse
import glob, gzip, operator, os, random, re, sys, time, urllib2, zipfile
from datetime import date

import requests #For using Cytoscape's cyREST API
import json
import subprocess

import pandas as pd
from pandas import pivot_table

import py2neo
from py2neo import authenticate, Graph, Node, Relationship

from tqdm import *

#Constants and Options
directories = ["output","databases", "og_info"]

upid_match = re.compile('[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}')

nowstring = (date.today()).isoformat()

#Classes

#Functions

def build_bgid_conv_file(filename):
	#Given the name of a BioGRID ID file, produces a new file
	#specific to converting UniprotAC IDs to BioGRID IDs.
	
	outfilename = 'BIOGRID-conversion.txt'
	
	i = 0
	with open(filename) as bgidfile:
		filelen = sum(1 for line in bgidfile) -1
		bgidfile.seek(0)
		pbar = tqdm(total=filelen)
		with open(outfilename, "w+b") as outfile:
			outfile.write("BIOGRID_ID\tUNIPROT_ACCESSION\n")
			for line in bgidfile:
				pbar.update(1)
				#Skip the header first
				if line[0:10] != "BIOGRID_ID":
					pass
				else:
					break
			for line in bgidfile:
				pbar.update(1)
				splitline = line.split("\t")
				if splitline[2] in ["UNIPROT-ACCESSION", "SWISS-PROT", "TREMBL"]:
					outline = "%s\t%s\n" % (splitline[0], splitline[1])
					outfile.write(outline)
	
	pbar.close()
	return outfilename
	
def build_bgid_conv_map(filename):
	#Given the name of a file produced by build_bgid_conv_file(),
	#return a dict with UniprotAC IDs as values and BioGRID IDs as keys.
	#Note that this is not 1:1 as one BioGRID ID may include multiple 
	#(usually redundant) Uniprot IDs.
	#These are retaining in the dictionary.
	
	conv_map = {}
	
	with open(filename) as bgidconvfile:
		bgidconvfile.readline() #Skip header
		for line in bgidconvfile:
			splitline = (line.rstrip()).split("\t")
			bgid = splitline[0]
			upid = splitline[1]
			if bgid not in conv_map:
				conv_map[bgid] = [upid]
			else:
				if upid not in conv_map[bgid]:
					conv_map[bgid].append(upid)
	
	return conv_map
			
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
	Or, if we have them already, locate them
	'''
	
	if name == "IntAct":
		baseURL = "ftp://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/"
		dbfilename = "intact.zip"
		
	if name == "BioGRID":
		baseURL = "https://thebiogrid.org/downloads/archives/Latest%20Release/"
		dbfilename = "BIOGRID-ALL-LATEST.mitab.zip"
		bgidfilename = "BIOGRID-IDENTIFIERS-LATEST.tab.zip" #BioGRID internal ID conversions

	dbfilepath = baseURL + dbfilename
	
	#Database names may change
	if name == "IntAct":
		outfilename = "intact.txt"
	elif name == "BioGRID":
		biogrid_db_list = glob.glob('BIOGRID-ALL-*.mitab.txt')
		if len(biogrid_db_list) > 0:
			outfilename = biogrid_db_list[0]
		else: 
			outfilename = "BIOGRID-ALL-LATEST.mitab.txt" 
			#This won't actually be the filename
		
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
		
	if name == "BioGRID": #May need the BioGRID ID conversion file
		biogrid_idfile_list = glob.glob('BIOGRID-IDENTIFIERS-*.tab.txt')
		if len(biogrid_idfile_list) > 0:
			print("Found BioGRID ID conversion file: %s" % biogrid_idfile_list[0])
			dl_bgidfile = False
		else:
			dl_bgidfile = True
		
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
						
	if name == "BioGRID" and dl_bgidfile:
		print("Downloading and decompressing BioGRID ID conversion file.")
		dbfilepath = baseURL + bgidfilename
		print("Downloading from %s" % dbfilepath)
		response = urllib2.urlopen(dbfilepath)
		
		compressed_file = open(os.path.basename(dbfilepath), "w+b")
		#Start local compressed file
		chunk = 1048576
		while 1:
			data = (response.read(chunk)) #Read one Mb at a time
			compressed_file.write(data)
			if not data:
				print("\n%s file download complete." % bgidfilename)
				compressed_file.close()
				break
			sys.stdout.flush()
			sys.stdout.write(".")
		
		try:
			with zipfile.ZipFile(bgidfilename, "r") as infile:
				for filename in infile.namelist():
					print(filename)
					outfile = open(filename, 'w+b')
					outfile.write(infile.read(filename))
					outfile.close()
					outfilename = filename
					break #Just want one file.
		except zipfile.BadZipfile as e:
			sys.exit("Something is wrong with the BioGRID ID file.\n"
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
			
	os.chdir("..")		
	
	return target_OG_maps, prot_OG_maps, prots_without_OG
	
def search_int_file(ids, filenames, db, target_ogs, all_og_map, bgidconvmap):
	'''
	Searches a PSI-MI TAB format set or sets of interactions (filenames) 
	for the Uniprot IDs provided (ids) and for and IDs with the same
	OGs as the provided IDs, using target_ogs and the og map.
	Indexes file first, searches IDs in the index, then returns
	interactions at the specified indices.
	Also returns interactions in short form, as proteins, corresponding
	OGs, and counts/IDs of publications each interaction was observed
	in.
	
	Have multiple files to deal with - 
	load IntAct, then BioGRID, as the latter uses more general IDs
	for interactors.
	Note that prot_match_count and og_match_count may include 
	double counts if the same interaction is present
	in both BioGRID and IntAct (this will be fixed soon).
	'''
	
	ids = set(ids) #ids should be unique, plus sets are more efficient
	
	all_int = {} #Just the interactors from each interaction.
					#Keys are line numbers
					#Values are tuples of interactors
					#Specific to each database.
					
	all_og_int = {} #Same as all_int but values are tuples of OGs
					#the proteins belong to.
					#Specific to each database.
	
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
							#[protA, protB, OG_A, OG_B, 
							#	pub_count, pub_ids,
							#	pred_type_A, pred_type_B]
							
	if db == "previous":
		os.chdir(directories[0])
	else:
		os.chdir(directories[1])
		
	#Load IntAct interactions or previous results
	with open(filenames[0]) as intactfile:
		#Count the lines in the file first to determine how many PPI
		filelen = sum(1 for line in intactfile) -1
		print("Searching %s interactions from IntAct." % filelen)
		
		#Progbar
		pbar = tqdm(total=filelen)
		
		all_int = {} 
		all_og_int = {} 
		
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
			pbar.update(1)
			
		pbar.close()
		
		#Now search interactor pairs for target interactors
		get_these_lines = set()
		
		#Search protein IDs first
		for line_num in all_int:
			for interactor in all_int[line_num]:
				if interactor in ids: #Option to search for match in both goes here
					get_these_lines.add(line_num)
					prot_match_count = prot_match_count +1
					break
		
		#Now search OG IDs
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
	
	if db == "full":
		
		all_int = {}
		all_og_int = {}
		
		with open(filenames[1]) as biogridfile:
			#Count the lines in the file first to determine how many PPI
			filelen = sum(1 for line in biogridfile) -1
			print("Searching %s interactions from BioGRID." % filelen)
			
			#Progbar
			pbar = tqdm(total=filelen)
			
			biogridfile.seek(0)
			biogridfile.readline() #Skip header
			
			i = 1
			#bgidconvmap
			#Load all interactor pairs as original IDs and as OGs
			for line in biogridfile:
				splitline = (line.rstrip()).split("\t")
				#Convert BioGRID ID to Uniprot
				if (splitline[2].split(":"))[0] == "biogrid" and \
					(splitline[3].split(":"))[0] == "biogrid":
					bgidA = (((splitline[2].split(":"))[1]).split("|"))[0]
					bgidB = (((splitline[3].split(":"))[1]).split("|"))[0]
					#Just using one of the corresponding UPIDs for simplicity
					#Not all BGIDs map to UPIDs.
					if bgidA in bgidconvmap:
						interactorA = bgidconvmap[bgidA][0]
					else:
						interactorA = "0"
						
					if bgidB in bgidconvmap:
						interactorB = bgidconvmap[bgidB][0]
					else:
						interactorB = "0"
						
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
				pbar.update(1)
				
			pbar.close()
			
			#Now search interactor pairs for target interactors
			get_these_lines = set()
			
			#Search protein IDs first
			for line_num in all_int:
				for interactor in all_int[line_num]:
					if interactor in ids: #Option to search for match in both goes here
						get_these_lines.add(line_num)
						prot_match_count = prot_match_count +1
						break
			
			#Now search OG IDs
			for line_num in all_og_int:
				for interactor in all_og_int[line_num]:
					if interactor in target_ogs:
						get_these_lines.add(line_num)
						og_match_count = og_match_count +1
						break
			
			biogridfile.seek(0)
			biogridfile.readline() #Skip header 	
			
			j = 1
			
			#Add all interactions from the file
			for line in biogridfile:
				if j in get_these_lines:
					outline = (line.rstrip()).split("\t")
					#Ensure that the first two items are uniprotkb:
					outline[0] = "uniprotkb:%s" % all_int[j][0]
					outline[1] = "uniprotkb:%s" % all_int[j][1]
					complete_interactions.append(outline)
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
		
		if prot_A in ids:
			pred_type_A = "Target"
		else:
			if og_A in target_ogs:
				pred_type_A = "Shared_OG"
			else:
				pred_type_A = "Other"
			
		if prot_B in ids:
			pred_type_B = "Target"
		else:
			if og_B in target_ogs:
				pred_type_B = "Shared_OG"
			else:
				pred_type_B = "Other"
		
		values = [prot_A, prot_B, og_A, og_B, pub_count, pub_ids, 
					pred_type_A, pred_type_B]
		
		short_interactions.append(values)
	
	#Recalculate total unique protein and OG interactors matched
	#based on the short_interactions content
	matched_prots = set()
	matched_ogs = set()
	for interaction in short_interactions:
		for item in interaction[:2]:
			matched_prots.add(item)
		for item in interaction[2:4]:
			matched_ogs.add(item)
	
	prot_match_count = len(matched_prots)
	og_match_count = len(matched_ogs)
						 
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
		header = "prot_A\tprot_B\tOG_A\tOG_B\tPubCount\tPublications\tType_A\tType_B"
	
	os.chdir(directories[0])
	
	with open(filename, 'wb') as outfile:
		#Write the header
		outfile.write("%s\n" % header)
		for interaction in ppi:
			flatline = "\t".join(interaction)
			outfile.write("%s\n" % flatline)
		
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
	
def show_graph(interactions):
	#Takes interactions in short format and visualizes with Cytoscape
	#Returns False if it fails, specifically if it can't produce the
	#graph
	port = '1234'
	ip = 'localhost'
	base = 'http://' + ip + ':' + port + '/v1/'
	headers = {'Content-Type': 'application/json'}
	
	pedges = ""
	oedges = ""
	
	try:
		res = requests.get(base)
		requests.delete(base + 'networks')
		
		#Assemble the input
		for ppi in interactions:
			new_pedge = "%s\t%s" % (ppi[0], ppi[1])
			pedges = pedges + new_pedge + "\n"
			new_oedge = "%s\t%s" % (ppi[2], ppi[3])
			oedges = oedges + new_oedge + "\n"
		
		print("Sending protein network to Cytoscape.")
		
		res = requests.post(base + 'networks?format=edgelist&collection=pining&title=Proteins', 
							data=pedges, headers=headers)
							
		p_res_dict = res.json()
		p_suid = p_res_dict['networkSUID']
		requests.get(base + 'apply/layouts/circular/' + str(p_suid))
		
		print("Sending OG network to Cytoscape.")
		
		res = requests.post(base + 'networks?format=edgelist&collection=pining&title=OGs', 
							data=oedges, headers=headers)
		
		o_res_dict = res.json()
		o_suid = o_res_dict['networkSUID']
		requests.get(base + 'apply/layouts/circular/' + str(o_suid))
		
	except requests.exceptions.ConnectionError:
		print("**Tried to output networks to Cytoscape but it " 
				"may not be running.")
		return False

def is_service_running(name):
	#Checks if a linux service is running.
	#See https://stackoverflow.com/questions/17541044/how-can-i-make-the-python-program-to-check-linux-services
    with open(os.devnull, 'wb') as hide_output:
        exit_code = subprocess.Popen(['service', name, 'status'], \
			stdout=hide_output, stderr=hide_output).wait()
        return exit_code == 0

def create_graphdb(interactions):
	#Takes interactions in short format and produces a graph database
	#using Neo4j through py2neo.
	#This database includes both OG and protein information, such that
	#"interacts with" relationships are between proteins,
	#but "is a member of" relationships are between proteins and OGs.
	
	graph_build_complete = False
	
	#Assemble the input
	prot_list = []
	og_list = []
	for ppi in interactions:
		for name in (ppi[0], ppi[1]):
			prot_list.append(name)
		for name in (ppi[2], ppi[3]):
			og_list.append(name)
	
	#Unique names only
	prot_list = list(set(prot_list))
	og_list = list(set(og_list))
	node_dict = {}
	
	#Start the service if needed
	if not is_service_running('neo4j'):
		print("Starting neo4j service....")
		os.system("sudo service neo4j start")
		time.sleep(5) #Take a few moments to let service start
	authenticate("localhost:7474", "neo4j", "pining")
	
	while not graph_build_complete:
		try:
			g = Graph("http://localhost:7474/db/data/")
			g.delete_all()
			tx = g.begin()
	
			#Assemble the nodes in the graph
			print("Assembling nodes.")
			for name in prot_list:
				this_node = Node("protein", name=name)
				node_dict[name] = this_node
				tx.create(this_node)
			for name in og_list:
				this_node = Node("OG", name=name)
				node_dict[name] = this_node
				tx.create(this_node)
	
			#Assemble relationships in the graph
			#Assumes each interaction involves one PPI and 2 OGs, 
			#for 3 total realationships at most
			print("Assembling relationships.")
			pbar = tqdm(total=len(interactions)*3) 
			for ppi in interactions:
				tx.create(Relationship(node_dict[ppi[0]], "interacts_with", node_dict[ppi[1]]))
				pbar.update(1)
				tx.create(Relationship(node_dict[ppi[0]], "member_of", node_dict[ppi[2]]))
				pbar.update(1)
				tx.create(Relationship(node_dict[ppi[1]], "member_of", node_dict[ppi[3]]))
				pbar.update(1)
			tx.commit()
			
			pbar.close()
			print("Access graph at http://localhost:7474")
			graph_build_complete = True
			
		except py2neo.packages.httpstream.http.SocketError as e:
			print("**Error accessing the Neo4j database: %s" % e)
			print("**Please try accessing the server at http://localhost:7474/")
		except (py2neo.database.status.Unauthorized,
				py2neo.database.status.GraphError) as e:
			#This may be a new database - the user needs to handle
			#Really wish there was a better way to handle this
			print("**If this is a new database, you may need to set a password.")
			print("**Please try accessing the server at http://localhost:7474/")
			print("**The default username is \"neo4j\" and the password is \"neo4j\".")
			print("**Please change the password to \"pining\".\n")
			raw_input("When done, press any key to continue.")
			#I know, not Python 3 compatible, thank you
	
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
	
	#Check for BioGRID original ID file or new conversion file
	bgidconvfilename = 'BIOGRID-conversion.txt'
	if os.path.isfile(bgidconvfilename):
		print("Found BioGRID ID conversion map.")
		
	else:
		bgidfile_list = glob.glob('BIOGRID-IDENTIFIERS-*.tab.txt')
		if len(bgidfile_list) > 1:
			sys.exit("Found more than one BioGRID ID file.\n"
						"Check for duplicates.")
						
		if len(bgidfile_list) == 0:
			sys.exit("No BioGRID ID file found. Please re-download database "
					"to ensure IDs match database version.")
	
		if len(bgidfile_list) == 1:
			print("Found BioGRID ID file. Using it to build ID conversion map...")
			bgidfile = bgidfile_list[0]
			bgidconvfilename = build_bgid_conv_file(bgidfile)
	
	biogridfile_loc = os.path.join(directories[1], biogridfilename)
	bgidconvfile_loc = os.path.join(directories[1], bgidconvfilename)
		
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
	
	#Set up the BioGRID Uniprot/internal ID conversion map
	print("Mapping Uniprot IDs to BioGRID IDs.")
	bgidconvmap = build_bgid_conv_map(bgidconvfile_loc)
	bgid_conv_count = len(bgidconvmap)
	all_upids = [item for idlist in bgidconvmap.values() for item in idlist]
	upid_conv_count = len(set(all_upids))
	print("Conversion map covers %s Uniprot protein IDs "
			"and %s unique BioGRID IDs." % (upid_conv_count, bgid_conv_count))
	
	'''
	Map proteins to OGs.
	og_map includes only target protein IDs as keys
	all_og_map includes ALL protein IDs.
	We use this second map to get OGs for interactors outside the target set.
	'''
	
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
		print("Searching IntAct and BioGRID interactions for provided protein IDs.")
		db = "full"
		complete_interactions, prot_match_count, og_match_count, \
			short_interactions \
			= search_int_file(prot_ids, [intactfilename, biogridfilename], 
								db, these_ogs, all_og_map, bgidconvmap)
	else:
		have_output = True
		print("Found existing interaction output file: %s" % ppioutfilename)
		print("Analyzing previous interaction output.")
		db = "previous"
		complete_interactions, prot_match_count, og_match_count, \
			short_interactions \
			= search_int_file(prot_ids, [ppioutfilename], db, these_ogs, 
								all_og_map, bgidconvmap)
	
	ppi_count = len(complete_interactions)
	short_count = len(short_interactions)
	
	print("\nFound %s total interactions involving target protein IDs "
			"or shared OGs." % ppi_count)
	print("%s unique proteins (including orthologs) are involved in "
			"these interactions." % prot_match_count)
	print("%s unique OGs are involved in these interactions." % og_match_count)
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
	
	print("See output files:\n")
	for output_type in outputs:
		print("%s: %s\n" % (output_type, outputs[output_type]))
	
	print("In the short interaction file, interactors may be \n"
			"Target if included in the input file, \n"
			"Shared_OG if in the same OG as a protein in the input,\n"
			"or Other if not in a shared OG.\n")
	 
	print("Saving graph database.")
	create_graphdb(short_interactions)
	 
	show_graph(short_interactions)

	print("Complete.")

if __name__ == "__main__":
	sys.exit(main())
