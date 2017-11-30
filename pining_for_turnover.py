#!/usr/bin/python
#pining_for_turnover.py
'''
A tool for working with protein-protein interaction networks and
comparing subgraphs on the basis of averaged protein abundance levels
and the temporal differences of those abundances (i.e., protein
turnover).

This tool works with Neo4j databases and output files produced by or in 
the same format of those produced by pining_for_new_data.py.
It tries to start the Neo4j server if it is not already running.

'''
__author__= "Harry Caufield"
__email__ = "j.harry.caufield@gmail.com"

import glob, os, subprocess, sys, time

import pandas as pd
from pandas import pivot_table

import py2neo
from py2neo import authenticate, Graph, Node, Relationship

#Constants and Options
directories = ["output","databases", "og_info"]

#Classes

#Functions

def access_graphdb():
	#Verifies that an existing Neo4j database exists and that it is
	#populated. Does not modify the database.
	
	if not is_service_running('neo4j'):
		print("Starting neo4j service....")
		os.system("sudo service neo4j start")
		time.sleep(5)
	authenticate("localhost:7474", "neo4j", "pining")
	
	try:
		#Just access graph and retrieve stats
		g = Graph("http://localhost:7474/db/data/")
		interactions = (g.data("MATCH p=()-[r:interacts_with]->() RETURN p"))
		member_ofs = (g.data("MATCH p=()-[r:member_of]->() RETURN p"))
		if not interactions or not member_ofs:
			sys.exit("Graph is empty or is missing expected relationships. "
						"Exiting.")
		else:
			print("Graph contains %s protein interactions and %s OG memberships."
					% (len(interactions), len(member_ofs)))
	
	except (py2neo.packages.httpstream.http.SocketError,
			py2neo.database.status.Unauthorized) as e:
		print("**Error accessing the Neo4j database: %s" % e)
		print("**Please try accessing the server at http://localhost:7474/")
		sys.exit()
		
def add_values_to_graphdb(data):
	#Adds values, provided in pandas dataframe, to Neo4j database.
	#Assumes input dataframe is indexed using OG IDs.
	#Requires unique OG IDs, so this must be resolved before passing
	#input to this method.
	#Just tests for now.
	
	#Melt the data frame so each observation is one row.
	flatdata = pd.melt(data.reset_index(), id_vars=['id'])
	
	#We don't know how many conditions there are -
	#only that we have an id, a value, and n conditions in between.
	flatdata["condition"] = flatdata.apply(lambda row: '|'.join(map(str, row.iloc[1:-1])), axis=1)
	flatdata = flatdata.drop(flatdata.columns[1:-2], axis=1)
	
	#Now add it to the graph database
	try:
		g = Graph("http://localhost:7474/db/data/")
		for index, row in flatdata.iterrows():
			condition_name = row['condition']
			cypher_string = "MATCH (n:OG {name:$id}) MERGE ({`%s`: $value})" % condition_name
			g.run(cypher_string, parameters = {'id': row['id'],
									'value': row['value']})
		#for index in set(data.index):
		#	node_id = index
		#	print(data.loc[index])
	
	except (py2neo.packages.httpstream.http.SocketError,
			py2neo.database.status.Unauthorized) as e:
		print("**Error accessing the Neo4j database: %s" % e)
		print("**Please try accessing the server at http://localhost:7474/")
		sys.exit()
	
def get_og_dict():
	#Uses Uniprot to OG maps to build dict.
	#Assumes the new data module has already created a map file.
	#Returns error otherwise.
	
	og_dict = {}
	
	os.chdir(directories[2])
	
	map_file_list = glob.glob('uniprot_og_maps_*')
	
	if len(map_file_list) > 1:
		sys.exit("More than one OG map file found. "
					"Only one is necessary.")
	if len(map_file_list) == 0:
		sys.exit("Could not find an OG map file.\n"
					"You may need to run the new data module first.")
	
	with open(map_file_list[0]) as map_file:
		for line in map_file:
			splitline = (line.rstrip()).split("\t")
			og_dict[splitline[0]] = splitline[1]
	
	size = len(og_dict)
	
	print("Loaded %s protein IDs and their corresponding OGs." % size)
	
	os.chdir("..")
	
	return og_dict

def get_wide_data():
	#Returns contents of the wide data file produced by 
	#pining_for_new_data
	
	os.chdir(directories[0])
	
	wide_file_list = glob.glob('*_wide.txt')
	
	if len(wide_file_list) > 1:
		sys.exit("More than one wide format data file found. "
					"Only one is necessary.")
	if len(wide_file_list) == 0:
		sys.exit("Could not find a wide format data file.\n"
					"You may need to run the new data module first.")
					
	wide_file_name = wide_file_list[0]
	print("Loading data from %s." % wide_file_name)
	
	data = pd.read_csv(wide_file_name, sep='\t', skipinitialspace=True, 
						header=[0,1,2], index_col=0)
	#print(data)
	shape = data.shape
	print("Data contains %s protein IDs and %s sets of values."
			% (shape[0], shape[1]))
	
	os.chdir("..")
	
	return data
	
def is_service_running(name):
	#Checks if a linux service is running.
	#See https://stackoverflow.com/questions/17541044/how-can-i-make-the-python-program-to-check-linux-services
    with open(os.devnull, 'wb') as hide_output:
        exit_code = subprocess.Popen(['service', name, 'status'], \
			stdout=hide_output, stderr=hide_output).wait()
        return exit_code == 0
        
def map_data_to_OGs(data_frame, og_dict):
	#Given a data frame of values and a protein-to-OG dictionary,
	#returns the same data frame with corresponding OG IDs
	#in the first column.
	#This assumes the OG inherits the values of the protein -
	#this may not always be true but allows values to extend
	#across species.
	#Note that this will likely produce duplicate index values!
	#This will be fixed soon. Trust me, please.
	
	og_data = data_frame
	
	og_data.index = og_data.index.to_series().map(og_dict)
	
	print(og_data)
	
	return og_data
		
#Main
def main():
	print("** pining - turnover module **")
	
	print("Searching for wide format data file.")
	data_to_map = get_wide_data()
	
	print("Retreiving OG IDs.")
	og_dict = get_og_dict()
	
	print("Mapping data to OGs.")
	og_data = map_data_to_OGs(data_to_map, og_dict)
	print("Mapped.")
	
	print("Locating interaction database...")
	access_graphdb()
	
	print("Mapping data to interaction database...")
	add_values_to_graphdb(og_data)
	
	print("Complete.")
	
if __name__ == "__main__":
	sys.exit(main())
