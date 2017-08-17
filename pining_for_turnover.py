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
	
def get_wide_data():
	#Returns the wide data file produced by pining_for_new_data
	
	os.chdir(directories[0])
	
	wide_file_list = glob.glob('*_wide.txt')
	
	if len(wide_file_list) > 1:
		sys.exit("More than one wide format data file found. "
					"Only one is necessary.")
	if len(wide_file_list) == 0:
		sys.exit("Could not find a wide format data file.")
	
	wide_file_name = wide_file_list[0]
	
	return wide_file_name
	
def is_service_running(name):
	#Checks if a linux service is running.
	#See https://stackoverflow.com/questions/17541044/how-can-i-make-the-python-program-to-check-linux-services
    with open(os.devnull, 'wb') as hide_output:
        exit_code = subprocess.Popen(['service', name, 'status'], \
			stdout=hide_output, stderr=hide_output).wait()
        return exit_code == 0

#Main
def main():
	print("** pining - turnover module **")
	
	print("Searching for wide format data file.")
	wide_file_name = get_wide_data()
	print("Found file: %s" % wide_file_name)
	
	print("Locating interaction database...")
	access_graphdb()
	
if __name__ == "__main__":
	sys.exit(main())
