
import numpy as np
import matplotlib.pyplot as plt
import scipy
import scipy.optimize
import matplotlib.pyplot as plt
import matplotlib.lines as lines
import itertools
import csv
import sys
import math

#upload master file of id's and names
master_file = open("master_pitchid.csv")
master_file.readline() # read past the header line


ids_unformatted = csv.reader(master_file, delimiter=',', quotechar='"') 

ids = []

#file is in format with first and last names in different columns, 
#this combines names and pairs them with corresponding id number, [name,id] 
for line in ids_unformatted:
	names_and_ids = []
	#adds full names to list (first + last)
	names_and_ids.append(line[1]+" "+line[0])

	# adds the ids to the full names list, for each name
	try: # do this within a try/except
		names_and_ids.append(float(line[4].strip()))
		ids.append(names_and_ids)
	except ValueError: # if there was a value error, the line was bad
		continue


id_numbers=[i[1] for i in ids] #list of just names from master list
id_names=[i[0] for i in ids] #list of just id #'s from master list, in same order as name list

bad_names=['Zachary Neal','Carl Edwards Jr','J.C. Ramirez','Joshua Osich','Nate Karns','Jorge De La Rosa','Rubby De La Rosa','Robbie Ross','Vincent Velasquez'] 
replacements=['Zach Neal','Carl Edwards Jr.','JC Ramirez','Josh Osich','Nathan Karns','Jorge de la Rosa','Rubby de la Rosa','Robbie Ross Jr.','Vince Velasquez']

def id_name_match(filename):
	
	deception_file = open(filename) #open deception score file
	deception_file.readline() # read past the header line

	deception_list = csv.reader(deception_file, delimiter=',', quotechar='"')
	
	#create list of just id #'s from deception score file
	d_ids = []
	for line in deception_list:
		d_ids+=[float(line[0])]	

	d_names=[] #initializing deception score file name list
	
	#if id in deception file is in the master id list, add the corresponding name to the name list
	#this keeps the order of the name list inline with the order of the deception score file id list
	for i in d_ids:
		d_names+=[id_names[id_numbers.index(i)]]

	for i in range(len(bad_names)):
		try: 
			d_names[d_names.index(bad_names[i])]=replacements[i]
		except:
			continue

	return d_names
	


#list of deception score files and created name list file for each year
#and function that prints out exported file

#this can be used to export a name file for each deception score file
#if if deception score file is in the same folder as the program


#files=[["16_nq_deception_short.csv",'16_id_names_short.csv'],["15_nq_deception_short.csv",'15_id_names_short.csv'],["14_nq_deception_short.csv",'14_id_names_short.csv'],["13_nq_deception_short.csv",'13_id_names_short.csv'],["12_nq_deception_short.csv",'12_id_names_short.csv'],["11_nq_deception_short.csv",'11_id_names_short.csv'],["10_nq_deception_short.csv",'10_id_names_short.csv'],["09_nq_deception_short.csv",'09_id_names_short.csv'],["08_nq_deception_short.csv",'08_id_names_short.csv']]

#for i in files:
#	sys.stdout = open(i[1], 'w')
#	for n in id_name_match(i[0]):
#		print n





