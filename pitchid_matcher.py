
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
	

FB_CH=id_name_match("input_filename")

sys.stdout = open("output_filename", 'w')





