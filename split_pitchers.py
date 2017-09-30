import numpy as np 
import argparse
import sys
import os

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('input', help="Input file of processed and annotated pitch_table.", type=str)
	parser.add_argument('directory', help="Title of new directory where split pitchers will be stored.", type=str)
	parser.add_argument('year', help="Year of the file.", type=str)

	args = parser.parse_args()

	# create the new sub directory to store the created graphs
	path = os.path.dirname(os.path.abspath(__file__))
	path = path +'/' + args.directory
	if os.path.exists(path):
		print "A directory with the specified name already exists. Please delete this directory or choose a different name."
		sys.exit()
	os.makedirs(path)

	# open the input file
	infile = open(args.input)

	# store lines in a dictionary, with the id as a key
	pitcher_dict = {}
	for line in infile:
		split = line.split(',')

		# if this ID hasn't been seen
		if split[0] not in pitcher_dict:
			pitcher_dict[split[0]] = []
		pitcher_dict[split[0]].append(line)

	# create a new file for each pitcher ID that is seen
	for pid in pitcher_dict:
		outfile = open(args.directory + '/' + args.year + '_' + pid + '.csv', 'w')
		for line in pitcher_dict[pid]:
			outfile.write(line)
		outfile.close()


if __name__ == "__main__":
	main()