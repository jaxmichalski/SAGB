import numpy as np 
import sys
import math
import os
import argparse

# dimensions
global d
d = 4

# pitch type labels
global pitch_labels
pitch_labels = ['AB', 'CH', 'CU', 'EP', 'FA', 'FC', 'FF', 'FO', 'FS', 'FT', 'IN', 'KC', 'KN', 'SC', 'SI', 'SL', 'UN']

# finds the deception for a single pitch
# deception is the P of misclassification
def pitchDeception(pitch, classes, pis, mus, inv_meanSigma, det_meanSigma, label):
	# calculate the leading term
	k = 1. / ( (((2 * np.pi)**0.5)**d) * det_meanSigma)
	
	# caculate the probability for each class
	probs = {}
	for c in classes:
		t1 = pitch - mus[c]
		exp_term = -0.5 * np.matmul(np.matmul(np.transpose(t1), inv_meanSigma), t1)
		new_prob = pis[c] * k * math.exp(exp_term[0][0])
		probs[c] = new_prob

	# return the covnerse of the label P divided by the sum of the p's
	return 1 - (probs[label] / sum(probs.values()))

# takes in a dictioanry of pitch types
def pitcherDeception(classData):
	# calcualte num_samples
	num_samples = 0
	for c in classData:
		num_samples += len(classData)

	# calculate mus, pi_cs, and sigmas for each class
	classes = [] # list of classes
	pis = {} # rates of each pitch
	mus = {} # mus for each class
	sigmas = {} # covariance for each class
	kludge = 5*(10**(-6))

	for c in classData.keys():
		# gets the data for this class
		locs = classData[c]

		classes.append(c)
		# pi
		pis[c] = (float(locs.shape[0]) / num_samples)
		# mu
		mus[c] = (np.transpose(np.array([np.mean(locs, axis=0)]))) # takes the mean along each col
		# finds covariance
		cov = np.cov(np.transpose(locs))
		sigmas[c] = cov + kludge*np.identity(4)

	# calculate the meanSigma, then its inverse
	meanSigma = np.zeros((4,4))
	for label in classes:
		meanSigma = np.add(meanSigma, pis[label]*sigmas[label])

	inv_meanSigma = np.linalg.inv(meanSigma)

	det_meanSigma = np.linalg.det(meanSigma)

	# array to return; contains triplets of label, count, and deception
	ret = {}
	# calculate the average deception (aka misclassification P) for each class
	for label in classes:
		# the data points for the class
		locs = classData[label]
		class_deception = 0
		for pitch in locs:
			class_deception += pitchDeception(np.transpose(np.array([pitch])), classes, pis, mus, inv_meanSigma, det_meanSigma, label)
		
		# divide by count
		class_deception = class_deception / len(locs)
		ret[label]= [len(locs), class_deception]

	# add the average deception
	total_count = 0
	total_deception = 0.
	for r in ret:
		total_count += ret[r][0]
		total_deception += ret[r][1] * ret[r][0]
	total_deception *= 1. / total_count
	ret['total'] = [total_count, total_deception]
	return ret

def create_line(year, pid, pd):
	# line starts with the ID and the year
	new_line = pid + ',' + year
	# for each label in the global list
	for label in pitch_labels:
		if label not in pd: # pitcher doesn't have this pitch type, add blank fields
			new_line += ',' + ','
		else:
			new_line += ',' + (',').join([str(i) for i in pd[label]])


	new_line += ',' + (',').join([str(i) for i in pd['total']])

	return new_line + '\n'

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('directory', help="Title of directory containing the split pitchers.", type=str)
	parser.add_argument('output', help="Name of output file.", type=str)
	parser.add_argument('year', help="Year of the data being analyzed", type=str)
	args = parser.parse_args()

	# open the output file and write the header line
	outfile = open(args.output, 'w')
	header = "PID,year"
	for label in pitch_labels:
		header += ',' + label + '_count,' + label + '_deception'
	outfile.write(header + ',total_count,total_deception\n')

	print len(os.listdir(args.directory))
	counter = 0

	# for each split pither file in the directory
	for filename in os.listdir(args.directory):
		counter += 1
		print counter
		# open up the pitcher file
		infile = open(args.directory + '/' + filename)
		pid = filename.split('.')[0].split('_')[1]
		# print pid

		# create a dictionary of the data for this pitcher
		# keys are the pitch labels, entries are the four location values
		classData = {}
		for line in infile:
			line = line.strip().split(',')

			# ignore pitch outs
			if line[2] == 'PO':
				continue

			# if the label has not been seen before, initialize
			if line[2] not in classData:
				classData[line[2]] = np.array([[float(i) for i in line[3:-1]]])
			else:
				classData[line[2]] = np.vstack((classData[line[2]], np.array([[float(i) for i in line[3:-1]]])))

		# remove any classes with fewer than three pitches
		for c in classData.keys():
			if len(classData[c]) < 3:
				classData.pop(c)

		# if all classes have been removed, skip
		if len(classData) == 0:
			print "skipped pid " + pid
			continue

		# get the deception for this pitcher
		outfile.write(create_line(args.year, pid, pitcherDeception(classData)))


if __name__=="__main__":
	main()