import numpy as np 
import sys
import math

# dimensions
global d
d = 4

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

# takes in an array of arrays, data
# the first col is the labels (pitch types)
def pitcherDeception(data):
	# dictionary of classes; key is label, val is nparray of samples
	# also counts the number of samples
	classData = {}
	num_samples = 0
	for d in data:
		num_samples += 1
		# if the label has not been seen before, initialize
		if d[0] not in classData:
			classData[d[0]] = np.array([d[1:]]) 
		else: # otherwise, vstack
			classData[d[0]] = np.vstack((classData[d[0]], np.array([d[1:]])))

	# calculate mus, pi_cs, and sigmas for each class
	classes = [] # list of classes
	pis = {} # rates of each pitch
	mus = {} # mus for each class
	sigmas = {} # covariance for each class

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
		sigmas[c] = cov

	# calculate the meanSigma, then its inverse
	meanSigma = np.zeros((4,4))
	for label in classes:
		meanSigma = np.add(meanSigma, pis[label]*sigmas[label])
	inv_meanSigma = np.linalg.inv(meanSigma)
	det_meanSigma = np.linalg.det(meanSigma)

	# calculate the average deception (aka misclassification P) for each class
	for label in classes:
		# the data points for the class
		locs = classData[label]
		class_deception = 0
		for pitch in locs:
			class_deception += pitchDeception(np.transpose(np.array([pitch])), classes, pis, mus, inv_meanSigma, det_meanSigma, label)
		
		# divide by count
		class_deception = class_deception / len(locs)
		print label 
		print len(locs)
		print class_deception

def main():
	# read in processed data and store it as an array of arrays
	infile = open('bumgarner_processed.txt')
	data = []
	for line in infile:
		new_data = []
		line = line.split('\t')
		new_data = [line[0]] + [float(i) for i in line[1:]]
		data.append(new_data)

	pitcherDeception(data)

if __name__=="__main__":
	main()
