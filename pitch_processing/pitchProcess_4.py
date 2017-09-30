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
import argparse

# properties of the baseball
m=.145 #kg
r=.0364 #meters
p=1.23 #kg/meters^3
A=np.pi*(r**2) #meters^2

global B
B=(1./(2*m))*p*A
global dragCoeff
dragCoeff = 0.3
global g
g = 9.81
global DIST_STEP
DIST_STEP = 0.1524 # half a foot in meters

# performs conversions on a coordinate array
# input is list in the following format:
# x0,z0,vx0,vy0,vz0,w,theta
# returns a list in the same format
def metricConvert(coords):
	coords[0] *= 0.3048
	coords[1] *= 0.3048
	coords[2] *= 0.3048
	coords[3] *= 0.3048
	coords[4] *= 0.3048
	coords[5] *= 0.104719755
	coords[6] *= np.pi /  180
	return coords

# caculates the spin factor from velocity and omega
def spinFactor(v,w): #v in m/s and w in rad/s
    S = (r*w) / v
    return S

# calculates the lift coefficient of the magnus force
def liftCoeff(v,w):
    Cl = ((17 * spinFactor(v,w)) / 30.) + (7. / 75)
    return Cl

# calculates the magnus components given velocity, spin, and angle
def magnusComponents(v_x, v_y, v_z, w, theta):
	# cross products for each component
	cross = [-v_y*np.sin(theta), v_z*np.cos(theta) - v_x*np.sin(theta), v_y*np.cos(theta)]
	cross = [w*i for i in cross]

	# calculate denominator of the magnus components
	magnus_d = (((cross[0]**2) + (cross[1]**2) + (cross[2]**2))**.5)
	magnus = [i / magnus_d for i in cross]

	return magnus

# calculates the acceleration of each component
# in_coords is np.array with format:
# x_pos, y_pos, z_pos, x_vel, y_vel, z_vel
def acceleration(in_coords, w, theta):
	v_x = in_coords[3]
	v_y = in_coords[4]
	v_z = in_coords[5]

	# calculates absolue velocity
	v = ((v_x**2) + (v_y**2) + (v_z**2))**0.5

	# gets the magnus components
	magnus = magnusComponents(v_x, v_y, v_z, w, theta)

	# precalculate whatever possible
	Cl = liftCoeff(v,w)*v**2

	a_x = B*(Cl*magnus[0] - dragCoeff*v_x*abs(v_x))
	a_y = B*(Cl*magnus[1] + dragCoeff*v_y**2)
	a_z = B*(Cl*magnus[2] + dragCoeff*v_z**2) - g


	return np.array([v_x, v_y, v_z, a_x, a_y, a_z])

# finds the next set of coordinates through numerical integration
# takes in a np.array in the following format:
# x_pos, y_pos, z_pos, x_vel, y_vel, z_vel
def nextCoordinates(in_coords, dt, w, theta):
	k0 = dt*acceleration(in_coords, w, theta)
	k1 = dt*acceleration(in_coords + k0, w, theta)
	next_coords = in_coords + 0.5*(k0 + k1)


	return next_coords

# while loop construction for iteratively finding
# position closest to target y
def iterativeIntegrate(start_coords, w, theta, dt, target_y):
	dif = abs(start_coords[1] - target_y)
	old_coords = start_coords

	time = 0.0

	while True: # iterate until return
		# calculate new coordinates and update the dif
		new_coords = nextCoordinates(old_coords, dt, w, theta)
		new_dif = abs(new_coords[1] - target_y)

		# if the dif has gotten smaller, continue
		if new_dif < dif:
			dif = new_dif
			old_coords = new_coords
			time += dt
		else: # otherwise, return the former coordinates, as they were closest
			print time
			print old_coords
			sys.exit()
			return old_coords

# receives a start_state in following format:
# x0,z0,vx0,vy0,vz0,theta,w
def pitchProcess(start_state):
	# converts the start state to metric
	start_state = metricConvert(start_state)
	y = 15.24 # starting y in meters
	theta = start_state[6] # angle
	w = start_state[5] # spin rate

	# extracts the elements of the start state needed for coordinates
	start_coords = np.array([start_state[0], y] + start_state[1:5])

	# setting dt
	dt = 0.5 / 1000

	# append the start position as the first result (50')
	results = []
	results.append(start_coords)
	# print results
	while y > 0: # gets new values until y passes 0 (home plate)
		y = y - DIST_STEP # decrement y by half a foot
		results.append(iterativeIntegrate(results[-1], w, theta, dt, 0))
		# print results[-1]
		# raw_input()

	out = []
	for r in results:
		out.append(r[0])
		out.append(r[2])

	# include the spin theta in the final result
	out.append(start_state[6])

	return out
	
# processes a baseball savant pitcher file into an array of arrays
# each element is an array with the relavent values
def fileProcess(filename):
	infile = open(filename)
	infile.readline() # read past the header line

	reader = csv.reader(infile, delimiter=',', quotechar='"')

	# pitcher_id, pitcher_handedness, pitch_type
	str_indices = [18,20,57]
	# x0, z0, vx0, vy0, vz0, spin_rate, spin_dir
	q_indices = [45,47,48,49,50,62,61]

	data = []

	for line in reader:
		new_data = []
		# add the str_indices to the front of new_data
		for s in str_indices:
			new_data.append(line[s])

		# add the q_indices as floats
		try: # do this within a try/except
			for i in q_indices:
				new_data.append(float(line[i].strip()))
			data.append(new_data)
		except ValueError: # if there was a value error, the line was bad
			continue

	return data

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('input', help="Input Baseball Savant file.", type=str)
	parser.add_argument('output', help="Output file.", type=str)
	args = parser.parse_args()

	print args.input
	data = fileProcess(args.input)

	# array of arrays
	# pitch_type, then pairs of x/z at 25', then 23'8"
	results = []
	counter = 0
	for d in data:
		counter += 1
		if counter % 1000 == 0:
			print counter
		new_result = d[:3]
		new_result += pitchProcess(d[3:])
		results.append(new_result)

	# file to write out
	outfile = open(args.output, 'w')
	# create header line
	header = '##ID,HAND,PITCH_TYPE'
	for i in sorted(range(101), reverse=True):
		header += ',' + str(float(i)/2) + ','
	header += ',SPIN_DIR\n'
	outfile.write(header)

	for r in results:
		outfile.write(','.join([str(i) for i in r]) + '\n')

	outfile.close()

if __name__=="__main__":
 	main()

# print pitchProcess([-1.983, 5.924, 2.69, -118.296, -0.12, 2093.328, 65.298], [7.62, 7.215])
# print pitchProcess([2.489, 5.889, -10.06, -133.987, -5.45, 152.584, 1783.096], [7.62, 7.215])
