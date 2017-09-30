
import numpy as np
import matplotlib.pyplot as plt
import scipy
import scipy.optimize
import matplotlib.pyplot as plt
import matplotlib.lines as lines
from mpl_toolkits.mplot3d import Axes3D
import itertools
import csv
import sys
import timeit



#variables used in acceleration equation
pi=np.pi
m=.145 #kg
r=.0364 #meters
p=1.23 #kg/meters^3
A=pi*(r**2) #meters^2
B=(1./(2*m))*p*A
g=9.81


#to make typing sin and cos easier to use
def sin(x):
    return np.sin(x)
def cos(x):
    return np.cos(x)


def DecPoi(x0,z0,vx0,vy0,vz0,w,theta):
	#Coverting from feet to meters, ft/sec to meters/sec, degrees to radians, and rev/min to radians/sec
	x0*=0.3048
	y0=15.24
	z0*=0.3048
	vx0*=0.3048
	vy0*=0.3048
	vz0*=0.3048
	w*=0.104719755
	theta*=(pi/180.)

	#Setting the initial state of the ball
	initial_state = np.array([x0,y0,z0,vx0,vy0,vz0]) #inital state (x0,y0,z0,vx0,vy0,vz0, w, theta))

	    
	#function that gives spin factor
	def S(v,w): #v in m/s and w in rad/s
	    S=(r*w)/v
	    return S
	#function that gives the lift coefficient for the Magnus force, using the spin factory
	def Cl(v,w):
	    Cl=((17*S(v,w))/30.)+(7./75)
	    return Cl

	Cd=0.3 #known value for drag coefficient of baseball


	#function that gives direction of the Magnus force, which is in the direction of w cross v
	def Magnus_components(v_x, v_y, v_z, w, theta):
	    cross=[-w*v_y*sin(theta),(w*v_z*cos(theta)-w*v_x*sin(theta)), w*v_y*cos(theta)]
	    Magnusx=(cross[0])/((cross[0]**2+cross[1]**2+cross[2]**2)**.5)
	    Magnusy=(cross[1])/((cross[0]**2+cross[1]**2+cross[2]**2)**.5)
	    Magnusz=(cross[2])/((cross[0]**2+cross[1]**2+cross[2]**2)**.5)
	    return [Magnusx, Magnusy, Magnusz] #unit vector for Magnus force direction

	def accel(state,time): #input is state (x,y,z,v_x,v_y,v_z) and time, returns array with (v_x,v_y,v_z,a_x,a_y,a_z) from force equation for ball in air
	    global a_x
	    global a_y
	    global a_z
	    v_x=state[3]
	    v_y=state[4]
	    v_z=state[5]
	    v=((v_x**2)+(v_y**2)+(v_z**2))**.5
	    a_x=B*((Cl(v,w)*v**2*Magnus_components(v_x,v_y,v_z,w,theta)[0])-(Cd*v_x*abs(v_x)))
	    a_y=B*((Cl(v,w)*v**2*Magnus_components(v_x,v_y,v_z,w,theta)[1])+(Cd*v_y**2))
	    a_z=B*((Cl(v,w)*v**2*Magnus_components(v_x,v_y,v_z,w,theta)[2])+(Cd*v_z**2))-g
	    return np.array([v_x,v_y,v_z,a_x,a_y, a_z])

	#numerical integration to find x,y,z and vx,vy,vz after the initial position
	def rk2(y, time, dt, derivs): 
	    k0 = dt*derivs(y, time)
	    k1 = dt*derivs(y+k0, time+dt) 
	    y_next = y+0.5*(k0+k1)
	    return y_next        
	    

	time=1 #total time integrating over
	N= 3000#number of steps
	dt = time/float(N) #dt of 0.003 

	answerRK = np.zeros([N,6]) #intializing output of rk2, six columns of zeros, N entries long
	answerRK[0,:] = initial_state #first column is initial state
	x1=answerRK[:,0] #x values are first column
	y1=answerRK[:,1] #y values are second column
	z1=answerRK[:,2] #y values are third column
	vx1=answerRK[:,3] #vx values are fourth column
	vy1=answerRK[:,4] #vy values are fifth column
	vz1=answerRK[:,5] #vz values are 6th column

	#going through rk N times to get through entire wanted_times
	for j in range (N-1): 
	    answerRK[j+1] = rk2(answerRK[j], 0, dt , accel)
	 
	#finding index number for when y=1.417 ft for pitch, this will tell us when the ball has reached the front of homeplate
	def yf_index(y):
	    # for i in y[1:]:
	    #     if i<=0.437 and i>=.421: #when 0<y<0.03, the ball is at homeplate
	    #         index_num=y.tolist().index(i)
	    # return index_num 
	    upper_bound = len(y) - 1 
	    lower_bound = 0
	    i = len(y)/2
	    while y[i] > 0.437 or y[i] < .421:
	    	if y[i] > 0.437:
	    		lower_bound = i
	    		i = (i + upper_bound) / 2
	    	else: 
	    		upper_bound = i
	    		i = i / 2
	    return i 

	 

	#converting after the initial point data to feet
	y_index = yf_index(y1)
	x_data1=3.28084*x1[1:y_index]
	y_data1=3.28084*y1[1:y_index]
	z_data1=3.28084*z1[1:y_index]



	#creating list of x,y,z coordinates
	coordinates=np.zeros([len(y_data1),3]) 
	for i in range(0,len(y_data1)):
	    coordinates[i]=[x_data1[i], y_data1[i], z_data1[i]]

	def ydecision_index(y):
	    # for i in y[1:]:
	    #     if i<=24 and i>=23.7:
	    #         index_num_dec=y.tolist().index(i)
	    # return index_num_dec
	    upper_bound = len(y) - 1 
	    lower_bound = 0
	    i = len(y)/2
	    while y[i] > 24 or y[i] < 23.7:
	    	if y[i] > 24:
	    		lower_bound = i
	    		i = (i + upper_bound) / 2
	    	else: 
	    		upper_bound = i
	    		i = i / 2
	    return i 

	def second_index(y):
	    # for i in y[1:]:
	    #     if i<=24 and i>=23.7:
	    #         index_num_dec=y.tolist().index(i)
	    # return index_num_dec
	    upper_bound = len(y) - 1 
	    lower_bound = 0
	    i = len(y)/2
	    while y[i] > 25.1 or y[i] < 24.8:
	    	if y[i] > 25.1:
	    		lower_bound = i
	    		i = (i + upper_bound) / 2
	    	else: 
	    		upper_bound = i
	    		i = i / 2
	    return i 

	decision_point=coordinates[ydecision_index(y_data1)]
	second_z_point = coordinates[second_index(y_data1)]
	final_coords=coordinates[len(y_data1)-1]
	two_valued_array = np.append(decision_point,second_z_point)

	return two_valued_array


def DPexporter(PitcherData,ExportFile):
	x0=[]
	z0=[]
	vx0=[]
	vy0=[]
	vz0=[]
	w=[]
	theta=[]

	with open(PitcherData, 'rb') as DP:
		Data = csv.reader(DP, delimiter=',', quotechar='|')
		for i in Data:
			x0+=[float(i[0])]
			z0+=[float(i[1])]
			vx0+=[float(i[2])]
			vy0+=[float(i[3])]
			vz0+=[float(i[4])]
			w+=[float(i[5])]
			theta+=[float(i[6])]

	#Decision Point and Second point data
	DPx=np.zeros(len(x0))
	DPy=np.zeros(len(x0))
	DPz=np.zeros(len(x0))
	SPx=np.zeros(len(x0))
	SPy=np.zeros(len(x0))
	SPz=np.zeros(len(x0))

	for i in range(0,len(x0)):	
		x = DecPoi(x0[i],z0[i],vx0[i],vy0[i],vz0[i],w[i],theta[i])
		DPx[i]+=[x[0]]
		DPy[i]+=[x[1]]
		DPz[i]+=[x[2]]
		SPx[i]+=[x[3]]
		SPy[i]+=[x[4]]
		SPz[i]+=[x[5]]



	sys.stdout = open(ExportFile, 'w')
	for i in range(len(DPx)):
		print DPx[i],',', DPy[i],',',DPz[i]
		print SPx[i],',', SPy[i],',',SPz[i]
	return 

DPexporter('KershawCurve.csv', 'KershawCurveClusterData.csv')


