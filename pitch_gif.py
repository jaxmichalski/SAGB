import numpy as np
import matplotlib.pyplot as plt
import scipy
import scipy.optimize
import scipy.stats
from scipy.stats import linregress
from scipy.stats import zscore
import matplotlib.lines as lines
import itertools
import csv
import sys
import math
from mpl_toolkits.mplot3d import Axes3D
from heapq import nsmallest
from heapq import nlargest
from operator import *
import pandas as pd
import glob
import os


#variables used in acceleration equation
pi=np.pi
m=.145 #kg
r=.0364 #meters
p=1.23 #kg/meters^3
A=pi*(r**2) #meters^2
B=(1./(2*m))*p*A #coefficient in acceleraton eq's
g=9.81

#y values of home plate, the pitcher's mound, and the decision points
mound=16.76 #in meters
home=.431 #in meters
dp1= 7.62 #25 feet from home
dp2=7.215 #23.8 feet from home


#to make typing sin and cos easier to use
def sin(x):
    return np.sin(x)
def cos(x):
    return np.cos(x)




####Chossing Pitchers######

#input:
#x0,z0,vx0,vy0,vz0,w,theta
#x0(horizontally across plate in meters), z0(vertically across plate), vx0 (in meters/sec),
#vy0, vz0, spin rate (radians/sec), spin direction (in radians)

# Curveball
p_ckershaw = [0.957,6.545,-0.62,-108.712,2.57,1849.671,347.094]

# Fastball
p_asanchez = [-1.219,6.076,1.45,-141.313,-7.45,2024.318,209.778]

# Slider
p_crodon = [1.467,6.131,-6.21,-130.198,-6.5, 1190.121,266.36]

# Slider
p_csale = [2.559, 5.679, -9.6, -114.459, -1.28, 1041.628, 302.356]

# Slider
p_kherrera = [-1.983, 5.924, 2.69, -118.296, -0.12, 2093.328, 65.298]

# Splitter
p_twalker = [2.122,5.992,5.25,-135.013,-6.97, 2163.977, 240.932]

#Curveball
p_rhill = [2.832,4.558, -8.73,-107.952,2.33, 1486.03, 289.566]

# Fastball
p_nsyndergaard = [-0.44, 6.231, 4.23,-145.797, -9.7, 2630.391, 198.865]

# Changeup
p_cwhitely = [-0.407, 5.892, 4.71, -123.094, -6.59, 1910.055, 244.976]

p_jverlander = [-2.1720000000000002, 6.3436000000000003, 3.5764, -119.59, 1.165, 2998.0, 170]

pitch_names = ['kershaw_curveball', 'sanchez_fastball', 'rodon_slider', 'sale_slider', \
'herrera_slider' , 'walker_splitter', 'hill_curveball', 'syndergaard_fastball', \
'whitely_changeup', 'verlander_curveball']

pitches = [p_ckershaw, p_asanchez, p_crodon, p_csale, p_kherrera, p_twalker, p_rhill, p_nsyndergaard, p_cwhitely, p_jverlander]


print 'Choose Pitch'
option = 1
for i in pitch_names:
    print str(option) + '.' + i
    option += 1

pitch_index = raw_input('Enter Pitch Number:')

x0, z0, vx0, vy0, vz0 , w, theta = pitches[int(pitch_index) - 1]






####Generating Trajectory######



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
initial_state = [x0, y0, z0, vx0, vy0, vz0 , w, theta]

    
#function that gives spin factor
def spinfactor(v,w): #v in m/s and w in rad/s
    S=(r*w)/float(v)
    return S
#function that gives the lift coefficient for the Magnus force, using the spin factory
def liftcoefficient(v,w):
    Cl=((17*spinfactor(v,w))/30.)+(7./75)
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

    v=float(((v_x**2)+(v_y**2)+(v_z**2))**.5) #to use for magnus force

    magnus = Magnus_components(v_x, v_y, v_z, w, theta)

    Cl=liftcoefficient(v,w)*(v**2) 

    #acceleration equations dervied from forces, F=ma
    a_x=B*((Cl*magnus[0])-(Cd*v_x*v_x))
    a_y=B*((Cl*magnus[1])+(Cd*v_y**2))
    a_z=B*((Cl*magnus[2])+(Cd*v_z**2))-g
    
    return [v_x,v_y,v_z,a_x,a_y, a_z]

#numerical integration to find x,y,z and vx,vy,vz after the initial position
def rk2(state, time, dt, derivs): 
    state_constant = state[6:]   
    state = state[0:6]
    k0 = [dt * x for x in derivs(state, time)]
    
    
    k1 = [dt * x for x in derivs(map(add, state, k0), time + dt)] 
    
    change = [0.5 * x for x in map(add, k0, k1)]

    state_next = map(add, state, change) + state_constant

    return state_next       

#numerical integration to find x,y,z and vx,vy,vz before the initial position (backwards in time)
def rk2_negative(state, time, dt, derivs): 
    state_constant = state[6:]   
    state = state[0:6]

    k0 = [dt * x for x in derivs(state, time)]
    k1 = [dt * x for x in derivs((map(sub, state, k0)), time - dt)]
    
    change = [0.5 * x for x in (map(add, k0, k1))]
    
    state_next = map(sub, state, change) + state_constant            

    return state_next

#function that numerically integrates until at target y value, input is (initial state [x,y,z,vx,vy,vz], step size, target y position, either rk2 or rk2_negative)
def iterativeIntegrate(start_coords, dt, target_y,rk):

    dif = abs(start_coords[1] - target_y) #original difference in y_initial and y_target
    old_coords = start_coords #sets start coordinates as old coordinates to calculate new dif with new coordinated
    trajectory=[start_coords] #initializing list of pitch trajectory
    while True: # iterate until return
        # calculate new coordinates and update the dif
        new_coords = rk(old_coords, 1, dt, accel) #using Runge-Kutta to find next coordinate
        new_dif = abs(new_coords[1] - target_y) #calculate new y difference using new coordinate

        # if the dif has gotten smaller, continue
        if new_dif < dif:
            trajectory.append(new_coords) #add coordinate to trajectory list
            dif = new_dif #update old dif with new dif
            old_coords = new_coords #reset and do it all over until...
        else: # otherwise, return the former coordinates, as they were closest
            return trajectory

def iterativeIntegrate_point(start_coords, dt, target_y,rk):
    
    dif = abs(start_coords[1] - target_y) #original difference in y_initial and y_target
    old_coords = start_coords #sets start coordinates as old coordinates to calculate new dif with new coordinated
    trajectory=[start_coords] #initializing list of pitch trajectory
    while True: # iterate until return
        # calculate new coordinates and update the dif
        new_coords = rk(old_coords, 1, dt, accel) #using Runge-Kutta to find next coordinate
        new_dif = abs(new_coords[1] - target_y) #calculate new y difference using new coordinate

        # if the dif has gotten smaller, continue
        if new_dif < dif:
            trajectory.append(new_coords) #add coordinate to trajectory list
            dif = new_dif #update old dif with new dif
            old_coords = new_coords #reset and do it all over until...
        else: # otherwise, return the former coordinates, as they were closest
            return new_coords


time=1 #total time integrating over
N= 4000#number of steps
dt = time/float(N) #dt of 0.00025 

#coordinates in "positive time," converted from list of arrays to vertically stacked array, y_target is home plate 
coordinates_pos = iterativeIntegrate(initial_state, dt, home, rk2)

#coordinates in "negative time," converted from list of arrays to vertically stacked array, y_target is pitcher's mound
coordinates_neg = iterativeIntegrate(initial_state, dt, mound, rk2_negative)

#decision point coordinates, converted to feet  
#DP1=3.28084*iterativeIntegrate_point(initial_state, dt, dp1, rk2)
#DP2=3.28084*iterativeIntegrate_point(initial_state, dt, dp2, rk2)

#converting to feet from meters, extracting into new lists
x_pos = [3.28084 * x[0] for x in coordinates_pos] #x values are first column
y_pos = [3.28084 * x[1] for x in coordinates_pos] #y values are second column
z_pos = [3.28084 * x[2] for x in coordinates_pos] #y values are third column
vx_pos =[3.28084 * x[3] for x in coordinates_pos] #vx values are fourth column
vy_pos = [3.28084 * x[4] for x in coordinates_pos] #vy values are fifth column
vz_pos = [3.28084 * x[5] for x in coordinates_pos] #vz values are 6th column

#converting to feet from meters, flipping lists so that order is in "positive time," extracting into new lists
x_neg = [3.28084 * x[0] for x in coordinates_neg] #x values are first column
y_neg = [3.28084 * x[1] for x in coordinates_neg] #y values are second column
z_neg = [3.28084 * x[2] for x in coordinates_neg] #y values are third column
vx_neg =[3.28084 * x[3] for x in coordinates_neg] #vx values are fourth column
vy_neg = [3.28084 * x[4] for x in coordinates_neg] #vy values are fifth column
vz_neg = [3.28084 * x[5] for x in coordinates_neg] #vz values are 6th column

x_neg = x_neg[::-1]
y_neg = y_neg[::-1]
z_neg = z_neg[::-1]


#combining positive and negative time coordinates
x_total = x_neg + x_pos
y_total = y_neg + y_pos
z_total = z_neg + z_pos

#total time of the pitch is dt times the amount of points in the entire trajectory set
time_total = len(y_total) * dt

#organzing into list of ([x,y,z])

coords_for_plot = zip(x_total, y_total, z_total)

name = pitch_names[int(pitch_index)-1]
filename = name + '_xyz.csv' 
df_write = pd.DataFrame(data = coords_for_plot, columns=['X', 'Y', 'Z'])

print len(coords_for_plot)
df_write.to_csv(filename)









####Creating Trajectory Plots######
cur_path = os.getcwd()

Location = cur_path + "/"+filename
df_read = pd.read_csv(Location)

x_total = df_read['X']
y_total = df_read['Y']
z_total = df_read['Z']

print 'X:', len(list(x_total))
print 'Y:', len(list(y_total))
print 'Z:', len(list(z_total))



#data = [x for x in data if len(x) == 3]

#else:
#    for x in data[-1]:
#        if type(x) is not float:
#            data.pop()

os.remove(name+'_xyz.csv')

#x_f = float(x_total[-1])
#y_f = float(y_total[-1])
#z_f = float(z_total[-1])
#35 total frames rather than 2000+

x_total = [x_total[i] for i in range(len(x_total)) if i%45==0] + [list(x_total)[-1]]
y_total = [y_total[i] for i in range(len(y_total)) if i%45==0] + [list(y_total)[-1]]
z_total = [z_total[i] for i in range(len(z_total)) if i%45==0] + [list(z_total)[-1]]

#x_total += [x_f]
#y_total += [y_f]
#z_total += [z_f]


print "Frames:",len(x_total)

#creating strikezone and pitchers mound visuals for plot
home = y_total[-1]
xl = -1
xr = 1
zb = 1.75
zt = 3.5
#creating strikezone and pitchers mound visuals for plot
x_strikezone_l = [xl for x in range(100)]
x_strikezone_r = [xr for x in range(100)]
x_strikezone = np.linspace(xl, xr, 100)

z_strikezone_bottom = [zb for z in range(100)]
z_strikezone_top = [zt for z in range(100)]
z_strikezone = np.linspace(zb, zt, 100)

y_strikezone = [home for y in range(100)]

y_mound = [60.5 for y in range(100)]
x_mound = np.linspace(-2, 2, 100)




def animation(angle1, angle2, axes, angle3):
    count = 0
    x_graph = []
    y_graph = []
    z_graph = []

    for i in range(len(x_total)):
     
    # Make the plot
        
        x_graph.append(x_total[i])
        y_graph.append(y_total[i])
        z_graph.append(z_total[i])

        fraction = count/float(len(x_total))
        count += 1
        
        ymax = int(round(y_graph[-1] + 5))
        yticks = range(ymax - 10, ymax, 1)
        
        print axes, count

        if axes == 'all':
            fig = plt.figure(figsize=(30,14))
            #fig.subplots_adjust(wspace=0.4)
            for n in range(1,4):
                ax = fig.add_subplot(1, 3, n, projection= '3d')
                PT=ax.plot(x_graph, y_graph, z_graph, lw=3, c="b")
                EP1=ax.scatter(x_graph[-1], y_graph[-1], z_graph[-1], s=60, c='k', marker="x")

                SZ_bottom=ax.plot(x_strikezone, y_strikezone, z_strikezone_bottom, c='k', lw=1)
                SZ_top=ax.plot(x_strikezone, y_strikezone, z_strikezone_top, c='k', lw=1)
                SZ_left=ax.plot(x_strikezone_l, y_strikezone, z_strikezone, c='k', lw=1)
                SZ_right=ax.plot(x_strikezone_r, y_strikezone, z_strikezone, c='k', lw=1)

                if n == 1:
                   
                    ax.set_xlabel('')

                    z_progress = list(np.linspace(0, fraction * 9, 100))
                    x_progress = [4 for _ in range(100)]
                    y_progress = [y_graph[-1] - 5 for _ in range(100)]
                    ax.plot(x_progress, y_progress, z_progress, lw=10, c="r", alpha=.3)
                    ax.set_xlabel('Horizontal Movement')
                    ax.set_zlabel('Height', rotation=angle3[n-1])
                if n == 2:
                	ax.set_ylabel('')
                  	if x_graph[0] > 0:
	                    factor = -1
	                else:
	                    factor = 1
	          
	                x_progress = list(np.linspace(-5*factor, -5*factor + (fraction*10)*factor, 100))
	                y_progress = [y_graph[-1]-60 for _ in range(100)]
	                z_progress = [-1 for _ in range(100)]
	                ax.plot(x_progress, y_progress, z_progress, lw=10, c="r", alpha=.3)
	             	ax.set_xlabel('Horizontal Movement')
	             	ax.set_zlabel('Height', rotation=angle3[n-1])
                   

                if n == 3:
                    x_progress = [-4 for _ in range(100)]
                    y_progress = list(np.linspace(y_graph[-1] + 5, y_graph[-1] + 5 - (fraction*10), 100))
                    z_progress = [0 for _ in range(100)]
                    ax.plot(x_progress, y_progress, z_progress, lw=10, c="r", alpha=.3)
                    ax.set_ylabel('Distance from Home Plate')
                    ax.set_xlabel('Horizontal Movement')

                
                ax.zaxis.set_rotate_label(False)
                
                ax.set_xlim3d(-5,5)
                ax.set_ylim3d(y_graph[-1] - 5, y_graph[-1] + 5)
                ax.set_zlim3d(-1,9)
                ax.set_yticks([y_graph[-1]])
                ax.set_zticks(range(0, 10))
                ax.set_frame_on(False)
             
                     
                ax.view_init(angle1[n-1] ,angle2[n-1])
                ax.dist = 7
            plt.subplots_adjust(wspace=.01, hspace=.01)
            filename = cur_path + '/' + str(count) + '.png'
            if i == len(x_total)-1:
            	for _ in range(10):
            		plt.savefig(filename, dpi=96, bbox_inches='tight',  pad_inches = 0, transparent = True)
            		count += 1
            		filename = cur_path + '/' + str(count) + '.png'
            else:
            	plt.savefig(filename, dpi=96, bbox_inches='tight',  pad_inches = 0, transparent = True)
            plt.close()



                
        
        else:
            fig = plt.figure(figsize=(8,10))
            #fig = plt.figure()
            ax = Axes3D(fig)

            PT=ax.plot(x_graph, y_graph, z_graph, lw=3, c="b") #pitch trajectory
            #RP=ax.scatter(x_total[0], y_total[0], z_total[0], s=50, c='y', marker="o") #release point
            EP1=ax.scatter(x_graph[-1], y_graph[-1], z_graph[-1], s=60, c='k', marker="x") #end point 
            #EP2=ax.scatter(x_total[-1], y_total[-1], z_total[-1], s=20, c='r', marker="o") #end point
            #DP1=ax.scatter(DP1[0], DP1[1], DP1[2], s=40, c='g', marker='o')


            #strike zone and mound
            SZ_bottom=ax.plot(x_strikezone, y_strikezone, z_strikezone_bottom, c='k', lw=1)
            SZ_top=ax.plot(x_strikezone, y_strikezone, z_strikezone_top, c='k', lw=1)
            SZ_left=ax.plot(x_strikezone_l, y_strikezone, z_strikezone, c='k', lw=1)
            SZ_right=ax.plot(x_strikezone_r, y_strikezone, z_strikezone, c='k', lw=1)
            #Mound=ax.plot(x_mound, y_mound, zs=0,c='silver', lw=5)

            if axes == 'xz':
                ax.set_ylabel('')

                if x_graph[0] > 0:
                    factor = -1
                else:
                    factor = 1
          
                x_progress = list(np.linspace(-5*factor, -5*factor + (fraction*10)*factor, 100))
                y_progress = [y_graph[-1]-60 for _ in range(100)]
                z_progress = [-1 for _ in range(100)]
                ax.plot(x_progress, y_progress, z_progress, lw=10, c="r", alpha=.3)
             
            else:
                ax.set_ylabel('Distance from Home Plate')

            if axes == 'yz':
                ax.set_xlabel('')

                z_progress = list(np.linspace(0, fraction * 9, 100))
                x_progress = [4 for _ in range(100)]
                y_progress = [y_graph[-1] - 5 for _ in range(100)]
                ax.plot(x_progress, y_progress, z_progress, lw=10, c="r", alpha=.3)
               
            
            else:
                ax.set_xlabel('Horizontal Movement')

            if axes == 'yx':
                x_progress = [-4 for _ in range(100)]
                y_progress = list(np.linspace(y_graph[-1] + 5, y_graph[-1] + 5 - (fraction*10), 100))
                z_progress = [0 for _ in range(100)]
                ax.plot(x_progress, y_progress, z_progress, lw=10, c="r", alpha=.3)


            ax.zaxis.set_rotate_label(False)
            ax.set_zlabel('Height', rotation=angle3)
            ax.set_xlim3d(-5,5)
            ax.set_ylim3d(y_graph[-1] - 5, y_graph[-1] + 5)
            ax.set_zlim3d(-1,9)
            ax.set_yticks([y_graph[-1]])
            ax.set_zticks(range(0, 10))
            ax.set_frame_on(False)
         
                 
            ax.view_init(angle1 ,angle2)

         
        # Save it
            filename = cur_path + '/' + str(count) + '.png'
            if i == len(x_total)-1:
            	for _ in range(10):
            		plt.savefig(filename, dpi=96, bbox_inches='tight',  pad_inches = 0, transparent = True)
            		count += 1
            		filename = cur_path + '/' + str(count) + '.png'
            else:
            	plt.savefig(filename, dpi=96, bbox_inches='tight',  pad_inches = 0, transparent = True)

            plt.close()

    file_list = glob.glob('*png')

    list.sort(file_list, key=lambda x: int(x.split('.png')[0]))
    with open('image_list.txt', 'w') as file1:
        for item in file_list:
            file1.write("%s\n" % item)
    

    outfile = cur_path + '/' + name +'/' + name + str(axes) + '.gif'
    os.system('convert -dispose previous -delay 10 -loop 0 @image_list.txt ' + outfile)
   
    os.remove('image_list.txt')


    print name + str(axes) + '.gif'
    for n in range(1, len(x_total) + 10):
       os.remove(str(n) + '.png')


path = cur_path + '/' + name

if not os.path.exists(path):
    os.mkdir(path)



animation(0, 90, 'xz', 90)
animation(0, 0, 'yz', 90)
animation(90, 180, 'yx', 0)
animation((0,0,90), (0,90,180), 'all', (90,90,0))






