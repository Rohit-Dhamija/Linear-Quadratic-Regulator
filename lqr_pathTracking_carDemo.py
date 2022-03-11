#!/usr/bin/env python
import rospy
from prius_msgs.msg import Control    
from nav_msgs.msg import Odometry 
from nav_msgs.msg import Path
import numpy as np    
import math   
from scipy.spatial import KDTree
from scipy.linalg import solve_discrete_are                                                                
from scipy.signal import cont2discrete, lti, dlti, dstep
from tf.transformations import euler_from_quaternion 

import matplotlib.pyplot as plt																								

pi = math.pi

global x,y,xm,ym,V,Vx,Vy,theta,throttle,steer    
global total_path_points                                                                                                                                    
total_path_points = 0                                                                                                                                                                            
global path                                                                                                                                       

V_ref = 10	
kp,ki,kd = 5,3,0.0001
delta_T = 0.1

Q1 = 1500 
Q2 = 1500
Q3 = 1500
Q4 = 1500

R = 100


error_allowed = 0.1

steer_constant = 40*pi/180                                                          # steer_constant = steer/steer_input
m = 1380                                                                            # mass of vehicle in kg
l_car = 3
mf = mr = m/2 																		                                                                                          
lf = l_car*(1-mf/m)
lr = l_car*(1-mr/m)
cf = cr = 188225                                                                                          

I_z = lf*lr*(mf + mr)





def pid():
	
	global throttle
	e_prev,sum_e = 0,0
	
	e =  V_ref - V
	sum_e = sum_e + e
	throttle = kp*e + ki*sum_e + kd*(e - e_prev)
	
	throttle = min(throttle,1)
	throttle = max(throttle,-1)
 
	
def pathfunc(Path):                                                                                                                                       
	
	global total_path_points,path
	if total_path_points == 0:
		
		total_path_points = len(Path.poses)
		path = np.zeros((total_path_points,2))													

		for i in range(0,total_path_points):                                                                                                                
			path[i][0] = Path.poses[i].pose.position.x
			path[i][1] = Path.poses[i].pose.position.y   

def odomfunc(odom):
	
	global x,y,xm,ym,V,Vx,Vy,theta
	x = odom.pose.pose.position.x                                                                                   ## odom of the rear wheel
	y = odom.pose.pose.position.y 
	Vx = odom.twist.twist.linear.x
	Vy = odom.twist.twist.linear.y
	V = math.sqrt(Vx**2 + Vy**2)

	quaternions =  odom.pose.pose.orientation                                                       
	quaternions_list = [quaternions.x,quaternions.y,quaternions.z,quaternions.w]
	roll,pitch,yaw = euler_from_quaternion(quaternions_list)
	theta = yaw

	xm = x + (l_car/2)*math.cos(theta)                                                                                   ## odom of centre
	ym = y + (l_car/2)*math.sin(theta)




def my_mainfunc():
	rospy.init_node('mpc_multipleShooting_pathTracking_carDemo', anonymous=True)
	rospy.Subscriber('/base_pose_ground_truth' , Odometry, odomfunc)            
	rospy.Subscriber('/astroid_path', Path, pathfunc)

	instance = rospy.Publisher('prius', Control, queue_size=10)

	rate = rospy.Rate(10)
	rate.sleep()  

	msg = Control()


	ecg_previous,theta_e_previous = 0,0
	
	while ( math.sqrt( (xm-path[total_path_points-1][0])**2 + (ym-path[total_path_points-1][1])**2 ) )  > error_allowed:		
	

		A = np.array([ (0 , 1                              , 0                   , 0                                   ),
					   (0 , -(cf+cr)/(m*(V+1e-9))          , (cf+cr)/m           , (lr*cr-lf*cf)/(m*(V+1e-9))          ),  
					   (0 , 0                              , 0                   , 1                                   ),
					   (0 , (lr*cr-lf*cf)/(I_z*(V+1e-9))   , (lf*cf-lr*cr)/I_z   , -(lf**2*cf+lr**2*cr)/(I_z*(V+1e-9)) ) ])

			
		B = np.array([ [0]         ,
		               [cf/m]      ,
					   [0]         ,
					   [lf*cf/I_z]  ])                

		
		C = np.array([(1., 0,0,0)])																				####### C and D are waste matrices bcz the scipy.cnt2discrete takes 4 matrices
		D = np.array([[0.]])																					#######

		
	
		l_system = lti(A, B, C, D)


		
		d_system = cont2discrete((A, B, C, D), delta_T, method="bilinear")
		s, x_d = dstep(d_system)
		

		A_disc = d_system[0]
		B_disc = d_system[1]

		print ("VELOCITY (in m/s) =", round(V,2))
				
		Q = np.eye(4, dtype = 'f')
		Q[0,0] = Q1
		Q[1,1] = Q2
		Q[2,2] = Q3
		Q[3,3] = Q4
		
		P = solve_discrete_are(A_disc,B_disc,Q,R)

		F = -(   np.dot(    np.linalg.inv( R +  np.dot( B.T,np.dot(P,B) ) )  ,    np.dot( B.T,np.dot(P,A) )      )    )                                                                                                                    
		
		ecg,close_index = KDTree(path).query(np.array([xm,ym])) 					 

																				
		if close_index == total_path_points:
			theta_path = math.atan( (path[close_index][1]-path[close_index-1][1]) / (path[close_index][0]-path[close_index-1][0]+ 1e-10) )   
		else:
			theta_path = math.atan( (path[close_index+1][1]-path[close_index][1]) / (path[close_index+1][0]-path[close_index][0]+ 1e-10) )     

		theta_e = theta_path - theta


		X = np.array([ (ecg)                      ,
		               (ecg-ecg_previous)         , 
		               (theta_e)                  , 
		               (theta_e-theta_e_previous)  ])



		ecg_previous = ecg
		theta_e_previous = theta_e

		

		steer = (np.dot(-F,X))[0]
	
		steer = min(steer,40*pi/180)
		steer = max(steer,-40*pi/180)
		steer_input = steer/steer_constant

		pid()
		global throttle
	   
		
		msg.throttle = throttle                                  
		msg.brake = 0.0 
		msg.steer = steer_input
		msg.shift_gears =2
		if throttle < 0:
			msg.shift_gears =3                                              # reverse gear
			throttle = -throttle        

		instance.publish(msg)
		rate.sleep()


	print ("PATH TRACKED")
	msg.throttle = 0                                  
	msg.brake = 1 
	msg.steer = 0
	msg.shift_gears =1                                                                  # nuetral gear    
	instance.publish(msg)



if __name__ == '__main__':
   
	try:
		my_mainfunc()
	except rospy.ROSInterruptException:
		pass
