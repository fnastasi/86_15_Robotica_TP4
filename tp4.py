#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 15:28:29 2019

@author: franco
"""

# In[]
"""
INFO IMPORTANTE
    roll: http://softmc.servotronix.com/wiki/SCARA_robot

"""
# In[]

from numpy import *
import matplotlib.pyplot as plt
from tp1_modificado import *

# In[]

a1 = 200
a2 = 200
a = 200 # Como a1 = a2 se define a = a1 = a2

# In[]
# Defino una función que me devuelve la matriz homogenea que representa 
# la rototraslación a partir del criterio Denavit-Hartenberg 

# ángulos deben estar en radianes
def DH_hom_mat(theta,d,a,alpha):
    
    # Se ponen en rango [-180,180)los ángulos
    #theta = get_in_range(theta)
    #alpha = get_in_range(alpha)
    
    
    R = array(
        [[cos(theta), -sin(theta)*cos(alpha), sin(theta)*sin(alpha)],
        [sin(theta), cos(theta)*cos(alpha), -cos(theta)*sin(alpha)],
        [0,sin(alpha) , cos(alpha)]])
    P = array([[a*cos(theta),
               a*sin(theta),
               d]]).T
    A = concatenate((R,P),axis = 1) # Junta R y P -> [R | P]
    
    # Junta Ry P con el vector [0001] -> [  R    | P ]
    #                                    [ 0 0 0 | 1 ]
    
    A = concatenate((A,array([[0,0,0,1]])), axis = 0) 
    A[absolute(A)<tol] = 0
    return A

# In[]
def prob_dir_scara(t1,t2,d3,t4):
    
    [t1,t2,t4] = array([t1,t2,t4])*pi/180
    
    A_1_0  = DH_hom_mat (t1,0,a1,0)
        
    A_2_1  = DH_hom_mat (t2,0,a2, 0)
    
    A_3_2  = DH_hom_mat (0,d3,0, 0)
        
    A_4_3  = DH_hom_mat (t4,0,0, 0)
        
    A = linalg.multi_dot([A_1_0, A_2_1, A_3_2, A_4_3])

    g = 1 if sign(sin(t2)) >= 0 else -1
    
    return [A,g]


# In[]
    
def prob_inv_scara(A,g):
    
    px,py,pz = A[0:3,3]
    d3 = pz
    c2 = (px**2 + py**2 - a1**2 - a2**2 )/(2*a1*a2) 
    if abs(c2) <= 1:
        s2 = g*sqrt(1-c2**2)
        t2 = arctan2(s2,c2)
    else:
        print("Error: Posición no alcanzable")
        return [0,0,0,0]
        
    P = array([px, py])
    M = array([[a2*c2+a1, -a2*s2],[a2*s2,a2*c2+a1]])
    c1,s1 = linalg.solve(M,P)
    t1 = arctan2(s1,c1)
    
    #R_1_0  = DH_hom_mat(t1,0,a1,0)[0:3,0:3]
    #R_2_1  = DH_hom_mat(t2,0,a2, 0)[0:3,0:3]
    #R_2_0 = dot(R_1_0,R_2_1)
    #R_4_0 = A[0:3,0:3]
    #R_4_3 = dot(R_2_0.T,R_4_0)
    #t4 = arctan2(R_4_3[1,0],R_4_3[1,1])
    
    t124 = arctan2(A[1,0],A[0,0])
    t4 = t124 - t1- t2
    
    t1 = t1*180/pi
    t2 = t2*180/pi
    t4 = t4*180/pi
    
    tope_mec = t2>= 150 or t2<= -150 or d3>=-50 or d3<=-250
    if tope_mec:
        print("Atención: Para alcanzar la POSE, se sobrepasan los topes mecánicos por lo que no es alcanzable")
    return t1,t2,d3,t4

# In[]
    
def resolv_inv_scara(x,y,z,roll,g):
    R = array ([[cos(roll),-sin(roll), 0],
                [sin(roll),cos(roll),0],
                [0,0,1]])
    P = array([[x,y,z]])
    A = concatenate( (R,P.T),axis=1)
    A = concatenate((A,array([[0,0,0,1]])))
    print(A)
    return prob_inv_scara(A,g)


# In[]
    
# Prueba problema indirecto:
    x=200
    y=200
    z=200
    roll= 0
    g=1
    t1,t2,d3,t4 = resolv_inv_scara(x,y,z,roll,g)


# In[]
    
    tacc = 200
    v1max = 90
    v2max = 180
    v3max = 1000
    v4max = 360
    vmax = array([v1max,v2max,v3max,v4max])
    dt = 1e-4   # diferencial de tiempo para hacer los gráficos.
                # Es conveniente que divida a tacc
    
# In[]

def gen_tray(POSE_vec,td_vec):
    
    
    
    
    
    
    
    
    
    
    
    
    
# In[]
    
def gen_tray_POSE(var_art_pos,td):
    
    
    """
    Debe devolver:  - La matriz con los valores de las variables articulares
                    - La matriz con los valores de las derivadas de las variables articulares
                    - vector de los valores de T en cada tramos (esto es para graficar después)
    """
    
    # var_art_pos es vector de los valores a donde quiere ir las variables articuladas
    #   var_art_pos[0] es la posición de la que sale. Debe tener minimo 2 valores   
    # td tiempo deseado 
    
    
    cant_tray = len(var_art_pos[0,:]) -1 # Cantidad de trayectorias
    
    # Inicialización de las variables articulares
    var_art = array([],
                    [],
                    [],
                    [])
    
    # Inicialización de matrix A, B y C
    A = zeros(4,cant_tray)
    B = zeros(4,cant_tray)
    C = zeros(4,cant_tray)
    A[:,0] = var_art_pos[:,0]
    
    for i in arange(cant_tray):
        
        # Seteo valores de B, C y delta A y delta C en cada trayectoria
        B = var_art_pos[:,i]
        C = var_art_pos[:,i+1]
        DA = A-B
        DC = C-B
        
        # Calculo del tiempo en que se mueven los ejes
        T = amax([2*tacc,DC/vmax,td[i]])
        
        # Calculo de las variables articulares en la zona 1
        var_art_zoneI = zeros(4, )
    
    
    
    
    
    
