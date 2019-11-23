#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 21:58:07 2019

@author: franco
"""

# In[]
def main():
    POSE1 = array([200,200,-100,0,1])
    POSE2 = array([200,200,-200,0,-1])
    POSE_vec = array([POSE1,POSE2,])
    td = array([1])
    dt = 0.1 # en mm
    gen_tray(POSE_vec,td,dt)


# In[]
from tp4 import *   

if __name__ == "__main__":
    main()
    
