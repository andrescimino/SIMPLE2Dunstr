#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 14 17:37:26 2024

@author: sebastian
"""

class Point(object):
    def __init__(self, x, y):
        self.x = x
        self.y = y

class Line(object):
    def __init__(self, p1:Point, p2:Point, x,y):
        self.p1 = p1
        self.p2 = p2
        self.p3=Point(x, y)
        
    def slope(self):
        return (self.p2.y - self.p1.y) / (self.p2.x - self.p1.x) #which value is being accessed with self.p2.y
    
    def intercept(self):
        return (self.p2.y + self.p1.y) /2 + (self.p2.y - self.p1.y) / (self.p2.x - self.p1.x)*(self.p2.x + self.p1.x) /2 #which value is being accessed with self.p2.y  
point1=Point(1,1)

point2=Point(2,2)

line1=Line(p1=point1,p2= point2, x=1,y=2)
print("slope=", line1.slope())
print("intercept=", line1.intercept())