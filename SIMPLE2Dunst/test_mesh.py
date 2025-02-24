#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 14 19:16:41 2024

@author: andres cimino

prueba simple de las clases vertex, face y element 
para una malla simple de 1 cuadrilatero y un triangulo
"""

#import numpy as np
import matplotlib.pyplot as plt
import data_struc as ds

v1=ds.vertex(coor_x=0,coor_y=0, index=0)
v2=ds.vertex(coor_x=1,coor_y=0, index=1)
v3=ds.vertex(coor_x=2,coor_y=0, index=2)
v4=ds.vertex(coor_x=1,coor_y=1, index=3)
v5=ds.vertex(coor_x=0,coor_y=1, index=4)
vertices=[v1,v2,v3,v4,v5]


## faces
face1=ds.face(v1,v2, 1,-1, index=0)
face2=ds.face(v2,v3, 2,-1, index=1)
face3=ds.face(v3,v4, 2,-1, index=2)
face4=ds.face(v2,v4, 1,2, index=3)
face5=ds.face(v4,v5, 1,-1, index=4)
face6=ds.face(v5,v1, 1,-1, index=5)


faces= [face1, face2, face3, face4, face5, face6]


n1=face1.normal(); n2=face2.normal(); n3=face3.normal();
n4=face4.normal(); n5=face5.normal(); n6=face6.normal();
normals=[n1,n2,n3,n4,n5, n6]



c1=face1.center(); c2=face2.center(); c3=face3.center();
c4=face4.center(); c5=face5.center(); c6=face6.center();
facecenters=[c1,c2,c3,c4,c5, c6]

# elementos
connect_faces1=[1, 4, 5, 6]
connect_faces2=[2,3, 4]

connect_cells1=[2] ; connect_cells2=[1]

cell1=ds.cell(connect_faces1, connect_cells1, index=1 )
cell2=ds.cell(connect_faces2, connect_cells2, index=2 )

centroid1=cell1.centroid([face1, face4, face5, face6])
centroid2=cell2.centroid([face2, face3, face4])

centroids=[centroid1, centroid2]
""" graficos"""
plt.figure(1)

for face_i in faces:
#     # grafico de lineas entre los vertices que componen cada cara (para ver que esten todas las caras)
    plt.plot([face_i.v1.x, face_i.v2.x],  [face_i.v1.y , face_i.v2.y] ) 

for vertex_i in vertices:
    #grafico de posiciones de los vertices (con circulos) 
    plt.plot( vertex_i.x, vertex_i.y,   'o'   )

for i in range(len(normals)):
    plt.quiver(facecenters[i][ 0],facecenters[i][1], normals[i][0], normals[i][1] )

for centroid_i in centroids:
    plt.plot(centroid_i[0], centroid_i[1], 'xr')
#plt.legend()
plt.axis("equal")