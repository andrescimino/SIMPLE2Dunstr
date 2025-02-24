#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 29 19:07:31 2024

@author: acimino
"""

import numpy as np
import matplotlib.pyplot as plt
import data_struc as ds

vertex_file="points2.dat"
face_file="faces3.dat"
owners_file="owners3.dat"
neighbours_file="neighbours3.dat"
boundaries_file="boundaries3.dat"


the_mesh=ds.mesh(vertex_file, face_file, owners_file, neighbours_file, boundaries_file)
# cálculos de propiedades de áreas
[face.get_center() for face in the_mesh.faces]
[face.get_normal() for face in the_mesh.faces]
[face.get_area() for face in the_mesh.faces]

# cálculos de propiedades de celdas
[cell.get_faces_signs([the_mesh.faces[j] for j in cell.ifaces]) for cell in the_mesh.cells]
[cell.get_centroid([the_mesh.faces[j] for j in cell.ifaces]) for cell in the_mesh.cells]
#calcula los volúmenes (áreas) de las celdas
[cell.get_volume2([the_mesh.faces[j] for j in cell.ifaces]) for cell in the_mesh.cells]
# encuentra los vecinos a cada celda
[cell.get_neighbours([the_mesh.faces[j] for j in cell.ifaces]) for cell in the_mesh.cells]



""" inicializamos un campo de presiones lineal en x e y"""

init_pres=lambda x, y  :10-x
init_vel=lambda x, y  :np.array([4*10* (y/2-y**2/4), 0*x])

#prueba con funciones tipo lambda
pressure_field_t0=ds.field(the_Type="scalar", the_locale="cells", fun=init_pres, the_mesh=the_mesh, step=0)
velocity_field_t0=ds.field(the_Type="vector", the_locale="cells", fun=init_vel, the_mesh=the_mesh, step=0)

#prueba con un valor constante
pressure_field_t1=ds.field(the_Type="scalar", the_locale="cells", the_values=1, the_mesh=the_mesh, step=1)
velocity_field_t1=ds.field(the_Type="vector", the_locale="cells", the_values=1, the_mesh=the_mesh, step=1)

""" graficos """

plt.figure(1)
# geometria, normales, , indices y conectividades 
for face_i in the_mesh.faces:
#     # grafico de lineas entre los vertices que componen cada cara (para ver que esten todas las caras)
    plt.quiver(face_i.v1.x, face_i.v1.y,  face_i.v2.x -face_i.v1.x , face_i.v2.y-face_i.v1.y , width=5e-3) 
    plt.plot([face_i.v1.x ,face_i.v2.x],  [face_i.v1.y ,face_i.v2.y])
# for i in range(the_mesh.nvertices):
#     #grafico de posiciones de los vertices (con circulos) 
#     plt.plot( the_mesh.vertices[i].x, the_mesh.vertices[i].y,   'og'   )
#     plt.text(  the_mesh.vertices[i].x-0.1,  the_mesh.vertices[i].y+0.1,   str(i) , color="green" )
    
for vertex in the_mesh.vertices:
    #grafico de posiciones de los vertices (con circulos) 
    plt.plot( vertex.x, vertex.y,   'og'   )
    plt.text(  vertex.x-0.1,  vertex.y+0.1,   str(vertex.index) , color="green" )
# gráfico de normales y números de caras
for face in the_mesh.faces:
    plt.quiver(face.center[ 0],face.center[1],face.normal[0], face.normal[1] )
    plt.text(face.center[ 0],face.center[1]+.1, str(face.index),  color= 'black')

for cell in the_mesh.cells:
      plt.plot(cell.centroid[0], cell.centroid[1], 'xr')
      plt.text(cell.centroid[0], cell.centroid[1]+.1, str(cell.index),  color= 'red')
#plt.legend()
plt.axis("equal")
#plt.xlim([-0,3])
#plt.ylim([-0,2])

# #datos numéricos para verificación
volumes=[cell.volume for cell in the_mesh.cells]
neighbours=[cell.icells for cell in the_mesh.cells]

plt.figure(2)

#volumenes de celdas
plt.plot(range(the_mesh.ncells), volumes, 'o')


from mpl_toolkits import  mplot3d
plt.figure(3)

#ax=plt.axes(projection="3d")
#distribucion de presones en el instante inicial
x=np.array([cell.centroid[0] for cell in the_mesh.cells])
y=np.array([cell.centroid[1] for cell in the_mesh.cells])

X,Y= np.meshgrid(x,y)
z=pressure_field_t0.values[0:the_mesh.ncells]
u=velocity_field_t0.values[0:the_mesh.ncells,0]
v=velocity_field_t0.values[0:the_mesh.ncells, 1]

#ax.plot3D(x,y,z, 'xr')
f, ax = plt.subplots(1,2, sharex=True, sharey=True)
#ax[0].tripcolor(x,y,z)
ax[1].tricontourf(x,y,z, 20) # choose 20 contour levels, just to show how good its interpolation is
#ax[1].plot(x,y, 'ko ')
#ax[0].plot(x,y, 'ko ')
ax[0].quiver(x,y, u,v)

plt.xlim([-0,3])
plt.ylim([-0,2])
