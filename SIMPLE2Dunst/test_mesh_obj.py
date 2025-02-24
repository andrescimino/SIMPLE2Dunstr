
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 19 19:28:22 2024

@author: acimino
"""
import numpy as np
import matplotlib.pyplot as plt
import data_struc as ds


#declara los archivos de entrada con coordenadas y conectividades.
# en formato analogo a openfoam
vertex_file="points2.dat"
face_file="faces3.dat"
owners_file="owners3.dat"
neighbours_file="neighbours3.dat"
boundaries_file="boundaries3.dat"

# parseo de coordenadas de vertices
fpoints=open(vertex_file, 'r')
nvertices=int(fpoints.readline())
vertices=fpoints.readlines()
vertices.pop(0) #elimina la primera fila que contiene solo un parentesis
vertices.pop(-1) #elimina la ultima fila que contiene solo un parentesis
vertices=[vertex.split(" ") for vertex in vertices] #separa los puntos por espacio
vertices=[[num.strip("(|)\n") for num in vertex]  for vertex in vertices] #elimina los parentesis al inicio y final de cada fila
#castea a punto flotante
vertices= [[float(num) for num in vertex]  for vertex in vertices]

fpoints.close()

vertices_coord=np.array(vertices)#transforma la lista de vertices a un array de numpy

#### lista con vertices

vertex_list=[ds.vertex(coor_x=vertices_coord[i, 0],coor_y=vertices_coord[i, 1], index=i) for i in range(nvertices)]


""" caras"""

#abre el archivo que contiene las conectividades de las caras
ffaces=open(face_file, 'r')
nfaces=int(ffaces.readline())
faces_connect=ffaces.readlines()
faces_connect.pop(0) #elimina la primera fila que contiene solo un parentesis
faces_connect.pop(-1) #elimina la ultima fila que contiene solo un parentesis
 #elimina el primer 2 (ya que en toda malla 2D las caras tienen 2 puntos) y los parentesis, 
 #y separa por espacio
faces_connect=[face[1:-1].strip("(|)\n").split(" ") for face in faces_connect]

ffaces.close()
faces_connect= np.array([[int(num) for num in vertex]  for vertex in faces_connect])

#mesh_test=ds.mesh(vertex_file, face_file, owners_file, neighbors_file)

## asignación de owners 
fowners=open(owners_file, 'r')
#verifica que la cantidad de owners sea igual a la cant de faces
if int(fowners.readline())!=nfaces:
    print("error en los archivos de entrada")

owners=fowners.readlines()
owners.pop(0) #elimina la primera fila que contiene solo un parentesis
owners.pop(-1) #elimina la ultima fila que contiene solo un parentesis
owners=[int(i) for i in owners]
fowners.close()

#el numero de celdas es el indice mas grande del archivo owners mas 1 (ya que cuenta desde cero)
ncells=np.max(owners)+1

fneighbours=open(neighbours_file, 'r')

nneighbours=int(fneighbours.readline())
neighbours=fneighbours.readlines()
neighbours.pop(0) #elimina la primera fila que contiene solo un parentesis
neighbours.pop(-1) #elimina la ultima fila que contiene solo un parentesis
neighbours=[int(i) for i in neighbours]
n_boundary_faces=nfaces-nneighbours

fneighbours.close()



#armado de lista de caras

#arma una lista con instancias de face para cada cara. inicializa a todos los vecinos con -1 (como si fueran todas caras del contorno)
faces_list=[ds.face(v1=vertex_list[faces_connect[i,0]], v2=vertex_list[faces_connect[i,1]], index=i, owner=owners[i], neighbour=-1) for i in range(nfaces)]

# agrega las celdas vecinas asignadas a cada cara (si las tienen)
for i in range( nfaces-n_boundary_faces):
    #faces_list[i].neighbour=neighbours[i-n_boundary_faces]
    faces_list[i].neighbour=neighbours[i]

#calcula las normales a cada una de las caras
[face.get_normal() for  face in faces_list]
[face.get_center() for  face in faces_list]
[face.get_area() for  face in faces_list]



"""  celdas  """
#### armado de conectividades de celdas

cell_faces_connect=[[] for i in range(ncells)]
for i in range(nfaces):
    cell_faces_connect[owners[i]].append(i)
    
for i in range(nneighbours):
    cell_faces_connect[neighbours[i]].append(i)
    
    
    
cell_list=[ ds.cell(connect_faces=cell_faces_connect[i][:], index=i) for i in range(ncells)  ]
[cell_i.get_faces_signs([faces_list[j] for j in cell_i.ifaces]) for cell_i in cell_list]
#calcula centroides de celdas
[cell_i.get_centroid([faces_list[j] for j in cell_i.ifaces]) for cell_i in cell_list]
#calcula los volúmenes (áreas) de las celdas
[cell_i.get_volume2([faces_list[j] for j in cell_i.ifaces]) for cell_i in cell_list]
# encuentra los vecinos a cada celda
[cell_i.get_neighbours([faces_list[j] for j in cell_i.ifaces]) for cell_i in cell_list]

#factores de interpolación de las caras


[face.get_interp_factor(owner=cell_list[face.owner], neighbour=cell_list[face.neighbour]) for  face in faces_list]

""" objeto de malla"""

mesh_object=ds.mesh(vertex_file, face_file, owners_file, neighbours_file, boundaries_file)


""" graficos """

plt.figure(1)

for face_i in faces_list:
#     # grafico de lineas entre los vertices que componen cada cara (para ver que esten todas las caras)
    plt.quiver(face_i.v1.x, face_i.v1.y,  face_i.v2.x -face_i.v1.x , face_i.v2.y-face_i.v1.y , width=5e-3) 
    plt.plot([face_i.v1.x ,face_i.v2.x],  [face_i.v1.y ,face_i.v2.y])
for i in range(nvertices):
    #grafico de posiciones de los vertices (con circulos) 
    plt.plot( vertex_list[i].x, vertex_list[i].y,   'og'   )
    plt.text( vertex_list[i].x-0.1, vertex_list[i].y+0.1,   str(i) , color="green" )
# gráfico de normales y números de caras
for i in range(nfaces):
    plt.quiver(faces_list[i].center[ 0],faces_list[i].center[1],faces_list[i].normal[0], faces_list[i].normal[1] )
    plt.text(faces_list[i].center[ 0],faces_list[i].center[1]+.1, str(i),  color= 'black')

for i in range(ncells):
     plt.plot(cell_list[i].centroid[0], cell_list[i].centroid[1], 'xr')
     plt.text(cell_list[i].centroid[0], cell_list[i].centroid[1]+.1, str(i),  color= 'red')
#plt.legend()
plt.axis("equal")
#plt.xlim([-0,3])
#plt.ylim([-0,2])

#datos numéricos para verificación
volumes=[cell_i.volume for cell_i in cell_list]
neighbours=[cell_i.icells for cell_i in cell_list]

plt.figure(2)
plt.plot(range(ncells), volumes, 'o')