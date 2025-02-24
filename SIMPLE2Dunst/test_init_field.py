#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 20:19:09 2024

@author: sebastian
"""
"""inicializacion de la malla
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



the_Type="pressure"
the_locale="cells"
the_file="100/p.dat"


pressure_field_t0=ds.field(the_Type=the_Type, the_locale=the_locale,step=0, the_mesh=the_mesh, the_file=the_file)

the_Type="pressure"
the_locale="boundary"
the_file="100/p.dat"

pressure_field_bound=[[] for i in range(the_mesh.n_boundaries)]

# for b_values in pressure_field_bound:
#     b_values = ds.field(the_Type=the_Type, the_locale="boundary",step=0, the_mesh=the_mesh, the_file=the_file, bound=the_mesh.boundaries[pressure_field_bound.index(b_values)])

for i in range(the_mesh.n_boundaries):
    pressure_field_bound[i] = ds.field(the_Type=the_Type, the_locale="boundary",step=0, the_mesh=the_mesh, the_file=the_file, bound=the_mesh.boundaries[i])



"""testeo fuera de estructura de datos"""
match the_Type:
    case "pressure":
        Type="pressure"
        
        #parseo de archivo con datos de presion
        
        ffield=open(the_file, 'r')
        field_data=ffield.readlines()
        ffield.close()
        field_data.pop(0)#elimina la primera fila
        field_data.pop(-1)#elimina la ultima fila
        field_data=[line.strip("\n|;|,|[|]").split(" ") for line in field_data] #elimina saltos de lineas y tabulaciones


        field_data=[ [item for item in line if len(item)>0] for line in field_data  ]# eliminacaracteres vacios dentro de las lineas (caracteres de long 0)
        field_data=[ line for line in field_data if line !=[] ]# elimina filas vacias (caracteres de long 0)
        field_data=[ item for item in field_data if item !=['{']  ] 
        field_data=[ item for item in field_data if item !=['}']  ] 
        field_data=[ item for item in field_data if item !=[')']  ] 
        field_data=[ item for item in field_data if item !=['(']  ] 

        match the_locale:
            case "cells":
                #asigna valores en centroides de celdas
                #self.nvalues=the_mesh.ncells+the_mesh.n_boundary_faces
                # nvalues=the_mesh.ncells
                # values=np.zeros(nvalues)
               
                pressure_0_int=np.zeros([the_mesh.ncells])
                
                for line in field_data:
                    for item in line:
                        if item=="internalField":
                        
                            domain_dist_type= field_data[field_data.index(line)][line.index(item)+1]
                            domain_dist_indx= [field_data.index(line),line.index(item)+1]
               
                if domain_dist_type=="uniform":
                    pressure_0_int[:]=float(field_data[domain_dist_indx[0]][domain_dist_indx[1]+1])
                elif domain_dist_type=="nonuniform":
                    
                    num_values=int(field_data[domain_dist_indx[0]+1][0])
                    for i in range(num_values):  
                        pressure_0_int[i]=float(field_data[domain_dist_indx[0]+2+i][0])
                #asigna valores en caras del contorno
             
            case "boundaries":
                 #asigna valores en centroides de celdas
                 #self.nvalues=the_mesh.ncells+the_mesh.n_boundary_faces
                 # nvalues=the_mesh.n_boundary_faces
                 # values=np.zeros(nvalues)
                 #pass
                 
                 for line in field_data:
                     for item in line:
                         # if item=="boundaryField":
                        
                         #     boundaries_init_indx= field_data.index(line)
                             
                         boundary_field_name=np.zeros(the_mesh.n_boundaries)
                         for i in range(the_mesh.n_boundaries):               
                             if the_mesh.boundaries[i].name==line[0]:
                                 boundary_field_name[i]=line[0]
                                 boundary_field_values=np.zeros([the_mesh.boundaries[i].nfaces])
                   
                 #asigna valores en caras del contorno
            case "faces":
                # self.nvalues=the_mesh.nfaces
                # self.values=np.zeros(self.nvalues)
                pass
            case "boundaries":
                pass
            case "vertices":
                pass
            case _:
                print("unknown locale line 134") 
            
    case "velocity":
        Type="velocity"
        match the_locale:
            case "cells":
                nvalues=the_mesh.ncells+the_mesh.n_boundary_faces
                values=np.zeros((nvalues, 2)) #asigna el espacio en la memoria
            case "faces":
                pass
            case "boundaries":
                pass
            case "vertices":
                pass
            case _:
                print("unknown locale")
    case _:
        print("unknown field type")