#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 14 16:42:14 2024

@author: andres M Cimino

estructura de datos para una malla 2D para un codigo de volumenes 
finitos para mallas no estructuradas

"""
import numpy as np
"""puntos generados de la malla"""
class vertex(object):
    def __init__ (self, coor_x,coor_y, index):
        self.index=index
        self.x=coor_x
        self.y=coor_y
""" caras de interfaces, formadas por dos vertices"""        
class face(object):
    def __init__ (self, v1: vertex, v2 : vertex, owner, neighbour, index):
        self.v1=v1 #primer vertice
        self.v2=v2 #segundo vertice
        self.index=index
        self.owner=owner #celda que posee la cara
        self.neighbour=neighbour
        
    def get_normal(self):
        """calcula los vectores normales (unitarios) a cada cara. Como se trata de un caso 2D
        la normal es simplemente el vector entre v2 y v1 rotado 90 grados y normalizado
        los guarda como atributo del objeto face que se  pasa como argumento"""
        l_face=np.sqrt((self.v2.x-self.v1.x)**2+ (self.v2.y-self.v1.y)**2)
        n1_x=(self.v2.y-self.v1.y)/l_face; n1_y=-(self.v2.x-self.v1.x)/l_face 
        
        # devuelve las normales como array de numpy
        self.normal= np.array([n1_x, n1_y ])
    
    def get_normal2(self):
        """calcula los vectores normales (unitarios) a cada cara. Como se trata de un caso 2D
        la normal es simplemente el vector entre v2 y v1 rotado 90 grados y normalizado
        los guarda como atributo del objeto face que se  pasa como argumento
        """
        l_face=np.sqrt((self.v2.x-self.v1.x)**2+ (self.v2.y-self.v1.y)**2)
        n1_x=(self.v2.y-self.v1.y)/l_face; n2_x=-(self.v2.x-self.v1.x)/l_face 
        #devuelve las componentes de la normal como escalares. 
        self.n_x=n1_x;   self.n_y=n2_x; 
        
    
    def get_center(self):
        """calcula las coordenadas del centro de la cara. El promedio aritmetico de las coords
        de los vertices que lo componen. los guarda como atributo del objeto face que se  pasa como argumento"""
        self.center=np.array([(self.v2.x+self.v1.x)/2, (self.v2.y+self.v1.y)/2  ])
    
    def get_area(self):
        """calcula el area de la cara: long del segmento por ancho unitario"""
        self.area= np.sqrt((self.v2.x-self.v1.x)**2+ (self.v2.y-self.v1.y)**2 )
        
    def get_interp_factor(self, owner, neighbour):
        """ calcula el factor de interpolacion g_f haciendo productos escalares entre el vector que une el centroide 
        de la celda owner y el centroide de la cara con la normal a la cara y entre el centroide de la celda vecina
        
        gf=(dCf. n_f)/(dCf. n_f+dfF. n_f)))
        donde dCf= vector desde el centroide de C al centro de la cara f
        donde dfF= vector desde el centro de la cara f al centroide de la celda vecina F
        """
        if self.neighbour==-1: # caso de cara ubicada en el contorno
            g_f=1
        else:
            dCf=self.center-owner.centroid
            dfF=neighbour.centroid-self.center
        
            g_f=np.dot(dCf, self.normal)/(np.dot(dCf, self.normal)+np.dot(dfF, self.normal))
        self.g_f=g_f
        
    def calc_face_value(self, owner_value, neigh_value):
        """ cacula el valor de la variable en el centro de la cara interpolando 
        con valores en las celdas owner y neighbour
        devuelve un valor de salida face_value"""
        face_value=self.g_f*neigh_value+ (1-self.g_f)*owner_value
        return face_value
        
""" 
elementos o celdas
""" 

class cell(object):
    def __init__ (self, connect_faces,  index):
        """ inicializa el objeto, requiere un arreglo con las conectividades ( indices) de las caras que 
        encierran la celda y el indice (un entero) que se le asignara a la misma)"""
        self.index=index
        self.numfaces=len(connect_faces)
        #self.numneigh=len(connect_cells)
        self.ifaces=np.array(connect_faces)
        #self.icells=np.array(connect_cells)
        # self.connect_faces=connect_faces
        # self.connect_cells=connect_cells
        
    def get_faces_signs(self, faces):
        """ asigna signo positivo (normal hacia afuera) si la cara es propiedad de la celda considerada
            y signo negativo si pertenece a un vecino
            se pasa como argumento de entrada una lista con instancias de la clase face
            asigna un arreglo con 1 o -1 segun el signo correspondiente a cada cara
        """
        faces_sign=np.ones(len(faces))
        for i in range (len(faces)):
            if faces[i].owner != self.index :
                faces_sign[i]=-1
        self.faces_sign=faces_sign
        
    def get_centroid(self, faces):
        """ calcula el centroide de la celda,
        si la celda es triangular, es el promedio de las coordenadas de los vertices (centro geom)
        Si la celda no es triangular, se divide en tantos triangs como lados tenga. 
        el tercer vértice es el promedio de las coordenadas de los vértices
        el centroide es el promedio ponderado de los centroides de los triángulos (con sus áreas como factor de ponderación)
        asigna el valor como float como atributo del objeto
        """
        xg=0; yg=0;
        nfaces=self.numfaces
        ind_cell=self.index
        if nfaces<4:
            # si la celda es triangular el centroide del triangulo coincide con el promedio de las coordenadas de los 
            # vertices
            
            for i in range(nfaces):
                if faces[i].owner ==ind_cell:
                    xg=xg+faces[i].v1.x
                    yg=yg+faces[i].v1.y
                else:
                    xg=xg+faces[i].v2.x
                    yg=yg+faces[i].v2.y
            xg=xg/3; yg=yg/3;
        else:
            x_geom=0; y_geom=0;
            for i in range(nfaces):
                # calcula un centro geometrico como promedio aritmetico de las coords de los 
                #vertices
                if faces[i].owner ==ind_cell:
                    x_geom=x_geom+faces[i].v1.x
                    y_geom=y_geom+faces[i].v1.y
                else:
                    x_geom=x_geom+faces[i].v2.x
                    y_geom=y_geom+faces[i].v2.y
            x_geom=x_geom/nfaces; y_geom=y_geom/nfaces;
            #con este centro geom arma tantos triangulos como caras tenga la celda, con los 
            # vertices de  la cara y el centro geometrico
            for i in range(nfaces):
                # suma el centroide de cada uno de los i triangulos internos 
                xg=xg+(faces[i].v1.x+faces[i].v2.x+x_geom)/3
                yg=yg+(faces[i].v1.y+faces[i].v2.y+y_geom)/3
            xg=xg/nfaces; yg=yg/nfaces;
        self.centroid=np.array([xg, yg])
        
    def get_volume(self, faces):
        """ calcula el volumen de la celda. En este caso es igual al área de la misma por un
        ancho unitario. SI la celda tiene tres caras se calcula con la expr del producto
        vectorial para el área de un triángulo. Si tiene más de 3 caras se descompone en 
        n triángulos
        """
        nfaces=self.numfaces
        ind_cell=self.index
        if nfaces<4:
            # si la celda es triangular 
            if self.faces_sign[0] ==self.faces_sign[1]:
                x1=faces[0].v1.x
                y1=faces[0].v1.y
                x2=faces[0].v2.x
                y2=faces[0].v2.y
                x3=faces[1].v2.x
                y3=faces[1].v2.y
            else:
                x1=faces[0].v2.x
                y1=faces[0].v2.y
                x2=faces[0].v1.x
                y2=faces[0].v1.y
                x3=faces[1].v2.x
                y3=faces[1].v2.y

            self.volume=0.5*(x1*( y2-y3)+x2*( y3- y1)+x3*( y1-y2))
                
        else:
            x_geom=0; y_geom=0;
            Vol=0
            for i in range(nfaces):
                # calcula un centro geometrico como promedio aritmetico de las coords de los 
                #vertices
                if faces[i].owner ==ind_cell:
                    x_geom=x_geom+faces[i].v1.x
                    y_geom=y_geom+faces[i].v1.y
                else:
                    x_geom=x_geom+faces[i].v2.x
                    y_geom=y_geom+faces[i].v2.y
            x_geom=x_geom/nfaces; y_geom=y_geom/nfaces;
            #con este centro geom arma tantos triangulos como caras tenga la celda, con los 
            # vertices de  la cara y el centro geometrico
                
            for i in range(nfaces):
                #
                # if self.faces_sign[i] ==1:
                #     x1=faces[i].v1.x
                #     y1=faces[i].v1.y
                #     x2=faces[i].v2.x
                #     y2=faces[i].v2.y
                #     x3=x_geom
                #     y3=y_geom
                #     Vol=Vol+0.5*(x1*( y2-y3)+x2*( y3- y1)+x3*( y1-y2))
                # else:
                #     x1=faces[i].v2.x
                #     y1=faces[i].v2.y
                #     x2=faces[i].v1.x
                #     y2=faces[i].v1.y
                #     x3=x_geom
                #     y3=y_geom
                #     Vol=Vol+0.5*(x1*( y2-y3)+x2*( y3- y1)+x3*( y1-y2))
                x1=faces[i].v1.x
                y1=faces[i].v1.y
                x2=faces[i].v2.x
                y2=faces[i].v2.y
                x3=x_geom
                y3=y_geom
                Vol=Vol+np.abs(0.5*(x1*( y2-y3)+x2*( y3- y1)+x3*( y1-y2)))
                    
            self.volume=Vol
            
    def get_volume2(self, faces):
        """ calcula el volumen de la celda, asumiendo que se conoce el centroide de la misma
        Se calcula como la sumatoria de las areas de los n triángulos para evitar usar
        condiciones lógicas para cada cara y su orientacón
        """
        nfaces=self.numfaces
        x_geom, y_geom=self.centroid[0], self.centroid[1]
        Vol=0
        for i in range(nfaces):

            x1=faces[i].v1.x
            y1=faces[i].v1.y
            x2=faces[i].v2.x
            y2=faces[i].v2.y
            x3=x_geom
            y3=y_geom
            Vol=Vol+np.abs(0.5*(x1*( y2-y3)+x2*( y3- y1)+x3*( y1-y2)))
        self.volume=Vol
    
    def get_neighbours(self, faces):
        """ arma las conectividades de las celdas en base a la lista de vecinos de las caras"""
        
        nfaces=len(faces)
        ineighbours=[]        
        for i in range(nfaces):
            if self.faces_sign[i]==-1:
                ineighbours.append(faces[i].owner)
            else:
                ineighbours.append(faces[i].neighbour)
        self.icells=np.array(ineighbours)
        
    def calc_gradient(self, faces, face_value):
        """calcula el gradiente en el centro de la celda usando el método de green gauss"""
        nfaces=len(faces)
        Gradient=np.zeros(2)
        for i in range(nfaces):
            Gradient=Gradient+face_value[i]*faces[i].normal*faces[i].area
        Gradient=Gradient/self.volume
        return Gradient


class boundary(object):
    
    def __init__(self, Type: str,  ifaces: np.array, name:str, index):
        """se inicializa con el  tipo de condicion de contorno y con una lista de caras que la componen
        tipos:
            vel_inlet
            p_inlet
            p_outlet
            wall
        
        """
        self.Type=Type
        self.name=name
        self.nfaces=len(ifaces)
        self.ifaces=ifaces
        self.index=index
        

    def set_value(self, the_field, the_value):
        """establece valores en los centros de las caras dependiendo del tipo de condicion de contorno que se trate"""   
        match the_field.Type:
            case "scalar":
                match self.Type:
                    case "vel_inlet":
                        pass                        
                    case "p_inlet":
                        field.values=the_value
                    case "p_outlet":
                        field.values=the_value
                    case "wall":
                        pass 
                    case _:
                        return("error in BC ")
                        
                    
            case "vector":
                match self.Type:
                    case "vel_inlet":
                        pass                        
                    case "p_inlet":
                        field.values=the_value
                    case "p_outlet":
                        field.values=the_value
                    case "wall":
                        pass
                    case _:
                        return("error in BC ")
        
""" malla
 incorpora elementos de las clases vertex, face y cell como partes constituyentes
 """       

class mesh(object):
    def __init__(self, vertex_file, faces_file, owners_file, neighbours_file, boundaries_file):
        #abre el archivo de datos de vértices en modo solo lectura
        #usa formato de openfoam, pero para una malla 2D
        fpoints=open(vertex_file, 'r')#handle para archivo de datos
        self.nvertices=int(fpoints.readline()) #la primera linea del archivo tiene el num de vertices
        #determinacion de coordenadas del vertices
        vertices=fpoints.readlines()
        vertices.pop(0) #elimina la primera fila que contiene solo un parentesis
        vertices.pop(-1) #elimina la ultima fila que contiene solo un parentesis
        vertices=[vertex.split(" ") for vertex in vertices] #separa los puntos por espacio
        vertices=[[num.strip("(|)\n") for num in vertex]  for vertex in vertices] #elimina los parentesis al inicio y final de cada fila
        vertices= np.array([[float(num) for num in vertex]  for vertex in vertices])
        
        fpoints.close()
        self.vertices=[ vertex(coor_x=vertices[i,0], coor_y=vertices[i,1], index=i) for i in range(self.nvertices)]
                                  
        # y lo agregamos al objeto malla
        
        
        ffaces=open(faces_file, 'r')
        self.nfaces=int(ffaces.readline())
        faces_connect=ffaces.readlines()
        faces_connect.pop(0) #elimina la primera fila que contiene solo un parentesis
        faces_connect.pop(-1) #elimina la ultima fila que contiene solo un parentesis
         #elimina el primer 2 (ya que en toda malla 2D las caras tienen 2 puntos) y los parentesis, 
         #y separa por espacio
        faces_connect=[face[1:-1].strip("(|)\n").split(" ") for face in faces_connect]

        ffaces.close()
        faces_connect=np.array([[int(num) for num in vertex]  for vertex in faces_connect])
        #arma una lista para alojar los objetos face con los datos de las caras
        #self.faces=[i for i in range(self.nfaces)]
        
        
        
        
        ### asignación de owners 
        fowners=open(owners_file, 'r')
        #verifica que la cantidad de owners sea igual a la cant de faces
        if int(fowners.readline())!=self.nfaces:
            print("error en los archivos de entrada")
        
        owners=fowners.readlines()
        owners.pop(0) #elimina la primera fila que contiene solo un parentesis
        owners.pop(-1) #elimina la ultima fila que contiene solo un parentesis
        owners=[int(i) for i in owners]
        
        fowners.close()
        # obtiene el número de elementos del archivo de owners
        self.ncells=np.max(owners)+1
        
        # asignación de neighbours a las caras
        fneighbours=open(neighbours_file, 'r')

        nneighbours=int(fneighbours.readline())
        neighbours=fneighbours.readlines()
        neighbours.pop(0) #elimina la primera fila que contiene solo un parentesis
        neighbours.pop(-1) #elimina la ultima fila que contiene solo un parentesis
        neighbours=[int(i) for i in neighbours]
        self.n_boundary_faces=self.nfaces-nneighbours

        fneighbours.close()
        
        #lectura de archivo de condiciones de contorno
        fboundaries=open(boundaries_file, 'r')
        boundaries=fboundaries.readlines()

        fboundaries.close()

        boundaries=[data.strip(" |\n").split("\t") for data in boundaries]

        boundaries=[ [i for i in item if len(i)!=0] for item in boundaries ]
        boundaries=[ item for item in boundaries if item !=['{']  ] 
        boundaries=[ item for item in boundaries if item !=['}']  ] 



        self.n_boundaries=len(boundaries)//4
        bound_list=[]
        
        #armado de arreglo con boundaries

        for i in range(self.n_boundaries):
            name=boundaries[i*4][0]
            startface=int(boundaries[3+i*4][1])
            nfaces=int(boundaries[2+i*4][1])
            ifaces=np.arange(startface, startface+nfaces)
            Type=boundaries[1+4*i][1]
            bound_list.append(boundary(Type, ifaces, name, i))
            
        self.boundaries=bound_list
    
        
        """ armado de arreglo con faces
        
        """
        # for i in range(self.nfaces):
        #     #self.face.append(face(v1=faces_connect[0], v2=faces_connect[1], owner=1,neighbour=-1, index=i))
        #     self.faces[i]=face(v1=faces_connect[0], v2=faces_connect[1], owner=1,neighbour=-1, index=i)
        self.faces=[face(v1=self.vertices[faces_connect[i,0]], v2=self.vertices[faces_connect[i,1]], owner=owners[i],neighbour=-1, index=i) for i in range(self.nfaces)]
        for i in range( self.nfaces-self.n_boundary_faces):
            #faces_list[i].neighbour=neighbours[i-n_boundary_faces]
            self.faces[i].neighbour=neighbours[i]
            
        """ armado de arreglo con celdas  """
        #### armado de conectividades de celdas

        cell_faces_connect=[[] for i in range(self.ncells)]
        for i in range(self.nfaces):
            cell_faces_connect[owners[i]].append(i)
            
        for i in range(nneighbours):
            cell_faces_connect[neighbours[i]].append(i)

        self.cells=[ cell(connect_faces=cell_faces_connect[i][:], index=i) for i in range(self.ncells)  ]
        
        
        
class field(object):
    def __init__(self, the_Type: str, the_locale: str, step,  the_mesh: mesh, fun=None, the_values=None, the_file=None, bound=None):
        """ inicializa un campo escalar o vectorial en puntos conocidos de la malla
        - centroides de celdas, si locale= cells
        -centros de caras, si locale=faces
        -nodos o vertices si locale= vertices
        
        si es un campo escalar, asigna un solo valor a cada punto. almacena los datos como un arreglo de numpy
        de una dimension con n elementos (n es el numero de celdas, caras, etc, segun la locale).
        Si es un campo vectorial almacena en un arreglo de 2 *n de numpy
        
        puede asignarse los valores (por ej asignando de una iteracion anterior) o una funcion anonima 
        hay que asignar una opcion como si fuera un kwarg
        
        fun es una funcion anonima cualquiera de x e y para asignar el valor de la variable en las coordenadas de los puntos 
        
        the_file es un archivo de distribucion de presiones o velocidades de openfoam  en un instante de tiempo dado (por ej: p.dat en la carpeta 0)
        bound es un objeto de clase boundary, donde esta el nombre, numero de caras y conectividades de las caras que la componen. Por default esta como none para no interferir con
        la inicializacion del dominio
        
        """
        
        
        match the_Type:
            case "pressure":
                self.Type="pressure"
                
                #parseo de archivo con datos de presion
                if the_file != None:
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
                        self.nvalues=the_mesh.ncells
                        self.values=np.zeros(self.nvalues)
                        self.locale=the_locale
                        if fun != None:
                            #asignacion de valores pasando una funcion de x, y como argumento
                            self.values[0:the_mesh.ncells]=np.array([fun(cell.centroid[0], cell.centroid[1]) for cell in the_mesh.cells])
                        elif the_file != None:
                            # lectura de datos desde un archivo (como en openfoam)

                            for line in field_data:
                                for item in line:
                                    if item=="internalField":
        
                                        self.dist_type= field_data[field_data.index(line)][line.index(item)+1]
                                        domain_dist_indx= [field_data.index(line),line.index(item)+1]

                            if self.dist_type=="uniform":
                                self.values[:]=float(field_data[domain_dist_indx[0]][domain_dist_indx[1]+1])
                            elif self.dist_type=="nonuniform":
    
                                num_values=int(field_data[domain_dist_indx[0]+1][0])
                                for i in range(num_values):  
                                    self.values[i]=float(field_data[domain_dist_indx[0]+2+i][0])
                        else:
                            self.values[0:the_mesh.ncells]=the_values
                        #asigna valores en caras del contorno
                     
                    case "boundary":
                         #asigna valores en centroides de celdas
                         #self.nvalues=the_mesh.ncells+the_mesh.n_boundary_faces
                         self.locale=bound.name
                         self.nvalues=bound.nfaces
                         self.values=np.zeros(self.nvalues)
                         

                             
                         
                         if fun != None:
                             #asignacion de valores pasando una funcion de x, y como argumento
                             #self.values[0:self.nvalues]=np.array([fun(face.center[0], face.center[1]) for face in the_mesh.faces])
                             pass
                         elif the_file != None:
                             
                             
                             for line in field_data:
                                 for item in line:
                                      if item==self.locale:
                                          bound_indx=field_data.index(line)
                                          self.bound_type=field_data[bound_indx+1][1]
                                          if field_data[bound_indx+2][0]=="value":
                                              self.dist_type=field_data[bound_indx+2][1]
                                              if self.dist_type=="uniform":
                                                  self.values[:] =float(field_data[bound_indx+2][2])
                                              elif self.dist_type=="nonuniform":
                  
                                                  num_values=int(field_data[bound_indx+3][0])
                                                  for i in range(num_values):  
                                                      self.values[i]=float(field_data[bound_indx+4+i][0])
                                         #     boundaries_init_indx= field_data.index(line)
                                     
                                        
                         else:
                             
                             self.values[0:the_mesh.ncells]=the_values
                         #asigna valores en caras del contorno
                    case "faces":
                        self.nvalues=the_mesh.nfaces
                        self.values=np.zeros(self.nvalues)
                    case "vertices":
                        pass
                        
                    
            case "velocity":
                self.Type="velocity"
                match the_locale:
                    case "cells":
                        self.nvalues=the_mesh.ncells+the_mesh.n_boundary_faces
                        self.values=np.zeros((self.nvalues, 2)) #asigna el espacio en la memoria
                        if fun != None:
                            self.values[0:the_mesh.ncells]=np.array([fun(cell.centroid[0], cell.centroid[1]) for cell in the_mesh.cells])
                        else:
                            self.values[0:the_mesh.ncells]=the_values
                    case "faces":
                        pass
                    case "vertices":
                        pass
                    case _:
                        return("unknown locale")
            case _:
                return("unknown field type")
                