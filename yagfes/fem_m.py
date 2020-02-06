#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: hcorzopola

    YAGFES: Yet Another Groundwater Flow Equation Solver
    Copyright (C) 2019  H.A. Corzo Pola

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
    
#########################FINITE ELEMENT METHOD MODULE##########################

This module includes all object classes and functions required to use the
Finite Element Method (FEM) for Groundwater (GW) modelling.

"""
###############################################################################
"""
##########################LIBRARY IMPORTATION SECTION##########################
"""
from yagfes import aux_f #Import 'aux_f' module. This module contains a bunch of auxiliar functions.
import numpy as np #Import 'numpy' for numerical operations

"""
##########################LIBRARY IMPORTATION SECTION##########################
"""
###############################################################################
"""
##############################MESH CLASS SECTION###############################
"""
class mesh:
    """
    n_n = # of nodes            | n_e = # of elements
    m_p = Point Matrix          | m_c = Connectivity Matrix
    a_e = Element's Area        | bc = Boundary Conditions 
    """ 
    def __init__(self,file):
        #READ MESH FILE########################################################
        __f=open(aux_f.path+'/'+file) #Open the given file in 'read' mode

        ##READ NODES INFO######################################################
        for i in range(4):
            __f.readline() #Skip 4 lines
        ###Read number of nodes################################################
        self.n_n=int(__f.readline().rstrip('\n')) #Read number of nodes.
        #'rstrip('\n')' removes the 'new line' command at the end of the line.
        #'int()' converts the string to a number.
        #######################################################################
        ###Create an array ('Point Matrix') for storing the node coordinates###
        self.m_p=[0]*self.n_n
        for i in range(self.n_n):
            self.m_p[i]=[0]*3
        #The 'Point Matrix' has dimensions of 'n_n x 3'. This allow us to store
        # up to three coordinates per node.
        #######################################################################
        ###Retrieve node coordinates from the mesh file########################
        for i in range(self.n_n):
            __aux_n=__f.readline().split() #Store the line in an auxiliar array
            for j in range(3):
                self.m_p[i][j]=float(__aux_n[j+1]) #Store coordinates in the 'Point Matrix'
        #######################################################################
        ##READ NODES INFO######################################################
        
        ##READ ELEMENTS INFO###################################################
        for i in range(2):
            __f.readline() #Skip 2 lines
        ###Read number of elements#############################################
        self.n_e=int(__f.readline().rstrip('\n')) #Read number of elements.
        #'rstrip('\n')' removes the 'new line' command at the end of the line.
        #'int()' converts the string to a number.
        #This number accounts for two-node elements. We will update this later,
        # since the final number of elements might not be the same.
        #######################################################################
        ###Create an array ('Connectivity Matrix') for storing the nodes per element
        self.m_c=list()
        #The 'Connectivity Matrix' has dimensions of 'n_e x 3', since we are
        # working with triangular elements. However, we won't specify its
        # dimensions in the code, since we are uncertain of 'n_e' final value.
        #######################################################################
        ###Create an array for storing Boundary Conditions (BC's)##############
        self.bc=list()
        #This array will store:
        #1) type of BC (0=Dirichlet, 1=Neumann),
        #2) boundary ID (since a boundary can be made of multiple lines) and
        #3) nodes composing the 'line' on the boundary.
        #######################################################################
        ###Retrieve elements from the mesh file################################
        for i in range(self.n_e):
            __aux_e=__f.readline().split() #Store the line in an auxiliar array
            if __aux_e[1]=='1':
                #Store in 'bc'
                self.bc.append( list( map(int,__aux_e[3:]) ) )
                #'map()' applies the 'int()' function to the whole array 
                # '__aux_e[3:]'. 'list()' just formats the 'map()' output to
                # a list.
            else:
                #Store in 'm_c'
                self.m_c.append( list( map(int,__aux_e[5:]) ) )
                #'map()' applies the 'int()' function to the whole array 
                # '__aux_e[5:]'. 'list()' just formats the 'map()' output to
                # a list
        ###Fix node index######################################################
        for i in range(len(self.m_c)):
            for j in range(3):
                self.m_c[i][j]=self.m_c[i][j]-1
        #######################################################################
        ###Update 'n_e' value##################################################
        self.n_e=len(self.m_c)
        #######################################################################
        ##READ ELEMENTS INFO###################################################
        __f.close() #Close the file
        #READ MESH FILE########################################################
        
        #CALCULATE ELEMENT AREA################################################
        self.a_e=[0]*self.n_e #The array has a size equal to 'n_e'.
        for i in range(self.n_e):
            a=[self.m_p[self.m_c[i][1]][0]-self.m_p[self.m_c[i][0]][0],\
               self.m_p[self.m_c[i][1]][1]-self.m_p[self.m_c[i][0]][1],\
               self.m_p[self.m_c[i][1]][2]-self.m_p[self.m_c[i][0]][2]]
            b=[self.m_p[self.m_c[i][2]][0]-self.m_p[self.m_c[i][0]][0],\
               self.m_p[self.m_c[i][2]][1]-self.m_p[self.m_c[i][0]][1],\
               self.m_p[self.m_c[i][2]][2]-self.m_p[self.m_c[i][0]][2]]
            self.a_e[i]=0.5*((a[0]*b[1]-a[1]*b[0])**2+\
                    (a[0]*b[2]-a[2]*b[0])**2+\
                    (a[1]*b[2]-a[2]*b[1])**2)**0.5
        #CALCULATE ELEMENT AREA################################################
        
    def write_bc_file(self,file):
        #WRITE A BC FILE BASED ON MESH INFO####################################
        __f=open(aux_f.path+'/'+file,'w') #Open the given file in 'write' mode
        __nbc=len(self.bc)
        __f.write('#BC FILE\n') #Write document title

        ##DIRICHLET BC's#######################################################
        __f.write('#DIRICHLET BC\n') #Write section title
        ###Identify lines######################################################
        ####Create auxiliar arrays for storing line ID's and values############
        __lineID=list()
        __lineV=list()
        #######################################################################
        ###Loop for identifying lines and assigning them flow values###########
        for i in range(__nbc):
            #Check BC flag; 0=Dirichlet BC
            if self.bc[i][0]==0:
                #Check if 'Line ID' has already been stored
                if self.bc[i][1] not in __lineID:
                    __lineID.append(self.bc[i][1])
                    __type='0'
                    while __type!='1' and __type!='2':
                        __type=input('You have Dirichlet BC at line '+\
                                     str(self.bc[i][1])+\
                                     '. Please input the type of prescribed head line that you prefer: '+\
                                     '\n 1) Constant head\n 2) Gradient\n')
                        if __type!='1' and __type!='2':
                            print("That is not a valid option. Please, try again.")
                    if __type=='1': #Constant head
                        __v=input('Please, input the Dirichlet BC at line '+\
                              str(self.bc[i][1])+'.\n') #Read Dirichlet BC value
                        __lineV.append([int(__type),float(__v)])
                    else: #Gradient or __type==2
                        __v=input('Please, input the Dirichlet BC at the initial node (as defined on mesh file).\n') #Read Dirichlet BC value at initial node
                        __v2=input('Please, input the Dirichlet BC at the last node (as defined on mesh file).\n') #Read Dirichlet BC value at last node
                        __lineV.append([int(__type),float(__v),float(__v2)])
        ###Loop for identifying lines and assigning them flow values###########
        #######################################################################
        ##Assign head value to lines in BC file################################
        ###Create an auxiliar array for storing the Dirichlet nodes############
        __aux_d=list()
        __aux_dv=list()
        for i in range(len(__lineID)):
            if __lineV[i][0]==1: #Constant head
                ###Look for all 'elements' with the current line ID############
                for j in range(__nbc):
                    #Check if 'element' line ID equals current line ID
                    if self.bc[j][1]==__lineID[i]:
                        #Check the first node
                        if self.bc[j][2] not in __aux_d: #Check if node has already been added
                            __aux_d.append(self.bc[j][2])
                            __aux_dv.append(__lineV[i][1])
                        #Check the second node
                        if self.bc[j][3] not in __aux_d: #Check if node has already been added
                            __aux_d.append(self.bc[j][3])
                            __aux_dv.append(__lineV[i][1])
                ###############################################################
            else: #Gradient
                __l=0 #Auxiliar variable to count number of nodes on the current line
                ###Look for all 'elements' with the current line ID############
                for j in range(__nbc):
                    if self.bc[j][1]==__lineID[i]:
                        #Check the first node
                        if self.bc[j][2] not in __aux_d: #Check if node has already been added
                            __aux_d.append(self.bc[j][2])
                            __l+=1
                        #Check the second node
                        if self.bc[j][3] not in __aux_d: #Check if node has already been added
                            __aux_d.append(self.bc[j][3])
                            __l+=1
                __delta=(__lineV[i][2]-__lineV[i][1])/(__l-1)
                for j in range(__nbc):
                    for k in range(__l):
                        __aux_dv.append(__lineV[i][1]+k*__delta)
        ###WRITING LOOP########################################################
        __l=len(__aux_d)
        __f.write(str(__l)+'\n') #Write # of Dirichlet Nodes
        for i in range(__l):
            __f.write(str(__aux_d[i]-1)+','+str(__aux_dv[i])+'\n')
        #######################################################################
        #######################################################################
        
        ##NEUMANN BC's#########################################################
        __f.write('#NEUMANN BC\n') #Write the section title
        ###Identify lines######################################################
        ####Create auxiliar arrays for storing line ID's and values############
        __lineID=list()
        __lineV=list()
        #######################################################################
        ###Loop for identifying lines and assigning them flow values###########
        for i in range(__nbc):
            #Check BC flag; 1=Neumann BC
            if self.bc[i][0]==1:
                #Check if 'Line ID' has already been stored
                if self.bc[i][1] not in __lineID:
                    __lineID.append(self.bc[i][1])
                    __v=input('Please, input the Neumann BC at line '+\
                              str(self.bc[i][1])+'.\n') #Read Neumann BC value
                    __lineV.append(float(__v))
        ###Loop for identifying lines and assigning them flow values###########
        #######################################################################
        ##Assign flow value to each line in BC file############################
        __l=0 #Auxiliar variable to store the number of Neumann BC faces
        __aux_n=list()
        for i in range(len(__lineID)):
            ###Look for all 'elements' with the curre line ID##################
            for j in range(__nbc):
                #Check if 'element' line ID equals current line ID
                if self.bc[j][1]==__lineID[i]:
                    __aux_n.append(str(self.bc[j][2]-1)+','+\
                                   str(self.bc[j][3]-1)+\
                              ','+str(__lineV[i])+'\n')
                    __l+=1
            ###################################################################
        __f.write(str(__l)+'\n') #Write # of Neumann Faces
        ###Actually write flow values in file##################################
        for i in range(__l):
            __f.write(__aux_n[i])
        #######################################################################
        #######################################################################

        ##WELL AS NEUMANN BC's#################################################
        __f.write('#WELL AS NEUMANN BC\n') #Write the section title
        ###Loop for identifying lines and assigning them flow values###########
        __lineID=list() #auxiliar array to store line IDs
        __lineV=list() #auxiliar array to store pumping intervals
        __well=list() #auxiliar array to identify wells
        for i in range(__nbc):
            #Check BC flag; 2=Well expressed as a Neumann BC
            if self.bc[i][0]==2:
                #Check if 'Line ID' has already been stored
                if self.bc[i][1] not in __lineID:
                    __lineID.append(self.bc[i][1])
                    __file=input('Please input the pumping schedule file of the well at line '+\
                              str(self.bc[i][1])+'.\n') #Read well's pumping/ rate
                    ##RETRIEVE WELL INFORMATION FROM FILE######################
                    __f2=aux_f.read_dict(__file)
                    __lineV.append(__f2['rate'])
                    __well.append(__file)
                    ##RETRIEVE WELL INFORMATION FROM FILE######################
        ###Loop for identifying lines and assigning them flow values###########
        ###Loop for getting the perimeter of each Neumann BC line##############
        __lineP=list() #Array for storing the perimeter values
        for i in range(len(__lineID)):
            __p=0 #Auxiliar variable for storing the perimeter of each line
            ###Look for all 'elements' with the curre line ID##################
            for j in range(__nbc):
                if self.bc[j][1] in __lineID:
                    #Get well ID
                    __ID=__lineID.index(self.bc[j][1])
                    #Check if 'element' line ID equals current line ID
                    if __well[__ID]==__well[i]:
                        __l=((self.m_p[self.bc[j][3]-1][0]-self.m_p[self.bc[j][2]-1][0])**2+\
                             (self.m_p[self.bc[j][3]-1][1]-self.m_p[self.bc[j][2]-1][1])**2)**0.5
                        __p+=__l
            __lineP.append(__p)
        ###Loop for getting the perimeter of each Neumann BC line##############
        #######################################################################
        ##Assign flow value to each line in BC file############################
        __l=0 #Auxiliar variable to store the number of Neumann BC faces
        __aux_n=list()
        for i in range(len(__lineID)):
            ###Look for all 'elements' with the curre line ID##################
            for j in range(__nbc):
                #Check if 'element' line ID equals current line ID
                if self.bc[j][1]==__lineID[i]:
                    __line=str(self.bc[j][2]-1)+','+\
                    str(self.bc[j][3]-1) #Node on Neumann Face
                    if isinstance(__lineV[i],list): #Transient-state
                        __n=len(__lineV[i])
                        __line+=','+str(__n)
                        for k in range(__n):
                            __line+=','+str(__lineV[i][k][0])+','+\
                            str(__lineV[i][k][1])+','+\
                            str(__lineV[i][k][2]/__lineP[i])
                    else: #Steady-state
                        __n=1
                        __line+=','+str(__n)+','+str(__lineV[i]/__lineP[i])
                    __line+='\n' #New line character
                    __aux_n.append(__line)
                    __l+=1
            ###################################################################
        __f.write(str(__l)+'\n') #Write # of Neumann Faces
        ###Actually write flow values in file##################################
        for i in range(__l):
            __f.write(__aux_n[i])
        #######################################################################
        #######################################################################
        
        __f.close() #Close the file
        #WRITE A BC FILE BASED ON MESH INFO####################################
        
    def coord2node(self,x,y):
        #RETURN THE NEAREST NODE TO THE GIVEN COORDINATES######################
        ##Find the nearest node to the given coordinates by brute force########
        __min=100 #Auxiliar variable to store the previous minimum distance
        for i in range(len(self.m_p)):
            __d=((self.m_p[i][0]-x)**2+\
                 (self.m_p[i][1]-y)**2)**0.5
            #If '__d' is lower than '__min', replace '__min' by '__d' and
            # assume 'i' as the closest node.
            if __d<__min:
                __min=__d #Replace '__min' by '__d'
                __node=i #Assume 'i' is the closes node to '(x,y)'
            #If '__d' is equal to 0, then 'i' is the closest node to '(x,y)'
            elif __d==0: #Check if the node is equal to 0
                __node=i
                break
        #######################################################################
        return __node
        #RETURN THE NEAREST NODE TO THE GIVEN COORDINATES######################
        
    def coord2element(self,x,y):
        #RETURN THE NEAREST NODE TO THE GIVEN COORDINATES######################
        ##Find the nearest node to the given coordinates by brute force########
        __min=1000 #Auxiliar variable to store the previous minimum distance
        for i in range(len(self.m_c)):
            __n1=self.m_c[i][0]
            __n2=self.m_c[i][1]
            __n3=self.m_c[i][2]
            __d=((self.m_p[__n1][0]-x)**2+(self.m_p[__n1][1]-y)**2)**0.5+\
            ((self.m_p[__n2][0]-x)**2+(self.m_p[__n2][1]-y)**2)**0.5+\
            ((self.m_p[__n3][0]-x)**2+(self.m_p[__n3][1]-y)**2)**0.5
            #If '__d' is lower than '__min', replace '__min' by '__d' and
            # assume 'i' as the closest node.
            if __d<__min:
                __min=__d #Replace '__min' by '__d'
                __element=i #Assume 'i' is the closes node to '(x,y)'
        #######################################################################
        return __element
        #RETURN THE NEAREST NODE TO THE GIVEN COORDINATES######################
"""
##############################MESH CLASS SECTION###############################
"""
###############################################################################
"""
#########################PHYSICAL MEDIUM CLASS SECTION#########################
"""
class phymed:
    """
    n_n = # of nodes                | n_e = # of elements
    k_e = Hydraulic Conductivity    | ss_e = Specific Storage
    bc_d = Dirichlet BC's           | bc_n = Neumann BC's
    hh_n = Hydraulic Head           | q = Pumping/Recharge Well Rate
    """
    def __init__(self,config,n_e,n_n):
        #RETRIEVE MEDIUM INFO FROM THE CONFIG FILE#############################
        self.steady=config['steady']
        ##IF WE HAVE A TRANSIENT MEDIUM (STEADY=FALSE), READ TIME PARAMETERS###
        if self.steady==False:
            self.dt=config['dt']
            self.time_steps=int(config['total_time']/config['dt'])+1
        #######################################################################
        #######################################################################
        #INITIALIZE ALL THE ARRAYS USED FOR STORING PHYSICAL PROPERTIES########
        ##INITIALIZE HYDRAULIC HEAD ARRAY######################################
        self.hh_n=[0]*n_n
        #######################################################################
        ##INITIALIZE HYDRAULIC CONDUCTIVITY ARRAY##############################
        self.k_e=[1]*n_e #Heterogeneity
        for i in range(n_e): #Anisotropy
            self.k_e[i]=[float(config['Kx']),float(config['Ky'])] #Hydraulic Conductivity Array (n_e x 2)
        #######################################################################
        ##INITIALIZE SPECIFIC STORAGE ARRAY####################################
        self.ss_e=[float(config['Ss'])]*n_e
        #######################################################################
        ##INITIALIZE WELL CONDITIONS ARRAY#####################################
        self.q=[0]*n_n
        #######################################################################
        ##ADD TEMPORAL AXIS IF MEDIUM IS IN TRANSIENT-STATE####################
        if self.steady==False: #Check steady-state flag
            ###TEMPORAL AXIS FOR NODE PROPERTIES
            for i in range(n_n):
                self.q[i]=[0]*self.time_steps
            #############################
            self.hh_n=[0]*self.time_steps
            for i in range(self.time_steps):
                self.hh_n[i]=[0]*n_n
        #######################################################################
        #INITIALIZE ALL THE ARRAYS USED FOR STORING BOUNDARY CONDITIONS########
        self.bc_d=list()
        self.bc_n=list()
        self.well_asbcn=list()
        #######################################################################
        #INITIALIZE DICTIONARY FOR STORING LAYER/ZONE PROPERTIES###############
        ##INITIALIZE AUXILIAR ARRAY FOR STORING NODE IDS
        __elements=[0]*n_e
        for i in range(n_e):
            __elements[i]=i
        ##DEFINE DICTIONARY WITH ZONE PROPERTIES
        self.zones={0:{'Ss':float(config['Ss']),
                       'Kx':float(config['Kx']),
                       'Ky':float(config['Ky']),
                       'elements':__elements}}
        #######################################################################
        del __elements,n_n,n_e,config
        
    #READ BOUNDARY CONDITIONS##################################################
    def read_bc_file(self,file):
        #READ BC FILE##########################################################
        __f=open(aux_f.path+'/'+file) #Open the given file in 'read' mode
        
        #DIRICHLET BC's########################################################
        for i in range(2):
            __f.readline() #Skip 2 lines
        __n=int(__f.readline().rstrip('\n')) #Read number of Dirichlet nodes
        ##READ DIRICHLET BC's FROM FILE########################################
        for i in range(__n):
            __aux=__f.readline().rstrip('\n').split(',') #Read and split at ','
            self.bc_d.append([int(__aux[0]),float(__aux[1])]) #Append to 'bc_d'
        #######################################################################
        #######################################################################
        
        #NEUMANN BC's##########################################################
        __f.readline() #Skip 1 line
        __n=int(__f.readline().rstrip('\n')) #Read number of Neumann faces
        ##READ NEUMANN BC's FROM FILE##########################################
        for i in range(__n):
            __aux=__f.readline().rstrip('\n').split(',')
            self.bc_n.append([int(__aux[0]),int(__aux[1]),float(__aux[2])])
        #######################################################################
        #######################################################################
        
        #WELL AS NEUMANN BC's##################################################
        __f.readline() #Skip 1 line
        __n=int(__f.readline().rstrip('\n')) #No. of nodes on the BC
        if self.steady==True: #Steady-state
            for i in range(__n):
                __aux=__f.readline().rstrip('\n').split(',')
                self.well_asbcn.append([int(__aux[0]),int(__aux[1]),float(__aux[3])])
        else: #Transient-state
            for i in range(__n):
                __aux=__f.readline().rstrip('\n').split(',')
                #__aux[2] has the number of pumping intervals
                __v=list()
                for j in range(int(__aux[2])):
                    __v.append([float(__aux[3+j*3]),float(__aux[4+j*3]),float(__aux[5+j*3])])
                self.well_asbcn.append([int(__aux[0]),int(__aux[1]),int(__aux[2]),__v])
        #WELL AS NEUMANN BC's##################################################
        
        __f.close()
        #READ BC FILE##########################################################
        del __n,__aux
    #READ BOUNDARY CONDITIONS##################################################
        
    #READ INITIAL CONDITIONS###################################################
    def read_hh0(self,file):
        #READS AN HHD FILE AND ASSIGNS IT TO HH_0##############################
        __f=open(aux_f.path+'/'+file)
        
        __n=int(__f.readline()) #Read number of nodes
        for i in range(__n):
            __aux=__f.readline().split()
            self.hh_n[0][i]=float(__aux[1])
        
        __f.close()
        #READS AN HHD FILE AND ASSIGNS IT TO HH_0##############################
        del __n,__aux
    #READ INITIAL CONDITIONS###################################################
"""
#########################PHYSICAL MEDIUM CLASS SECTION#########################
"""
###############################################################################
"""
###########################FEM OBJECT CLASS SECTION############################
"""
class fem:
    """
    mesh = Mesh Object          | phymed = Physical Medium Object
    m_mass = Mass Matrix        | m_stiffness = Stiffnes Matrix
    v_load = Load Vector        |
    """  
    def __init__(self,init):
        self.type='FEM'
        #READ INIT FILE########################################################
        self.config=aux_f.read_dict(init)
        #READ INIT FILE########################################################
        #CREATE MESH AND PHYMED OBJECTS########################################
        self.mesh=mesh(self.config['mesh_file'])
        self.phymed=phymed(self.config,self.mesh.n_e,self.mesh.n_n)
        #CREATE MESH AND PHYMED OBJECTS########################################
        #INITIALIZE MATRICES ARRAYS
        self.m_mass=list() #INITIALIZE ARRAY FOR MASS MATRIX
        self.m_stiffness=list() #INITIALIZE ARRAY FOR STIFFNESS MATRIX
        #INITIALIZE ARRAY FOR LOAD VECTOR
        self.v_load=list()
        #######################################################################
        
###############################################################################
#####ASSEMBLE MATRICES AND VECTORS#############################################
###############################################################################
    ##ASSEMBLE MASS MATRIX#####################################################
    def mass_m_assembly(self):
        #INITIALIZE ARRAY FOR MASS MATRIX
        self.m_mass=[0]*self.mesh.n_n
        for i in range(self.mesh.n_n):
            self.m_mass[i]=[0]*self.mesh.n_n
        #############################
        for k in range(self.mesh.n_e):
            ##CALCULATE LOCAL ELEMENT MASS MATRIX##############################
            for i in range(3): #Local Element Mass Matrix Index 'i'
                for j in range(3): #Local Element Mass Matrix Index 'j'
                    if i==j: #Main diagonal
                        __lmm_ij=self.mesh.a_e[k]*self.phymed.ss_e[k]/6
                    else:
                        __lmm_ij=self.mesh.a_e[k]*self.phymed.ss_e[k]/12
                    (global_i,global_j)=self.local2global(k,i,j)
                    self.m_mass[global_i][global_j]+=__lmm_ij #Add to Global Mass Matrix
            ###################################################################
    ##ASSEMBLE MASS MATRIX#####################################################
    
    ##ASSEMBLE STIFFNESS MATRIX################################################
    def stiff_m_assembly(self):
        #INITIALIZE ARRAY FOR STIFFNESS MATRIX
        self.m_stiffness=[0]*self.mesh.n_n
        for i in range(self.mesh.n_n):
            self.m_stiffness[i]=[0]*self.mesh.n_n
        ##############################
        for k in range(self.mesh.n_e):
            ##CALCULATE LOCAL ELEMENT STIFFNESS MATRIX#########################
            ###GET ELEMENT THE HAT FUNCTION FOR ELEMENT K######################
            (__a,__b,__c)=self.hat_phi(k)
            ###################################################################
            for i in range(3): #Local Element Stiffness Matrix Index 'i'
                for j in range(3): #Local Element Stiffness Matrix Index 'j'
                    __lsm_ij=(self.phymed.k_e[k][0]*(__b[i]*__b[j])+\
                    self.phymed.k_e[k][1]*(__c[i]*__c[j]))*self.mesh.a_e[k]
                    (global_i,global_j)=self.local2global(k,i,j)
                    self.m_stiffness[global_i][global_j]+=__lsm_ij #Add to Global Stiffness Matrix
            ###################################################################
    ##ASSEMBLE STIFFNESS MATRIX################################################
    
    ##ASSEMBLE LOAD VECTOR#####################################################
    def load_v_assembly(self):
        #Steady-state##########################################################
        if self.phymed.steady==True:
            #INITIALIZE ARRAY FOR LOAD VECTOR
            self.v_load=[0]*self.mesh.n_n
            #ASSIGN VALUES
            for i in range(self.mesh.n_n):
                self.v_load[i]-=self.phymed.q[i]
        #Steady-state##########################################################
        
        #Transient-state#######################################################
        else:
            #INITIALIZE ARRAY FOR LOAD VECTOR
            self.v_load=[0]*self.phymed.time_steps
            for i in range(self.phymed.time_steps):
                self.v_load[i]=[0]*self.mesh.n_n
            #ASSIGN VALUES
            for i in range(self.mesh.n_n):
                for k in range(self.phymed.time_steps):
                    self.v_load[k][i]-=self.phymed.q[i][k]
        #Transient-state#######################################################
    ##ASSEMBLE LOAD VECTOR#####################################################
    
    ###########################################################################
    #MISC FUNCTIONS############################################################
    ###########################################################################
    
    ##LOCAL2GLOBAL#############################################################
    def local2global(self,k,i,j):
        #'k'=element index, 'i'=local element index i, 'j'=local element index j
        global_i=self.mesh.m_c[k][i]
        global_j=self.mesh.m_c[k][j]
        return (global_i,global_j)
    ##LOCAL2GLOBAL#############################################################
    
    ##COMPUTE HAT FUNCTION OF ELEMENT 'K'######################################
    def hat_phi(self,k):
        __n1=self.mesh.m_c[k][0] #Node 1
        __n2=self.mesh.m_c[k][1] #Node 2
        __n3=self.mesh.m_c[k][2] #Node 3
        __j=2*self.mesh.a_e[k] #Jacobian
        a=[(self.mesh.m_p[__n2][0]*self.mesh.m_p[__n3][1]-self.mesh.m_p[__n3][0]*self.mesh.m_p[__n2][1])/__j,\
           (self.mesh.m_p[__n3][0]*self.mesh.m_p[__n1][1]-self.mesh.m_p[__n1][0]*self.mesh.m_p[__n3][1])/__j,\
           (self.mesh.m_p[__n1][0]*self.mesh.m_p[__n2][1]-self.mesh.m_p[__n2][0]*self.mesh.m_p[__n1][1])/__j]
        b=[(self.mesh.m_p[__n2][1]-self.mesh.m_p[__n3][1])/__j,\
           (self.mesh.m_p[__n3][1]-self.mesh.m_p[__n1][1])/__j,\
           (self.mesh.m_p[__n1][1]-self.mesh.m_p[__n2][1])/__j]
        c=[(self.mesh.m_p[__n3][0]-self.mesh.m_p[__n2][0])/__j,\
           (self.mesh.m_p[__n1][0]-self.mesh.m_p[__n3][0])/__j,\
           (self.mesh.m_p[__n2][0]-self.mesh.m_p[__n1][0])/__j]
        return (a,b,c)
    ##COMPUTE HAT FUNCTION OF ELEMENT 'K'######################################
    
    ###########################################################################
    #MISC FUNCTIONS############################################################
    ###########################################################################    
###############################################################################
#####ASSEMBLE MATRICES AND VECTORS#############################################
###############################################################################
    
###############################################################################
#####APPLY BOUNDARY CONDITIONS#################################################
###############################################################################
    ##APPLY DIRICHLET BC's#####################################################
    def apply_bc_d(self):
        #UNIT CONVERSION#######################################################
        uF=1
        if self.config['length'].lower()!=self.config['head'].lower(): #Check if conversion is required
            uF=aux_f.unit_dict[self.config['length'].lower()]/aux_f.unit_dict[self.config['head'].lower()] #Get conversion factor
        #UNIT CONVERSION#######################################################
        #Sort the 'bc_d' array for an easier time 'popping' nodes.
        self.phymed.bc_d.sort()
        __counter=0 #Auxiliar counter used for keeping track of how many nodes have been removed
        #Steady-state##########################################################
        if self.phymed.steady==True:
            for k in range(len(self.phymed.bc_d)):
                __node=self.phymed.bc_d[k][0]-__counter
                ##REMOVE NODE FROM HYDRAULIC HEAD VECTOR#######################
                del self.phymed.hh_n[__node]
                ##REMOVE NODE FROM LOAD VECTOR#################################
                del self.v_load[__node]
                ##REMOVE ROW FROM STIFFNESS MATRIX#############################
                del self.m_stiffness[__node]
                ##ADDITIONAL OPERATIONS########################################
                for i in range(len(self.v_load)):
                    self.v_load[i]-=self.m_stiffness[i][__node]*\
                    self.phymed.bc_d[k][1]*uF
                    ###REMOVE COLUMN FROM STIFFNESS MATRIX#####################
                    del self.m_stiffness[i][__node]
                __counter+=1
        #Steady-state##########################################################
        #Transient-state#######################################################
        if self.phymed.steady==False:
            for k in range(len(self.phymed.bc_d)):
                __node=self.phymed.bc_d[k][0]-__counter
                ##REMOVE NODE FROM HYDRAULIC HEAD VECTOR#######################
                for i in range(self.phymed.time_steps):
                    del self.phymed.hh_n[i][__node]
                ##REMOVE NODE FROM LOAD VECTOR#################################
                for i in range(self.phymed.time_steps):
                    del self.v_load[i][__node]
                ##REMOVE ROW FROM STIFFNESS MATRIX#############################
                del self.m_stiffness[__node]
                ##REMOVE ROW FROM MASS MATRIX##################################
                del self.m_mass[__node]
                ##ADDITIONAL OPERATIONS########################################
                for i in range(len(self.v_load[0])):
                    for j in range(self.phymed.time_steps):
                        self.v_load[j][i]-=self.m_stiffness[i][__node]*\
                        self.phymed.bc_d[k][1]*uF
                    ###REMOVE COLUMN FROM STIFFNESS MATRIX#####################
                    del self.m_stiffness[i][__node]
                    ##REMOVE COLUMN FROM MASS MATRIX###########################
                    del self.m_mass[i][__node]
                __counter+=1
        #Transient-state#######################################################
    ##APPLY DIRICHLET BC's#####################################################
    
    ##APPLY NEUMANN BC's#######################################################
    def apply_bc_n(self):
        #UNIT CONVERSION#######################################################
        uF=1
        if self.config['length'].lower()!=self.config['head'].lower(): #Check if conversion is required
            uF=(aux_f.unit_dict[self.config['length'].lower()])/(aux_f.unit_dict[self.config['head'].lower()]) #Get conversion factor            
        #UNIT CONVERSION#######################################################
        #Steady-state##########################################################
        if self.phymed.steady==True:
            for k in range(len(self.phymed.bc_n)):
                __l=((self.mesh.m_p[self.phymed.bc_n[k][1]][0]-self.mesh.m_p[self.phymed.bc_n[k][0]][0])**2+\
                     (self.mesh.m_p[self.phymed.bc_n[k][1]][1]-self.mesh.m_p[self.phymed.bc_n[k][0]][1])**2)**0.5
                self.v_load[self.phymed.bc_n[k][0]]+=self.phymed.bc_n[k][2]*__l*self.config['aq_thickness']*uF/2
                self.v_load[self.phymed.bc_n[k][1]]+=self.phymed.bc_n[k][2]*__l*self.config['aq_thickness']*uF/2
        #Steady-state##########################################################
        #Transient-state#######################################################
        else: #Transient-state
            for k in range(len(self.phymed.bc_n)):
                __l=((self.mesh.m_p[self.phymed.bc_n[k][1]][0]-self.mesh.m_p[self.phymed.bc_n[k][0]][0])**2+\
                     (self.mesh.m_p[self.phymed.bc_n[k][1]][1]-self.mesh.m_p[self.phymed.bc_n[k][0]][1])**2)**0.5
                for i in range(self.phymed.time_steps):
                    self.v_load[i][self.phymed.bc_n[k][0]]+=self.phymed.bc_n[k][2]*__l*self.config['aq_thickness']*uF/2 #RECORDATORIO: AGREGAR PASO TEMPORAL UNA VEZ QUE EL MÉTODO FUNCIONE
                    self.v_load[i][self.phymed.bc_n[k][1]]+=self.phymed.bc_n[k][2]*__l*self.config['aq_thickness']*uF/2 #RECORDATORIO: AGREGAR PASO TEMPORAL UNA VEZ QUE EL MÉTODO FUNCIONE
        #Transient-state#######################################################
    ##APPLY NEUMANN BC's#######################################################
    
    ##REBUILD HH ARRAY#########################################################
    def rebuild_hh(self):
        #UNIT CONVERSION#######################################################
        if self.config['length'].lower()!=self.config['head'].lower(): #Check if conversion is required
            uF=aux_f.unit_dict[self.config['head'].lower()]/aux_f.unit_dict[self.config['length'].lower()] #Get conversion factor
            for i in range(len(self.phymed.hh_n)):
                self.phymed.hh_n[i]*=uF
        #UNIT CONVERSION#######################################################
        #Sort the 'bc_d' array for an easier time 'restoring' nodes.
        self.phymed.bc_d.sort()
        ##ADD BACK THE DIRICHLET BC'S##########################################
        if self.phymed.steady==True: #Steady-state
            for i in range(len(self.phymed.bc_d)):
                self.phymed.hh_n.insert(self.phymed.bc_d[i][0],self.phymed.bc_d[i][1])
        else: #Transient-state
            for k in range(self.phymed.time_steps):
                for i in range(len(self.phymed.bc_d)):
                    self.phymed.hh_n[k].insert(self.phymed.bc_d[i][0],self.phymed.bc_d[i][1])
        ##ADD BACK THE DIRICHLET BC'S##########################################
    ##REBUILD HH ARRAY#########################################################  
###############################################################################
#####APPLY BOUNDARY CONDITIONS#################################################
###############################################################################

###############################################################################
#####MODIFY PHYSICAL MEDIUM####################################################
###############################################################################
    ##ADD WELLS################################################################
    def add_well(self,file): #Add well as constant function over an element
        #UNIT CONVERSION#######################################################
        uF=1
        if self.config['length'].lower()!=self.config['head'].lower(): #Check if conversion is required
            uF=(aux_f.unit_dict[self.config['length'].lower()]**3)/(aux_f.unit_dict[self.config['head'].lower()]**3) #Get conversion factor            
        #UNIT CONVERSION#######################################################
        ##RETRIEVE WELL INFORMATION FROM FILE##################################
        __f=aux_f.read_dict(file)
        ##RETRIEVE WELL INFORMATION FROM FILE##################################
        __element=self.mesh.coord2element(__f['x'],__f['y'])
        ##ASSIGN FLOW VALUE TO 'Q' ARRAY#######################################
        if self.phymed.steady==True: #Steady-state
            for i in range(3):
                self.phymed.q[self.mesh.m_c[__element][i]]=-__f['rate']*uF/\
                (3*self.config['aq_thickness']) #RECORDATORIO:PUEDE QUE HAYA QUE DIVIDIR ENTRE EL ESPESOR DEL ACUÍFERO
        else: #Transient-state
            __n=len(__f['rate']) #NUMBER OF PUMPING INTERVALS
            for i in range(__n):
                #IDENTIFY TIME-STEPS WITHIN THE PUMPING INTERVAL###############
                #IDENTIFY 'k_min'
                for k in range(self.phymed.time_steps):
                    if k*self.phymed.dt>=__f['rate'][i][0]:
                        t_min=k
                        break
                #IDENTIFY 'k_max'
                for k in range(t_min,self.phymed.time_steps):
                    if k*self.phymed.dt>=__f['rate'][i][1]:
                        t_max=k+1
                        break
                for j in range(3):
                    for k in range(t_min,t_max):
                        self.phymed.q[self.mesh.m_c[__element][j]][k]=-__f['rate'][i][2]*uF/(3*self.config['aq_thickness']) #RECORDATORIO:PUEDE QUE HAYA QUE DIVIDIR ENTRE EL ESPESOR DEL ACUÍFERO
                        
    def add_well2(self,file): #Add well as a point source over a node
        #UNIT CONVERSION#######################################################
        uF=1
        if self.config['length'].lower()!=self.config['head'].lower(): #Check if conversion is required
            uF=(aux_f.unit_dict[self.config['length'].lower()]**3)/(aux_f.unit_dict[self.config['head'].lower()]**3) #Get conversion factor            
        #UNIT CONVERSION#######################################################
        ##RETRIEVE WELL INFORMATION FROM FILE##################################
        __f=aux_f.read_dict(file)
        ##RETRIEVE WELL INFORMATION FROM FILE##################################
        __node=self.mesh.coord2node(__f['x'],__f['y'])
        ##ASSIGN FLOW VALUE TO 'Q' ARRAY#######################################
        if self.phymed.steady==True: #Steady-state
            for i in range(3):
                self.phymed.q[__node]=-__f['rate']*uF/\
                (self.config['aq_thickness']) #RECORDATORIO:PUEDE QUE HAYA QUE DIVIDIR ENTRE EL ESPESOR DEL ACUÍFERO
        else: #Transient-state
            __n=len(__f['rate']) #NUMBER OF PUMPING INTERVALS
            for i in range(__n):
                #IDENTIFY TIME-STEPS WITHIN THE PUMPING INTERVAL###############
                #IDENTIFY 'k_min'
                for k in range(self.phymed.time_steps):
                    if k*self.phymed.dt>=__f['rate'][i][0]:
                        t_min=k
                        break
                #IDENTIFY 'k_max'
                for k in range(t_min,self.phymed.time_steps):
                    if k*self.phymed.dt>=__f['rate'][i][1]:
                        t_max=k+1
                        break
                for j in range(3):
                    for k in range(t_min,t_max):
                        self.phymed.q[__node][k]=-__f['rate'][i][2]*uF/(self.config['aq_thickness']) #RECORDATORIO:PUEDE QUE HAYA QUE DIVIDIR ENTRE EL ESPESOR DEL ACUÍFERO
                        
    def add_well3(self,file): #Add well as point source over an element
        #UNIT CONVERSION#######################################################
        uF=1
        if self.config['length'].lower()!=self.config['head'].lower(): #Check if conversion is required
            uF=(aux_f.unit_dict[self.config['length'].lower()]**3)/(aux_f.unit_dict[self.config['head'].lower()]**3) #Get conversion factor            
        #UNIT CONVERSION#######################################################
        ##RETRIEVE WELL INFORMATION FROM FILE##################################
        __f=aux_f.read_dict(file)
        ##RETRIEVE WELL INFORMATION FROM FILE##################################
        __element=self.mesh.coord2element(__f['x'],__f['y'])
        #Retrieve nodes for interpolation
        __n1=self.mesh.m_c[__element][0] #Node 1 index
        __n2=self.mesh.m_c[__element][1] #Node 2 index
        __n3=self.mesh.m_c[__element][2] #Node 3 index
        ##GET INTERPOLATION ELEMENT########################
        ##GET NATURAL SYSTEM VALUES########################
        __x21=self.mesh.m_p[__n2][0]-self.mesh.m_p[__n1][0]
        __y21=self.mesh.m_p[__n2][1]-self.mesh.m_p[__n1][1]
        __x31=self.mesh.m_p[__n3][0]-self.mesh.m_p[__n1][0]
        __y31=self.mesh.m_p[__n3][1]-self.mesh.m_p[__n1][1]
        __xi=(__y31*(__f['x']-self.mesh.m_p[__n1][0])-__x31*(__f['y']-self.mesh.m_p[__n1][1]))/\
        (__x21*__y31-__x31*__y21)
        __eta=(-__y21*(__f['x']-self.mesh.m_p[__n1][0])+__x21*(__f['y']-self.mesh.m_p[__n1][1]))/\
        (__x21*__y31-__x31*__y21)
        __N=[1-__xi-__eta,__xi,__eta]
        ##GET NATURAL SYSTEM VALUES########################
        ##COMPUTE DRAWDOWN VALUES##########################
        ##ASSIGN FLOW VALUE TO 'Q' ARRAY#######################################
        if self.phymed.steady==True: #Steady-state
            for i in range(3):
                self.phymed.q[self.mesh.m_c[__element][i]]=-__f['rate']*__N[i]*uF/\
                (self.config['aq_thickness']) #RECORDATORIO:PUEDE QUE HAYA QUE DIVIDIR ENTRE EL ESPESOR DEL ACUÍFERO
        else: #Transient-state
            __n=len(__f['rate']) #NUMBER OF PUMPING INTERVALS
            for i in range(__n):
                #IDENTIFY TIME-STEPS WITHIN THE PUMPING INTERVAL###############
                #IDENTIFY 'k_min'
                for k in range(self.phymed.time_steps):
                    if k*self.phymed.dt>=__f['rate'][i][0]:
                        t_min=k
                        break
                #IDENTIFY 'k_max'
                for k in range(t_min,self.phymed.time_steps):
                    if k*self.phymed.dt>=__f['rate'][i][1]:
                        t_max=k+1
                        break
                for j in range(3):
                    for k in range(t_min,t_max):
                        self.phymed.q[self.mesh.m_c[__element][j]][k]=-__f['rate'][i][2]*__N[j]*uF/\
                        (self.config['aq_thickness']) #RECORDATORIO:PUEDE QUE HAYA QUE DIVIDIR ENTRE EL ESPESOR DEL ACUÍFERO
    
    def add_well_asnbc(self):
        #UNIT CONVERSION#######################################################
        uF=1
        if self.config['length'].lower()!=self.config['head'].lower(): #Check if conversion is required
            uF=(aux_f.unit_dict[self.config['length'].lower()]**3)/(aux_f.unit_dict[self.config['head'].lower()]**3) #Get conversion factor            
        #UNIT CONVERSION#######################################################
        if self.phymed.steady==True: #Steady-state
            for k in range(len(self.phymed.well_asbcn)):
                __l=((self.mesh.m_p[self.phymed.well_asbcn[k][1]][0]-self.mesh.m_p[self.phymed.well_asbcn[k][0]][0])**2+\
                     (self.mesh.m_p[self.phymed.well_asbcn[k][1]][1]-self.mesh.m_p[self.phymed.well_asbcn[k][0]][1])**2)**0.5
                self.phymed.q[self.phymed.well_asbcn[k][0]]-=self.phymed.well_asbcn[k][2]*__l*uF/(2*self.config['aq_thickness'])
                self.phymed.q[self.phymed.well_asbcn[k][1]]-=self.phymed.well_asbcn[k][2]*__l*uF/(2*self.config['aq_thickness'])
        else: #Transient-state
            for k in range(len(self.phymed.well_asbcn)):
                __l=((self.mesh.m_p[self.phymed.well_asbcn[k][1]][0]-self.mesh.m_p[self.phymed.well_asbcn[k][0]][0])**2+\
                     (self.mesh.m_p[self.phymed.well_asbcn[k][1]][1]-self.mesh.m_p[self.phymed.well_asbcn[k][0]][1])**2)**0.5
                __n=self.phymed.well_asbcn[k][2]
                __v=self.phymed.well_asbcn[k][3]
                #IDENTIFY TIME-STEPS WITHIN THE PUMPING INTERVAL###############
                for i in range(__n):
                    #IDENTIFY 'k_min'
                    for j in range(self.phymed.time_steps):
                        if j*self.phymed.dt>=__v[i][0]:
                            t_min=j
                            break
                    #IDENTIFY 'k_max'
                    for j in range(t_min,self.phymed.time_steps):
                        if j*self.phymed.dt>=__v[i][1]:
                            t_max=j+1
                            break
                    for t in range(t_min,t_max):
                        self.phymed.q[self.phymed.well_asbcn[k][0]][t]-=__v[i][2]*__l*uF/(2*self.config['aq_thickness'])
                        self.phymed.q[self.phymed.well_asbcn[k][1]][t]-=__v[i][2]*__l*uF/(2*self.config['aq_thickness'])
    ##ADD WELLS################################################################
    
    #ZONE RELATED FUNCTIONS####################################################
    ##DEFINE NEW ZONE##########################################################
    def new_zone(self,zone_id,Ss:float,Kx:float,Ky:float,x_min:float,x_max:float,y_min:float,y_max:float): #Defines a new zone in the domain as a rectangle
        if zone_id not in self.phymed.zones: #CHECK IF THE ZONE EXISTS
            #IDENTIFY ALL ELEMENTS INSIDE THE ZONE#############################
            __ID=list() #Auxiliar list for storing the ID values of the elements inside the zone
            for i in range(self.mesh.n_e):
                #Get element centroid
                __cx=(self.mesh.m_p[self.mesh.m_c[i][0]][0]+\
                      self.mesh.m_p[self.mesh.m_c[i][1]][0]+\
                      self.mesh.m_p[self.mesh.m_c[i][2]][0])/3
                __cy=(self.mesh.m_p[self.mesh.m_c[i][0]][1]+\
                      self.mesh.m_p[self.mesh.m_c[i][1]][1]+\
                      self.mesh.m_p[self.mesh.m_c[i][2]][1])/3
                if x_min<=__cx<=x_max and y_min<=__cy<=y_max:
                    __ID.append(i)
            #IDENTIFY ALL ELEMENTS INSIDE THE ZONE#############################
            #ZONE ASSIGNATION LOOP#############################################
            for i in __ID:
                self.phymed.k_e[i]=[Kx,Ky]
                self.phymed.ss_e[i]=Ss
                #Look for the element in all other zones and remove it from their lists
                for k in self.phymed.zones:
                    if i in self.phymed.zones[k]['elements']: #Check if node is in a zone 'k'
                        __e_ID=self.phymed.zones[k]['elements'].index(i) #Get ID within 'nodes' list
                        del self.phymed.zones[k]['elements'][__e_ID] #Remove node from zone 'k'
            self.phymed.zones[zone_id]={'Ss':Ss,
                             'Kx':Kx,
                             'Ky':Ky,
                             'elements':__ID}
            #ZONE ASSIGNATION LOOP#############################################
            del zone_id,Ss,Kx,Ky,x_min,x_max,y_min,y_max,__cx,__cy,__e_ID,__ID
        else:
            print('Zone {0} already exists'.format(zone_id))
            del zone_id,Ss,Kx,Ky,x_min,x_max,y_min,y_max
    ##DEFINE NEW ZONE##########################################################
    
    ##APPEND NEW RECTANGLE TO EXISTING ZONE####################################
    def expand_zone(self,zone_id,x_min:float,x_max:float,y_min:float,y_max:float): #Adds a new rectangle to an existing zone
        if zone_id in self.phymed.zones: #CHECK IF THE ZONE EXISTS
            #IDENTIFY ALL ELEMENTS INSIDE THE ZONE#############################
            __ID=list() #Auxiliar list for storing the ID values of the elements inside the zone
            for i in range(self.mesh.n_e):
                #Get element centroid
                __cx=(self.mesh.m_p[self.mesh.m_c[i][0]][0]+\
                      self.mesh.m_p[self.mesh.m_c[i][1]][0]+\
                      self.mesh.m_p[self.mesh.m_c[i][2]][0])/3
                __cy=(self.mesh.m_p[self.mesh.m_c[i][0]][1]+\
                      self.mesh.m_p[self.mesh.m_c[i][1]][1]+\
                      self.mesh.m_p[self.mesh.m_c[i][2]][1])/3
                if x_min<=__cx<=x_max and y_min<=__cy<=y_max:
                    __ID.append(i)
            #IDENTIFY ALL ELEMENTS INSIDE THE ZONE#############################
            #ZONE ASSIGNATION LOOP#############################################
            for i in __ID:
                #Look for the element in all other zones and remove it from their lists
                for k in self.phymed.zones:
                    if i in self.phymed.zones[k]['elements']: #Check if node is in a zone 'k'
                        __e_ID=self.phymed.zones[k]['elements'].index(i) #Get ID within 'nodes' list
                        del self.phymed.zones[k]['elements'][__e_ID] #Remove node from zone 'k'
                #Add elements to zone
                if i not in self.phymed.zones[zone_id]['elements']:
                    self.phymed.zones[zone_id]['elements'].append(i)
            self.update_zone(zone_id,self.phymed.zones[zone_id]['Ss'],self.phymed.zones[zone_id]['Kx'],self.phymed.zones[zone_id]['Ky'])
            #ZONE ASSIGNATION LOOP#############################################
            del zone_id,x_min,x_max,y_min,y_max,__cx,__cy,__e_ID,__ID
        else:
            print("Zone {0} doesn't exist".format(zone_id))
            del zone_id,x_min,x_max,y_min,y_max
    ##APPEND NEW RECTANGLE TO EXISTING ZONE####################################
    
    ##UPDATE ZONE PARAMETERS###################################################
    def update_zone(self,zone_id,Ss:float,Kx:float,Ky:float): #Update zone properties
        if zone_id in self.phymed.zones: #CHECK IF THE ZONE EXISTS
            self.phymed.zones[zone_id]['Ss']=Ss #Update Specific Storage value
            self.phymed.zones[zone_id]['Kx']=Kx #Update Kx value
            self.phymed.zones[zone_id]['Ky']=Ky #Update Ky value
            for i in self.phymed.zones[zone_id]['elements']:
                self.phymed.k_e[i]=[Kx,Ky] #Update element value
                self.phymed.ss_e[i]=Ss #Update element value
            del zone_id,Ss,Kx,Ky
        else:
            print("Zone {0} doesn't exist".format(zone_id))
            del zone_id,Ss,Kx,Ky
    ##UPDATE ZONE PARAMETERS###################################################
    
    ##DELETE ZONE##############################################################
    def delete_zone(self,zone_id):
        if zone_id in self.phymed.zones: #CHECK IF THE ZONE EXISTS
            for i in self.phymed.zones[zone_id]['elements']:
                self.phymed.zones[0]['elements'].append(i) #Add element to default zone
                self.phymed.k_e[i]=[self.phymed.zones[0]['Kx'],self.phymed.zones[0]['Ky']] #Update element value
                self.phymed.ss_e[i]=self.phymed.zones[0]['Ss'] #Update element value
            del self.phymed.zones[zone_id],zone_id
        else:
            print("Zone {0} doesn't exist".format(zone_id))
            del zone_id
    ##DELETE ZONE##############################################################
    #ZONE RELATED FUNCTIONS####################################################
###############################################################################
#####MODIFY PHYSICAL MEDIUM####################################################
###############################################################################

###############################################################################
#####SYSTEM SOLVING############################################################
###############################################################################
    #SOLVE LINEAR SYSTEM#######################################################
    def solve(self):
        if self.phymed.steady==True: #Steady-state
            print('ASSEMBLING GLOBAL MATRICES...')
            ##ASSEMBLY GLOBAL MATRICES#########################################
            self.stiff_m_assembly() #ASSEMBLY STIFFNESS MATRIX
            self.load_v_assembly() #ASSEMBLY LOAD VECTOR
            ###################################################################
            
            print('APPLYING BCs...')
            ##APPLY BC'S#######################################################
            self.apply_bc_n()
            self.apply_bc_d()
            ##APPLY BC'S#######################################################
            
            print('COMPUTING HYDRAULIC HEAD DISTRIBUTION...')
            ##GET HYDRAULIC HEAD DISTRIBUTION##################################    
            __aux_hh=np.matmul(np.linalg.inv(self.m_stiffness),self.v_load)
            self.phymed.hh_n=__aux_hh.tolist()
            self.rebuild_hh()
            ##GET HYDRAULIC HEAD DISTRIBUTION##################################
        else: #Transient-state
            print('ASSEMBLING GLOBAL MATRICES...')
            ##ASSEMBLY GLOBAL MATRICES#########################################
            self.mass_m_assembly() #ASSEMBLY MASS MATRIX
            self.stiff_m_assembly() #ASSEMBLY STIFFNESS MATRIX
            self.load_v_assembly() #ASSEMBLY LOAD VECTOR
            ###################################################################
            
            print('APPLYING BCs...')
            ##APPLY BC'S#######################################################
            self.apply_bc_n()
            self.apply_bc_d()
            ##APPLY BC'S#######################################################
            
            print('ASSEMBLING AUXILIAR MATRICES...')
            ##GET AUXILIAR MATRICES############################################
            ##__A=MASS+STIFFNESS*dt
            ###INITIALIZE ARRAY
            __n_n=len(self.m_mass)
            __A=[0]*__n_n
            for i in range(__n_n):
                __A[i]=[0]*__n_n
            ###ASSIGN VALUES
            #DIVIDE MASS MATRIX BY DT
            for i in range(__n_n):
                for j in range(__n_n):
                    self.m_mass[i][j]/=self.phymed.dt
            for i in range(__n_n):
                for j in range(__n_n):
                    __A[i][j]=self.m_mass[i][j]+self.m_stiffness[i][j]#*self.phymed.dt
            __Ainv=np.linalg.inv(__A)
            ##__V=LOAD*dt+MASS*hh_t-1
            ###INITIALIZE ARRAY
            __V=[0]*self.phymed.time_steps
            for i in range(self.phymed.time_steps):
                __V[i]=[0]*__n_n
            ###ASSIGN VALUES
            for i in range(self.phymed.time_steps):
                for j in range(__n_n):
                        __V[i][j]=self.v_load[i][j]#*self.phymed.dt
            ##GET AUXILIAR MATRICES############################################
            
            print('SOLVING FOR EACH TIME STEP...')
            for k in range(self.phymed.time_steps-1): #TIME-LOOP
                __step=k+1
                print(str(100*__step/self.phymed.time_steps)+'%...')
                
                ##UPDATE AUXILIAR MATRICES#####################################
                __M=np.matmul(self.m_mass,self.phymed.hh_n[__step-1])
                for i in range(__n_n):
                    __V[__step][i]+=__M[i]
                ##UPDATE AUXILIAR MATRICES#####################################
                
                ##GET HYDRAULIC HEAD DISTRIBUTION##############################
                ###CORREGIR
                print('COMPUTING HYDRAULIC HEAD DISTRIBUTION...')
                __aux_hh=np.matmul(__Ainv,__V[__step])
                self.phymed.hh_n[__step]=__aux_hh.tolist()
                ##GET HYDRAULIC HEAD DISTRIBUTION##############################
            ###HH VALUES#######################################################
            self.rebuild_hh()
            print('100%')
            ###HH VALUES#######################################################
    #SOLVE LINEAR SYSTEM#######################################################
###############################################################################
#####SYSTEM SOLVING############################################################
###############################################################################
    
###############################################################################
#####OUPUTS####################################################################
###############################################################################
    ##EXPORT RESULTS TO GMSH'S .POS FORMAT#####################################
    def export_results(self,file,step=None):
        #WRITE AN OUTPUT FILE CONTAINING THE HYDRAULIC HEAD DISTRIBUTION#######
        __f=open(aux_f.path+'/output/'+file+'_hhd.msh','w') #Open the given file in 'write' mode
        __fOrig=open(aux_f.path+'/'+self.config['mesh_file']) #Original mesh file
        
        #WRITING THE FILE HEADERS##############################################
        #Copy mesh file to the new file
        for line in __fOrig:
            __f.write(line)
        #WRITING THE FILE HEADERS##############################################
        #LOOP FOR WRITING VALUES ON THE FILE###################################
        __f.write('$NodeData\n')
        __f.write('{0}\n'.format(1)) #One string-tag:
        __f.write('"Hydraulic Head Distribution"\n') #Name of view
        __f.write('{0}\n'.format(1)) #One real-tag:
        if(step==None):
            __f.write('{0:f}\n'.format(0)) #Simulation time = 0.0
            __f.write('{0}\n'.format(3)) #Three integer tags:
            __f.write('{0}\n'.format(0)) #Current step
            __f.write('{0}\n'.format(1)) #Scalar field
            __f.write('{0}\n'.format(self.mesh.n_n)) #Number of nodes
            for i in range(self.mesh.n_n):
                __f.write('{0} {1}\n'.format(i+1,self.phymed.hh_n[i]) )
        else:
            __f.write('{0:16f}\n'.format(step*self.phymed.dt)) #Simulation time = step * dt
            __f.write('{0}\n'.format(3)) #Three integer tags:
            __f.write('{0}\n'.format(step)) #Current step
            __f.write('{0}\n'.format(1)) #Scalar field
            __f.write('{0}\n'.format(self.mesh.n_n)) #Number of nodes
            for i in range(self.mesh.n_n):
                __f.write('{0} {1}\n'.format(i+1,self.phymed.hh_n[step][i]) )
        __f.write('$EndNodeData')
        #LOOP FOR WRITING VALUES ON THE FILE###################################
        __fOrig.close()
        __f.close()
    ##EXPORT RESULTS TO GMSH'S .POS FORMAT#####################################
    
    ##DRAWDOWN AT A GIVEN POINT################################################
    def drawdown_at(self,x,y,file=None):
        #COMPUTE DRAWDOWN VALUE AT X,Y#########################################
        ##GET INTERPOLATION ELEMENT########################
        __n=self.mesh.coord2element(x,y) #Get element ID
        #Retrieve nodes for interpolation
        __n1=self.mesh.m_c[__n][0] #Node 1 index
        __n2=self.mesh.m_c[__n][1] #Node 2 index
        __n3=self.mesh.m_c[__n][2] #Node 3 index
        ##GET INTERPOLATION ELEMENT########################
        ##GET NATURAL SYSTEM VALUES########################
        __x21=self.mesh.m_p[__n2][0]-self.mesh.m_p[__n1][0]
        __y21=self.mesh.m_p[__n2][1]-self.mesh.m_p[__n1][1]
        __x31=self.mesh.m_p[__n3][0]-self.mesh.m_p[__n1][0]
        __y31=self.mesh.m_p[__n3][1]-self.mesh.m_p[__n1][1]
        __xi=(__y31*(x-self.mesh.m_p[__n1][0])-__x31*(y-self.mesh.m_p[__n1][1]))/\
        (__x21*__y31-__x31*__y21)
        __eta=(-__y21*(x-self.mesh.m_p[__n1][0])+__x21*(y-self.mesh.m_p[__n1][1]))/\
        (__x21*__y31-__x31*__y21)
        __N1=1-__xi-__eta
        __N2=__xi
        __N3=__eta
        ##GET NATURAL SYSTEM VALUES########################
        ##COMPUTE DRAWDOWN VALUES##########################
        __h=[0]*self.phymed.time_steps #Auxiliar array to store head values
        drawdown=[0]*self.phymed.time_steps #Initialize array
        for k in range(self.phymed.time_steps):
            __h[k]=self.phymed.hh_n[k][__n1]*__N1+\
            self.phymed.hh_n[k][__n2]*__N2+\
            self.phymed.hh_n[k][__n3]*__N3
            drawdown[k]=__h[0]-__h[k]
        ##COMPUTE DRAWDOWN VALUES##########################
        #COMPUTE DRAWDOWN VALUE AT X,Y#########################################
        if file!=None:
            __f=open(aux_f.path+'/output/'+file+'.dwn','w')
            __f.write('#FILE TYPE\n')
            __f.write('{0}\n'.format(0)) #Write file type: 0=Time
            __f.write('#UNITS\n')
            __f.write(self.config['time']+','+self.config['head']+'\n') #Write model time and head units
            __f.write('#NUMBER OF TIME-STEPS\n')
            __f.write('{0}\n'.format(self.phymed.time_steps)) #Write number of time-steps
            __f.write('#OBSERVATION POINT COORDINATES\n')
            __f.write('{0:10f} {1:10f}\n'.format(x,y)) #Write observation point intended coordinates
            __f.write('#TIME-STEP VALUES\n')
            for i in range(self.phymed.time_steps):
                __f.write('{0:.7E} '.format(i*self.phymed.dt))
            __f.write('\n')
            __f.write('#DRADOWN VALUES\n')
            for i in range(self.phymed.time_steps):
                __f.write('{0:.7E} '.format(drawdown[i]))
            __f.close()
        return drawdown
    ##DRAWDOWN AT A GIVEN POINT################################################
    
    ##DRAWDOWN PROFILE#########################################################
    def drawdown_profile(self,initial_point,final_point,npts,step,file=None):
        ##COMPUTE DRAWDOWN OVER A PROFILE######################################
        ###GET OBSERVATION POINT COORDINATES###################################
        obs_pts=[0]*npts #Auxiliar array to store the coords of obs. points
        dx=(final_point[0]-initial_point[0])/(npts-1)
        dy=(final_point[1]-initial_point[1])/(npts-1)
        for i in range(npts):
            xi=initial_point[0]+i*dx #'x'
            yi=initial_point[1]+i*dy #'y'
            di=i*(dx**2+dy**2)**0.5 #Distance to the initial point (over profile)
            obs_pts[i]=[xi,yi,di]
        ###GET OBSERVATION POINT COORDINATES###################################
        #COMPUTE DRAWDOWN VALUE AT EACH OBSERVATION POINT######################
        drawdown=[0]*npts #Initialize array
        for i in range(npts):
            ##GET INTERPOLATION ELEMENT########################
            x=obs_pts[i][0]
            y=obs_pts[i][1]
            __n=self.mesh.coord2element(x,y) #Get element ID
            #Retrieve nodes for interpolation
            __n1=self.mesh.m_c[__n][0] #Node 1 index
            __n2=self.mesh.m_c[__n][1] #Node 2 index
            __n3=self.mesh.m_c[__n][2] #Node 3 index
            ##GET INTERPOLATION ELEMENT########################
            ##GET NATURAL SYSTEM VALUES########################
            __x21=self.mesh.m_p[__n2][0]-self.mesh.m_p[__n1][0]
            __y21=self.mesh.m_p[__n2][1]-self.mesh.m_p[__n1][1]
            __x31=self.mesh.m_p[__n3][0]-self.mesh.m_p[__n1][0]
            __y31=self.mesh.m_p[__n3][1]-self.mesh.m_p[__n1][1]
            __xi=(__y31*(x-self.mesh.m_p[__n1][0])-__x31*(y-self.mesh.m_p[__n1][1]))/\
            (__x21*__y31-__x31*__y21)
            __eta=(-__y21*(x-self.mesh.m_p[__n1][0])+__x21*(y-self.mesh.m_p[__n1][1]))/\
            (__x21*__y31-__x31*__y21)
            __N1=1-__xi-__eta
            __N2=__xi
            __N3=__eta
            ##GET NATURAL SYSTEM VALUES########################
            ##COMPUTE DRAWDOWN VALUES##########################
            __h0=self.phymed.hh_n[0][__n1]*__N1+\
            self.phymed.hh_n[0][__n2]*__N2+\
            self.phymed.hh_n[0][__n3]*__N3
            __hi=self.phymed.hh_n[step][__n1]*__N1+\
            self.phymed.hh_n[step][__n2]*__N2+\
            self.phymed.hh_n[step][__n3]*__N3
            drawdown[i]=__h0-__hi
            ##COMPUTE DRAWDOWN VALUES##########################
        #COMPUTE DRAWDOWN VALUE AT EACH OBSERVATION POINT######################
        if file!=None:
            __f=open(aux_f.path+'/output/'+file+'.dwn','w')
            __f.write('#FILE TYPE\n')
            __f.write('{0}\n'.format(1)) #Write file type: 1=Profile
            __f.write('#UNITS\n')
            __f.write(self.config['length']+','+self.config['head']+'\n') #Write model length and head units
            __f.write('#NUMBER OF OBSERVATION POINTS\n')
            __f.write('{0}\n'.format(npts)) #Write number of observation points
            __f.write('#INITIAL POINT COORDINATES\n')
            __f.write('{0:10f} {1:10f}\n'.format(initial_point[0],initial_point[1])) #Write initial point coordinates
            __f.write('#OBSERVATION POINTS DISTANCE FROM INITIAL POINT\n')
            for i in range(npts):
                __f.write('{0:.7E} '.format(i*(dx**2+dy**2)**0.5))
            __f.write('\n')
            __f.write('#DRADOWN VALUES\n')
            for i in range(npts):
                __f.write('{0:.7E} '.format(drawdown[i]))
            __f.close()
        return drawdown
    ##DRAWDOWN PROFILE#########################################################
    
    ##AUXILIAR EXPORTATION OF INITIAL CONDITIONS###############################
    def export_hh0(self,file):
        #WRITE AN OUTPUT FILE CONTAINING THE HYDRAULIC HEAD DISTRIBUTION#######
        __f=open(aux_f.path+'/'+file+'.hh0','w') #Open the given file in 'write' mode
        
        __f.write(str(self.mesh.n_n)+'\n')
        for i in range(self.mesh.n_n):
            __f.write('{0} {1:.7E}\n'.format(i,self.phymed.hh_n[i]))
        __f.close()
    ##AUXILIAR EXPORTATION OF INITIAL CONDITIONS###############################
###############################################################################
#####OUPUTS####################################################################
###############################################################################
"""
##############################MESH CLASS SECTION###############################
"""