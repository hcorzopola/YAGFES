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
    nNodes = # of nodes                 | nElements = # of elements
    pointMatrix = Point Matrix          | cMatrix = Connectivity Matrix
    elementArea = Element's Area        | bc = Boundary Conditions 
    """ 
    def __init__(self,file):
        #READ MESH FILE########################################################
        __f=open(aux_f.path+'/'+file) #Open the given file in 'read' mode

        ##READ NODES INFO######################################################
        for i in range(4):
            __f.readline() #Skip 4 lines
        ###Read number of nodes################################################
        self.nNodes=int(__f.readline().rstrip('\n')) #Read number of nodes.
        #######################################################################
        ###Create an array ('Point Matrix') for storing the node coordinates###
        self.pointMatrix=[0]*self.nNodes #The 'Point Matrix' has dimensions of 'nNodes x 3'.
        #This allow us to store up to three coordinates per node.
        #######################################################################
        ###Retrieve node coordinates from the mesh file########################
        for i in range(self.nNodes):
            __aux_n=__f.readline().split()[1:] #Store the line in an auxiliar array. Skip node number.
            self.pointMatrix[i]=[float(x) for x in __aux_n]
        #######################################################################
        ##READ NODES INFO######################################################
        
        ##READ ELEMENTS INFO###################################################
        for i in range(2):
            __f.readline() #Skip 2 lines
        ###Read number of elements#############################################
        self.nElements=int(__f.readline().rstrip('\n')) #Read number of elements.
        #This number accounts for two-node elements. We will update this later,
        #since the final number of elements might not be the same.
        #######################################################################
        ###Create an array ('Connectivity Matrix') for storing the nodes per element
        self.cMatrix=list()
        #The 'Connectivity Matrix' has dimensions of 'nElements x 3', since we are
        # working with triangular elements. However, we won't specify its
        # dimensions in the code, since we are uncertain of 'nElements' final value.
        #######################################################################
        ###Create an array for storing Boundary Conditions (BC's)##############
        self.bc=list()
        #This array will store:
        #1) type of BC (0=Dirichlet, 1=Neumann),
        #2) boundary ID (since a boundary can be made of multiple lines) and
        #3) nodes composing the 'line' on the boundary.
        #######################################################################
        ###Retrieve elements from the mesh file################################
        for i in range(self.nElements):
            __aux_e=__f.readline().split() #Store the line in an auxiliar array
            if __aux_e[1]=='1': #Linear elements
                self.bc.append( list( map(int,__aux_e[3:]) ) ) #Store element in 'bc'
            else: #Triangular elements
                self.cMatrix.append( list( map(int,__aux_e[5:]) ) ) #Store element in 'cMatrix'
        ###Fix node index######################################################
        for i in range(len(self.cMatrix)):
            for j in range(3):
                self.cMatrix[i][j]=self.cMatrix[i][j]-1
        #######################################################################
        ###Update 'nElements' value##################################################
        self.nElements=len(self.cMatrix)
        #######################################################################
        ##READ ELEMENTS INFO###################################################
        __f.close() #Close the file
        #READ MESH FILE########################################################
        
        #CALCULATE ELEMENT AREA################################################
        self.elementArea=[0]*self.nElements #The array has a size equal to 'nElements'.
        for i in range(self.nElements):
            a=[self.pointMatrix[self.cMatrix[i][1]][0]-self.pointMatrix[self.cMatrix[i][0]][0],\
               self.pointMatrix[self.cMatrix[i][1]][1]-self.pointMatrix[self.cMatrix[i][0]][1],\
               self.pointMatrix[self.cMatrix[i][1]][2]-self.pointMatrix[self.cMatrix[i][0]][2]]
            b=[self.pointMatrix[self.cMatrix[i][2]][0]-self.pointMatrix[self.cMatrix[i][0]][0],\
               self.pointMatrix[self.cMatrix[i][2]][1]-self.pointMatrix[self.cMatrix[i][0]][1],\
               self.pointMatrix[self.cMatrix[i][2]][2]-self.pointMatrix[self.cMatrix[i][0]][2]]
            self.elementArea[i]=0.5*((a[0]*b[1]-a[1]*b[0])**2+\
                    (a[0]*b[2]-a[2]*b[0])**2+\
                    (a[1]*b[2]-a[2]*b[1])**2)**0.5
        #CALCULATE ELEMENT AREA################################################
        
    def writeBCfile(self,file):
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
        __lineV=list() #auxiliar array to store pumping intervals (stress periods)
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
                    __f2=aux_f.readDict(__file)
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
                    if __ID==i:
                        __l=((self.pointMatrix[self.bc[j][3]-1][0]-self.pointMatrix[self.bc[j][2]-1][0])**2+\
                             (self.pointMatrix[self.bc[j][3]-1][1]-self.pointMatrix[self.bc[j][2]-1][1])**2)**0.5
                        __p+=__l
            __lineP.append(__p)
        ###Loop for getting the perimeter of each Neumann BC line##############
        #######################################################################
        ##Assign flow value to each line in BC file############################
        __l=0 #Auxiliar variable to store the number of Neumann BC faces
        __aux_n=list()
        for i in range(len(__lineID)):
            ###Look for all 'elements' with the current line ID################
            for j in range(__nbc):
                #Check if 'element' line ID equals current line ID
                if self.bc[j][1]==__lineID[i]:
                    __line=str(self.bc[j][2]-1)+','+\
                    str(self.bc[j][3]-1) #Node on Neumann Face
                    if isinstance(__lineV[i],list): #Transient-state
                        __n=len(__lineV[i]) #Number of stress periods
                        __line+=','+str(__n)
                        for k in range(__n):
                            __line+=','+str(__lineV[i][k][0])+','+\
                            str(__lineV[i][k][1])+','+\
                            str(__lineV[i][k][2]/__lineP[i])
                    else: #Steady-state
                        __n=1 #Single stress period
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
        for i in range(len(self.pointMatrix)):
            __d=((self.pointMatrix[i][0]-x)**2+\
                 (self.pointMatrix[i][1]-y)**2)**0.5
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
        for i in range(len(self.cMatrix)):
            __n1=self.cMatrix[i][0]
            __n2=self.cMatrix[i][1]
            __n3=self.cMatrix[i][2]
            __d=((self.pointMatrix[__n1][0]-x)**2+(self.pointMatrix[__n1][1]-y)**2)**0.5+\
            ((self.pointMatrix[__n2][0]-x)**2+(self.pointMatrix[__n2][1]-y)**2)**0.5+\
            ((self.pointMatrix[__n3][0]-x)**2+(self.pointMatrix[__n3][1]-y)**2)**0.5
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
    nNodes = # of nodes             | nElements = # of elements
    k_e = Hydraulic Conductivity    | ss_e = Specific Storage
    bc_d = Dirichlet BC's           | bc_n = Neumann BC's
    h_n = Hydraulic Head            | q = Pumping/Recharge Well Rate
    """
    def __init__(self,config,nElements,nNodes):
        #INITIALIZE ALL THE ARRAYS USED FOR STORING PHYSICAL PROPERTIES########
        ##INITIALIZE HYDRAULIC CONDUCTIVITY ARRAY##############################
        self.k_e=[1]*nElements #Heterogeneity
        for i in range(nElements): #Anisotropy
            self.k_e[i]=[float(config['Kx']),float(config['Ky'])] #Hydraulic Conductivity Array (nElements x 2)
        #######################################################################
        ##INITIALIZE SPECIFIC STORAGE ARRAY####################################
        self.ss_e=[float(config['Ss'])]*nElements
        #######################################################################
        ##INITIALIZE WELL CONDITIONS ARRAY#####################################
        self.q=dict()
        #######################################################################
        ##INITIALIZE HYDRAULIC HEAD ARRAY######################################
        self.steady=config['steady'] #Check if we are working with transient-state or steady-state flow
        if self.steady==True: #Steady-state flow
            self.h_n=[0]*nNodes
        else: #Transient-state flow
            ##"ADD" TEMPORAL AXIS IF MEDIUM IS IN TRANSIENT-STATE##############
            ##READ TIME PARAMETERS#############################################
            self.dt=config['dt']
            self.timeSteps=int(config['total_time']/config['dt'])+1
            ##MODIFY HYDRAULIC HEAD ARRAY FOR TRANSIENT-STATE FLOW#############
            self.h_n=[0]*2 #To store 'previous' and 'current' states
            for i in range(2):
                self.h_n[i]=[0]*nNodes
        #######################################################################
        #INITIALIZE ALL THE ARRAYS USED FOR STORING BOUNDARY CONDITIONS########
        self.bc_d=list()
        self.bc_n=list()
        self.well_asBC=list()
        #######################################################################
        #INITIALIZE DICTIONARY FOR STORING LAYER/ZONE PROPERTIES###############
        ##INITIALIZE AUXILIAR ARRAY FOR STORING NODE IDS
        __elements=[0]*nElements
        for i in range(nElements):
            __elements[i]=i
        ##DEFINE DICTIONARY WITH ZONE PROPERTIES
        self.zones={0:{'Ss':float(config['Ss']),
                       'Kx':float(config['Kx']),
                       'Ky':float(config['Ky']),
                       'elements':__elements}}
        #######################################################################
        
    #READ BOUNDARY CONDITIONS##################################################
    def readBCfile(self,file):
        #READ BC FILE##########################################################
        __f=open(aux_f.path+'/'+file) #Open the given file in 'read' mode
        
        #DIRICHLET BC's########################################################
        for i in range(2):
            __f.readline() #Skip 2 lines
        __n=int(__f.readline().rstrip('\n')) #Read number of Dirichlet nodes
        self.bc_d=[0]*__n
        ##READ DIRICHLET BC's FROM FILE########################################
        for i in range(__n):
            __aux=__f.readline().rstrip('\n').split(',') #Read and split at ','
            self.bc_d[i]=[int(__aux[0]),float(__aux[1])]
        #######################################################################
        #######################################################################
        
        #NEUMANN BC's##########################################################
        __f.readline() #Skip 1 line
        __n=int(__f.readline().rstrip('\n')) #Read number of Neumann faces
        self.bc_n=[0]*__n
        ##READ NEUMANN BC's FROM FILE##########################################
        for i in range(__n):
            __aux=__f.readline().rstrip('\n').split(',')
            self.bc_n[i]=[int(__aux[0]),int(__aux[1]),float(__aux[2])]
        #######################################################################
        #######################################################################
        
        #WELL AS NEUMANN BC's##################################################
        __f.readline() #Skip 1 line
        __n=int(__f.readline().rstrip('\n')) #No. of nodes on the BC
        self.well_asBC=[0]*__n
        if self.steady==True: #Steady-state
            for i in range(__n):
                __aux=__f.readline().rstrip('\n').split(',')
                self.well_asBC[i]=[int(__aux[0]),int(__aux[1]),float(__aux[3])]
        else: #Transient-state
            for i in range(__n):
                __aux=__f.readline().rstrip('\n').split(',')
                #__aux[2] has the number of pumping intervals
                __v=[0]*__aux
                for j in range(int(__aux[2])):
                    __v[i]=[float(__aux[3+j*3]),float(__aux[4+j*3]),float(__aux[5+j*3])]
                self.well_asBC[i]=[int(__aux[0]),int(__aux[1]),int(__aux[2]),__v]
        #WELL AS NEUMANN BC's##################################################
        
        __f.close()
        #READ BC FILE##########################################################
    #READ BOUNDARY CONDITIONS##################################################
        
    #READ INITIAL CONDITIONS###################################################
    def read_h0(self,file):
        #READS AN HHD FILE AND ASSIGNS IT TO HH_0##############################
        __f=open(aux_f.path+'/'+file)
        
        __n=int(__f.readline()) #Read number of nodes
        for i in range(__n):
            __aux=__f.readline().split()
            self.h_n[0][i]=float(__aux[1])
        
        __f.close()
        #READS AN HHD FILE AND ASSIGNS IT TO HH_0##############################
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
    mesh = Mesh Object          | phyMed = Physical Medium Object
    massMatrix = Mass Matrix    | stiffnessMatrix = Stiffnes Matrix
    loadVector = Load Vector    |
    """  
    def __init__(self,init):
        self.type='FEM'
        #READ INIT FILE########################################################
        self.config=aux_f.readDict(init)
        #READ INIT FILE########################################################
        #CREATE MESH AND PHYMED OBJECTS########################################
        self.mesh=mesh(self.config['mesh_file'])
        self.phyMed=phymed(self.config,self.mesh.nElements,self.mesh.nNodes)
        #CREATE MESH AND PHYMED OBJECTS########################################
        #INITIALIZE MATRICES ARRAYS
        self.massMatrix=list() #INITIALIZE ARRAY FOR MASS MATRIX
        self.stiffnessMatrix=list() #INITIALIZE ARRAY FOR STIFFNESS MATRIX
        #INITIALIZE ARRAY FOR LOAD VECTOR
        self.loadVector=list()
        #######################################################################
        
###############################################################################
#####ASSEMBLE MATRICES AND VECTORS#############################################
###############################################################################
    ##ASSEMBLE MASS MATRIX#####################################################
    def __assembleMassMatrix(self):
        #INITIALIZE ARRAY FOR MASS MATRIX
        self.massMatrix=[0]*self.mesh.nNodes
        for i in range(self.mesh.nNodes):
            self.massMatrix[i]=[0]*self.mesh.nNodes
        #############################
        for k in range(self.mesh.nElements):
            ##CALCULATE LOCAL ELEMENT MASS MATRIX##############################
            for i in range(3): #Local Element Mass Matrix Index 'i'
                for j in range(3): #Local Element Mass Matrix Index 'j'
                    if i==j: #Main diagonal
                        __lmm_ij=self.mesh.elementArea[k]*self.phyMed.ss_e[k]/6
                    else:
                        __lmm_ij=self.mesh.elementArea[k]*self.phyMed.ss_e[k]/12
                    (global_i,global_j)=self.__local2global(k,i,j)
                    self.massMatrix[global_i][global_j]+=__lmm_ij #Add to Global Mass Matrix
            ###################################################################
    ##ASSEMBLE MASS MATRIX#####################################################
    
    ##ASSEMBLE STIFFNESS MATRIX################################################
    def __assembleStiffnessMatrix(self):
        #INITIALIZE ARRAY FOR STIFFNESS MATRIX
        self.stiffnessMatrix=[0]*self.mesh.nNodes
        for i in range(self.mesh.nNodes):
            self.stiffnessMatrix[i]=[0]*self.mesh.nNodes
        ##############################
        for k in range(self.mesh.nElements):
            ##CALCULATE LOCAL ELEMENT STIFFNESS MATRIX#########################
            ###GET ELEMENT THE HAT FUNCTION FOR ELEMENT K######################
            (__a,__b,__c)=self.hat_phi(k)
            ###################################################################
            for i in range(3): #Local Element Stiffness Matrix Index 'i'
                for j in range(3): #Local Element Stiffness Matrix Index 'j'
                    __lsm_ij=(self.phyMed.k_e[k][0]*(__b[i]*__b[j])+\
                    self.phyMed.k_e[k][1]*(__c[i]*__c[j]))*self.mesh.elementArea[k]
                    (global_i,global_j)=self.__local2global(k,i,j)
                    self.stiffnessMatrix[global_i][global_j]+=__lsm_ij #Add to Global Stiffness Matrix
            ###################################################################
    ##ASSEMBLE STIFFNESS MATRIX################################################
    
    ##ASSEMBLE LOAD VECTOR#####################################################
    def __assembleLoadVector(self):
        #Steady-state##########################################################
        if self.phyMed.steady==True:
            #INITIALIZE ARRAY FOR LOAD VECTOR
            self.loadVector=[0]*self.mesh.nNodes
            #ASSIGN VALUES
            for i in self.phyMed.q:
                self.loadVector[i]=-self.phyMed.q[i]
        #Steady-state##########################################################
        #Transient-state#######################################################
        else:
            #INITIALIZE ARRAY FOR LOAD VECTOR
            self.loadVector=[0]*self.phyMed.timeSteps
            for i in range(self.phyMed.timeSteps):
                self.loadVector[i]=[0]*self.mesh.nNodes
            #ASSIGN VALUES
            for i in self.phyMed.q:
                nStressPeriods=self.phyMed.q[i][0] #number of stress periods
                for j in range(nStressPeriods):
                    t0=self.phyMed.q[i][1][j][0] #Start of stress period 'j'
                    tf=self.phyMed.q[i][1][j][1] #End of stress period 'j'
                    v=self.phyMed.q[i][1][j][2] #Pumping rate of stress period 'j'
                    for k in range(t0,tf):
                        self.loadVector[k][i]=-v
        #Transient-state#######################################################
    ##ASSEMBLE LOAD VECTOR#####################################################
    
    ###########################################################################
    #MISC FUNCTIONS############################################################
    ###########################################################################
    
    ##LOCAL2GLOBAL#############################################################
    def __local2global(self,k,i,j):
        #'k'=element index, 'i'=local element index i, 'j'=local element index j
        global_i=self.mesh.cMatrix[k][i]
        global_j=self.mesh.cMatrix[k][j]
        return (global_i,global_j)
    ##LOCAL2GLOBAL#############################################################
    
    ##COMPUTE HAT FUNCTION OF ELEMENT 'K'######################################
    def hat_phi(self,k):
        __n1=self.mesh.cMatrix[k][0] #Node 1
        __n2=self.mesh.cMatrix[k][1] #Node 2
        __n3=self.mesh.cMatrix[k][2] #Node 3
        __j=2*self.mesh.elementArea[k] #Jacobian
        a=[(self.mesh.pointMatrix[__n2][0]*self.mesh.pointMatrix[__n3][1]-self.mesh.pointMatrix[__n3][0]*self.mesh.pointMatrix[__n2][1])/__j,\
           (self.mesh.pointMatrix[__n3][0]*self.mesh.pointMatrix[__n1][1]-self.mesh.pointMatrix[__n1][0]*self.mesh.pointMatrix[__n3][1])/__j,\
           (self.mesh.pointMatrix[__n1][0]*self.mesh.pointMatrix[__n2][1]-self.mesh.pointMatrix[__n2][0]*self.mesh.pointMatrix[__n1][1])/__j]
        b=[(self.mesh.pointMatrix[__n2][1]-self.mesh.pointMatrix[__n3][1])/__j,\
           (self.mesh.pointMatrix[__n3][1]-self.mesh.pointMatrix[__n1][1])/__j,\
           (self.mesh.pointMatrix[__n1][1]-self.mesh.pointMatrix[__n2][1])/__j]
        c=[(self.mesh.pointMatrix[__n3][0]-self.mesh.pointMatrix[__n2][0])/__j,\
           (self.mesh.pointMatrix[__n1][0]-self.mesh.pointMatrix[__n3][0])/__j,\
           (self.mesh.pointMatrix[__n2][0]-self.mesh.pointMatrix[__n1][0])/__j]
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
    def __applyDirichletBC(self):
        #UNIT CONVERSION#######################################################
        uF=1
        if self.config['length'].lower()!=self.config['head'].lower(): #Check if conversion is required
            uF=aux_f.unitDict[self.config['length'].lower()]/aux_f.unitDict[self.config['head'].lower()] #Get conversion factor
        #UNIT CONVERSION#######################################################
        #Sort the 'bc_d' array for an easier time 'popping' nodes.
        self.phyMed.bc_d.sort()
        __counter=0 #Auxiliar counter used for keeping track of how many nodes have been removed
        #Steady-state##########################################################
        if self.phyMed.steady==True:
            for k in range(len(self.phyMed.bc_d)):
                __node=self.phyMed.bc_d[k][0]-__counter
                ##REMOVE NODE FROM HYDRAULIC HEAD VECTOR#######################
                del self.phyMed.h_n[__node]
                ##REMOVE NODE FROM LOAD VECTOR#################################
                del self.loadVector[__node]
                ##REMOVE ROW FROM STIFFNESS MATRIX#############################
                del self.stiffnessMatrix[__node]
                ##ADDITIONAL OPERATIONS########################################
                for i in range(len(self.loadVector)):
                    self.loadVector[i]-=self.stiffnessMatrix[i][__node]*\
                    self.phyMed.bc_d[k][1]*uF
                    ###REMOVE COLUMN FROM STIFFNESS MATRIX#####################
                    del self.stiffnessMatrix[i][__node]
                __counter+=1
        #Steady-state##########################################################
        #Transient-state#######################################################
        if self.phyMed.steady==False:
            for k in range(len(self.phyMed.bc_d)):
                __node=self.phyMed.bc_d[k][0]-__counter
                ##REMOVE NODE FROM HYDRAULIC HEAD VECTOR#######################
                for i in range(self.phyMed.timeSteps):
                    del self.phyMed.h_n[i][__node]
                ##REMOVE NODE FROM LOAD VECTOR#################################
                for i in range(self.phyMed.timeSteps):
                    del self.loadVector[i][__node]
                ##REMOVE ROW FROM STIFFNESS MATRIX#############################
                del self.stiffnessMatrix[__node]
                ##REMOVE ROW FROM MASS MATRIX##################################
                del self.massMatrix[__node]
                ##ADDITIONAL OPERATIONS########################################
                for i in range(len(self.loadVector[0])):
                    for j in range(self.phyMed.timeSteps):
                        self.loadVector[j][i]-=self.stiffnessMatrix[i][__node]*\
                        self.phyMed.bc_d[k][1]*uF
                    ###REMOVE COLUMN FROM STIFFNESS MATRIX#####################
                    del self.stiffnessMatrix[i][__node]
                    ##REMOVE COLUMN FROM MASS MATRIX###########################
                    del self.massMatrix[i][__node]
                __counter+=1
        #Transient-state#######################################################
    ##APPLY DIRICHLET BC's#####################################################
    
    ##APPLY NEUMANN BC's#######################################################
    def __applyNeumannBC(self):
        #UNIT CONVERSION#######################################################
        uF=1
        if self.config['length'].lower()!=self.config['head'].lower(): #Check if conversion is required
            uF=(aux_f.unitDict[self.config['length'].lower()])/(aux_f.unitDict[self.config['head'].lower()]) #Get conversion factor            
        #UNIT CONVERSION#######################################################
        #Steady-state##########################################################
        if self.phyMed.steady==True:
            for k in range(len(self.phyMed.bc_n)):
                __l=((self.mesh.pointMatrix[self.phyMed.bc_n[k][1]][0]-self.mesh.pointMatrix[self.phyMed.bc_n[k][0]][0])**2+\
                     (self.mesh.pointMatrix[self.phyMed.bc_n[k][1]][1]-self.mesh.pointMatrix[self.phyMed.bc_n[k][0]][1])**2)**0.5
                self.loadVector[self.phyMed.bc_n[k][0]]+=self.phyMed.bc_n[k][2]*__l*self.config['aq_thickness']*uF/2
                self.loadVector[self.phyMed.bc_n[k][1]]+=self.phyMed.bc_n[k][2]*__l*self.config['aq_thickness']*uF/2
        #Steady-state##########################################################
        #Transient-state#######################################################
        else: #Transient-state
            for k in range(len(self.phyMed.bc_n)):
                __l=((self.mesh.pointMatrix[self.phyMed.bc_n[k][1]][0]-self.mesh.pointMatrix[self.phyMed.bc_n[k][0]][0])**2+\
                     (self.mesh.pointMatrix[self.phyMed.bc_n[k][1]][1]-self.mesh.pointMatrix[self.phyMed.bc_n[k][0]][1])**2)**0.5
                for i in range(self.phyMed.timeSteps):
                    self.loadVector[i][self.phyMed.bc_n[k][0]]+=self.phyMed.bc_n[k][2]*__l*self.config['aq_thickness']*uF/2 #RECORDATORIO: AGREGAR PASO TEMPORAL UNA VEZ QUE EL MÉTODO FUNCIONE
                    self.loadVector[i][self.phyMed.bc_n[k][1]]+=self.phyMed.bc_n[k][2]*__l*self.config['aq_thickness']*uF/2 #RECORDATORIO: AGREGAR PASO TEMPORAL UNA VEZ QUE EL MÉTODO FUNCIONE
        #Transient-state#######################################################
    ##APPLY NEUMANN BC's#######################################################
    
    ##REBUILD HH ARRAY#########################################################
    def __rebuild_h(self):
        #UNIT CONVERSION#######################################################
        if self.config['length'].lower()!=self.config['head'].lower(): #Check if conversion is required
            uF=aux_f.unitDict[self.config['head'].lower()]/aux_f.unitDict[self.config['length'].lower()] #Get conversion factor
            for i in range(len(self.phyMed.h_n)):
                self.phyMed.h_n[i]*=uF
        #UNIT CONVERSION#######################################################
        #Sort the 'bc_d' array for an easier time 'restoring' nodes.
        self.phyMed.bc_d.sort()
        ##ADD BACK THE DIRICHLET BC'S##########################################
        if self.phyMed.steady==True: #Steady-state
            for i in range(len(self.phyMed.bc_d)):
                self.phyMed.h_n.insert(self.phyMed.bc_d[i][0],self.phyMed.bc_d[i][1])
        else: #Transient-state
            for k in range(self.phyMed.timeSteps):
                for i in range(len(self.phyMed.bc_d)):
                    self.phyMed.h_n[k].insert(self.phyMed.bc_d[i][0],self.phyMed.bc_d[i][1])
        ##ADD BACK THE DIRICHLET BC'S##########################################
    ##REBUILD HH ARRAY#########################################################  
###############################################################################
#####APPLY BOUNDARY CONDITIONS#################################################
###############################################################################

###############################################################################
#####MODIFY PHYSICAL MEDIUM####################################################
###############################################################################
    ##ADD WELLS################################################################
    def addWell(self,file): #Add well as point source over an element
        #UNIT CONVERSION#######################################################
        uF=1
        if self.config['length'].lower()!=self.config['head'].lower(): #Check if conversion is required
            uF=(aux_f.unitDict[self.config['length'].lower()]**3)/(aux_f.unitDict[self.config['head'].lower()]**3) #Get conversion factor            
        #UNIT CONVERSION#######################################################
        ##RETRIEVE WELL INFORMATION FROM FILE##################################
        __f=aux_f.readDict(file)
        ##RETRIEVE WELL INFORMATION FROM FILE##################################
        __element=self.mesh.coord2element(__f['x'],__f['y'])
        #Retrieve nodes for interpolation
        __n1=self.mesh.cMatrix[__element][0] #Node 1 index
        __n2=self.mesh.cMatrix[__element][1] #Node 2 index
        __n3=self.mesh.cMatrix[__element][2] #Node 3 index
        ##GET INTERPOLATION ELEMENT########################
        ##GET NATURAL SYSTEM VALUES########################
        __x21=self.mesh.pointMatrix[__n2][0]-self.mesh.pointMatrix[__n1][0]
        __y21=self.mesh.pointMatrix[__n2][1]-self.mesh.pointMatrix[__n1][1]
        __x31=self.mesh.pointMatrix[__n3][0]-self.mesh.pointMatrix[__n1][0]
        __y31=self.mesh.pointMatrix[__n3][1]-self.mesh.pointMatrix[__n1][1]
        __xi=(__y31*(__f['x']-self.mesh.pointMatrix[__n1][0])-__x31*(__f['y']-self.mesh.pointMatrix[__n1][1]))/\
        (__x21*__y31-__x31*__y21)
        __eta=(-__y21*(__f['x']-self.mesh.pointMatrix[__n1][0])+__x21*(__f['y']-self.mesh.pointMatrix[__n1][1]))/\
        (__x21*__y31-__x31*__y21)
        __N=[1-__xi-__eta,__xi,__eta]
        ##GET NATURAL SYSTEM VALUES########################
        ##COMPUTE DRAWDOWN VALUES##########################
        ##ASSIGN FLOW VALUE TO 'Q' ARRAY#######################################
        if self.phyMed.steady==True: #Steady-state
            for i in range(3):
                __node=self.mesh.cMatrix[__element][i]
                v=-__f['rate']*__N[i]*uF/\
                (self.config['aq_thickness'])
                #Check if the node already exists in the 'q' array
                if __node in self.phyMed.q: #Node already has a well on it
                    self.phyMed.q[__node]+=v
                else: #New well
                    self.phyMed.q[__node]=v
        else: #Transient-state
            __n=len(__f['rate']) #NUMBER OF PUMPING INTERVALS
            for i in range(__n):
                #IDENTIFY TIME-STEPS WITHIN THE PUMPING INTERVAL###############
                #IDENTIFY 'k_min'
                for k in range(self.phyMed.timeSteps):
                    if k*self.phyMed.dt>=__f['rate'][i][0]:
                        t_min=k
                        break
                #IDENTIFY 'k_max'
                for k in range(t_min,self.phyMed.timeSteps):
                    if k*self.phyMed.dt>=__f['rate'][i][1]:
                        t_max=k+1
                        break
                #ASSIGN FLOW VALUE AT INTERVAL 'i'#############################
                for j in range(3):
                    __node=self.mesh.cMatrix[__element][j]
                    offset=0
                    #Check if the node already exists in the 'q' array
                    if __node in self.phyMed.q: #Node already has a well on it
                        offset=self.phyMed.q[__node][0] #Get current number of stress periods
                        self.phyMed.q[__node][0]+=__n #Add number of stress periods to the current counter
                        for i in range(__n):
                            self.phyMed.q[__node][1].append(0) #Add placeholders to store the new stress periods
                    else: #New well
                        self.phyMed.q[__node]=[__n,[0]*__n]
                    v=-__f['rate'][i][2]*__N[j]*uF/\
                        (self.config['aq_thickness'])
                    self.phyMed.q[__node][1][i+offset][0]=t_min #Start of stress period 'i'
                    self.phyMed.q[__node][1][i+offset][1]=t_max #End of stress period 'i'
                    self.phyMed.q[__node][1][i+offset][2]=v #Pumping rate of stress period 'i'

    def addWell_onElement(self,file): #Add well as constant function over an element
        #UNIT CONVERSION#######################################################
        uF=1
        if self.config['length'].lower()!=self.config['head'].lower(): #Check if conversion is required
            uF=(aux_f.unitDict[self.config['length'].lower()]**3)/(aux_f.unitDict[self.config['head'].lower()]**3) #Get conversion factor            
        #UNIT CONVERSION#######################################################
        ##RETRIEVE WELL INFORMATION FROM FILE##################################
        __f=aux_f.readDict(file)
        ##RETRIEVE WELL INFORMATION FROM FILE##################################
        __element=self.mesh.coord2element(__f['x'],__f['y'])
        ##ASSIGN FLOW VALUE TO 'Q' ARRAY#######################################
        if self.phyMed.steady==True: #Steady-state
            v=-__f['rate']*uF/\
            (3*self.config['aq_thickness'])
            for i in range(3):
                __node=self.mesh.cMatrix[__element][i]
                #Check if the node already exists in the 'q' array
                if __node in self.phyMed.q: #Node already has a well on it
                    self.phyMed.q[__node]+=v
                else: #New well
                    self.phyMed.q[__node]=v
        else: #Transient-state
            __n=len(__f['rate']) #NUMBER OF PUMPING INTERVALS
            for i in range(__n):
                #IDENTIFY TIME-STEPS WITHIN THE PUMPING INTERVAL###############
                #IDENTIFY 'k_min'
                for k in range(self.phyMed.timeSteps):
                    if k*self.phyMed.dt>=__f['rate'][i][0]:
                        t_min=k
                        break
                #IDENTIFY 'k_max'
                for k in range(t_min,self.phyMed.timeSteps):
                    if k*self.phyMed.dt>=__f['rate'][i][1]:
                        t_max=k+1
                        break
                #ASSIGN FLOW VALUE AT INTERVAL 'i'#############################
                v=-__f['rate'][i][2]*uF/\
                (3*self.config['aq_thickness'])
                for j in range(3):
                    __node=self.mesh.cMatrix[__element][j]
                    offset=0
                    #Check if the node already exists in the 'q' array
                    if __node in self.phyMed.q: #Node already has a well on it
                        offset=self.phyMed.q[__node][0] #Get current number of stress periods
                        self.phyMed.q[__node][0]+=__n #Add number of stress periods to the current counter
                        for i in range(__n):
                            self.phyMed.q[__node][1].append(0) #Add placeholders to store the new stress periods
                    else: #New well
                        self.phyMed.q[__node]=[__n,[0]*__n]
                    self.phyMed.q[__node][1][i+offset][0]=t_min #Start of stress period 'i'
                    self.phyMed.q[__node][1][i+offset][1]=t_max #End of stress period 'i'
                    self.phyMed.q[__node][1][i+offset][2]=v #Pumping rate of stress period 'i'
                        
    def addWell_onNode(self,file): #Add well as a point source over a node
        #UNIT CONVERSION#######################################################
        uF=1
        if self.config['length'].lower()!=self.config['head'].lower(): #Check if conversion is required
            uF=(aux_f.unitDict[self.config['length'].lower()]**3)/(aux_f.unitDict[self.config['head'].lower()]**3) #Get conversion factor            
        #UNIT CONVERSION#######################################################
        ##RETRIEVE WELL INFORMATION FROM FILE##################################
        __f=aux_f.readDict(file)
        ##RETRIEVE WELL INFORMATION FROM FILE##################################
        __node=self.mesh.coord2node(__f['x'],__f['y'])
        ##ASSIGN FLOW VALUE TO 'Q' ARRAY#######################################
        if self.phyMed.steady==True: #Steady-state
            v=-__f['rate']*uF/\
            (self.config['aq_thickness'])
            #Check if the node already exists in the 'q' array
            if __node in self.phyMed.q: #Node already has a well on it
                self.phyMed.q[__node]+=v
            else: #New well
                self.phyMed.q[__node]=v
        else: #Transient-state
            __n=len(__f['rate']) #NUMBER OF PUMPING INTERVALS
            offset=0
            #Check if the node already exists in the 'q' array
            if __node in self.phyMed.q: #Node already has a well on it
                offset=self.phyMed.q[__node][0] #Get current number of stress periods
                self.phyMed.q[__node][0]+=__n #Add number of stress periods to the current counter
                for i in range(__n):
                    self.phyMed.q[__node][1].append(0) #Add placeholders to store the new stress periods
            else: #New well
                self.phyMed.q[__node]=[__n,[0]*__n]
            for i in range(__n):
                #IDENTIFY TIME-STEPS WITHIN THE PUMPING INTERVAL###############
                #IDENTIFY 'k_min'
                for k in range(self.phyMed.timeSteps):
                    if k*self.phyMed.dt>=__f['rate'][i][0]:
                        t_min=k
                        break
                #IDENTIFY 'k_max'
                for k in range(t_min,self.phyMed.timeSteps):
                    if k*self.phyMed.dt>=__f['rate'][i][1]:
                        t_max=k+1
                        break
                #ASSIGN FLOW VALUE AT INTERVAL 'i'#############################
                v=-__f['rate'][i][2]*uF/\
                (self.config['aq_thickness'])
                self.phyMed.q[__node][1][i+offset][0]=t_min #Start of stress period 'i'
                self.phyMed.q[__node][1][i+offset][1]=t_max #End of stress period 'i'
                self.phyMed.q[__node][1][i+offset][2]=v #Pumping rate of stress period 'i'
    
    def addWell_asBC(self):
        #UNIT CONVERSION#######################################################
        uF=1
        if self.config['length'].lower()!=self.config['head'].lower(): #Check if conversion is required
            uF=(aux_f.unitDict[self.config['length'].lower()]**3)/(aux_f.unitDict[self.config['head'].lower()]**3) #Get conversion factor            
        #UNIT CONVERSION#######################################################
        if self.phyMed.steady==True: #Steady-state
            for k in range(len(self.phyMed.well_asBC)):
                __l=((self.mesh.pointMatrix[self.phyMed.well_asBC[k][1]][0]-self.mesh.pointMatrix[self.phyMed.well_asBC[k][0]][0])**2+\
                     (self.mesh.pointMatrix[self.phyMed.well_asBC[k][1]][1]-self.mesh.pointMatrix[self.phyMed.well_asBC[k][0]][1])**2)**0.5
                v=-self.phyMed.well_asBC[k][2]*__l*uF/\
                (2*self.config['aq_thickness'])
                for i in range(2):
                    __node=self.phyMed.well_asBC[k][i]
                    #Check if the node already exists in the 'q' array
                    if __node in self.phyMed.q: #Node already has a well on it
                        self.phyMed.q[__node]+=v
                    else: #New well
                        self.phyMed.q[__node]=v
        else: #Transient-state
            for k in range(len(self.phyMed.well_asBC)):
                __l=((self.mesh.pointMatrix[self.phyMed.well_asBC[k][1]][0]-self.mesh.pointMatrix[self.phyMed.well_asBC[k][0]][0])**2+\
                     (self.mesh.pointMatrix[self.phyMed.well_asBC[k][1]][1]-self.mesh.pointMatrix[self.phyMed.well_asBC[k][0]][1])**2)**0.5
                __n=self.phyMed.well_asBC[k][2]
                __v=self.phyMed.well_asBC[k][3]
                #IDENTIFY TIME-STEPS WITHIN THE PUMPING INTERVAL###############
                for i in range(__n):
                    #IDENTIFY 'k_min'
                    for j in range(self.phyMed.timeSteps):
                        if j*self.phyMed.dt>=__v[i][0]:
                            t_min=j
                            break
                    #IDENTIFY 'k_max'
                    for j in range(t_min,self.phyMed.timeSteps):
                        if j*self.phyMed.dt>=__v[i][1]:
                            t_max=j+1
                            break
                    for j in range(2):
                        __node=self.phyMed.well_asBC[k][j]
                        offset=0
                        #Check if the node already exists in the 'q' array
                        if __node in self.phyMed.q: #Node already has a well on it
                            offset=self.phyMed.q[__node][0] #Get current number of stress periods
                            self.phyMed.q[__node][0]+=__n #Add number of stress periods to the current counter
                            for i in range(__n):
                                self.phyMed.q[__node][1].append(0) #Add placeholders to store the new stress periods
                        else: #New well
                            self.phyMed.q[__node]=[__n,[0]*__n]
                        v=-__v[i][2]*__l*uF/\
                        (2*self.config['aq_thickness'])
                        self.phyMed.q[__node][1][i+offset][0]=t_min #Start of stress period 'i'
                        self.phyMed.q[__node][1][i+offset][1]=t_max #End of stress period 'i'
                        self.phyMed.q[__node][1][i+offset][2]=v #Pumping rate of stress period 'i'
    ##ADD WELLS################################################################
    
    #ZONE RELATED FUNCTIONS####################################################
    ##DEFINE NEW ZONE##########################################################
    def newZone(self,zone_id,Ss:float,Kx:float,Ky:float,x_min:float,x_max:float,y_min:float,y_max:float): #Defines a new zone in the domain as a rectangle
        if zone_id not in self.phyMed.zones: #CHECK IF THE ZONE EXISTS
            #IDENTIFY ALL ELEMENTS INSIDE THE ZONE#############################
            __ID=list() #Auxiliar list for storing the ID values of the elements inside the zone
            for i in range(self.mesh.nElements):
                #Get element centroid
                __cx=(self.mesh.pointMatrix[self.mesh.cMatrix[i][0]][0]+\
                      self.mesh.pointMatrix[self.mesh.cMatrix[i][1]][0]+\
                      self.mesh.pointMatrix[self.mesh.cMatrix[i][2]][0])/3
                __cy=(self.mesh.pointMatrix[self.mesh.cMatrix[i][0]][1]+\
                      self.mesh.pointMatrix[self.mesh.cMatrix[i][1]][1]+\
                      self.mesh.pointMatrix[self.mesh.cMatrix[i][2]][1])/3
                if x_min<=__cx<=x_max and y_min<=__cy<=y_max:
                    __ID.append(i)
            #IDENTIFY ALL ELEMENTS INSIDE THE ZONE#############################
            #ZONE ASSIGNATION LOOP#############################################
            for i in __ID:
                #Assign zone values
                self.phyMed.k_e[i]=[Kx,Ky]
                self.phyMed.ss_e[i]=Ss
                #Look for the element in all other zones and remove it from their lists
                for k in self.phyMed.zones:
                    if i in self.phyMed.zones[k]['elements']: #Check if node is in a zone 'k'
                        __e_ID=self.phyMed.zones[k]['elements'].index(i) #Get ID within 'nodes' list
                        del self.phyMed.zones[k]['elements'][__e_ID] #Remove node from zone 'k'
            self.phyMed.zones[zone_id]={'Ss':Ss,
                             'Kx':Kx,
                             'Ky':Ky,
                             'elements':__ID}
            #ZONE ASSIGNATION LOOP#############################################
        else:
            print('Zone {0} already exists'.format(zone_id))
    ##DEFINE NEW ZONE##########################################################
    
    ##APPEND NEW RECTANGLE TO EXISTING ZONE####################################
    def expandZone(self,zone_id,x_min:float,x_max:float,y_min:float,y_max:float): #Adds a new rectangle to an existing zone
        if zone_id in self.phyMed.zones: #CHECK IF THE ZONE EXISTS
            #IDENTIFY ALL ELEMENTS INSIDE THE ZONE#############################
            __ID=list() #Auxiliar list for storing the ID values of the elements inside the zone
            for i in range(self.mesh.nElements):
                #Get element centroid
                __cx=(self.mesh.pointMatrix[self.mesh.cMatrix[i][0]][0]+\
                      self.mesh.pointMatrix[self.mesh.cMatrix[i][1]][0]+\
                      self.mesh.pointMatrix[self.mesh.cMatrix[i][2]][0])/3
                __cy=(self.mesh.pointMatrix[self.mesh.cMatrix[i][0]][1]+\
                      self.mesh.pointMatrix[self.mesh.cMatrix[i][1]][1]+\
                      self.mesh.pointMatrix[self.mesh.cMatrix[i][2]][1])/3
                if x_min<=__cx<=x_max and y_min<=__cy<=y_max:
                    __ID.append(i)
            #IDENTIFY ALL ELEMENTS INSIDE THE ZONE#############################
            #ZONE ASSIGNATION LOOP#############################################
            for i in __ID:
                #Assign zone values
                self.phyMed.k_e[i]=[self.phyMed.zones[zone_id]['Kx'],self.phyMed.zones[zone_id]['Ky']]
                self.phyMed.ss_e[i]=self.phyMed.zones[zone_id]['Ss']
                #Look for the element in all other zones and remove it from their lists
                for k in self.phyMed.zones:
                    if i in self.phyMed.zones[k]['elements']: #Check if node is in a zone 'k'
                        __e_ID=self.phyMed.zones[k]['elements'].index(i) #Get ID within 'nodes' list
                        del self.phyMed.zones[k]['elements'][__e_ID] #Remove node from zone 'k'
                #Add elements to zone
                if i not in self.phyMed.zones[zone_id]['elements']:
                    self.phyMed.zones[zone_id]['elements'].append(i)
            #ZONE ASSIGNATION LOOP#############################################
        else:
            print("Zone {0} doesn't exist".format(zone_id))
    ##APPEND NEW RECTANGLE TO EXISTING ZONE####################################
    
    ##UPDATE ZONE PARAMETERS###################################################
    def updateZone(self,zone_id,Ss:float,Kx:float,Ky:float): #Update zone properties
        if zone_id in self.phyMed.zones: #CHECK IF THE ZONE EXISTS
            self.phyMed.zones[zone_id]['Ss']=Ss #Update Specific Storage value
            self.phyMed.zones[zone_id]['Kx']=Kx #Update Kx value
            self.phyMed.zones[zone_id]['Ky']=Ky #Update Ky value
            for i in self.phyMed.zones[zone_id]['elements']:
                self.phyMed.k_e[i]=[Kx,Ky] #Update element value
                self.phyMed.ss_e[i]=Ss #Update element value
        else:
            print("Zone {0} doesn't exist".format(zone_id))
    ##UPDATE ZONE PARAMETERS###################################################
    
    ##DELETE ZONE##############################################################
    def deleteZone(self,zone_id):
        if zone_id in self.phyMed.zones: #CHECK IF THE ZONE EXISTS
            for i in self.phyMed.zones[zone_id]['elements']:
                self.phyMed.zones[0]['elements'].append(i) #Add element to default zone
                self.phyMed.k_e[i]=[self.phyMed.zones[0]['Kx'],self.phyMed.zones[0]['Ky']] #Update element value
                self.phyMed.ss_e[i]=self.phyMed.zones[0]['Ss'] #Update element value
            del self.phyMed.zones[zone_id]
        else:
            print("Zone {0} doesn't exist".format(zone_id))
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
        if self.phyMed.steady==True: #Steady-state
            print('ASSEMBLING GLOBAL MATRICES...')
            ##ASSEMBLE GLOBAL MATRICES#########################################
            self.__assembleStiffnessMatrix() #ASSEMBLE STIFFNESS MATRIX
            self.__assembleLoadVector() #ASSEMBLE LOAD VECTOR
            ###################################################################
            
            print('APPLYING BCs...')
            ##APPLY BC'S#######################################################
            self.__applyNeumannBC()
            self.__applyDirichletBC()
            ##APPLY BC'S#######################################################
            
            print('COMPUTING HYDRAULIC HEAD DISTRIBUTION...')
            ##GET HYDRAULIC HEAD DISTRIBUTION##################################    
            __aux_hh=np.matmul(np.linalg.inv(self.stiffnessMatrix),self.loadVector)
            self.phyMed.h_n=__aux_hh.tolist()
            self.__rebuild_h()
            ##GET HYDRAULIC HEAD DISTRIBUTION##################################
        else: #Transient-state
            print('ASSEMBLING GLOBAL MATRICES...')
            ##ASSEMBLE GLOBAL MATRICES#########################################
            self.__assembleMassMatrix() #ASSEMBLE MASS MATRIX
            self.__assembleStiffnessMatrix() #ASSEMBLE STIFFNESS MATRIX
            self.__assembleLoadVector() #ASSEMBLE LOAD VECTOR
            ###################################################################
            
            print('APPLYING BCs...')
            ##APPLY BC'S#######################################################
            self.__applyNeumannBC()
            self.__applyDirichletBC()
            ##APPLY BC'S#######################################################
            
            print('ASSEMBLING AUXILIAR MATRICES...')
            ##GET AUXILIAR MATRICES############################################
            ##__A=MASS+STIFFNESS*dt
            ###INITIALIZE ARRAY
            __nNodes=len(self.massMatrix)
            __A=[0]*__nNodes
            for i in range(__nNodes):
                __A[i]=[0]*__nNodes
            ###ASSIGN VALUES
            #DIVIDE MASS MATRIX BY DT
            for i in range(__nNodes):
                for j in range(__nNodes):
                    self.massMatrix[i][j]/=self.phyMed.dt
            for i in range(__nNodes):
                for j in range(__nNodes):
                    __A[i][j]=self.massMatrix[i][j]+self.stiffnessMatrix[i][j]#*self.phyMed.dt
            __Ainv=np.linalg.inv(__A)
            ##__V=LOAD*dt+MASS*hh_t-1
            ###INITIALIZE ARRAY
            __V=[0]*self.phyMed.timeSteps
            for i in range(self.phyMed.timeSteps):
                __V[i]=[0]*__nNodes
            ###ASSIGN VALUES
            for i in range(self.phyMed.timeSteps):
                for j in range(__nNodes):
                        __V[i][j]=self.loadVector[i][j]#*self.phyMed.dt
            ##GET AUXILIAR MATRICES############################################
            
            print('SOLVING FOR EACH TIME STEP...')
            for k in range(self.phyMed.timeSteps-1): #TIME-LOOP
                __step=k+1
                print(str(100*__step/self.phyMed.timeSteps)+'%...')
                
                ##UPDATE AUXILIAR MATRICES#####################################
                __M=np.matmul(self.massMatrix,self.phyMed.h_n[__step-1])
                for i in range(__nNodes):
                    __V[__step][i]+=__M[i]
                ##UPDATE AUXILIAR MATRICES#####################################
                
                ##GET HYDRAULIC HEAD DISTRIBUTION##############################
                ###CORREGIR
                print('COMPUTING HYDRAULIC HEAD DISTRIBUTION...')
                __aux_hh=np.matmul(__Ainv,__V[__step])
                self.phyMed.h_n[__step]=__aux_hh.tolist()
                ##GET HYDRAULIC HEAD DISTRIBUTION##############################
            ###HH VALUES#######################################################
            self.__rebuild_h()
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
    def exportResults(self,file,step=None):
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
            __f.write('{0}\n'.format(self.mesh.nNodes)) #Number of nodes
            for i in range(self.mesh.nNodes):
                __f.write('{0} {1}\n'.format(i+1,self.phyMed.h_n[i]) )
        else:
            __f.write('{0:16f}\n'.format(step*self.phyMed.dt)) #Simulation time = step * dt
            __f.write('{0}\n'.format(3)) #Three integer tags:
            __f.write('{0}\n'.format(step)) #Current step
            __f.write('{0}\n'.format(1)) #Scalar field
            __f.write('{0}\n'.format(self.mesh.nNodes)) #Number of nodes
            for i in range(self.mesh.nNodes):
                __f.write('{0} {1}\n'.format(i+1,self.phyMed.h_n[step][i]) )
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
        __n1=self.mesh.cMatrix[__n][0] #Node 1 index
        __n2=self.mesh.cMatrix[__n][1] #Node 2 index
        __n3=self.mesh.cMatrix[__n][2] #Node 3 index
        ##GET INTERPOLATION ELEMENT########################
        ##GET NATURAL SYSTEM VALUES########################
        __x21=self.mesh.pointMatrix[__n2][0]-self.mesh.pointMatrix[__n1][0]
        __y21=self.mesh.pointMatrix[__n2][1]-self.mesh.pointMatrix[__n1][1]
        __x31=self.mesh.pointMatrix[__n3][0]-self.mesh.pointMatrix[__n1][0]
        __y31=self.mesh.pointMatrix[__n3][1]-self.mesh.pointMatrix[__n1][1]
        __xi=(__y31*(x-self.mesh.pointMatrix[__n1][0])-__x31*(y-self.mesh.pointMatrix[__n1][1]))/\
        (__x21*__y31-__x31*__y21)
        __eta=(-__y21*(x-self.mesh.pointMatrix[__n1][0])+__x21*(y-self.mesh.pointMatrix[__n1][1]))/\
        (__x21*__y31-__x31*__y21)
        __N1=1-__xi-__eta
        __N2=__xi
        __N3=__eta
        ##GET NATURAL SYSTEM VALUES########################
        ##COMPUTE DRAWDOWN VALUES##########################
        __h=[0]*self.phyMed.timeSteps #Auxiliar array to store head values
        drawdown=[0]*self.phyMed.timeSteps #Initialize array
        for k in range(self.phyMed.timeSteps):
            __h[k]=self.phyMed.h_n[k][__n1]*__N1+\
            self.phyMed.h_n[k][__n2]*__N2+\
            self.phyMed.h_n[k][__n3]*__N3
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
            __f.write('{0}\n'.format(self.phyMed.timeSteps)) #Write number of time-steps
            __f.write('#OBSERVATION POINT COORDINATES\n')
            __f.write('{0:10f} {1:10f}\n'.format(x,y)) #Write observation point intended coordinates
            __f.write('#TIME-STEP VALUES\n')
            for i in range(self.phyMed.timeSteps):
                __f.write('{0:.7E} '.format(i*self.phyMed.dt))
            __f.write('\n')
            __f.write('#DRADOWN VALUES\n')
            for i in range(self.phyMed.timeSteps):
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
            __n1=self.mesh.cMatrix[__n][0] #Node 1 index
            __n2=self.mesh.cMatrix[__n][1] #Node 2 index
            __n3=self.mesh.cMatrix[__n][2] #Node 3 index
            ##GET INTERPOLATION ELEMENT########################
            ##GET NATURAL SYSTEM VALUES########################
            __x21=self.mesh.pointMatrix[__n2][0]-self.mesh.pointMatrix[__n1][0]
            __y21=self.mesh.pointMatrix[__n2][1]-self.mesh.pointMatrix[__n1][1]
            __x31=self.mesh.pointMatrix[__n3][0]-self.mesh.pointMatrix[__n1][0]
            __y31=self.mesh.pointMatrix[__n3][1]-self.mesh.pointMatrix[__n1][1]
            __xi=(__y31*(x-self.mesh.pointMatrix[__n1][0])-__x31*(y-self.mesh.pointMatrix[__n1][1]))/\
            (__x21*__y31-__x31*__y21)
            __eta=(-__y21*(x-self.mesh.pointMatrix[__n1][0])+__x21*(y-self.mesh.pointMatrix[__n1][1]))/\
            (__x21*__y31-__x31*__y21)
            __N1=1-__xi-__eta
            __N2=__xi
            __N3=__eta
            ##GET NATURAL SYSTEM VALUES########################
            ##COMPUTE DRAWDOWN VALUES##########################
            __h0=self.phyMed.h_n[0][__n1]*__N1+\
            self.phyMed.h_n[0][__n2]*__N2+\
            self.phyMed.h_n[0][__n3]*__N3
            __hi=self.phyMed.h_n[step][__n1]*__N1+\
            self.phyMed.h_n[step][__n2]*__N2+\
            self.phyMed.h_n[step][__n3]*__N3
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
    def export_h0(self,file):
        #WRITE AN OUTPUT FILE CONTAINING THE HYDRAULIC HEAD DISTRIBUTION#######
        __f=open(aux_f.path+'/'+file+'.hh0','w') #Open the given file in 'write' mode
        
        __f.write(str(self.mesh.nNodes)+'\n')
        for i in range(self.mesh.nNodes):
            __f.write('{0} {1:.7E}\n'.format(i,self.phyMed.h_n[i]))
        __f.close()
    ##AUXILIAR EXPORTATION OF INITIAL CONDITIONS###############################
###############################################################################
#####OUPUTS####################################################################
###############################################################################
"""
##############################MESH CLASS SECTION###############################
"""