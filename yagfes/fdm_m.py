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
    
#######################FINITE DIFFERENCES METHOD MODULE########################

This module includes all object classes and functions required to use the
Finite Differences Method (FDM) for Groundwater (GW) modelling.

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
    nx = # of horizontal nodes      | ny = # of vertical nodes
    dx = 'x' step size              | dy = 'y' step size
    pointMatrix = Point Matrix
    """
    def __init__(self,file):
        #READ MESH FILE########################################################
        __f=aux_f.readDict(file) #Open the given file
        
        ##BUILD MESH###########################################################
        ###MESH CHARACTERISTICS################################################
        self.dx=__f['dx']
        self.nx=round((__f['max_x']-__f['min_x'])/self.dx)+1
        self.dy=__f['dy']
        self.ny=round((__f['max_y']-__f['min_y'])/self.dy)+1
        ###MESH COORDINATES####################################################
        self.nNodes=self.nx*self.ny
        self.pointMatrix=[0]*self.nNodes #Initialize the array containing the node coordinates
        ########################
        for i in range(self.nx):
            for j in range(self.ny):
                __node=i+j*self.nx
                self.pointMatrix[__node]=[i*self.dx+__f['min_x'],j*self.dy+__f['min_y']]
        ##BUILD MESH###########################################################
        ##INITIALIZE AUXILIAR ARRAY TO KEEP TRACK OF NODE ID's#################
        ##CREATE AUXILIAR LIST TO KEEP TRACK OF ACTUAL NODE INDEXES############
        self.ID=[0]*self.nNodes
        for i in range(self.nNodes):
            self.ID[i]=i
        #READ MESH FILE########################################################
        
    def writeBCfile(self,file):
        #WRITE A BC FILE BASED ON MESH INFO####################################
        __f=open(aux_f.path+'/'+file,'w') #Open the given file in 'write' mode
        
        __f.write('#BC FILE\n') #Write document title
        ##MESH LEFT
        __f.write('#MESH LEFT\n')
        __f.write('x_min=[')
        __bc_type=input('Please input the boundary type at the left boundary of the mesh (x='+str(self.pointMatrix[0][0])+'):\n 1) Dirichlet\n 2) Neumann\n')
        __f.write(str(int(__bc_type)-1))
        if int(__bc_type)==1: #Dirichlet BC
            __type='0'
            while __type!='1' and __type!='2':
                __type=input('You chose a Dirichlet BC. '+\
                             'Please input the type of prescribed head line that you prefer: '+\
                             '\n 1) Constant head\n 2) Gradient\n')
                if __type!='1' and __type!='2':
                    print("That is not a valid option. Please, try again.")
            if __type=='1': #Constant head
                __v=input('Please, input the Dirichlet BC at the left boundary.\n') #Read Dirichlet BC value
                for i in range(self.ny):
                    __f.write(','+str(__v))
            else: #Gradient head
                __v=input('Please, input the Dirichlet BC at the first node of the left boundary (y='+str(self.pointMatrix[0][1])+').\n') #Read Dirichlet BC value
                __v2=input('Please, input the Dirichlet BC at the last node of the left boundary (y='+str(self.pointMatrix[self.nx*(self.ny-1)][1])+').\n') #Read Dirichlet BC value
                __delta=(float(__v2)-float(__v))/(self.ny-1)
                for i in range(self.ny):
                    __f.write(','+str(float(__v)+i*__delta))
            __f.write(']\n')
        else: #Neumann BC
            __v=input('Please, input the Neumann BC.\n')
            __f.write(','+str(__v)+']\n')
        ##MESH RIGHT
        __f.write('#MESH RIGHT\n')
        __f.write('x_max=[')
        __bc_type=input('Please input the boundary type at the right boundary of the mesh (x='+str(self.pointMatrix[self.nx-1][0])+'):\n 1) Dirichlet\n 2) Neumann\n')
        __f.write(str(int(__bc_type)-1))
        if int(__bc_type)==1: #Dirichlet BC
            __type='0'
            while __type!='1' and __type!='2':
                __type=input('You chose a Dirichlet BC. '+\
                             'Please input the type of prescribed head line that you prefer: '+\
                             '\n 1) Constant head\n 2) Gradient\n')
                if __type!='1' and __type!='2':
                    print("That is not a valid option. Please, try again.")
            if __type=='1': #Constant head
                __v=input('Please, input the Dirichlet BC at the right boundary.\n') #Read Dirichlet BC value
                for i in range(self.ny):
                    __f.write(','+str(__v))
            else: #Gradient head
                __v=input('Please, input the Dirichlet BC at the first node of the right boundary (y='+str(self.pointMatrix[0][1])+').\n') #Read Dirichlet BC value
                __v2=input('Please, input the Dirichlet BC at the last node of the right boundary (y='+str(self.pointMatrix[self.nx*(self.ny-1)][1])+').\n') #Read Dirichlet BC value
                __delta=(float(__v2)-float(__v))/(self.ny-1)
                for i in range(self.ny):
                    __f.write(','+str(float(__v)+i*__delta))
            __f.write(']\n')
        else: #Neumann BC
            __v=input('Please, input the Neumann BC.\n')
            __f.write(','+str(__v)+']\n')
        ##MESH TOP
        __f.write('#MESH TOP\n')
        __f.write('y_min=[')
        __bc_type=input('Please input the boundary type at the top of the mesh (y='+str(self.pointMatrix[0][1])+'):\n 1) Dirichlet\n 2) Neumann\n')
        __f.write(str(int(__bc_type)-1))
        if int(__bc_type)==1: #Dirichlet BC
            __type='0'
            while __type!='1' and __type!='2':
                __type=input('You chose a Dirichlet BC. '+\
                             'Please input the type of prescribed head line that you prefer: '+\
                             '\n 1) Constant head\n 2) Gradient\n')
                if __type!='1' and __type!='2':
                    print("That is not a valid option. Please, try again.")
            if __type=='1': #Constant head
                __v=input('Please, input the Dirichlet BC at the top boundary.\n') #Read Dirichlet BC value
                for i in range(self.nx):
                    __f.write(','+str(__v))
            else: #Gradient head
                __v=input('Please, input the Dirichlet BC at the first node of the top boundary (x='+str(self.pointMatrix[0][0])+').\n') #Read Dirichlet BC value
                __v2=input('Please, input the Dirichlet BC at the last node of the top boundary (x='+str(self.pointMatrix[self.nx-1][0])+').\n') #Read Dirichlet BC value
                __delta=(float(__v2)-float(__v))/(self.nx-1)
                for i in range(self.nx):
                    __f.write(','+str(float(__v)+i*__delta))
            __f.write(']\n')
        else: #Neumann BC
            __v=input('Please, input the Neumann BC.\n')
            __f.write(','+str(__v)+']\n')
        ##MESH BOTTOM
        __f.write('#MESH BOTTOM\n')
        __f.write('y_max=[')
        __bc_type=input('Please input the boundary type at the bottom of the mesh (y='+str(self.pointMatrix[self.nx*(self.ny-1)][1])+'):\n 1) Dirichlet\n 2) Neumann\n')
        __f.write(str(int(__bc_type)-1))
        if int(__bc_type)==1: #Dirichlet BC
            __type='0'
            while __type!='1' and __type!='2':
                __type=input('You chose a Dirichlet BC. '+\
                             'Please input the type of prescribed head line that you prefer: '+\
                             '\n 1) Constant head\n 2) Gradient\n')
                if __type!='1' and __type!='2':
                    print("That is not a valid option. Please, try again.")
            if __type=='1': #Constant head
                __v=input('Please, input the Dirichlet BC at the bottom boundary.\n') #Read Dirichlet BC value
                for i in range(self.nx):
                    __f.write(','+str(__v))
            else: #Gradient head
                __v=input('Please, input the Dirichlet BC at the first node of the bottom boundary (x='+str(self.pointMatrix[0][0])+').\n') #Read Dirichlet BC value
                __v2=input('Please, input the Dirichlet BC at the last node of the bottom boundary (x='+str(self.pointMatrix[self.nx-1][0])+').\n') #Read Dirichlet BC value
                __delta=(float(__v2)-float(__v))/(self.nx-1)
                for i in range(self.nx):
                    __f.write(','+str(float(__v)+i*__delta))
            __f.write(']\n')
        else: #Neumann BC
            __v=input('Please, input the Neumann BC.\n')
            __f.write(','+str(__v)+']\n')
        #WRITE A BC FILE BASED ON MESH INFO####################################
        __f.close()
        
    def coord2node(self,x,y):
        #RETURN THE NEAREST NODE TO THE GIVEN COORDINATES######################
        ##Find the nearest node to the given coordinates by brute force########
        __min=100 #Auxiliar variable to store the previous minimum distance
        for i in range(self.nx):
            for j in range(self.ny):
                __node_ID=i+j*self.nx
                __d=((self.pointMatrix[__node_ID][0]-x)**2+\
                     (self.pointMatrix[__node_ID][1]-y)**2)**0.5
                #If '__d' is lower than '__min', replace '__min' by '__d' and
                # assume '__node_ID' as the closest node to '(x,y)'.
                if __d<__min:
                    __min=__d #Replace '__min' by '__d'
                    node_ID=__node_ID #Assume '__node_ID' is the closes node to '(x,y)'
                #If '__d' is equal to 0, then '__node_ID' is @ '(x,y)'
                elif __d==0: #Check if the node is equal to 0
                    node_ID=__node_ID
                    break
        #######################################################################
        return node_ID
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
    k_n = Hydraulic Conductivity    | ss_n = Specific Storage
    bc_d = Dirichlet BC's           | bc_n = Neumann BC's
    h_n = Hydraulic Head            | q = Pumping/Recharge Well Rate
    """
    def __init__(self,config,nNodes):
        #INITIALIZE ALL THE ARRAYS USED FOR STORING PHYSICAL PROPERTIES########
        ##INITIALIZE HYDRAULIC CONDUCTIVITY ARRAY##############################
        self.k_n=[1]*(nNodes)
        for i in range(nNodes):
            self.k_n[i]=[float(config['Kx']),float(config['Ky'])] #Hydraulic Conductivity Array ((nx x ny) x 2)
        #######################################################################
        ##INITIALIZE SPECIFIC STORAGE ARRAY####################################
        self.ss_n=[float(config['Ss'])]*(nNodes)
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
        self.bc_d=dict()
        self.bc_n=dict()
        #######################################################################
        #INITIALIZE DICTIONARY FOR STORING LAYER/ZONE PROPERTIES###############
        ##INITIALIZE AUXILIAR ARRAY FOR STORING NODE IDS
        __nodes=[0]*nNodes
        for i in range(nNodes):
            __nodes[i]=i
        ##DEFINE DICTIONARY WITH ZONE PROPERTIES
        self.zones={0:{'Ss':float(config['Ss']),
                       'Kx':float(config['Kx']),
                       'Ky':float(config['Ky']),
                       'nodes':__nodes}}
        #######################################################################
    
    #READ BOUNDARY CONDITIONS##################################################
    def readBCfile(self,file):
        #READ BC FILE##########################################################
        __f=aux_f.readDict(file) #Open the given file in 'read' mode
        
        #MESH LEFT
        if __f['x_min'][0]==0: #Dirichlet BC
            self.bc_d['x_min']=[0]*(len(__f['x_min'])-1)
            for i in range(len(__f['x_min'])-1):
                self.bc_d['x_min'][i]=__f['x_min'][i+1]
        else: #Neumann BC
            self.bc_n['x_min']=__f['x_min'][1]
        #MESH RIGHT
        if __f['x_max'][0]==0: #Dirichlet BC
            self.bc_d['x_max']=[0]*(len(__f['x_max'])-1)
            for i in range(len(__f['x_max'])-1):
                self.bc_d['x_max'][i]=__f['x_max'][i+1]
        else: #Neumann BC
            self.bc_n['x_max']=__f['x_max'][1]
        #MESH TOP
        if __f['y_min'][0]==0: #Dirichlet BC
            self.bc_d['y_min']=[0]*(len(__f['y_min'])-1)
            for i in range(len(__f['y_min'])-1):
                self.bc_d['y_min'][i]=__f['y_min'][i+1]
        else: #Neumann BC
            self.bc_n['y_min']=__f['y_min'][1]
        #MESH BOTTOM
        if __f['y_max'][0]==0: #Dirichlet BC
            self.bc_d['y_max']=[0]*(len(__f['y_max'])-1)
            for i in range(len(__f['y_max'])-1):
                self.bc_d['y_max'][i]=__f['y_max'][i+1]
        else: #Neumann BC
            self.bc_n['y_max']=__f['y_max'][1]
        #READ BC FILE##########################################################
        __f.close()
    #READ BOUNDARY CONDITIONS##################################################
        
    #READ INITIAL CONDITIONS###################################################
    def read_h0(self,file):
        #READS AN HHD FILE AND ASSIGNS IT TO HH_0##############################
        __f=open(aux_f.path+'/'+file)
        
        __f.readline() #Skip time value
        __nrow=int(__f.readline().rstrip('\n').split('=')[1]) #Read # of rows
        __ncol=int(__f.readline().rstrip('\n').split('=')[1]) #Read # of columns
        
        __f.readline() #Skip x_boundaries
        __f.readline() #Skip y_boundaries
        
        for j in range(__nrow):
            __aux=__f.readline().split()
            for i in range(__ncol):
                self.h_n[0][i+j*__ncol]=float(__aux[i])
        
        __f.close()
        #READS AN HHD FILE AND ASSIGNS IT TO HH_0##############################
    #READ INITIAL CONDITIONS###################################################
"""
#########################PHYSICAL MEDIUM CLASS SECTION#########################
"""
###############################################################################
"""
###########################FDM OBJECT CLASS SECTION############################
"""
class fdm:
    """
    mesh = Mesh Object                  | phyMed = Physical Medium Object
    coeffMatrix = Coefficient Matrix    | loadVector = Load Vector
    """  
    def __init__(self,init):
        self.type='FDM'
        #READ INIT FILE########################################################
        self.config=aux_f.readDict(init)
        #READ INIT FILE########################################################
        #CREATE MESH AND PHYMED OBJECTS########################################
        self.mesh=mesh(self.config['mesh_file'])
        self.phyMed=phymed(self.config,self.mesh.nNodes)
        #CREATE MESH AND PHYMED OBJECTS########################################
        #INITIALIZE MATRICES ARRAYS
        self.coeffMatrix=list() #INITIALIZE ARRAY FOR COEFFICIENTS MATRIX
        self.loadVector=list() #INITIALIZE ARRAY FOR LOAD VECTOR
        #######################################################################
    
###############################################################################
#####ASSEMBLE MATRICES AND VECTORS#############################################
###############################################################################
    ##ASSEMBLE COEFFICIENTS MATRIX#############################################
    def __assembleCoefficientMatrix(self):
        #INITIALIZE ARRAY FOR COEFFICIENTS MATRIX
        self.coeffMatrix=[0]*self.mesh.nNodes
        for i in range(self.mesh.nNodes):
            self.coeffMatrix[i]=[0]*self.mesh.nNodes
        #########################################
        for j in range(self.mesh.ny):
            for i in range(self.mesh.nx):
                __node=i+j*self.mesh.nx
                ##COMPUTE AND ASSIGN X-COEFFICIENTS############################
                if i==0: #x_min
                    ##GET 'EQUIVALENT CONDUCTIVIES'############################
                    __Kxb=self.phyMed.k_n[__node][0]
                    __Kxf=2/(1/self.phyMed.k_n[__node][0]+1/self.phyMed.k_n[__node+1][0])
                    ##GET 'EQUIVALENT CONDUCTIVIES'############################
                    ##GET 'CONDUCTANCE' COEFFICIENTS###########################
                    __Cxb=__Kxb/(self.mesh.dx**2)
                    __Cxf=__Kxf/(self.mesh.dx**2)
                    ##GET 'CONDUCTANCE' COEFFICIENTS###########################
                    ##ASSIGN VALUES TO MATRIX##################################
                    self.coeffMatrix[__node][__node+1]=(__Cxb+__Cxf)/self.phyMed.ss_n[__node] #X+ FLOW
                    ##ASSIGN VALUES TO MATRIX##################################
                elif i==self.mesh.nx-1: #x_max
                    ##GET 'EQUIVALENT CONDUCTIVIES'############################
                    __Kxb=2/(1/self.phyMed.k_n[__node][0]+1/self.phyMed.k_n[__node-1][0])
                    __Kxf=self.phyMed.k_n[__node][0]
                    ##GET 'EQUIVALENT CONDUCTIVIES'############################
                    ##GET 'CONDUCTANCE' COEFFICIENTS###########################
                    __Cxb=__Kxb/(self.mesh.dx**2)
                    __Cxf=__Kxf/(self.mesh.dx**2)
                    ##GET 'CONDUCTANCE' COEFFICIENTS###########################
                    ##ASSIGN VALUES TO MATRIX##################################
                    self.coeffMatrix[__node][__node-1]=(__Cxb+__Cxf)/self.phyMed.ss_n[__node] #X- FLOW
                    ##ASSIGN VALUES TO MATRIX##################################
                else:
                    ##GET 'EQUIVALENT CONDUCTIVIES'############################
                    __Kxb=2/(1/self.phyMed.k_n[__node][0]+1/self.phyMed.k_n[__node-1][0])
                    __Kxf=2/(1/self.phyMed.k_n[__node][0]+1/self.phyMed.k_n[__node+1][0])
                    ##GET 'EQUIVALENT CONDUCTIVIES'############################
                    ##GET 'CONDUCTANCE' COEFFICIENTS###########################
                    __Cxb=__Kxb/(self.mesh.dx**2)
                    __Cxf=__Kxf/(self.mesh.dx**2)
                    ##GET 'CONDUCTANCE' COEFFICIENTS###########################
                    ##ASSIGN VALUES TO MATRIX##################################
                    self.coeffMatrix[__node][__node-1]=__Cxb/self.phyMed.ss_n[__node] #X- FLOW
                    self.coeffMatrix[__node][__node+1]=__Cxf/self.phyMed.ss_n[__node] #X+ FLOW
                    ##ASSIGN VALUES TO MATRIX##################################
                ##COMPUTE AND ASSIGN X-COEFFICIENTS############################
                
                ##COMPUTE Y-COEFFICIENTS#######################################
                if j==0: #y_min
                    ##GET 'EQUIVALENT CONDUCTIVIES'############################
                    __Kyb=self.phyMed.k_n[__node][1]
                    __Kyf=2/(1/self.phyMed.k_n[__node][1]+1/self.phyMed.k_n[__node+self.mesh.nx][1])
                    ##GET 'EQUIVALENT CONDUCTIVIES'############################
                    ##GET 'CONDUCTANCE' COEFFICIENTS###########################
                    __Cyb=__Kyb/(self.mesh.dy**2)
                    __Cyf=__Kyf/(self.mesh.dy**2)
                    ##GET 'CONDUCTANCE' COEFFICIENTS###########################
                    ##ASSIGN VALUES TO MATRIX##################################
                    self.coeffMatrix[__node][__node+self.mesh.nx]=(__Cyb+__Cyf)/self.phyMed.ss_n[__node] #Y+ FLOW
                    ##ASSIGN VALUES TO MATRIX##################################
                elif j==self.mesh.ny-1: #y_max
                    ##GET 'EQUIVALENT CONDUCTIVIES'############################
                    __Kyb=2/(1/self.phyMed.k_n[__node][1]+1/self.phyMed.k_n[__node-self.mesh.nx][1])
                    __Kyf=self.phyMed.k_n[__node][1]
                    ##GET 'EQUIVALENT CONDUCTIVIES'############################
                    ##GET 'CONDUCTANCE' COEFFICIENTS###########################
                    __Cyb=__Kyb/(self.mesh.dy**2)
                    __Cyf=__Kyf/(self.mesh.dy**2)
                    ##GET 'CONDUCTANCE' COEFFICIENTS###########################
                    ##ASSIGN VALUES TO MATRIX##################################
                    self.coeffMatrix[__node][__node-self.mesh.nx]=(__Cyb+__Cyf)/self.phyMed.ss_n[__node] #Y- FLOW
                    ##ASSIGN VALUES TO MATRIX##################################
                else:
                    ##GET 'EQUIVALENT CONDUCTIVIES'############################
                    __Kyb=2/(1/self.phyMed.k_n[__node][1]+1/self.phyMed.k_n[__node-self.mesh.nx][1])
                    __Kyf=2/(1/self.phyMed.k_n[__node][1]+1/self.phyMed.k_n[__node+self.mesh.nx][1])
                    ##GET 'EQUIVALENT CONDUCTIVIES'############################
                    ##GET 'CONDUCTANCE' COEFFICIENTS###########################
                    __Cyb=__Kyb/(self.mesh.dy**2)
                    __Cyf=__Kyf/(self.mesh.dy**2)
                    ##GET 'CONDUCTANCE' COEFFICIENTS###########################
                    ##ASSIGN VALUES TO MATRIX##################################
                    self.coeffMatrix[__node][__node-self.mesh.nx]=__Cyb/self.phyMed.ss_n[__node] #Y- FLOW
                    self.coeffMatrix[__node][__node+self.mesh.nx]=__Cyf/self.phyMed.ss_n[__node] #Y+ FLOW
                    ##ASSIGN VALUES TO MATRIX##################################
                ##COMPUTE Y-COEFFICIENTS#######################################
                
                ##CHECK SYSTEM STATE###########################################
                if self.phyMed.steady==True: #Steady-state
                    __Ct=0
                else: #Transient-state
                    __Ct=1/self.phyMed.dt
                ##CHECK SYSTEM STATE###########################################
                
                ##ASSIGN VALUES TO CENTRAL NODE################################
                self.coeffMatrix[__node][__node]=-((__Cxb+__Cxf+__Cyb+__Cyf)/self.phyMed.ss_n[__node]+__Ct)
                ##ASSIGN VALUES TO CENTRAL NODE################################
    ##ASSEMBLE COEFFICIENTS MATRIX#############################################
    
    ##ASSEMBLE LOAD VECTOR#####################################################
    def __assembleLoadVector(self):
        #Steady-state##########################################################
        if self.phyMed.steady==True:
            #INITIALIZE ARRAY FOR LOAD VECTOR
            self.loadVector=[0]*self.mesh.nNodes
            #ASSIGN VALUES
            for i in self.phyMed.q:
                self.loadVector[i]=self.phyMed.q[i]
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
                        self.loadVector[k][i]=v
        #Transient-state#######################################################
    ##ASSEMBLE LOAD VECTOR#####################################################
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
        #Steady-State##########################################################
        if self.phyMed.steady==True:
            ##MESH LEFT
            if 'x_min' in self.phyMed.bc_d:
                for i in range(self.mesh.ny):
                    __node=i*self.mesh.nx
                    if __node in self.mesh.ID:
                        __rID=self.mesh.ID.index(__node) #Retrieves the real ID for __node
                        del self.mesh.ID[__rID] #Remove ID from the ID list
                        ##REMOVE NODE FROM HYDRAULIC HEAD VECTOR###################
                        del self.phyMed.h_n[__rID]
                        ##REMOVE NODE FROM LOAD VECTOR#############################
                        del self.loadVector[__rID]
                        ###REMOVE ROW FROM COEFFICIENT MATRIX######################
                        del self.coeffMatrix[__rID]
                        ##ADDITIONAL OPERATIONS####################################
                        for j in range(len(self.loadVector)):
                            self.loadVector[j]-=self.coeffMatrix[j][__rID]*self.phyMed.bc_d['x_min'][i]*uF
                            ##REMOVE COLUMN FROM COEFFICIENT MATRIX################
                            del self.coeffMatrix[j][__rID]
            ##MESH RIGHT
            if 'x_max' in self.phyMed.bc_d:
                for i in range(self.mesh.ny):
                    __node=(i+1)*self.mesh.nx-1
                    if __node in self.mesh.ID:
                        __rID=self.mesh.ID.index(__node) #Retrieves the real ID for __node
                        del self.mesh.ID[__rID] #Remove ID from the ID list
                        ##REMOVE NODE FROM HYDRAULIC HEAD VECTOR###################
                        del self.phyMed.h_n[__rID]
                        ##REMOVE NODE FROM LOAD VECTOR#############################
                        del self.loadVector[__rID]
                        ###REMOVE ROW FROM COEFFICIENT MATRIX######################
                        del self.coeffMatrix[__rID]
                        ##ADDITIONAL OPERATIONS####################################
                        for j in range(len(self.loadVector)):
                            self.loadVector[j]-=self.coeffMatrix[j][__rID]*self.phyMed.bc_d['x_max'][i]*uF
                            ##REMOVE COLUMN FROM COEFFICIENT MATRIX################
                            del self.coeffMatrix[j][__rID]
            ##MESH TOP
            if 'y_min' in self.phyMed.bc_d:
                for i in range(self.mesh.nx):
                    __node=i
                    if __node in self.mesh.ID:
                        __rID=self.mesh.ID.index(__node) #Retrieves the real ID for __node
                        del self.mesh.ID[__rID] #Remove ID from the ID list
                        ##REMOVE NODE FROM HYDRAULIC HEAD VECTOR###################
                        del self.phyMed.h_n[__rID]
                        ##REMOVE NODE FROM LOAD VECTOR#############################
                        del self.loadVector[__rID]
                        ###REMOVE ROW FROM COEFFICIENT MATRIX######################
                        del self.coeffMatrix[__rID]
                        ##ADDITIONAL OPERATIONS####################################
                        for j in range(len(self.loadVector)):
                            self.loadVector[j]-=self.coeffMatrix[j][__rID]*self.phyMed.bc_d['y_min'][i]*uF
                            ##REMOVE COLUMN FROM COEFFICIENT MATRIX################
                            del self.coeffMatrix[j][__rID]
            ##MESH BOTTOM
            if 'y_max' in self.phyMed.bc_d:
                for i in range(self.mesh.nx):
                    __node=i+self.mesh.nx*(self.mesh.ny-1)
                    if __node in self.mesh.ID:
                        __rID=self.mesh.ID.index(__node) #Retrieves the real ID for __node
                        del self.mesh.ID[__rID] #Remove ID from the ID list
                        ##REMOVE NODE FROM HYDRAULIC HEAD VECTOR###################
                        del self.phyMed.h_n[__rID]
                        ##REMOVE NODE FROM LOAD VECTOR#############################
                        del self.loadVector[__rID]
                        ###REMOVE ROW FROM COEFFICIENT MATRIX######################
                        del self.coeffMatrix[__rID]
                        ##ADDITIONAL OPERATIONS####################################
                        for j in range(len(self.loadVector)):
                            self.loadVector[j]-=self.coeffMatrix[j][__rID]*self.phyMed.bc_d['y_max'][i]*uF
                            ##REMOVE COLUMN FROM COEFFICIENT MATRIX################
                            del self.coeffMatrix[j][__rID]
        #Steady-State##########################################################        
        #Transient-state#######################################################
        else:
            ##MESH LEFT
            if 'x_min' in self.phyMed.bc_d:
                for i in range(self.mesh.ny):
                    __node=i*self.mesh.nx
                    if __node in self.mesh.ID:
                        __rID=self.mesh.ID.index(__node) #Retrieves the real ID for __node
                        del self.mesh.ID[__rID] #Remove ID from the ID list
                        ##REMOVE NODE FROM HYDRAULIC HEAD VECTOR###################
                        for j in range(self.phyMed.timeSteps):
                            del self.phyMed.h_n[j][__rID]
                        ##REMOVE NODE FROM LOAD VECTOR#############################
                        for j in range(self.phyMed.timeSteps):
                            del self.loadVector[j][__rID]
                        ###REMOVE ROW FROM COEFFICIENT MATRIX######################
                        del self.coeffMatrix[__rID]
                        ##ADDITIONAL OPERATIONS####################################
                        for j in range(len(self.loadVector[0])):
                            for k in range(self.phyMed.timeSteps):
                                self.loadVector[k][j]-=self.coeffMatrix[j][__rID]*self.phyMed.bc_d['x_min'][i]*uF
                            ##REMOVE COLUMN FROM COEFFICIENT MATRIX################
                            del self.coeffMatrix[j][__rID]
            ##MESH RIGHT
            if 'x_max' in self.phyMed.bc_d:
                for i in range(self.mesh.ny):
                    __node=(i+1)*self.mesh.nx-1
                    if __node in self.mesh.ID:
                        __rID=self.mesh.ID.index(__node) #Retrieves the real ID for __node
                        del self.mesh.ID[__rID] #Remove ID from the ID list
                        ##REMOVE NODE FROM HYDRAULIC HEAD VECTOR###################
                        for j in range(self.phyMed.timeSteps):
                            del self.phyMed.h_n[j][__rID]
                        ##REMOVE NODE FROM LOAD VECTOR#############################
                        for j in range(self.phyMed.timeSteps):
                            del self.loadVector[j][__rID]
                        ###REMOVE ROW FROM COEFFICIENT MATRIX######################
                        del self.coeffMatrix[__rID]
                        ##ADDITIONAL OPERATIONS####################################
                        for j in range(len(self.loadVector[0])):
                            for k in range(self.phyMed.timeSteps):
                                self.loadVector[k][j]-=self.coeffMatrix[j][__rID]*self.phyMed.bc_d['x_max'][i]*uF
                            ##REMOVE COLUMN FROM COEFFICIENT MATRIX################
                            del self.coeffMatrix[j][__rID]
            ##MESH TOP
            if 'y_min' in self.phyMed.bc_d:
                for i in range(self.mesh.nx):
                    __node=i
                    if __node in self.mesh.ID:
                        __rID=self.mesh.ID.index(__node) #Retrieves the real ID for __node
                        del self.mesh.ID[__rID] #Remove ID from the ID list
                        ##REMOVE NODE FROM HYDRAULIC HEAD VECTOR###################
                        for j in range(self.phyMed.timeSteps):
                            del self.phyMed.h_n[j][__rID]
                        ##REMOVE NODE FROM LOAD VECTOR#############################
                        for j in range(self.phyMed.timeSteps):
                            del self.loadVector[j][__rID]
                        ###REMOVE ROW FROM COEFFICIENT MATRIX######################
                        del self.coeffMatrix[__rID]
                        ##ADDITIONAL OPERATIONS####################################
                        for j in range(len(self.loadVector[0])):
                            for k in range(self.phyMed.timeSteps):
                                self.loadVector[k][j]-=self.coeffMatrix[j][__rID]*self.phyMed.bc_d['y_min'][i]*uF
                            ##REMOVE COLUMN FROM COEFFICIENT MATRIX################
                            del self.coeffMatrix[j][__rID]
            ##MESH BOTTOM
            if 'y_max' in self.phyMed.bc_d:
                for i in range(self.mesh.nx):
                    __node=i+self.mesh.nx*(self.mesh.ny-1)
                    if __node in self.mesh.ID:
                        __rID=self.mesh.ID.index(__node) #Retrieves the real ID for __node
                        del self.mesh.ID[__rID] #Remove ID from the ID list
                        ##REMOVE NODE FROM HYDRAULIC HEAD VECTOR###################
                        for j in range(self.phyMed.timeSteps):
                            del self.phyMed.h_n[j][__rID]
                        ##REMOVE NODE FROM LOAD VECTOR#############################
                        for j in range(self.phyMed.timeSteps):
                            del self.loadVector[j][__rID]
                        ###REMOVE ROW FROM COEFFICIENT MATRIX######################
                        del self.coeffMatrix[__rID]
                        ##ADDITIONAL OPERATIONS####################################
                        for j in range(len(self.loadVector[0])):
                            for k in range(self.phyMed.timeSteps):
                                self.loadVector[k][j]-=self.coeffMatrix[j][__rID]*self.phyMed.bc_d['y_max'][i]*uF
                            ##REMOVE COLUMN FROM COEFFICIENT MATRIX################
                            del self.coeffMatrix[j][__rID]
        #Transient-state#######################################################
    ##APPLY DIRICHLET BC's#####################################################
    
    ##APPLY NEUMANN BC's#######################################################
    def __applyNeumannBC(self):
        #UNIT CONVERSION#######################################################
        uF=1
        if self.config['length'].lower()!=self.config['head'].lower(): #Check if conversion is required
            uF=(aux_f.unitDict[self.config['length'].lower()])/(aux_f.unitDict[self.config['head'].lower()]) #Get conversion factor            
        #UNIT CONVERSION#######################################################
        ##INCREASE MESH SIZE###################################################
        #ASSIGN FICTIONAL HEAD VALUES TO MATCH FLOW CONDITION##################
        if self.phyMed.steady==True: #Steady-state
            ###MESH LEFT
            if 'x_min' in self.phyMed.bc_n:
                for i in range(self.mesh.ny):
                    __node=i*self.mesh.nx
                    __f=2*self.phyMed.bc_n['x_min']*self.config['aq_thickness']/self.mesh.dx
                    self.loadVector[__node]-=__f*uF/self.phyMed.ss_n[__node]
            ###MESH RIGHT
            if 'x_max' in self.phyMed.bc_n:
                for i in range(self.mesh.ny):
                    __node=(i+1)*self.mesh.nx-1
                    __f=2*self.phyMed.bc_n['x_max']*self.config['aq_thickness']/self.mesh.dx
                    self.loadVector[__node]-=__f*uF/self.phyMed.ss_n[__node]
            ###MESH TOP
            if 'y_min' in self.phyMed.bc_n:
                for i in range(self.mesh.nx):
                    __node=i
                    __f=2*self.phyMed.bc_n['y_min']*self.config['aq_thickness']/self.mesh.dy
                    self.loadVector[__node]-=__f*uF/self.phyMed.ss_n[__node]
            ###MESH BOTTOM
            if 'y_max' in self.phyMed.bc_n:
                for i in range(self.mesh.nx):
                    __node=i+self.mesh.nx*(self.mesh.ny-1)
                    __f=2*self.phyMed.bc_n['y_max']*self.config['aq_thickness']/self.mesh.dy
                    self.loadVector[__node]-=__f*uF/self.phyMed.ss_n[__node]
        else: #Transient-state
            ###MESH LEFT
            if 'x_min' in self.phyMed.bc_n:
                for i in range(self.mesh.ny):
                    __node=i*self.mesh.nx
                    for j in range(self.phyMed.timeSteps):
                        __f=2*self.phyMed.bc_n['x_min']*self.config['aq_thickness']/self.mesh.dx
                        self.loadVector[j][__node]-=__f*uF/self.phyMed.ss_n[__node]
            ###MESH RIGHT
            if 'x_max' in self.phyMed.bc_n:
                for i in range(self.mesh.ny):
                    __node=(i+1)*self.mesh.nx-1
                    for j in range(self.phyMed.timeSteps):
                        __f=2*self.phyMed.bc_n['x_max']*self.config['aq_thickness']/self.mesh.dx
                        self.loadVector[j][__node]-=__f*uF/self.phyMed.ss_n[__node]
            ###MESH TOP
            if 'y_min' in self.phyMed.bc_n:
                for i in range(self.mesh.nx):
                    __node=i
                    for j in range(self.phyMed.timeSteps):
                        __f=2*self.phyMed.bc_n['y_min']*self.config['aq_thickness']/self.mesh.dy
                        self.loadVector[j][__node]-=__f*uF/self.phyMed.ss_n[__node]
            ###MESH BOTTOM
            if 'y_max' in self.phyMed.bc_n:
                for i in range(self.mesh.nx):
                    __node=i+self.mesh.nx*(self.mesh.ny-1)
                    for j in range(self.phyMed.timeSteps):
                        __f=2*self.phyMed.bc_n['y_max']*self.config['aq_thickness']/self.mesh.dy
                        self.loadVector[j][__node]-=__f*uF/self.phyMed.ss_n[__node]
        #ASSIGN FICTIONAL HEAD VALUES TO MATCH FLOW CONDITION##################
    ##APPLY NEUMANN BC's#######################################################
    
    ##REBUILD HH ARRAY#########################################################
    def __rebuild_h(self):
        #UNIT CONVERSION#######################################################
        if self.config['length'].lower()!=self.config['head'].lower(): #Check if conversion is required
            uF=aux_f.unitDict[self.config['head'].lower()]/aux_f.unitDict[self.config['length'].lower()] #Get conversion factor
            for i in range(len(self.phyMed.h_n)):
                self.phyMed.h_n[i]*=uF
        #UNIT CONVERSION#######################################################
        #Steady-State##########################################################
        if self.phyMed.steady==True:
            ##Re-initialize hydraulic head array###############################
            __n=len(self.mesh.ID)
            __hh=[0]*__n
            ###Pass values to new array########################################
            for i in range(__n):
                __hh[i]=self.phyMed.h_n[i]
            ###Reset hydraulic head array and pass values######################
            self.phyMed.h_n=[0]*self.mesh.nNodes
            for i in range(__n):
                self.phyMed.h_n[self.mesh.ID[i]]=__hh[i]
            ##Re-initialize hydraulic head array###############################
            ##MESH LEFT
            if 'x_min' in self.phyMed.bc_d:
                for i in range(self.mesh.ny):
                    __node=i*self.mesh.nx
                    self.phyMed.h_n[__node]=self.phyMed.bc_d['x_min'][i]
            ##MESH RIGHT
            if 'x_max' in self.phyMed.bc_d:
                for i in range(self.mesh.ny):
                    __node=(i+1)*self.mesh.nx-1
                    self.phyMed.h_n[__node]=self.phyMed.bc_d['x_max'][i]
            ##MESH TOP
            if 'y_min' in self.phyMed.bc_d:
                for i in range(self.mesh.nx):
                    __node=i
                    self.phyMed.h_n[__node]=self.phyMed.bc_d['y_min'][i]
            ##MESH BOTTOM
            if 'y_max' in self.phyMed.bc_d:
                for i in range(self.mesh.nx):
                    __node=i+self.mesh.nx*(self.mesh.ny-1)
                    self.phyMed.h_n[__node]=self.phyMed.bc_d['y_max'][i]
        #Steady-State##########################################################
        #Transient-State#######################################################
        else:
            __n=len(self.mesh.ID)
            __hh=[0]*__n
            for j in range(self.phyMed.timeSteps):
                ##Re-initialize hydraulic head array###########################
                ###Pass values to new array####################################
                for i in range(__n):
                    __hh[i]=self.phyMed.h_n[j][i]
                ###Reset hydraulic head array and pass values##################
                self.phyMed.h_n[j]=[0]*self.mesh.nNodes
                for i in range(__n):
                    self.phyMed.h_n[j][self.mesh.ID[i]]=__hh[i]
                ##Re-initialize hydraulic head array###########################
                ##MESH LEFT
                if 'x_min' in self.phyMed.bc_d:
                    for i in range(self.mesh.ny):
                        __node=i*self.mesh.nx
                        self.phyMed.h_n[j][__node]=self.phyMed.bc_d['x_min'][i]
                ##MESH RIGHT
                if 'x_max' in self.phyMed.bc_d:
                    for i in range(self.mesh.ny):
                        __node=(i+1)*self.mesh.nx-1
                        self.phyMed.h_n[j][__node]=self.phyMed.bc_d['x_max'][i]
                ##MESH TOP
                if 'y_min' in self.phyMed.bc_d:
                    for i in range(self.mesh.nx):
                        __node=i
                        self.phyMed.h_n[j][__node]=self.phyMed.bc_d['y_min'][i]
                ##MESH BOTTOM
                if 'y_max' in self.phyMed.bc_d:
                    for i in range(self.mesh.nx):
                        __node=i+self.mesh.nx*(self.mesh.ny-1)
                        self.phyMed.h_n[j][__node]=self.phyMed.bc_d['y_max'][i]
        #Transient-State#######################################################
    ##REBUILD HH ARRAY#########################################################
###############################################################################
#####APPLY BOUNDARY CONDITIONS#################################################
###############################################################################

###############################################################################
#####MODIFY PHYSICAL MEDIUM####################################################
###############################################################################
    ##ADD WELLS################################################################
    def addWell(self,file):
        #UNIT CONVERSION#######################################################
        uF=1
        if self.config['length'].lower()!=self.config['head'].lower(): #Check if conversion is required
            uF=(aux_f.unitDict[self.config['length'].lower()]**3)/(aux_f.unitDict[self.config['head'].lower()]**3) #Get conversion factor            
        #UNIT CONVERSION#######################################################
        ##RETRIEVE WELL INFORMATION FROM FILE##################################
        __f=aux_f.readDict(file)
        ##RETRIEVE WELL INFORMATION FROM FILE##################################
        __node=self.mesh.coord2node(__f['x'],__f['y']) #Retrieve node indexes
        ##ASSIGN FLOW VALUE TO 'Q' ARRAY#######################################
        if self.phyMed.steady==True: #Steady-state
            v=-__f['rate']*uF/\
            (self.mesh.dx*self.mesh.dy*self.config['aq_thickness'])/self.phyMed.ss_n[__node]
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
                    (self.mesh.dx*self.mesh.dy*self.config['aq_thickness'])/self.phyMed.ss_n[__node]
                self.phyMed.q[__node][1][i+offset][0]=t_min #Start of stress period 'i'
                self.phyMed.q[__node][1][i+offset][1]=t_max #End of stress period 'i'
                self.phyMed.q[__node][1][i+offset][2]=v #Pumping rate of stress period 'i'
    ##ADD WELLS################################################################
    
    #ZONE RELATED FUNCTIONS####################################################
    ##DEFINE NEW ZONE##########################################################
    def newZone(self,zone_id,Ss:float,Kx:float,Ky:float,x_min:float,x_max:float,y_min:float,y_max:float): #Defines a new zone in the domain as a rectangle
        if zone_id not in self.phyMed.zones: #CHECK IF THE ZONE EXISTS
            #IDENTIFY ALL NODES INSIDE THE ZONE################################
            __x_min=self.mesh.pointMatrix[0][0]
            __y_min=self.mesh.pointMatrix[0][1]
            #Given the way that nodes are ordered, we can easily identify all the
            #'i' and 'j' indexes inside the rectangle as follows:
            ##IDENTIFY 'i_min'
            for i in range(self.mesh.nx):
                if i*self.mesh.dx+__x_min>=x_min:
                    i_min=i
                    break
            ##IDENTIFY 'i_max'
            for i in range(i_min,self.mesh.nx):
                if i*self.mesh.dx+__x_min>=x_max:
                    i_max=i+1
                    break
            ##IDENTIFY 'j_min'
            for i in range(self.mesh.ny):
                if i*self.mesh.dy+__y_min>=y_min:
                    j_min=i
                    break
            ##IDENTIFY 'j_max'
            for i in range(j_min,self.mesh.ny):
                if i*self.mesh.dy+__y_min>=y_max:
                    j_max=i+1
                    break
            #IDENTIFY ALL NODES INSIDE THE ZONE################################  
            #ZONE ASSIGNATION LOOP#############################################
            __ID=list() #Auxiliar list for storing the ID values of the nodes inside the zone
            for j in range(j_min,j_max):
                for i in range(i_min,i_max):
                    __node=i+j*self.mesh.nx
                    #Assign zone values
                    self.phyMed.k_n[__node]=[Kx,Ky]
                    self.phyMed.ss_n[__node]=Ss
                    __ID.append(__node) #Append value to auxiliar list
                    #Look for the node in all other zones and remove it from their lists
                    for k in self.phyMed.zones:
                        if __node in self.phyMed.zones[k]['nodes']: #Check if node is in a zone 'k'
                            __node_ID=self.phyMed.zones[k]['nodes'].index(__node) #Get ID within 'nodes' list
                            del self.phyMed.zones[k]['nodes'][__node_ID] #Remove node from zone 'k'
            self.phyMed.zones[zone_id]={'Ss':Ss,
                             'Kx':Kx,
                             'Ky':Ky,
                             'nodes':__ID}
            #ZONE ASSIGNATION LOOP#############################################
        else:
            print('Zone {0} already exists'.format(zone_id))
    ##DEFINE NEW ZONE##########################################################
    
    ##APPEND NEW RECTANGLE TO EXISTING ZONE####################################
    def expandZone(self,zone_id,x_min:float,x_max:float,y_min:float,y_max:float): #Adds a new rectangle to an existing zone
        if zone_id in self.phyMed.zones: #CHECK IF THE ZONE EXISTS
            #IDENTIFY ALL NODES INSIDE THE ZONE################################
            __x_min=self.mesh.pointMatrix[0][0]
            __y_min=self.mesh.pointMatrix[0][1]
            #Given the way that nodes are ordered, we can easily identify all the
            #'i' and 'j' indexes inside the rectangle as follows:
            ##IDENTIFY 'i_min'
            for i in range(self.mesh.nx):
                if i*self.mesh.dx+__x_min>=x_min:
                    i_min=i
                    break
            ##IDENTIFY 'i_max'
            for i in range(i_min,self.mesh.nx):
                if i*self.mesh.dx+__x_min>=x_max:
                    i_max=i+1
                    break
            ##IDENTIFY 'j_min'
            for i in range(self.mesh.ny):
                if i*self.mesh.dy+__y_min>=y_min:
                    j_min=i
                    break
            ##IDENTIFY 'j_max'
            for i in range(j_min,self.mesh.ny):
                if i*self.mesh.dy+__y_min>=y_max:
                    j_max=i+1
                    break
            #IDENTIFY ALL NODES INSIDE THE ZONE################################  
            #ZONE ASSIGNATION LOOP#############################################
            __ID=list() #Auxiliar list for storing the ID values of the nodes inside the zone
            for j in range(j_min,j_max):
                for i in range(i_min,i_max):
                    __node=i+j*self.mesh.nx
                    #Assign zone values
                    self.phyMed.k_n[__node]=[self.phyMed.zones[zone_id]['Kx'],self.phyMed.zones[zone_id]['Ky']]
                    self.phyMed.ss_n[__node]=self.phyMed.zones[zone_id]['Ss']
                    __ID.append(__node) #Append value to auxiliar list
                    #Look for the node in all other zones and remove it from their lists
                    for k in self.phyMed.zones:
                        if __node in self.phyMed.zones[k]['nodes']: #Check if node is in a zone 'k'
                            __node_ID=self.phyMed.zones[k]['nodes'].index(__node) #Get ID within 'nodes' list
                            del self.phyMed.zones[k]['nodes'][__node_ID] #Remove node from zone 'k'
            #Add nodes to zone
            for i in __ID:
                if i not in self.phyMed.zones[zone_id]['nodes']:
                    self.phyMed.zones[zone_id]['nodes'].append(i)
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
            for i in self.phyMed.zones[zone_id]['nodes']:
                self.phyMed.k_n[i]=[Kx,Ky] #Update node value
                self.phyMed.ss_n[i]=Ss #Update node value
        else:
            print("Zone {0} doesn't exist".format(zone_id))
    ##UPDATE ZONE PARAMETERS###################################################
    
    ##DELETE ZONE##############################################################
    def deleteZone(self,zone_id):
        if zone_id in self.phyMed.zones: #CHECK IF THE ZONE EXISTS
            for i in self.phyMed.zones[zone_id]['nodes']:
                self.phyMed.zones[0]['nodes'].append(i) #Add node to default zone
                self.phyMed.k_n[i]=[self.phyMed.zones[0]['Kx'],self.phyMed.zones[0]['Ky']] #Update node value
                self.phyMed.ss_n[i]=self.phyMed.zones[0]['Ss'] #Update node value
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
    def solve(self,tol=0.001):
        if self.phyMed.steady==True: #Steady-state
            print('ASSEMBLING GLOBAL MATRICES...')
            ##ASSEMBLE GLOBAL MATRICES#########################################
            self.__assembleCoefficientMatrix() #ASSEMBLE COEFFICIENTS MATRIX
            self.__assembleLoadVector() #ASSEMBLE LOAD VECTOR
            ##ASSEMBLE GLOBAL MATRICES#########################################
            
            print('APPLYING BCs...')
            ##APPLY BC's#######################################################
            self.__applyNeumannBC()
            self.__applyDirichletBC()
            ##APPLY BC's#######################################################
            
            print('COMPUTING HYDRAULIC HEAD DISTRIBUTION...')
            ##GET HYDRAULIC HEAD DISTRIBUTION##################################
            __aux_hh=np.matmul(np.linalg.inv(self.coeffMatrix),self.loadVector)
            self.phyMed.h_n=__aux_hh.tolist()
            self.__rebuild_h()
            ##GET HYDRAULIC HEAD DISTRIBUTION################################## 
        else: #Transient-state
            print('ASSEMBLING GLOBAL MATRICES...')
            ##ASSEMBLE GLOBAL MATRICES#########################################
            self.__assembleCoefficientMatrix() #ASSEMBLE COEFFICIENTS MATRIX
            self.__assembleLoadVector() #ASSEMBLE LOAD VECTOR
            ##ASSEMBLE GLOBAL MATRICES#########################################
            
            print('APPLYING BCs...')
            ##APPLY BC's#######################################################
            self.__applyNeumannBC()
            self.__applyDirichletBC()
            ##APPLY BC's#######################################################
            
            print('ASSEMBLING AUXILIAR MATRICES...')
            ##GET AUXILIAR MATRICES############################################
            __Cinv=np.linalg.inv(self.coeffMatrix)
            __n=len(self.coeffMatrix)
            __V=[0]*__n
            
            print('SOLVING FOR EACH TIME STEP...')
            for k in range(self.phyMed.timeSteps-1): #TIME-LOOP
                __step=k+1
                print(str(100*__step/self.phyMed.timeSteps)+'%...')
                
                ##UPDATE AUXILIAR MATRICES#####################################
                for i in range(__n):
                    __V[i]=self.loadVector[k][i]-(1/self.phyMed.dt)*self.phyMed.h_n[__step-1][i]
                ##UPDATE AUXILIAR MATRICES#####################################
                
                ##GET HYDRAULIC HEAD DISTRIBUTION##############################
                ###CORREGIR
                print('COMPUTING HYDRAULIC HEAD DISTRIBUTION...')
                __aux_hh=np.matmul(__Cinv,__V)
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
    ##EXPORT RESULTS TO NATIVE .HHD FORMAT#####################################
    def exportResults(self,file,step=None):
        #WRITE AN OUTPUT FILE CONTAINING THE HYDRAULIC HEAD DISTRIBUTION#######
        __f=open(aux_f.path+'/output/'+file+'.hhd','w') #Open the given file in 'write' mode
        
        #WRITING THE FILE HEADERS##############################################
        ##CURRENT TIPE STEP
        __f.write('time=')
        if step==None: #Steady-state
            __f.write('NA\n')
        else:
            __f.write(str(step*self.phyMed.dt)+'\n') #Writes current time based on the time-step and the dt values
        ##NUMBER OF ROWS
        __f.write('no_rows='+str(self.mesh.ny)+'\n') #Writes the number of rows of the mesh
        ##NUMBER OF COLUMNS
        __f.write('no_colums='+str(self.mesh.nx)+'\n') #Writes the number of colums of the mesh
        ##COORDINATES
        __f.write('x_boundaries=('+str(self.mesh.pointMatrix[0][0])+','+str(self.mesh.pointMatrix[self.mesh.nx-1][0])+')\n')
        __f.write('y_boundaries=('+str(self.mesh.pointMatrix[0][1])+','+str(self.mesh.pointMatrix[self.mesh.nx*(self.mesh.ny-1)][1])+')\n')
        #WRITING THE FILE HEADERS##############################################
        
        #LOOP FOR WRITING VALUES ON THE FILE###################################
        if step==None: #Steady-state
            for j in range(self.mesh.ny):
                for i in range(self.mesh.nx):
                    __node=i+j*self.mesh.nx
                    __f.write("{0:.7E} ".format(self.phyMed.h_n[__node]))
                __f.write('\n')
        else:
            for j in range(self.mesh.ny):
                for i in range(self.mesh.nx):
                    __node=i+j*self.mesh.nx
                    __f.write("{0:.7E} ".format(self.phyMed.h_n[step][__node]))
                __f.write('\n')
        #LOOP FOR WRITING VALUES ON THE FILE###################################
        
        __f.close()
    ##EXPORT RESULTS TO NATIVE .HHD FORMAT#####################################
    
    ##DRAWDOWN AT A GIVEN POINT################################################
    def drawdown_at(self,x,y,file=None):
        #COMPUTE DRAWDOWN VALUE AT X,Y#########################################
        __n=self.mesh.coord2node(x,y) #Get closest node index
        ##INTERPOLATE HEAD VALUE AT X,Y#################
        __h=[0]*self.phyMed.timeSteps #Auxiliar array to store head values
        ###CHECK HOW WE SHOULD INTERPOLATE THE HEAD VALUE
        #Do not interpolate if the node is at (x,y)
        if self.mesh.pointMatrix[__n][0]==x and self.mesh.pointMatrix[__n][1]==y: #Do not interpolate
            for k in range(self.phyMed.timeSteps):
                __h[k]=self.phyMed.h_n[k][__n]
        
        #Use a linear interpolation if the node is on 'x' or 'y'
        elif self.mesh.pointMatrix[__n][1]==y: #Interpolate over 'x'
            if self.mesh.pointMatrix[__n][0]<x: #(x,y) is to the right of the node
                for k in range(self.phyMed.timeSteps):
                    __h[k]=(x-self.mesh.pointMatrix[__n][0])*\
                    (self.phyMed.h_n[k][__n+1]-self.phyMed.h_n[k][__n])/\
                    (self.mesh.pointMatrix[__n+1][0]-self.mesh.pointMatrix[__n][0])+\
                    self.phyMed.h_n[k][__n]
            else: #(x,y) is to the left of the node
                for k in range(self.phyMed.timeSteps):
                    __h[k]=self.phyMed.h_n[k][__n]-\
                    (self.mesh.pointMatrix[__n][0]-x)*\
                    (self.phyMed.h_n[k][__n]-self.phyMed.h_n[k][__n-1])/\
                    (self.mesh.pointMatrix[__n][0]-self.mesh.pointMatrix[__n-1][0])
        elif self.mesh.pointMatrix[__n][0]==x: #Interpolate over 'y'
            if self.mesh.pointMatrix[__n][0]<y: #(x,y) is below the node
                for k in range(self.phyMed.timeSteps):
                    __h[k]=(y-self.mesh.pointMatrix[__n][1])*\
                    (self.phyMed.h_n[k][__n+self.mesh.nx]-self.phyMed.h_n[k][__n])/\
                    (self.mesh.pointMatrix[__n+self.mesh.nx][1]-self.mesh.pointMatrix[__n][1])+\
                    self.phyMed.h_n[k][__n]
            else: #(x,y) is above the node
                for k in range(self.phyMed.timeSteps):
                    __h[k]=self.phyMed.h_n[k][__n]-\
                    (self.mesh.pointMatrix[__n][1]-y)*\
                    (self.phyMed.h_n[k][__n]-self.phyMed.h_n[k][__n-self.mesh.nx])/\
                    (self.mesh.pointMatrix[__n][1]-self.mesh.pointMatrix[__n-self.mesh.nx][1])
        else:
            if self.mesh.pointMatrix[__n][0]<x: #(x,y) is to the right of the node
                if self.mesh.pointMatrix[__n][0]<y: #(x,y) is below the node
                    for k in range(self.phyMed.timeSteps):
                        hy1=(x-self.mesh.pointMatrix[__n+self.mesh.nx][0])*\
                        (self.phyMed.h_n[k][__n+1+self.mesh.nx]-self.phyMed.h_n[k][__n+self.mesh.nx])/\
                        (self.mesh.pointMatrix[__n+1+self.mesh.nx][0]-self.mesh.pointMatrix[__n+self.mesh.nx][0])+\
                        self.phyMed.h_n[k][__n+self.mesh.nx]
                        hy2=(x-self.mesh.pointMatrix[__n][0])*\
                        (self.phyMed.h_n[k][__n+1]-self.phyMed.h_n[k][__n])/\
                        (self.mesh.pointMatrix[__n+1][0]-self.mesh.pointMatrix[__n][0])+\
                        self.phyMed.h_n[k][__n]
                        __h[k]=(y-self.mesh.pointMatrix[__n][1])*\
                        (hy2-hy1)/\
                        (self.mesh.pointMatrix[__n+self.mesh.nx][1]-self.mesh.pointMatrix[__n][1])+\
                        hy2
                else: #(x,y) is above the node
                    for k in range(self.phyMed.timeSteps):
                        hy1=(x-self.mesh.pointMatrix[__n-self.mesh.nx][0])*\
                        (self.phyMed.h_n[k][__n+1-self.mesh.nx]-self.phyMed.h_n[k][__n-self.mesh.nx])/\
                        (self.mesh.pointMatrix[__n+1-self.mesh.nx][0]-self.mesh.pointMatrix[__n-self.mesh.nx][0])+\
                        self.phyMed.h_n[k][__n-self.mesh.nx]
                        hy2=(x-self.mesh.pointMatrix[__n][0])*\
                        (self.phyMed.h_n[k][__n+1]-self.phyMed.h_n[k][__n])/\
                        (self.mesh.pointMatrix[__n+1][0]-self.mesh.pointMatrix[__n][0])+\
                        self.phyMed.h_n[k][__n]
                        __h[k]=hy2-\
                        (self.mesh.pointMatrix[__n][1]-y)*\
                        (hy2-hy1)/\
                        (self.mesh.pointMatrix[__n][1]-self.mesh.pointMatrix[__n-self.mesh.nx][1])
            else: #(x,y) is to the left of the node
                if self.mesh.pointMatrix[__n][0]<y: #(x,y) is below the node
                    for k in range(self.phyMed.timeSteps):
                        hy1=self.phyMed.h_n[k][__n+self.mesh.nx]-\
                        (self.mesh.pointMatrix[__n+self.mesh.nx][0]-x)*\
                        (self.phyMed.h_n[k][__n+self.mesh.nx]-self.phyMed.h_n[k][__n-1+self.mesh.nx])/\
                        (self.mesh.pointMatrix[__n+self.mesh.nx][0]-self.mesh.pointMatrix[__n-1+self.mesh.nx][0])
                        hy2=self.phyMed.h_n[k][__n]-\
                        (self.mesh.pointMatrix[__n][0]-x)*\
                        (self.phyMed.h_n[k][__n]-self.phyMed.h_n[k][__n-1])/\
                        (self.mesh.pointMatrix[__n][0]-self.mesh.pointMatrix[__n-1][0])
                        __h[k]=(y-self.mesh.pointMatrix[__n][1])*\
                        (hy2-hy1)/\
                        (self.mesh.pointMatrix[__n+self.mesh.nx][1]-self.mesh.pointMatrix[__n][1])+\
                        hy2
                else: #(x,y) is above the node
                    for k in range(self.phyMed.timeSteps):
                        hy1=self.phyMed.h_n[k][__n-self.mesh.nx]-\
                        (self.mesh.pointMatrix[__n-self.mesh.nx][0]-x)*\
                        (self.phyMed.h_n[k][__n-self.mesh.nx]-self.phyMed.h_n[k][__n-1-self.mesh.nx])/\
                        (self.mesh.pointMatrix[__n-self.mesh.nx][0]-self.mesh.pointMatrix[__n-1-self.mesh.nx][0])
                        hy2=self.phyMed.h_n[k][__n]-\
                        (self.mesh.pointMatrix[__n][0]-x)*\
                        (self.phyMed.h_n[k][__n]-self.phyMed.h_n[k][__n-1])/\
                        (self.mesh.pointMatrix[__n][0]-self.mesh.pointMatrix[__n-1][0])
                        __h[k]=hy2-\
                        (self.mesh.pointMatrix[__n][1]-y)*\
                        (hy2-hy1)/\
                        (self.mesh.pointMatrix[__n][1]-self.mesh.pointMatrix[__n-self.mesh.nx][1])
        ##INTERPOLATE HEAD VALUE AT X,Y#################
        ##COMPUTE DRAWDOWN VALUES##############################################
        drawdown=[0]*self.phyMed.timeSteps #Initialize array
        for k in range(self.phyMed.timeSteps):
            drawdown[k]=__h[0]-__h[k]
        ##COMPUTE DRAWDOWN VALUES##############################################
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
        k=step
        for i in range(npts):
            x=obs_pts[i][0]
            y=obs_pts[i][1]
            __n=self.mesh.coord2node(x,y) #Get closest node index
            ##INTERPOLATE HEAD VALUE AT X,Y#################
            ###CHECK HOW WE SHOULD INTERPOLATE THE HEAD VALUE
            #Do not interpolate if the node is at (x,y)
            if self.mesh.pointMatrix[__n][0]==x and self.mesh.pointMatrix[__n][1]==y: #Do not interpolate
                __h0=self.phyMed.h_n[0][__n]
                __hk=self.phyMed.h_n[k][__n]
            #Use a linear interpolation if the node is on 'x' or 'y'
            elif self.mesh.pointMatrix[__n][1]==y: #Interpolate over 'x'
                if self.mesh.pointMatrix[__n][0]<x: #(x,y) is to the right of the node
                    __h0=(x-self.mesh.pointMatrix[__n][0])*\
                    (self.phyMed.h_n[0][__n+1]-self.phyMed.h_n[0][__n])/\
                    (self.mesh.pointMatrix[__n+1][0]-self.mesh.pointMatrix[__n][0])+\
                    self.phyMed.h_n[0][__n]
                    __hk=(x-self.mesh.pointMatrix[__n][0])*\
                    (self.phyMed.h_n[k][__n+1]-self.phyMed.h_n[k][__n])/\
                    (self.mesh.pointMatrix[__n+1][0]-self.mesh.pointMatrix[__n][0])+\
                    self.phyMed.h_n[k][__n]
                else: #(x,y) is to the left of the node
                    __h0=self.phyMed.h_n[0][__n]-\
                    (self.mesh.pointMatrix[__n][0]-x)*\
                    (self.phyMed.h_n[0][__n]-self.phyMed.h_n[0][__n-1])/\
                    (self.mesh.pointMatrix[__n][0]-self.mesh.pointMatrix[__n-1][0])
                    __hk=self.phyMed.h_n[k][__n]-\
                    (self.mesh.pointMatrix[__n][0]-x)*\
                    (self.phyMed.h_n[k][__n]-self.phyMed.h_n[k][__n-1])/\
                    (self.mesh.pointMatrix[__n][0]-self.mesh.pointMatrix[__n-1][0])
            elif self.mesh.pointMatrix[__n][0]==x: #Interpolate over 'y'
                if self.mesh.pointMatrix[__n][0]<y: #(x,y) is below the node
                    __h0=(y-self.mesh.pointMatrix[__n][1])*\
                    (self.phyMed.h_n[0][__n+self.mesh.nx]-self.phyMed.h_n[0][__n])/\
                    (self.mesh.pointMatrix[__n+self.mesh.nx][1]-self.mesh.pointMatrix[__n][1])+\
                    self.phyMed.h_n[0][__n]
                    __hk=(y-self.mesh.pointMatrix[__n][1])*\
                    (self.phyMed.h_n[k][__n+self.mesh.nx]-self.phyMed.h_n[k][__n])/\
                    (self.mesh.pointMatrix[__n+self.mesh.nx][1]-self.mesh.pointMatrix[__n][1])+\
                    self.phyMed.h_n[k][__n]
                else: #(x,y) is above the node
                    __h0=self.phyMed.h_n[0][__n]-\
                    (self.mesh.pointMatrix[__n][1]-y)*\
                    (self.phyMed.h_n[0][__n]-self.phyMed.h_n[0][__n-self.mesh.nx])/\
                    (self.mesh.pointMatrix[__n][1]-self.mesh.pointMatrix[__n-self.mesh.nx][1])
                    __hk=self.phyMed.h_n[k][__n]-\
                    (self.mesh.pointMatrix[__n][1]-y)*\
                    (self.phyMed.h_n[k][__n]-self.phyMed.h_n[k][__n-self.mesh.nx])/\
                    (self.mesh.pointMatrix[__n][1]-self.mesh.pointMatrix[__n-self.mesh.nx][1])
            else:
                if self.mesh.pointMatrix[__n][0]<x: #(x,y) is to the right of the node
                    if self.mesh.pointMatrix[__n][0]<y: #(x,y) is below the node
                        #Time-step=0
                        hy1_0=(x-self.mesh.pointMatrix[__n+self.mesh.nx][0])*\
                        (self.phyMed.h_n[0][__n+1+self.mesh.nx]-self.phyMed.h_n[0][__n+self.mesh.nx])/\
                        (self.mesh.pointMatrix[__n+1+self.mesh.nx][0]-self.mesh.pointMatrix[__n+self.mesh.nx][0])+\
                        self.phyMed.h_n[0][__n+self.mesh.nx]
                        hy2_0=(x-self.mesh.pointMatrix[__n][0])*\
                        (self.phyMed.h_n[0][__n+1]-self.phyMed.h_n[0][__n])/\
                        (self.mesh.pointMatrix[__n+1][0]-self.mesh.pointMatrix[__n][0])+\
                        self.phyMed.h_n[0][__n]
                        __h0=(y-self.mesh.pointMatrix[__n][1])*\
                        (hy2_0-hy1_0)/\
                        (self.mesh.pointMatrix[__n+self.mesh.nx][1]-self.mesh.pointMatrix[__n][1])+\
                        hy2_0
                        #Time-step=k
                        hy1_k=(x-self.mesh.pointMatrix[__n+self.mesh.nx][0])*\
                        (self.phyMed.h_n[k][__n+1+self.mesh.nx]-self.phyMed.h_n[k][__n+self.mesh.nx])/\
                        (self.mesh.pointMatrix[__n+1+self.mesh.nx][0]-self.mesh.pointMatrix[__n+self.mesh.nx][0])+\
                        self.phyMed.h_n[k][__n+self.mesh.nx]
                        hy2_k=(x-self.mesh.pointMatrix[__n][0])*\
                        (self.phyMed.h_n[k][__n+1]-self.phyMed.h_n[k][__n])/\
                        (self.mesh.pointMatrix[__n+1][0]-self.mesh.pointMatrix[__n][0])+\
                        self.phyMed.h_n[k][__n]
                        __hk=(y-self.mesh.pointMatrix[__n][1])*\
                        (hy2_k-hy1_k)/\
                        (self.mesh.pointMatrix[__n+self.mesh.nx][1]-self.mesh.pointMatrix[__n][1])+\
                        hy2_k
                    else: #(x,y) is above the node
                        #Time-step=0
                        hy1_0=(x-self.mesh.pointMatrix[__n-self.mesh.nx][0])*\
                        (self.phyMed.h_n[0][__n+1-self.mesh.nx]-self.phyMed.h_n[0][__n-self.mesh.nx])/\
                        (self.mesh.pointMatrix[__n+1-self.mesh.nx][0]-self.mesh.pointMatrix[__n-self.mesh.nx][0])+\
                        self.phyMed.h_n[0][__n-self.mesh.nx]
                        hy2_0=(x-self.mesh.pointMatrix[__n][0])*\
                        (self.phyMed.h_n[0][__n+1]-self.phyMed.h_n[0][__n])/\
                        (self.mesh.pointMatrix[__n+1][0]-self.mesh.pointMatrix[__n][0])+\
                        self.phyMed.h_n[0][__n]
                        __h0=hy2_0-\
                        (self.mesh.pointMatrix[__n][1]-y)*\
                        (hy2_0-hy1_0)/\
                        (self.mesh.pointMatrix[__n][1]-self.mesh.pointMatrix[__n-self.mesh.nx][1])
                        #Time-step=k
                        hy1_k=(x-self.mesh.pointMatrix[__n-self.mesh.nx][0])*\
                        (self.phyMed.h_n[k][__n+1-self.mesh.nx]-self.phyMed.h_n[k][__n-self.mesh.nx])/\
                        (self.mesh.pointMatrix[__n+1-self.mesh.nx][0]-self.mesh.pointMatrix[__n-self.mesh.nx][0])+\
                        self.phyMed.h_n[k][__n-self.mesh.nx]
                        hy2_k=(x-self.mesh.pointMatrix[__n][0])*\
                        (self.phyMed.h_n[k][__n+1]-self.phyMed.h_n[k][__n])/\
                        (self.mesh.pointMatrix[__n+1][0]-self.mesh.pointMatrix[__n][0])+\
                        self.phyMed.h_n[k][__n]
                        __hk=hy2_k-\
                        (self.mesh.pointMatrix[__n][1]-y)*\
                        (hy2_k-hy1_k)/\
                        (self.mesh.pointMatrix[__n][1]-self.mesh.pointMatrix[__n-self.mesh.nx][1])
                else: #(x,y) is to the left of the node
                    if self.mesh.pointMatrix[__n][0]<y: #(x,y) is below the node
                        #Time-step=0
                        hy1_0=self.phyMed.h_n[0][__n+self.mesh.nx]-\
                        (self.mesh.pointMatrix[__n+self.mesh.nx][0]-x)*\
                        (self.phyMed.h_n[0][__n+self.mesh.nx]-self.phyMed.h_n[0][__n-1+self.mesh.nx])/\
                        (self.mesh.pointMatrix[__n+self.mesh.nx][0]-self.mesh.pointMatrix[__n-1+self.mesh.nx][0])
                        hy2_0=self.phyMed.h_n[0][__n]-\
                        (self.mesh.pointMatrix[__n][0]-x)*\
                        (self.phyMed.h_n[0][__n]-self.phyMed.h_n[0][__n-1])/\
                        (self.mesh.pointMatrix[__n][0]-self.mesh.pointMatrix[__n-1][0])
                        __h0=(y-self.mesh.pointMatrix[__n][1])*\
                        (hy2_0-hy1_0)/\
                        (self.mesh.pointMatrix[__n+self.mesh.nx][1]-self.mesh.pointMatrix[__n][1])+\
                        hy2_0
                        #Time-step=k
                        hy1_k=self.phyMed.h_n[k][__n+self.mesh.nx]-\
                        (self.mesh.pointMatrix[__n+self.mesh.nx][0]-x)*\
                        (self.phyMed.h_n[k][__n+self.mesh.nx]-self.phyMed.h_n[k][__n-1+self.mesh.nx])/\
                        (self.mesh.pointMatrix[__n+self.mesh.nx][0]-self.mesh.pointMatrix[__n-1+self.mesh.nx][0])
                        hy2_k=self.phyMed.h_n[k][__n]-\
                        (self.mesh.pointMatrix[__n][0]-x)*\
                        (self.phyMed.h_n[k][__n]-self.phyMed.h_n[k][__n-1])/\
                        (self.mesh.pointMatrix[__n][0]-self.mesh.pointMatrix[__n-1][0])
                        __hk=(y-self.mesh.pointMatrix[__n][1])*\
                        (hy2_k-hy1_k)/\
                        (self.mesh.pointMatrix[__n+self.mesh.nx][1]-self.mesh.pointMatrix[__n][1])+\
                        hy2_k
                    else: #(x,y) is above the node
                        #Time-step=0
                        hy1_0=self.phyMed.h_n[0][__n-self.mesh.nx]-\
                        (self.mesh.pointMatrix[__n-self.mesh.nx][0]-x)*\
                        (self.phyMed.h_n[0][__n-self.mesh.nx]-self.phyMed.h_n[0][__n-1-self.mesh.nx])/\
                        (self.mesh.pointMatrix[__n-self.mesh.nx][0]-self.mesh.pointMatrix[__n-1-self.mesh.nx][0])
                        hy2_0=self.phyMed.h_n[0][__n]-\
                        (self.mesh.pointMatrix[__n][0]-x)*\
                        (self.phyMed.h_n[0][__n]-self.phyMed.h_n[0][__n-1])/\
                        (self.mesh.pointMatrix[__n][0]-self.mesh.pointMatrix[__n-1][0])
                        __h0=hy2_0-\
                        (self.mesh.pointMatrix[__n][1]-y)*\
                        (hy2_0-hy1_0)/\
                        (self.mesh.pointMatrix[__n][1]-self.mesh.pointMatrix[__n-self.mesh.nx][1])
                        #Time-step=k
                        hy1_k=self.phyMed.h_n[k][__n-self.mesh.nx]-\
                        (self.mesh.pointMatrix[__n-self.mesh.nx][0]-x)*\
                        (self.phyMed.h_n[k][__n-self.mesh.nx]-self.phyMed.h_n[k][__n-1-self.mesh.nx])/\
                        (self.mesh.pointMatrix[__n-self.mesh.nx][0]-self.mesh.pointMatrix[__n-1-self.mesh.nx][0])
                        hy2_k=self.phyMed.h_n[k][__n]-\
                        (self.mesh.pointMatrix[__n][0]-x)*\
                        (self.phyMed.h_n[k][__n]-self.phyMed.h_n[k][__n-1])/\
                        (self.mesh.pointMatrix[__n][0]-self.mesh.pointMatrix[__n-1][0])
                        __hk=hy2_k-\
                        (self.mesh.pointMatrix[__n][1]-y)*\
                        (hy2_k-hy1_k)/\
                        (self.mesh.pointMatrix[__n][1]-self.mesh.pointMatrix[__n-self.mesh.nx][1])
            ##INTERPOLATE HEAD VALUE AT X,Y#################
            ##COMPUTE DRAWDOWN VALUES##########################################
            drawdown[i]=__h0-__hk
            ##COMPUTE DRAWDOWN VALUES##########################################
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
###############################################################################
#####OUPUTS####################################################################
###############################################################################
"""
###########################FDM OBJECT CLASS SECTION############################
"""