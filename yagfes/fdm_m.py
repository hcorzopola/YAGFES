#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: hcorzopola

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
    x_axis = array with 'x' coords  | y_axis = array with 'y' coords
    """
    def __init__(self,file):
        #READ MESH FILE########################################################
        __f=aux_f.read_dict(file) #Open the given file
        
        ##BUILD MESH###########################################################
        ###MESH CHARACTERISTICS################################################
        self.dx=__f['dx']
        self.nx=round((__f['max_x']-__f['min_x'])/self.dx)+1
        self.dy=__f['dy']
        self.ny=round((__f['max_y']-__f['min_y'])/self.dy)+1
        ###MESH COORDINATES####################################################
        __n_n=self.nx*self.ny
        self.m_p=[0]*__n_n #Initialize the array containing the node coordinates
        ########################
        for i in range(self.nx):
            for j in range(self.ny):
                __node=i+j*self.nx
                self.m_p[__node]=[i*self.dx+__f['min_x'],j*self.dy+__f['min_y']]
        ##BUILD MESH###########################################################
        ##INITIALIZE AUXILIAR ARRAY TO KEEP TRACK OF NODE ID's#################
        ##CREATE AUXILIAR LIST TO KEEP TRACK OF ACTUAL NODE INDEXES############
        self.ID=[0]*__n_n
        for i in range(__n_n):
            self.ID[i]=i
        #READ MESH FILE########################################################
        
    def write_bc_file(self,file):
        #WRITE A BC FILE BASED ON MESH INFO####################################
        __f=open(aux_f.path+'/'+file,'w') #Open the given file in 'write' mode
        
        __f.write('#BC FILE\n') #Write document title
        ##MESH LEFT
        __f.write('#MESH LEFT\n')
        __f.write('x_min=[')
        __bc_type=input('Please input the boundary type at the left boundary of the mesh (x='+str(self.m_p[0][0])+'):\n 1) Dirichlet\n 2) Neumann\n')
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
                __v=input('Please, input the Dirichlet BC at the first node of the left boundary (y='+str(self.m_p[0][1])+').\n') #Read Dirichlet BC value
                __v2=input('Please, input the Dirichlet BC at the last node of the left boundary (y='+str(self.m_p[self.nx*(self.ny-1)][1])+').\n') #Read Dirichlet BC value
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
        __bc_type=input('Please input the boundary type at the right boundary of the mesh (x='+str(self.m_p[self.nx-1][0])+'):\n 1) Dirichlet\n 2) Neumann\n')
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
                __v=input('Please, input the Dirichlet BC at the first node of the right boundary (y='+str(self.m_p[0][1])+').\n') #Read Dirichlet BC value
                __v2=input('Please, input the Dirichlet BC at the last node of the right boundary (y='+str(self.m_p[self.nx*(self.ny-1)][1])+').\n') #Read Dirichlet BC value
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
        __bc_type=input('Please input the boundary type at the top of the mesh (y='+str(self.m_p[0][1])+'):\n 1) Dirichlet\n 2) Neumann\n')
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
                __v=input('Please, input the Dirichlet BC at the first node of the top boundary (x='+str(self.m_p[0][0])+').\n') #Read Dirichlet BC value
                __v2=input('Please, input the Dirichlet BC at the last node of the top boundary (x='+str(self.m_p[self.nx-1][0])+').\n') #Read Dirichlet BC value
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
        __bc_type=input('Please input the boundary type at the bottom of the mesh (y='+str(self.m_p[self.nx*(self.ny-1)][1])+'):\n 1) Dirichlet\n 2) Neumann\n')
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
                __v=input('Please, input the Dirichlet BC at the first node of the bottom boundary (x='+str(self.m_p[0][0])+').\n') #Read Dirichlet BC value
                __v2=input('Please, input the Dirichlet BC at the last node of the bottom boundary (x='+str(self.m_p[self.nx-1][0])+').\n') #Read Dirichlet BC value
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
                __d=((self.m_p[__node_ID][0]-x)**2+\
                     (self.m_p[__node_ID][1]-y)**2)**0.5
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
    k = Hydraulic Conductivity    | ss = Specific Storage
    bc_d = Dirichlet BC's         | bc_n = Neumann BC's
    hh = Hydraulic Head           | q = Pumping/Recharge Well Rate
    """
    def __init__(self,config,n_n):
        #RETRIEVE MEDIUM INFO FROM THE CONFIG FILE#############################
        self.steady=config['steady']
        ##IF WE HAVE A TRANSIENT MEDIUM (STEADY=FALSE), READ TIME PARAMETERS###
        if self.steady==False:
            self.dt=config['dt']
            self.time_steps=int(round(config['total_time']/config['dt']))+1
        #######################################################################
        #######################################################################
        #INITIALIZE ALL THE ARRAYS USED FOR STORING PHYSICAL PROPERTIES########
        ##INITIALIZE HYDRAULIC HEAD ARRAY######################################
        self.hh_n=[0]*(n_n)
        #######################################################################
        ##INITIALIZE HYDRAULIC CONDUCTIVITY ARRAY##############################
        self.k=[1]*(n_n)
        for i in range(n_n):
            self.k[i]=[float(config['Kx']),float(config['Ky'])] #Hydraulic Conductivity Array ((nx x ny) x 2)
        #######################################################################
        ##INITIALIZE SPECIFIC STORAGE ARRAY####################################
        self.ss=[float(config['Ss'])]*(n_n)
        #######################################################################
        ##INITIALIZE WELL CONDITIONS ARRAY#####################################
        self.q=[0]*(n_n)
        #######################################################################
        ##ADD TEMPORAL AXIS IF MEDIUM IS IN TRANSIENT-STATE####################
        if self.steady==False: #Check steady-state flag
            ###TEMPORAL AXIS FOR NODE PROPERTIES
            for i in range(n_n):
                self.q[i]=[0]*self.time_steps
            #############################
            self.hh_n=[0]*self.time_steps
            for k in range(self.time_steps):
                self.hh_n[k]=[0]*n_n
        #######################################################################
        #INITIALIZE ALL THE ARRAYS USED FOR STORING BOUNDARY CONDITIONS########
        self.bc_d=dict()
        self.bc_n=dict()
        #######################################################################
        #INITIALIZE DICTIONARY FOR STORING LAYER/ZONE PROPERTIES###############
        ##INITIALIZE AUXILIAR ARRAY FOR STORING NODE IDS
        __nodes=[0]*n_n
        for i in range(n_n):
            __nodes[i]=i
        ##DEFINE DICTIONARY WITH ZONE PROPERTIES
        self.zones={0:{'Ss':float(config['Ss']),
                       'Kx':float(config['Kx']),
                       'Ky':float(config['Ky']),
                       'nodes':__nodes}}
        #######################################################################
        del __nodes,n_n,config
    
    #READ BOUNDARY CONDITIONS##################################################
    def read_bc_file(self,file):
        #READ BC FILE##########################################################
        __f=aux_f.read_dict(file) #Open the given file in 'read' mode
        
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
        del __f
    #READ BOUNDARY CONDITIONS##################################################
        
    #READ INITIAL CONDITIONS###################################################
    def read_hh0(self,file):
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
                self.hh_n[0][i+j*__ncol]=float(__aux[i])
        
        __f.close()
        #READS AN HHD FILE AND ASSIGNS IT TO HH_0##############################
        del __nrow,__ncol,__aux
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
    """
    def __init__(self,init):
        self.type='FDM'
        #READ INIT FILE########################################################
        self.config=aux_f.read_dict(init)
        #READ INIT FILE########################################################
        #CREATE MESH AND PHYMED OBJECTS########################################
        self.mesh=mesh(self.config['mesh_file'])
        self.phymed=phymed(self.config,self.mesh.nx*self.mesh.ny)
        #CREATE MESH AND PHYMED OBJECTS########################################
        #INITIALIZE MATRICES ARRAYS
        self.m_coeff=list() #INITIALIZE ARRAY FOR COEFFICIENTS MATRIX
        self.v_load=list() #INITIALIZE ARRAY FOR LOAD VECTOR
        #######################################################################
    
###############################################################################
#####ASSEMBLE MATRICES AND VECTORS#############################################
###############################################################################
    ##ASSEMBLE COEFFICIENTS MATRIX#############################################
    def coeff_m_assembly(self):
        #INITIALIZE ARRAY FOR COEFFICIENTS MATRIX
        __n_n=self.mesh.nx*self.mesh.ny
        self.m_coeff=[0]*__n_n
        for i in range(__n_n):
            self.m_coeff[i]=[0]*__n_n
        #########################################
        for j in range(self.mesh.ny):
            for i in range(self.mesh.nx):
                __node=i+j*self.mesh.nx
                ##COMPUTE AND ASSIGN X-COEFFICIENTS############################
                if i==0: #x_min
                    ##GET 'EQUIVALENT CONDUCTIVIES'############################
                    __Kxb=self.phymed.k[__node][0]
                    __Kxf=2/(1/self.phymed.k[__node][0]+1/self.phymed.k[__node+1][0])
                    ##GET 'EQUIVALENT CONDUCTIVIES'############################
                    ##GET 'CONDUCTANCE' COEFFICIENTS###########################
                    __Cxb=__Kxb/(self.mesh.dx**2)
                    __Cxf=__Kxf/(self.mesh.dx**2)
                    ##GET 'CONDUCTANCE' COEFFICIENTS###########################
                    ##ASSIGN VALUES TO MATRIX##################################
                    self.m_coeff[__node][__node+1]=(__Cxb+__Cxf)/self.phymed.ss[__node] #X+ FLOW
                    ##ASSIGN VALUES TO MATRIX##################################
                elif i==self.mesh.nx-1: #x_max
                    ##GET 'EQUIVALENT CONDUCTIVIES'############################
                    __Kxb=2/(1/self.phymed.k[__node][0]+1/self.phymed.k[__node-1][0])
                    __Kxf=self.phymed.k[__node][0]
                    ##GET 'EQUIVALENT CONDUCTIVIES'############################
                    ##GET 'CONDUCTANCE' COEFFICIENTS###########################
                    __Cxb=__Kxb/(self.mesh.dx**2)
                    __Cxf=__Kxf/(self.mesh.dx**2)
                    ##GET 'CONDUCTANCE' COEFFICIENTS###########################
                    ##ASSIGN VALUES TO MATRIX##################################
                    self.m_coeff[__node][__node-1]=(__Cxb+__Cxf)/self.phymed.ss[__node] #X- FLOW
                    ##ASSIGN VALUES TO MATRIX##################################
                else:
                    ##GET 'EQUIVALENT CONDUCTIVIES'############################
                    __Kxb=2/(1/self.phymed.k[__node][0]+1/self.phymed.k[__node-1][0])
                    __Kxf=2/(1/self.phymed.k[__node][0]+1/self.phymed.k[__node+1][0])
                    ##GET 'EQUIVALENT CONDUCTIVIES'############################
                    ##GET 'CONDUCTANCE' COEFFICIENTS###########################
                    __Cxb=__Kxb/(self.mesh.dx**2)
                    __Cxf=__Kxf/(self.mesh.dx**2)
                    ##GET 'CONDUCTANCE' COEFFICIENTS###########################
                    ##ASSIGN VALUES TO MATRIX##################################
                    self.m_coeff[__node][__node-1]=__Cxb/self.phymed.ss[__node] #X- FLOW
                    self.m_coeff[__node][__node+1]=__Cxf/self.phymed.ss[__node] #X+ FLOW
                    ##ASSIGN VALUES TO MATRIX##################################
                ##COMPUTE AND ASSIGN X-COEFFICIENTS############################
                
                ##COMPUTE Y-COEFFICIENTS#######################################
                if j==0: #y_min
                    ##GET 'EQUIVALENT CONDUCTIVIES'############################
                    __Kyb=self.phymed.k[__node][1]
                    __Kyf=2/(1/self.phymed.k[__node][1]+1/self.phymed.k[__node+self.mesh.nx][1])
                    ##GET 'EQUIVALENT CONDUCTIVIES'############################
                    ##GET 'CONDUCTANCE' COEFFICIENTS###########################
                    __Cyb=__Kyb/(self.mesh.dy**2)
                    __Cyf=__Kyf/(self.mesh.dy**2)
                    ##GET 'CONDUCTANCE' COEFFICIENTS###########################
                    ##ASSIGN VALUES TO MATRIX##################################
                    self.m_coeff[__node][__node+self.mesh.nx]=(__Cyb+__Cyf)/self.phymed.ss[__node] #Y+ FLOW
                    ##ASSIGN VALUES TO MATRIX##################################
                elif j==self.mesh.ny-1: #y_max
                    ##GET 'EQUIVALENT CONDUCTIVIES'############################
                    __Kyb=2/(1/self.phymed.k[__node][1]+1/self.phymed.k[__node-self.mesh.nx][1])
                    __Kyf=self.phymed.k[__node][1]
                    ##GET 'EQUIVALENT CONDUCTIVIES'############################
                    ##GET 'CONDUCTANCE' COEFFICIENTS###########################
                    __Cyb=__Kyb/(self.mesh.dy**2)
                    __Cyf=__Kyf/(self.mesh.dy**2)
                    ##GET 'CONDUCTANCE' COEFFICIENTS###########################
                    ##ASSIGN VALUES TO MATRIX##################################
                    self.m_coeff[__node][__node-self.mesh.nx]=(__Cyb+__Cyf)/self.phymed.ss[__node] #Y- FLOW
                    ##ASSIGN VALUES TO MATRIX##################################
                else:
                    ##GET 'EQUIVALENT CONDUCTIVIES'############################
                    __Kyb=2/(1/self.phymed.k[__node][1]+1/self.phymed.k[__node-self.mesh.nx][1])
                    __Kyf=2/(1/self.phymed.k[__node][1]+1/self.phymed.k[__node+self.mesh.nx][1])
                    ##GET 'EQUIVALENT CONDUCTIVIES'############################
                    ##GET 'CONDUCTANCE' COEFFICIENTS###########################
                    __Cyb=__Kyb/(self.mesh.dy**2)
                    __Cyf=__Kyf/(self.mesh.dy**2)
                    ##GET 'CONDUCTANCE' COEFFICIENTS###########################
                    ##ASSIGN VALUES TO MATRIX##################################
                    self.m_coeff[__node][__node-self.mesh.nx]=__Cyb/self.phymed.ss[__node] #Y- FLOW
                    self.m_coeff[__node][__node+self.mesh.nx]=__Cyf/self.phymed.ss[__node] #Y+ FLOW
                    ##ASSIGN VALUES TO MATRIX##################################
                ##COMPUTE Y-COEFFICIENTS#######################################
                
                ##CHECK SYSTEM STATE###########################################
                if self.phymed.steady==True: #Steady-state
                    __Ct=0
                else: #Transient-state
                    __Ct=1/self.phymed.dt
                ##CHECK SYSTEM STATE###########################################
                
                ##ASSIGN VALUES TO CENTRAL NODE################################
                self.m_coeff[__node][__node]=-((__Cxb+__Cxf+__Cyb+__Cyf)/self.phymed.ss[__node]+__Ct)
                ##ASSIGN VALUES TO CENTRAL NODE################################
    ##ASSEMBLE COEFFICIENTS MATRIX#############################################
    
    ##ASSEMBLE LOAD VECTOR#####################################################
    def load_v_assembly(self):
        __n_n=self.mesh.nx*self.mesh.ny
        #Steady-state##########################################################
        if self.phymed.steady==True:
            #INITIALIZE ARRAY FOR LOAD VECTOR
            self.v_load=[0]*__n_n
            #ASSIGN VALUES
            for j in range(self.mesh.ny):
                for i in range(self.mesh.nx):
                    __node=i+j*self.mesh.nx
                    self.v_load[__node]=self.phymed.q[__node]
        #Steady-state##########################################################        
        #Transient-state#######################################################
        else:
            #INITIALIZE ARRAY FOR LOAD VECTOR
            self.v_load=[0]*self.phymed.time_steps
            for i in range(self.phymed.time_steps):
                self.v_load[i]=[0]*__n_n
            #ASSIGN VALUES
            for j in range(self.mesh.ny):
                for i in range(self.mesh.nx):
                    __node=i+j*self.mesh.nx
                    for k in range(self.phymed.time_steps):
                        self.v_load[k][__node]=self.phymed.q[__node][k]
        #Transient-state#######################################################
    ##ASSEMBLE LOAD VECTOR#####################################################
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
        #Steady-State##########################################################
        if self.phymed.steady==True:
            ##MESH LEFT
            if 'x_min' in self.phymed.bc_d:
                for i in range(self.mesh.ny):
                    __node=i*self.mesh.nx
                    if __node in self.mesh.ID:
                        __rID=self.mesh.ID.index(__node) #Retrieves the real ID for __node
                        del self.mesh.ID[__rID] #Remove ID from the ID list
                        ##REMOVE NODE FROM HYDRAULIC HEAD VECTOR###################
                        del self.phymed.hh_n[__rID]
                        ##REMOVE NODE FROM LOAD VECTOR#############################
                        del self.v_load[__rID]
                        ###REMOVE ROW FROM COEFFICIENT MATRIX######################
                        del self.m_coeff[__rID]
                        ##ADDITIONAL OPERATIONS####################################
                        for j in range(len(self.v_load)):
                            self.v_load[j]-=self.m_coeff[j][__rID]*self.phymed.bc_d['x_min'][i]*uF
                            ##REMOVE COLUMN FROM COEFFICIENT MATRIX################
                            del self.m_coeff[j][__rID]
            ##MESH RIGHT
            if 'x_max' in self.phymed.bc_d:
                for i in range(self.mesh.ny):
                    __node=(i+1)*self.mesh.nx-1
                    if __node in self.mesh.ID:
                        __rID=self.mesh.ID.index(__node) #Retrieves the real ID for __node
                        del self.mesh.ID[__rID] #Remove ID from the ID list
                        ##REMOVE NODE FROM HYDRAULIC HEAD VECTOR###################
                        del self.phymed.hh_n[__rID]
                        ##REMOVE NODE FROM LOAD VECTOR#############################
                        del self.v_load[__rID]
                        ###REMOVE ROW FROM COEFFICIENT MATRIX######################
                        del self.m_coeff[__rID]
                        ##ADDITIONAL OPERATIONS####################################
                        for j in range(len(self.v_load)):
                            self.v_load[j]-=self.m_coeff[j][__rID]*self.phymed.bc_d['x_max'][i]*uF
                            ##REMOVE COLUMN FROM COEFFICIENT MATRIX################
                            del self.m_coeff[j][__rID]
            ##MESH TOP
            if 'y_min' in self.phymed.bc_d:
                for i in range(self.mesh.nx):
                    __node=i
                    if __node in self.mesh.ID:
                        __rID=self.mesh.ID.index(__node) #Retrieves the real ID for __node
                        del self.mesh.ID[__rID] #Remove ID from the ID list
                        ##REMOVE NODE FROM HYDRAULIC HEAD VECTOR###################
                        del self.phymed.hh_n[__rID]
                        ##REMOVE NODE FROM LOAD VECTOR#############################
                        del self.v_load[__rID]
                        ###REMOVE ROW FROM COEFFICIENT MATRIX######################
                        del self.m_coeff[__rID]
                        ##ADDITIONAL OPERATIONS####################################
                        for j in range(len(self.v_load)):
                            self.v_load[j]-=self.m_coeff[j][__rID]*self.phymed.bc_d['y_min'][i]*uF
                            ##REMOVE COLUMN FROM COEFFICIENT MATRIX################
                            del self.m_coeff[j][__rID]
            ##MESH BOTTOM
            if 'y_max' in self.phymed.bc_d:
                for i in range(self.mesh.nx):
                    __node=i+self.mesh.nx*(self.mesh.ny-1)
                    if __node in self.mesh.ID:
                        __rID=self.mesh.ID.index(__node) #Retrieves the real ID for __node
                        del self.mesh.ID[__rID] #Remove ID from the ID list
                        ##REMOVE NODE FROM HYDRAULIC HEAD VECTOR###################
                        del self.phymed.hh_n[__rID]
                        ##REMOVE NODE FROM LOAD VECTOR#############################
                        del self.v_load[__rID]
                        ###REMOVE ROW FROM COEFFICIENT MATRIX######################
                        del self.m_coeff[__rID]
                        ##ADDITIONAL OPERATIONS####################################
                        for j in range(len(self.v_load)):
                            self.v_load[j]-=self.m_coeff[j][__rID]*self.phymed.bc_d['y_max'][i]*uF
                            ##REMOVE COLUMN FROM COEFFICIENT MATRIX################
                            del self.m_coeff[j][__rID]
        #Steady-State##########################################################        
        #Transient-state#######################################################
        else:
            ##MESH LEFT
            if 'x_min' in self.phymed.bc_d:
                for i in range(self.mesh.ny):
                    __node=i*self.mesh.nx
                    if __node in self.mesh.ID:
                        __rID=self.mesh.ID.index(__node) #Retrieves the real ID for __node
                        del self.mesh.ID[__rID] #Remove ID from the ID list
                        ##REMOVE NODE FROM HYDRAULIC HEAD VECTOR###################
                        for j in range(self.phymed.time_steps):
                            del self.phymed.hh_n[j][__rID]
                        ##REMOVE NODE FROM LOAD VECTOR#############################
                        for j in range(self.phymed.time_steps):
                            del self.v_load[j][__rID]
                        ###REMOVE ROW FROM COEFFICIENT MATRIX######################
                        del self.m_coeff[__rID]
                        ##ADDITIONAL OPERATIONS####################################
                        for j in range(len(self.v_load[0])):
                            for k in range(self.phymed.time_steps):
                                self.v_load[k][j]-=self.m_coeff[j][__rID]*self.phymed.bc_d['x_min'][i]*uF
                            ##REMOVE COLUMN FROM COEFFICIENT MATRIX################
                            del self.m_coeff[j][__rID]
            ##MESH RIGHT
            if 'x_max' in self.phymed.bc_d:
                for i in range(self.mesh.ny):
                    __node=(i+1)*self.mesh.nx-1
                    if __node in self.mesh.ID:
                        __rID=self.mesh.ID.index(__node) #Retrieves the real ID for __node
                        del self.mesh.ID[__rID] #Remove ID from the ID list
                        ##REMOVE NODE FROM HYDRAULIC HEAD VECTOR###################
                        for j in range(self.phymed.time_steps):
                            del self.phymed.hh_n[j][__rID]
                        ##REMOVE NODE FROM LOAD VECTOR#############################
                        for j in range(self.phymed.time_steps):
                            del self.v_load[j][__rID]
                        ###REMOVE ROW FROM COEFFICIENT MATRIX######################
                        del self.m_coeff[__rID]
                        ##ADDITIONAL OPERATIONS####################################
                        for j in range(len(self.v_load[0])):
                            for k in range(self.phymed.time_steps):
                                self.v_load[k][j]-=self.m_coeff[j][__rID]*self.phymed.bc_d['x_max'][i]*uF
                            ##REMOVE COLUMN FROM COEFFICIENT MATRIX################
                            del self.m_coeff[j][__rID]
            ##MESH TOP
            if 'y_min' in self.phymed.bc_d:
                for i in range(self.mesh.nx):
                    __node=i
                    if __node in self.mesh.ID:
                        __rID=self.mesh.ID.index(__node) #Retrieves the real ID for __node
                        del self.mesh.ID[__rID] #Remove ID from the ID list
                        ##REMOVE NODE FROM HYDRAULIC HEAD VECTOR###################
                        for j in range(self.phymed.time_steps):
                            del self.phymed.hh_n[j][__rID]
                        ##REMOVE NODE FROM LOAD VECTOR#############################
                        for j in range(self.phymed.time_steps):
                            del self.v_load[j][__rID]
                        ###REMOVE ROW FROM COEFFICIENT MATRIX######################
                        del self.m_coeff[__rID]
                        ##ADDITIONAL OPERATIONS####################################
                        for j in range(len(self.v_load[0])):
                            for k in range(self.phymed.time_steps):
                                self.v_load[k][j]-=self.m_coeff[j][__rID]*self.phymed.bc_d['y_min'][i]*uF
                            ##REMOVE COLUMN FROM COEFFICIENT MATRIX################
                            del self.m_coeff[j][__rID]
            ##MESH BOTTOM
            if 'y_max' in self.phymed.bc_d:
                for i in range(self.mesh.nx):
                    __node=i+self.mesh.nx*(self.mesh.ny-1)
                    if __node in self.mesh.ID:
                        __rID=self.mesh.ID.index(__node) #Retrieves the real ID for __node
                        del self.mesh.ID[__rID] #Remove ID from the ID list
                        ##REMOVE NODE FROM HYDRAULIC HEAD VECTOR###################
                        for j in range(self.phymed.time_steps):
                            del self.phymed.hh_n[j][__rID]
                        ##REMOVE NODE FROM LOAD VECTOR#############################
                        for j in range(self.phymed.time_steps):
                            del self.v_load[j][__rID]
                        ###REMOVE ROW FROM COEFFICIENT MATRIX######################
                        del self.m_coeff[__rID]
                        ##ADDITIONAL OPERATIONS####################################
                        for j in range(len(self.v_load[0])):
                            for k in range(self.phymed.time_steps):
                                self.v_load[k][j]-=self.m_coeff[j][__rID]*self.phymed.bc_d['y_max'][i]*uF
                            ##REMOVE COLUMN FROM COEFFICIENT MATRIX################
                            del self.m_coeff[j][__rID]
        #Transient-state#######################################################
    ##APPLY DIRICHLET BC's#####################################################
    
    ##APPLY NEUMANN BC's#######################################################
    def apply_bc_n(self):
        #UNIT CONVERSION#######################################################
        uF=1
        if self.config['length'].lower()!=self.config['head'].lower(): #Check if conversion is required
            uF=(aux_f.unit_dict[self.config['length'].lower()])/(aux_f.unit_dict[self.config['head'].lower()]) #Get conversion factor            
        #UNIT CONVERSION#######################################################
        ##INCREASE MESH SIZE###################################################
        #ASSIGN FICTIONAL HEAD VALUES TO MATCH FLOW CONDITION##################
        if self.phymed.steady==True: #Steady-state
            ###MESH LEFT
            if 'x_min' in self.phymed.bc_n:
                for i in range(self.mesh.ny):
                    __node=i*self.mesh.nx
                    __f=2*self.phymed.bc_n['x_min']*self.config['aq_thickness']/self.mesh.dx
                    self.v_load[__node]-=__f*uF/self.phymed.ss[__node]
            ###MESH RIGHT
            if 'x_max' in self.phymed.bc_n:
                for i in range(self.mesh.ny):
                    __node=(i+1)*self.mesh.nx-1
                    __f=2*self.phymed.bc_n['x_max']*self.config['aq_thickness']/self.mesh.dx
                    self.v_load[__node]-=__f*uF/self.phymed.ss[__node]
            ###MESH TOP
            if 'y_min' in self.phymed.bc_n:
                for i in range(self.mesh.nx):
                    __node=i
                    __f=2*self.phymed.bc_n['y_min']*self.config['aq_thickness']/self.mesh.dy
                    self.v_load[__node]-=__f*uF/self.phymed.ss[__node]
            ###MESH BOTTOM
            if 'y_max' in self.phymed.bc_n:
                for i in range(self.mesh.nx):
                    __node=i+self.mesh.nx*(self.mesh.ny-1)
                    __f=2*self.phymed.bc_n['y_max']*self.config['aq_thickness']/self.mesh.dy
                    self.v_load[__node]-=__f*uF/self.phymed.ss[__node]
        else: #Transient-state
            ###MESH LEFT
            if 'x_min' in self.phymed.bc_n:
                for i in range(self.mesh.ny):
                    __node=i*self.mesh.nx
                    for j in range(self.phymed.time_steps):
                        __f=2*self.phymed.bc_n['x_min']*self.config['aq_thickness']/self.mesh.dx
                        self.v_load[j][__node]-=__f*uF/self.phymed.ss[__node]
            ###MESH RIGHT
            if 'x_max' in self.phymed.bc_n:
                for i in range(self.mesh.ny):
                    __node=(i+1)*self.mesh.nx-1
                    for j in range(self.phymed.time_steps):
                        __f=2*self.phymed.bc_n['x_max']*self.config['aq_thickness']/self.mesh.dx
                        self.v_load[j][__node]-=__f*uF/self.phymed.ss[__node]
            ###MESH TOP
            if 'y_min' in self.phymed.bc_n:
                for i in range(self.mesh.nx):
                    __node=i
                    for j in range(self.phymed.time_steps):
                        __f=2*self.phymed.bc_n['y_min']*self.config['aq_thickness']/self.mesh.dy
                        self.v_load[j][__node]-=__f*uF/self.phymed.ss[__node]
            ###MESH BOTTOM
            if 'y_max' in self.phymed.bc_n:
                for i in range(self.mesh.nx):
                    __node=i+self.mesh.nx*(self.mesh.ny-1)
                    for j in range(self.phymed.time_steps):
                        __f=2*self.phymed.bc_n['y_max']*self.config['aq_thickness']/self.mesh.dy
                        self.v_load[j][__node]-=__f*uF/self.phymed.ss[__node]
        #ASSIGN FICTIONAL HEAD VALUES TO MATCH FLOW CONDITION##################
    ##APPLY NEUMANN BC's#######################################################
    
    ##REBUILD HH ARRAY#########################################################
    def rebuild_hh(self):
        #UNIT CONVERSION#######################################################
        if self.config['length'].lower()!=self.config['head'].lower(): #Check if conversion is required
            uF=aux_f.unit_dict[self.config['head'].lower()]/aux_f.unit_dict[self.config['length'].lower()] #Get conversion factor
            for i in range(len(self.phymed.hh_n)):
                self.phymed.hh_n[i]*=uF
        #UNIT CONVERSION#######################################################
        #Steady-State##########################################################
        if self.phymed.steady==True:
            ##Re-initialize hydraulic head array###############################
            __n_n=self.mesh.nx*self.mesh.ny
            __n=len(self.mesh.ID)
            __hh=[0]*__n
            ###Pass values to new array########################################
            for i in range(__n):
                __hh[i]=self.phymed.hh_n[i]
            ###Reset hydraulic head array and pass values######################
            self.phymed.hh_n=[0]*__n_n
            for i in range(__n):
                self.phymed.hh_n[self.mesh.ID[i]]=__hh[i]
            ##Re-initialize hydraulic head array###############################
            ##MESH LEFT
            if 'x_min' in self.phymed.bc_d:
                for i in range(self.mesh.ny):
                    __node=i*self.mesh.nx
                    self.phymed.hh_n[__node]=self.phymed.bc_d['x_min'][i]
            ##MESH RIGHT
            if 'x_max' in self.phymed.bc_d:
                for i in range(self.mesh.ny):
                    __node=(i+1)*self.mesh.nx-1
                    self.phymed.hh_n[__node]=self.phymed.bc_d['x_max'][i]
            ##MESH TOP
            if 'y_min' in self.phymed.bc_d:
                for i in range(self.mesh.nx):
                    __node=i
                    self.phymed.hh_n[__node]=self.phymed.bc_d['y_min'][i]
            ##MESH BOTTOM
            if 'y_max' in self.phymed.bc_d:
                for i in range(self.mesh.nx):
                    __node=i+self.mesh.nx*(self.mesh.ny-1)
                    self.phymed.hh_n[__node]=self.phymed.bc_d['y_max'][i]
        #Steady-State##########################################################
        #Transient-State#######################################################
        else:
            __n_n=self.mesh.nx*self.mesh.ny
            __n=len(self.mesh.ID)
            __hh=[0]*__n
            for j in range(self.phymed.time_steps):
                ##Re-initialize hydraulic head array###########################
                ###Pass values to new array####################################
                for i in range(__n):
                    __hh[i]=self.phymed.hh_n[j][i]
                ###Reset hydraulic head array and pass values##################
                self.phymed.hh_n[j]=[0]*__n_n
                for i in range(__n):
                    self.phymed.hh_n[j][self.mesh.ID[i]]=__hh[i]
                ##Re-initialize hydraulic head array###########################
                ##MESH LEFT
                if 'x_min' in self.phymed.bc_d:
                    for i in range(self.mesh.ny):
                        __node=i*self.mesh.nx
                        self.phymed.hh_n[j][__node]=self.phymed.bc_d['x_min'][i]
                ##MESH RIGHT
                if 'x_max' in self.phymed.bc_d:
                    for i in range(self.mesh.ny):
                        __node=(i+1)*self.mesh.nx-1
                        self.phymed.hh_n[j][__node]=self.phymed.bc_d['x_max'][i]
                ##MESH TOP
                if 'y_min' in self.phymed.bc_d:
                    for i in range(self.mesh.nx):
                        __node=i
                        self.phymed.hh_n[j][__node]=self.phymed.bc_d['y_min'][i]
                ##MESH BOTTOM
                if 'y_max' in self.phymed.bc_d:
                    for i in range(self.mesh.nx):
                        __node=i+self.mesh.nx*(self.mesh.ny-1)
                        self.phymed.hh_n[j][__node]=self.phymed.bc_d['y_max'][i]
        #Transient-State#######################################################
    ##REBUILD HH ARRAY#########################################################
###############################################################################
#####APPLY BOUNDARY CONDITIONS#################################################
###############################################################################

###############################################################################
#####MODIFY PHYSICAL MEDIUM####################################################
###############################################################################
    ##ADD WELLS################################################################
    def add_well(self,file):
        #UNIT CONVERSION#######################################################
        uF=1
        if self.config['length'].lower()!=self.config['head'].lower(): #Check if conversion is required
            uF=(aux_f.unit_dict[self.config['length'].lower()]**3)/(aux_f.unit_dict[self.config['head'].lower()]**3) #Get conversion factor            
        #UNIT CONVERSION#######################################################
        ##RETRIEVE WELL INFORMATION FROM FILE##################################
        __f=aux_f.read_dict(file)
        ##RETRIEVE WELL INFORMATION FROM FILE##################################
        __node=self.mesh.coord2node(__f['x'],__f['y']) #Retrieve node indexes
        ##ASSIGN FLOW VALUE TO 'Q' ARRAY#######################################
        if self.phymed.steady==True: #Steady-state
            self.phymed.q[__node]=-__f['rate']*uF/\
            (self.mesh.dx*self.mesh.dy*self.config['aq_thickness'])/self.phymed.ss[__node]
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
                #ASSIGN FLOW VALUE AT INTERVAL 'i'#############################
                for k in range(t_min,t_max):
                    self.phymed.q[__node][k]=-__f['rate'][i][2]*uF/\
                    (self.mesh.dx*self.mesh.dy*self.config['aq_thickness'])/self.phymed.ss[__node]
    ##ADD WELLS################################################################
    
    #ZONE RELATED FUNCTIONS####################################################
    ##DEFINE NEW ZONE##########################################################
    def new_zone(self,zone_id,Ss:float,Kx:float,Ky:float,x_min:float,x_max:float,y_min:float,y_max:float): #Defines a new zone in the domain as a rectangle
        if zone_id not in self.phymed.zones: #CHECK IF THE ZONE EXISTS
            #IDENTIFY ALL NODES INSIDE THE ZONE################################
            __x_min=self.mesh.m_p[0][0]
            __y_min=self.mesh.m_p[0][1]
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
                    self.phymed.k[__node]=[Kx,Ky]
                    self.phymed.ss[__node]=Ss
                    __ID.append(__node) #Append value to auxiliar list
                    #Look for the node in all other zones and remove it from their lists
                    for k in self.phymed.zones:
                        if __node in self.phymed.zones[k]['nodes']: #Check if node is in a zone 'k'
                            __node_ID=self.phymed.zones[k]['nodes'].index(__node) #Get ID within 'nodes' list
                            del self.phymed.zones[k]['nodes'][__node_ID] #Remove node from zone 'k'
            self.phymed.zones[zone_id]={'Ss':Ss,
                             'Kx':Kx,
                             'Ky':Ky,
                             'nodes':__ID}
            #ZONE ASSIGNATION LOOP#############################################
            del zone_id,Ss,Kx,Ky,x_min,x_max,y_min,y_max,i_min,i_max,j_min,j_max,__node,__node_ID,__ID
        else:
            print('Zone {0} already exists'.format(zone_id))
            del zone_id,Ss,Kx,Ky,x_min,x_max,y_min,y_max
    ##DEFINE NEW ZONE##########################################################
    
    ##APPEND NEW RECTANGLE TO EXISTING ZONE####################################
    def expand_zone(self,zone_id,x_min:float,x_max:float,y_min:float,y_max:float): #Adds a new rectangle to an existing zone
        if zone_id in self.phymed.zones: #CHECK IF THE ZONE EXISTS
            #IDENTIFY ALL NODES INSIDE THE ZONE################################
            __x_min=self.mesh.m_p[0][0]
            __y_min=self.mesh.m_p[0][1]
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
                    __ID.append(__node) #Append value to auxiliar list
                    #Look for the node in all other zones and remove it from their lists
                    for k in self.phymed.zones:
                        if __node in self.phymed.zones[k]['nodes']: #Check if node is in a zone 'k'
                            __node_ID=self.phymed.zones[k]['nodes'].index(__node) #Get ID within 'nodes' list
                            del self.phymed.zones[k]['nodes'][__node_ID] #Remove node from zone 'k'
            #Add nodes to zone
            for i in __ID:
                if i not in self.phymed.zones[zone_id]['nodes']:
                    self.phymed.zones[zone_id]['nodes'].append(i)
            self.update_zone(zone_id,self.phymed.zones[zone_id]['Ss'],self.phymed.zones[zone_id]['Kx'],self.phymed.zones[zone_id]['Ky'])
            #ZONE ASSIGNATION LOOP#############################################
            del zone_id,x_min,x_max,y_min,y_max,i_min,i_max,j_min,j_max,__node,__node_ID,__ID
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
            for i in self.phymed.zones[zone_id]['nodes']:
                self.phymed.k[i]=[Kx,Ky] #Update node value
                self.phymed.ss[i]=Ss #Update node value
            del zone_id,Ss,Kx,Ky
        else:
            print("Zone {0} doesn't exist".format(zone_id))
            del zone_id,Ss,Kx,Ky
    ##UPDATE ZONE PARAMETERS###################################################
    
    ##DELETE ZONE##############################################################
    def delete_zone(self,zone_id):
        if zone_id in self.phymed.zones: #CHECK IF THE ZONE EXISTS
            for i in self.phymed.zones[zone_id]['nodes']:
                self.phymed.zones[0]['nodes'].append(i) #Add node to default zone
                self.phymed.k[i]=[self.phymed.zones[0]['Kx'],self.phymed.zones[0]['Ky']] #Update node value
                self.phymed.ss[i]=self.phymed.zones[0]['Ss'] #Update node value
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
    def solve(self,tol=0.001):
        if self.phymed.steady==True: #Steady-state
            print('ASSEMBLING GLOBAL MATRICES...')
            ##ASSEMBLY GLOBAL MATRICES#########################################
            self.coeff_m_assembly() #ASSEMBLY COEFFICIENTS MATRIX
            self.load_v_assembly() #ASSEMBLY LOAD VECTOR
            ##ASSEMBLY GLOBAL MATRICES#########################################
            
            print('APPLYING BCs...')
            ##APPLY BC's#######################################################
            self.apply_bc_n()
            self.apply_bc_d()
            ##APPLY BC's#######################################################
            
            print('COMPUTING HYDRAULIC HEAD DISTRIBUTION...')
            ##GET HYDRAULIC HEAD DISTRIBUTION##################################
            __aux_hh=np.matmul(np.linalg.inv(self.m_coeff),self.v_load)
            self.phymed.hh_n=__aux_hh.tolist()
            self.rebuild_hh()
            ##GET HYDRAULIC HEAD DISTRIBUTION################################## 
        else: #Transient-state
            print('ASSEMBLING GLOBAL MATRICES...')
            ##ASSEMBLY GLOBAL MATRICES#########################################
            self.coeff_m_assembly() #ASSEMBLY COEFFICIENTS MATRIX
            self.load_v_assembly() #ASSEMBLY LOAD VECTOR
            ##ASSEMBLY GLOBAL MATRICES#########################################
            
            print('APPLYING BCs...')
            ##APPLY BC's#######################################################
            self.apply_bc_n()
            self.apply_bc_d()
            ##APPLY BC's#######################################################
            
            print('ASSEMBLING AUXILIAR MATRICES...')
            ##GET AUXILIAR MATRICES############################################
            __Cinv=np.linalg.inv(self.m_coeff)
            __n=len(self.m_coeff)
            __V=[0]*__n
            
            print('SOLVING FOR EACH TIME STEP...')
            for k in range(self.phymed.time_steps-1): #TIME-LOOP
                __step=k+1
                print(str(100*__step/self.phymed.time_steps)+'%...')
                
                ##UPDATE AUXILIAR MATRICES#####################################
                for i in range(__n):
                    __V[i]=self.v_load[k][i]-(1/self.phymed.dt)*self.phymed.hh_n[__step-1][i]
                ##UPDATE AUXILIAR MATRICES#####################################
                
                ##GET HYDRAULIC HEAD DISTRIBUTION##############################
                ###CORREGIR
                print('COMPUTING HYDRAULIC HEAD DISTRIBUTION...')
                __aux_hh=np.matmul(__Cinv,__V)
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
    ##EXPORT RESULTS TO NATIVE .HHD FORMAT#####################################
    def export_results(self,file,step=None):
        #WRITE AN OUTPUT FILE CONTAINING THE HYDRAULIC HEAD DISTRIBUTION#######
        __f=open(aux_f.path+'/output/'+file+'.hhd','w') #Open the given file in 'write' mode
        
        #WRITING THE FILE HEADERS##############################################
        ##CURRENT TIPE STEP
        __f.write('time=')
        if step==None: #Steady-state
            __f.write('NA\n')
        else:
            __f.write(str(step*self.phymed.dt)+'\n') #Writes current time based on the time-step and the dt values
        ##NUMBER OF ROWS
        __f.write('no_rows='+str(self.mesh.ny)+'\n') #Writes the number of rows of the mesh
        ##NUMBER OF COLUMNS
        __f.write('no_colums='+str(self.mesh.nx)+'\n') #Writes the number of colums of the mesh
        ##COORDINATES
        __f.write('x_boundaries=('+str(self.mesh.m_p[0][0])+','+str(self.mesh.m_p[self.mesh.nx-1][0])+')\n')
        __f.write('y_boundaries=('+str(self.mesh.m_p[0][1])+','+str(self.mesh.m_p[self.mesh.nx*(self.mesh.ny-1)][1])+')\n')
        #WRITING THE FILE HEADERS##############################################
        
        #LOOP FOR WRITING VALUES ON THE FILE###################################
        if step==None: #Steady-state
            for j in range(self.mesh.ny):
                for i in range(self.mesh.nx):
                    __node=i+j*self.mesh.nx
                    __f.write("{0:.7E} ".format(self.phymed.hh_n[__node]))
                __f.write('\n')
        else:
            for j in range(self.mesh.ny):
                for i in range(self.mesh.nx):
                    __node=i+j*self.mesh.nx
                    __f.write("{0:.7E} ".format(self.phymed.hh_n[step][__node]))
                __f.write('\n')
        #LOOP FOR WRITING VALUES ON THE FILE###################################
        
        __f.close()
    ##EXPORT RESULTS TO NATIVE .HHD FORMAT#####################################
    
    ##DRAWDOWN AT A GIVEN POINT################################################
    def drawdown_at(self,x,y,file=None):
        #COMPUTE DRAWDOWN VALUE AT X,Y#########################################
        __n=self.mesh.coord2node(x,y) #Get closest node index
        ##INTERPOLATE HEAD VALUE AT X,Y#################
        __h=[0]*self.phymed.time_steps #Auxiliar array to store head values
        ###CHECK HOW WE SHOULD INTERPOLATE THE HEAD VALUE
        #Do not interpolate if the node is at (x,y)
        if self.mesh.m_p[__n][0]==x and self.mesh.m_p[__n][1]==y: #Do not interpolate
            for k in range(self.phymed.time_steps):
                __h[k]=self.phymed.hh_n[k][__n]
        #Use a linear interpolation if the node is on 'x' or 'y'
        elif self.mesh.m_p[__n][1]==y: #Interpolate over 'x'
            if self.mesh.m_p[__n][0]<x: #(x,y) is to the right of the node
                for k in range(self.phymed.time_steps):
                    __h[k]=(x-self.mesh.m_p[__n][0])*\
                    (self.phymed.hh_n[k][__n+1]-self.phymed.hh_n[k][__n])/\
                    (self.mesh.m_p[__n+1][0]-self.mesh.m_p[__n][0])+\
                    self.phymed.hh_n[k][__n]
            else: #(x,y) is to the left of the node
                for k in range(self.phymed.time_steps):
                    __h[k]=self.phymed.hh_n[k][__n]-\
                    (self.mesh.m_p[__n][0]-x)*\
                    (self.phymed.hh_n[k][__n]-self.phymed.hh_n[k][__n-1])/\
                    (self.mesh.m_p[__n][0]-self.mesh.m_p[__n-1][0])
        elif self.mesh.m_p[__n][0]==x: #Interpolate over 'y'
            if self.mesh.m_p[__n][0]<y: #(x,y) is below the node
                for k in range(self.phymed.time_steps):
                    __h[k]=(y-self.mesh.m_p[__n][1])*\
                    (self.phymed.hh_n[k][__n+self.mesh.nx]-self.phymed.hh_n[k][__n])/\
                    (self.mesh.m_p[__n+self.mesh.nx][1]-self.mesh.m_p[__n][1])+\
                    self.phymed.hh_n[k][__n]
            else: #(x,y) is above the node
                for k in range(self.phymed.time_steps):
                    __h[k]=self.phymed.hh_n[k][__n]-\
                    (self.mesh.m_p[__n][1]-y)*\
                    (self.phymed.hh_n[k][__n]-self.phymed.hh_n[k][__n-self.mesh.nx])/\
                    (self.mesh.m_p[__n][1]-self.mesh.m_p[__n-self.mesh.nx][1])
        else:
            if self.mesh.m_p[__n][0]<x: #(x,y) is to the right of the node
                if self.mesh.m_p[__n][0]<y: #(x,y) is below the node
                    for k in range(self.phymed.time_steps):
                        hy1=(x-self.mesh.m_p[__n+self.mesh.nx][0])*\
                        (self.phymed.hh_n[k][__n+1+self.mesh.nx]-self.phymed.hh_n[k][__n+self.mesh.nx])/\
                        (self.mesh.m_p[__n+1+self.mesh.nx][0]-self.mesh.m_p[__n+self.mesh.nx][0])+\
                        self.phymed.hh_n[k][__n+self.mesh.nx]
                        hy2=(x-self.mesh.m_p[__n][0])*\
                        (self.phymed.hh_n[k][__n+1]-self.phymed.hh_n[k][__n])/\
                        (self.mesh.m_p[__n+1][0]-self.mesh.m_p[__n][0])+\
                        self.phymed.hh_n[k][__n]
                        __h[k]=(y-self.mesh.m_p[__n][1])*\
                        (hy2-hy1)/\
                        (self.mesh.m_p[__n+self.mesh.nx][1]-self.mesh.m_p[__n][1])+\
                        hy2
                else: #(x,y) is above the node
                    for k in range(self.phymed.time_steps):
                        hy1=(x-self.mesh.m_p[__n-self.mesh.nx][0])*\
                        (self.phymed.hh_n[k][__n+1-self.mesh.nx]-self.phymed.hh_n[k][__n-self.mesh.nx])/\
                        (self.mesh.m_p[__n+1-self.mesh.nx][0]-self.mesh.m_p[__n-self.mesh.nx][0])+\
                        self.phymed.hh_n[k][__n-self.mesh.nx]
                        hy2=(x-self.mesh.m_p[__n][0])*\
                        (self.phymed.hh_n[k][__n+1]-self.phymed.hh_n[k][__n])/\
                        (self.mesh.m_p[__n+1][0]-self.mesh.m_p[__n][0])+\
                        self.phymed.hh_n[k][__n]
                        __h[k]=hy2-\
                        (self.mesh.m_p[__n][1]-y)*\
                        (hy2-hy1)/\
                        (self.mesh.m_p[__n][1]-self.mesh.m_p[__n-self.mesh.nx][1])
            else: #(x,y) is to the left of the node
                if self.mesh.m_p[__n][0]<y: #(x,y) is below the node
                    for k in range(self.phymed.time_steps):
                        hy1=self.phymed.hh_n[k][__n+self.mesh.nx]-\
                        (self.mesh.m_p[__n+self.mesh.nx][0]-x)*\
                        (self.phymed.hh_n[k][__n+self.mesh.nx]-self.phymed.hh_n[k][__n-1+self.mesh.nx])/\
                        (self.mesh.m_p[__n+self.mesh.nx][0]-self.mesh.m_p[__n-1+self.mesh.nx][0])
                        hy2=self.phymed.hh_n[k][__n]-\
                        (self.mesh.m_p[__n][0]-x)*\
                        (self.phymed.hh_n[k][__n]-self.phymed.hh_n[k][__n-1])/\
                        (self.mesh.m_p[__n][0]-self.mesh.m_p[__n-1][0])
                        __h[k]=(y-self.mesh.m_p[__n][1])*\
                        (hy2-hy1)/\
                        (self.mesh.m_p[__n+self.mesh.nx][1]-self.mesh.m_p[__n][1])+\
                        hy2
                else: #(x,y) is above the node
                    for k in range(self.phymed.time_steps):
                        hy1=self.phymed.hh_n[k][__n-self.mesh.nx]-\
                        (self.mesh.m_p[__n-self.mesh.nx][0]-x)*\
                        (self.phymed.hh_n[k][__n-self.mesh.nx]-self.phymed.hh_n[k][__n-1-self.mesh.nx])/\
                        (self.mesh.m_p[__n-self.mesh.nx][0]-self.mesh.m_p[__n-1-self.mesh.nx][0])
                        hy2=self.phymed.hh_n[k][__n]-\
                        (self.mesh.m_p[__n][0]-x)*\
                        (self.phymed.hh_n[k][__n]-self.phymed.hh_n[k][__n-1])/\
                        (self.mesh.m_p[__n][0]-self.mesh.m_p[__n-1][0])
                        __h[k]=hy2-\
                        (self.mesh.m_p[__n][1]-y)*\
                        (hy2-hy1)/\
                        (self.mesh.m_p[__n][1]-self.mesh.m_p[__n-self.mesh.nx][1])
        ##INTERPOLATE HEAD VALUE AT X,Y#################
        ##COMPUTE DRAWDOWN VALUES##############################################
        drawdown=[0]*self.phymed.time_steps #Initialize array
        for k in range(self.phymed.time_steps):
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
        k=step
        for i in range(npts):
            x=obs_pts[i][0]
            y=obs_pts[i][1]
            __n=self.mesh.coord2node(x,y) #Get closest node index
            ##INTERPOLATE HEAD VALUE AT X,Y#################
            ###CHECK HOW WE SHOULD INTERPOLATE THE HEAD VALUE
            #Do not interpolate if the node is at (x,y)
            if self.mesh.m_p[__n][0]==x and self.mesh.m_p[__n][1]==y: #Do not interpolate
                __h0=self.phymed.hh_n[0][__n]
                __hk=self.phymed.hh_n[k][__n]
            #Use a linear interpolation if the node is on 'x' or 'y'
            elif self.mesh.m_p[__n][1]==y: #Interpolate over 'x'
                if self.mesh.m_p[__n][0]<x: #(x,y) is to the right of the node
                    __h0=(x-self.mesh.m_p[__n][0])*\
                    (self.phymed.hh_n[0][__n+1]-self.phymed.hh_n[0][__n])/\
                    (self.mesh.m_p[__n+1][0]-self.mesh.m_p[__n][0])+\
                    self.phymed.hh_n[0][__n]
                    __hk=(x-self.mesh.m_p[__n][0])*\
                    (self.phymed.hh_n[k][__n+1]-self.phymed.hh_n[k][__n])/\
                    (self.mesh.m_p[__n+1][0]-self.mesh.m_p[__n][0])+\
                    self.phymed.hh_n[k][__n]
                else: #(x,y) is to the left of the node
                    __h0=self.phymed.hh_n[0][__n]-\
                    (self.mesh.m_p[__n][0]-x)*\
                    (self.phymed.hh_n[0][__n]-self.phymed.hh_n[0][__n-1])/\
                    (self.mesh.m_p[__n][0]-self.mesh.m_p[__n-1][0])
                    __hk=self.phymed.hh_n[k][__n]-\
                    (self.mesh.m_p[__n][0]-x)*\
                    (self.phymed.hh_n[k][__n]-self.phymed.hh_n[k][__n-1])/\
                    (self.mesh.m_p[__n][0]-self.mesh.m_p[__n-1][0])
            elif self.mesh.m_p[__n][0]==x: #Interpolate over 'y'
                if self.mesh.m_p[__n][0]<y: #(x,y) is below the node
                    __h0=(y-self.mesh.m_p[__n][1])*\
                    (self.phymed.hh_n[0][__n+self.mesh.nx]-self.phymed.hh_n[0][__n])/\
                    (self.mesh.m_p[__n+self.mesh.nx][1]-self.mesh.m_p[__n][1])+\
                    self.phymed.hh_n[0][__n]
                    __hk=(y-self.mesh.m_p[__n][1])*\
                    (self.phymed.hh_n[k][__n+self.mesh.nx]-self.phymed.hh_n[k][__n])/\
                    (self.mesh.m_p[__n+self.mesh.nx][1]-self.mesh.m_p[__n][1])+\
                    self.phymed.hh_n[k][__n]
                else: #(x,y) is above the node
                    __h0=self.phymed.hh_n[0][__n]-\
                    (self.mesh.m_p[__n][1]-y)*\
                    (self.phymed.hh_n[0][__n]-self.phymed.hh_n[0][__n-self.mesh.nx])/\
                    (self.mesh.m_p[__n][1]-self.mesh.m_p[__n-self.mesh.nx][1])
                    __hk=self.phymed.hh_n[k][__n]-\
                    (self.mesh.m_p[__n][1]-y)*\
                    (self.phymed.hh_n[k][__n]-self.phymed.hh_n[k][__n-self.mesh.nx])/\
                    (self.mesh.m_p[__n][1]-self.mesh.m_p[__n-self.mesh.nx][1])
            else:
                if self.mesh.m_p[__n][0]<x: #(x,y) is to the right of the node
                    if self.mesh.m_p[__n][0]<y: #(x,y) is below the node
                        #Time-step=0
                        hy1_0=(x-self.mesh.m_p[__n+self.mesh.nx][0])*\
                        (self.phymed.hh_n[0][__n+1+self.mesh.nx]-self.phymed.hh_n[0][__n+self.mesh.nx])/\
                        (self.mesh.m_p[__n+1+self.mesh.nx][0]-self.mesh.m_p[__n+self.mesh.nx][0])+\
                        self.phymed.hh_n[0][__n+self.mesh.nx]
                        hy2_0=(x-self.mesh.m_p[__n][0])*\
                        (self.phymed.hh_n[0][__n+1]-self.phymed.hh_n[0][__n])/\
                        (self.mesh.m_p[__n+1][0]-self.mesh.m_p[__n][0])+\
                        self.phymed.hh_n[0][__n]
                        __h0=(y-self.mesh.m_p[__n][1])*\
                        (hy2_0-hy1_0)/\
                        (self.mesh.m_p[__n+self.mesh.nx][1]-self.mesh.m_p[__n][1])+\
                        hy2_0
                        #Time-step=k
                        hy1_k=(x-self.mesh.m_p[__n+self.mesh.nx][0])*\
                        (self.phymed.hh_n[k][__n+1+self.mesh.nx]-self.phymed.hh_n[k][__n+self.mesh.nx])/\
                        (self.mesh.m_p[__n+1+self.mesh.nx][0]-self.mesh.m_p[__n+self.mesh.nx][0])+\
                        self.phymed.hh_n[k][__n+self.mesh.nx]
                        hy2_k=(x-self.mesh.m_p[__n][0])*\
                        (self.phymed.hh_n[k][__n+1]-self.phymed.hh_n[k][__n])/\
                        (self.mesh.m_p[__n+1][0]-self.mesh.m_p[__n][0])+\
                        self.phymed.hh_n[k][__n]
                        __hk=(y-self.mesh.m_p[__n][1])*\
                        (hy2_k-hy1_k)/\
                        (self.mesh.m_p[__n+self.mesh.nx][1]-self.mesh.m_p[__n][1])+\
                        hy2_k
                    else: #(x,y) is above the node
                        #Time-step=0
                        hy1_0=(x-self.mesh.m_p[__n-self.mesh.nx][0])*\
                        (self.phymed.hh_n[0][__n+1-self.mesh.nx]-self.phymed.hh_n[0][__n-self.mesh.nx])/\
                        (self.mesh.m_p[__n+1-self.mesh.nx][0]-self.mesh.m_p[__n-self.mesh.nx][0])+\
                        self.phymed.hh_n[0][__n-self.mesh.nx]
                        hy2_0=(x-self.mesh.m_p[__n][0])*\
                        (self.phymed.hh_n[0][__n+1]-self.phymed.hh_n[0][__n])/\
                        (self.mesh.m_p[__n+1][0]-self.mesh.m_p[__n][0])+\
                        self.phymed.hh_n[0][__n]
                        __h0=hy2_0-\
                        (self.mesh.m_p[__n][1]-y)*\
                        (hy2_0-hy1_0)/\
                        (self.mesh.m_p[__n][1]-self.mesh.m_p[__n-self.mesh.nx][1])
                        #Time-step=k
                        hy1_k=(x-self.mesh.m_p[__n-self.mesh.nx][0])*\
                        (self.phymed.hh_n[k][__n+1-self.mesh.nx]-self.phymed.hh_n[k][__n-self.mesh.nx])/\
                        (self.mesh.m_p[__n+1-self.mesh.nx][0]-self.mesh.m_p[__n-self.mesh.nx][0])+\
                        self.phymed.hh_n[k][__n-self.mesh.nx]
                        hy2_k=(x-self.mesh.m_p[__n][0])*\
                        (self.phymed.hh_n[k][__n+1]-self.phymed.hh_n[k][__n])/\
                        (self.mesh.m_p[__n+1][0]-self.mesh.m_p[__n][0])+\
                        self.phymed.hh_n[k][__n]
                        __hk=hy2_k-\
                        (self.mesh.m_p[__n][1]-y)*\
                        (hy2_k-hy1_k)/\
                        (self.mesh.m_p[__n][1]-self.mesh.m_p[__n-self.mesh.nx][1])
                else: #(x,y) is to the left of the node
                    if self.mesh.m_p[__n][0]<y: #(x,y) is below the node
                        #Time-step=0
                        hy1_0=self.phymed.hh_n[0][__n+self.mesh.nx]-\
                        (self.mesh.m_p[__n+self.mesh.nx][0]-x)*\
                        (self.phymed.hh_n[0][__n+self.mesh.nx]-self.phymed.hh_n[0][__n-1+self.mesh.nx])/\
                        (self.mesh.m_p[__n+self.mesh.nx][0]-self.mesh.m_p[__n-1+self.mesh.nx][0])
                        hy2_0=self.phymed.hh_n[0][__n]-\
                        (self.mesh.m_p[__n][0]-x)*\
                        (self.phymed.hh_n[0][__n]-self.phymed.hh_n[0][__n-1])/\
                        (self.mesh.m_p[__n][0]-self.mesh.m_p[__n-1][0])
                        __h0=(y-self.mesh.m_p[__n][1])*\
                        (hy2_0-hy1_0)/\
                        (self.mesh.m_p[__n+self.mesh.nx][1]-self.mesh.m_p[__n][1])+\
                        hy2_0
                        #Time-step=k
                        hy1_k=self.phymed.hh_n[k][__n+self.mesh.nx]-\
                        (self.mesh.m_p[__n+self.mesh.nx][0]-x)*\
                        (self.phymed.hh_n[k][__n+self.mesh.nx]-self.phymed.hh_n[k][__n-1+self.mesh.nx])/\
                        (self.mesh.m_p[__n+self.mesh.nx][0]-self.mesh.m_p[__n-1+self.mesh.nx][0])
                        hy2_k=self.phymed.hh_n[k][__n]-\
                        (self.mesh.m_p[__n][0]-x)*\
                        (self.phymed.hh_n[k][__n]-self.phymed.hh_n[k][__n-1])/\
                        (self.mesh.m_p[__n][0]-self.mesh.m_p[__n-1][0])
                        __hk=(y-self.mesh.m_p[__n][1])*\
                        (hy2_k-hy1_k)/\
                        (self.mesh.m_p[__n+self.mesh.nx][1]-self.mesh.m_p[__n][1])+\
                        hy2_k
                    else: #(x,y) is above the node
                        #Time-step=0
                        hy1_0=self.phymed.hh_n[0][__n-self.mesh.nx]-\
                        (self.mesh.m_p[__n-self.mesh.nx][0]-x)*\
                        (self.phymed.hh_n[0][__n-self.mesh.nx]-self.phymed.hh_n[0][__n-1-self.mesh.nx])/\
                        (self.mesh.m_p[__n-self.mesh.nx][0]-self.mesh.m_p[__n-1-self.mesh.nx][0])
                        hy2_0=self.phymed.hh_n[0][__n]-\
                        (self.mesh.m_p[__n][0]-x)*\
                        (self.phymed.hh_n[0][__n]-self.phymed.hh_n[0][__n-1])/\
                        (self.mesh.m_p[__n][0]-self.mesh.m_p[__n-1][0])
                        __h0=hy2_0-\
                        (self.mesh.m_p[__n][1]-y)*\
                        (hy2_0-hy1_0)/\
                        (self.mesh.m_p[__n][1]-self.mesh.m_p[__n-self.mesh.nx][1])
                        #Time-step=k
                        hy1_k=self.phymed.hh_n[k][__n-self.mesh.nx]-\
                        (self.mesh.m_p[__n-self.mesh.nx][0]-x)*\
                        (self.phymed.hh_n[k][__n-self.mesh.nx]-self.phymed.hh_n[k][__n-1-self.mesh.nx])/\
                        (self.mesh.m_p[__n-self.mesh.nx][0]-self.mesh.m_p[__n-1-self.mesh.nx][0])
                        hy2_k=self.phymed.hh_n[k][__n]-\
                        (self.mesh.m_p[__n][0]-x)*\
                        (self.phymed.hh_n[k][__n]-self.phymed.hh_n[k][__n-1])/\
                        (self.mesh.m_p[__n][0]-self.mesh.m_p[__n-1][0])
                        __hk=hy2_k-\
                        (self.mesh.m_p[__n][1]-y)*\
                        (hy2_k-hy1_k)/\
                        (self.mesh.m_p[__n][1]-self.mesh.m_p[__n-self.mesh.nx][1])
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