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
    
###########################AUXILIAR FUNCTIONS MODULE###########################

This module includes all object classes and functions required for 
Groundwater (GW) modelling that were not included in the FEM module
nor the FDM one.

"""
###############################################################################
"""
##########################LIBRARY IMPORTATION SECTION##########################
"""
import ast #Import 'ast' library for file reading.
import numpy as np #Import 'numpy' for numerical operations
import matplotlib.pyplot as plot #Import 'matplotlib' library for showing contour maps
import matplotlib.tri as tri #Import 'matplotlib' library for showing contour maps
import os #Import for getting current path

"""
##########################AUXILIAR "GLOBAL" VARIABLES##########################
"""
path=os.getcwd()
#Dictionary containing translations for axis names, etc.
lang_dict={'es':['Distancia','Tiempo','Abatimiento'],
           'en':['Distance','Time','Drawdown'],
           'fr':['Distance','Temps','Dépression']}
#Dictionary for unit conversion
unit_dict={'mm':1000,'cm':100,'dm':10,'m':1,'dam':0.1,'hm':0.01,'km':0.001,\
           'in':39.3701,'ft':3.28084,'yd':1.09361,'mi':6.21371E-04}
"""
##########################AUXILIAR "GLOBAL" VARIABLES##########################
"""

"""
##########################LIBRARY IMPORTATION SECTION##########################
"""
###############################################################################
"""
############################VIEW AND EXPORT RESULTS############################
"""

#PLOT HYDRAULIC HEAD DISTRIBUTION##############################################
def plot_hh(model,npts=21,superpose=False,step=None,lang='en'):
    """
    model = 'fem' or 'fdm' object
    npts = # of interpolation points
    """
    #CONTOUR MAP###############################################################
    if superpose==False:
        plot.figure()
        cflow='Red'
    else:
        cflow='Blue'
    ##MAKE DATA GRID###########################################################
    ###RETRIEVE DATA FROM MODEL
    if model.type=='FEM':
        x_list=np.asarray(list(list(zip(*model.mesh.m_p))[0]))
        y_list=np.asarray(list(list(zip(*model.mesh.m_p))[1]))
        if step==None:
            z_list=np.asarray(model.phymed.hh_n)
        else:
            z_list=np.asarray(model.phymed.hh_n[step])
    else:
        x_list=np.asarray(list(list(zip(*model.mesh.m_p))[0]))
        y_list=np.asarray(list(list(zip(*model.mesh.m_p))[1]))
        if step==None:
            z_list=np.asarray(model.phymed.hh_n)
        else:
            z_list=np.asarray(model.phymed.hh_n[step])
    ###INTERPOLATE DATA
    x_axis=np.linspace(min(x_list),max(x_list),npts)
    y_axis=np.linspace(min(y_list),max(y_list),npts)
    triang_sp=tri.Triangulation(x_list,y_list)
    interpolate=tri.LinearTriInterpolator(triang_sp,z_list)
    ###ACTURAL GRIDDING
    xi,yi=np.meshgrid(x_axis,y_axis)
    hi=interpolate(xi,yi)
    ##MAKE DATA GRID###########################################################
    if superpose==False:
        contour_map=plot.contourf(xi,yi,hi) #Plot contourmap
        #STYLE COMMANDS########################################################
        plot.gca().invert_yaxis() #Flip 'y' axis
        plot.title(model.config['model_name'],fontweight='bold')
        plot.xlabel(lang_dict[lang][0]+' ('+model.config['length']+')')
        plot.ylabel(lang_dict[lang][0]+' ('+model.config['length']+')')
        #STYLE COMMANDS########################################################
        color_bar=plot.colorbar(contour_map) #Plot colorbar
#        color_bar.set_label('h ('+model.config['head']+')')
    #CONTOUR MAP###############################################################
    #QUIVER MAP################################################################
    ##CALCULATE FLOW VECTORS
    __dx=(max(x_list)-min(x_list))/(npts-1)
    __dy=(max(y_list)-min(y_list))/(npts-1)
    ###INITIALIZE AUXILIAR VECTORS
    x_flow=[0]*(npts)
    y_flow=[0]*(npts)
#    flow=[0]*(npts)
    for i in range(npts):
        x_flow[i]=[0]*(npts)
        y_flow[i]=[0]*(npts)
#        flow[i]=[0]*(npts)
    ###COMPUTE FLOW VALUES
    for i in range(1,npts-1):
        for j in range(1,npts-1):
            x_flow[i][j]=-(hi[i][j+1]-hi[i][j-1])/(2*__dx)
            y_flow[i][j]=(hi[i+1][j]-hi[i-1][j])/(2*__dy)
#            flow[i][j]=(x_flow[i][j]**2+y_flow[i][j]**2)**0.5
    ###LIST TO ARRAY
    ui=np.asarray(x_flow)
    vi=np.asarray(y_flow)
    ###LIST TO ARRAY
    plot.quiver(xi[1:npts-1],yi[1:npts-1],ui[1:npts-1],vi[1:npts-1],color=cflow)
    #QUIVER MAP################################################################
    plot.show()
#PLOT HYDRAULIC HEAD DISTRIBUTION##############################################

#PLOT DRAWDOWN CURVES##########################################################
def plot_drawdown(files:list,lang='en'):
    ##INITIALIZE FIGURE########################################################
    fig=plot.figure()
    plot.title(lang_dict[lang][2].upper())#,fontweight='bold')
    ##INITIALIZE FIGURE########################################################
    
    ##RETRIEVE DATA FROM EACH FILE#############################################
    for i in files:
        __curve_name=i[:-4]
        __f=open(path+'/output/'+i,'r')
        
        #READ HEADERS##########################################################
        __f.readline() #Skip a line
        __type=int(__f.readline()) #Read file type
        __f.readline() #Skip a line
        __units=__f.readline().strip('\n').split(',')
        __f.readline() #Skip a line
        __n=int(__f.readline()) #Read number of steps
        __f.readline() #Skip a line
        __coord=__f.readline().split(' ') #Split line
        #READ HEADERS##########################################################
        ##READ STEP VALUES################################################
        __f.readline() #Skip header line
        __k=__f.readline().split(' ')
        for j in range(__n):
            __k[j]=float(__k[j])
        del __k[__n]
        ##READ TIME-STEP VALUES################################################
        ##READ DRAWDOWN VALUES#################################################
        __f.readline() #Skip header line
        __dwn=__f.readline().split(' ')
        for j in range(__n):
            __dwn[j]=float(__dwn[j])
        del __dwn[__n]
        ##READ DRAWDOWN VALUES#################################################
        __f.close()
        plot.ylabel(lang_dict[lang][2]+' ('+__units[1]+')')
        if __type==0:
            plot.xlabel(lang_dict[lang][1]+' ('+__units[0]+')')
        elif __type==1:
            plot.xlabel(lang_dict[lang][0]+' ('+__units[0]+')')
        plot.plot(__k,__dwn,label=__curve_name) #'o' as third argument to display as points
    ##RETRIEVE DATA FROM EACH FILE#############################################
    
    ##ADITIONAL COMMANDS FOR PLOT##############################################
    plot.legend()
    plot.show()
    ##ADITIONAL COMMANDS FOR PLOT##############################################
    
#PLOT DRAWDOWN CURVES##########################################################
"""
############################VIEW AND EXPORT RESULTS############################
"""
###############################################################################
"""
################################MISC FUNCTIONS#################################
"""
#DICTIONARY-LIKE FILE READING##################################################
def read_dict(file):    
    #READ FILE AS IF IT WERE A PYTHON DICTIONARY###############################
    __f=open(path+'/'+file,'r') #Open the given file in 'read' mode
    ##CREATE A DICTIONARY TO STORE THE FILE INFO###############################
    content=dict()
    ###########################################################################
    for __line in __f:
        if __line[0]!='#' and __line.rstrip('\n')!='':
            __aline=__line.rstrip('\n').split('=') #Remove 'new line' character and split line @ '='
            content[__aline[0]]=ast.literal_eval(__aline[1])
    __f.close()
    return content
    #READ FILE AS IF IT WERE A PYTHON DICTIONARY###############################
#DICTIONARY-LIKE FILE READING##################################################
"""
################################MISC FUNCTIONS#################################
"""