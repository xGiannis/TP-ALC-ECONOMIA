# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
from funciones import inversaLU, calcularLU

"""
Editor de Spyder

Este es un archivo temporal.
"""

"""Descargar de este link, en la secci´on Adjuntos, la MIP Latinoamericana 2011 (Prelimi-
nar). Este es un archivo de Excel, el cual deben leer con la librer´ıa pandas de Python.
En la Hoja ’LAT IOT 2011’ se encontrar´an los flujos entre 40 sectores de 18 paises, ex-
presados en millones de d´olares. La ´ultima columna del archivo, con el t´ıtulo ’Output’
es el total producido para ese sector."""

p1="SLV"
p2="PAN"


nombreFile = "matrizlatina2011_compressed_0.xlsx"

archivo = pd.read_excel(nombreFile, sheet_name=1)
#%%%%%%%%%%%%%%%%%%%%%%%%
panama = archivo[archivo["Country_iso3"]==p2]
iPP= panama.filter(regex='^PAN', axis=1)


outputP=panama["Output"]

output01P= outputP.replace(0,1)






#%%%%%%%%%%%%%%%%%%%%%%%%

salvador = archivo[archivo["Country_iso3"]==p1]


iSS= salvador.filter(regex='^SLV', axis=1)

outputS = salvador["Output"]

output01S=outputS.replace(0,1)


#%%

iPS= panama.filter(regex='^SLV', axis=1)


#%%

iSP= salvador.filter(regex='^PAN', axis=1)
#%%
#ASP - APP
#APS - ASS
def coeficientesTecnicos(Z,P):
    P= np.diag(P)
    L, U, Per = calcularLU(P)
    P_inv = inversaLU(L, U, Per)
    A =Z@P_inv
    
    return  A

APP = coeficientesTecnicos(iPP, output01P)

ASP=coeficientesTecnicos(iSP, output01P)

ASS= coeficientesTecnicos(iSS, output01S) 

APS= coeficientesTecnicos(iPS, output01S)
#%%
#rastreamos coeficientes mayores a 1

APP1= APP > 1
indices = APP1.stack().index[APP1.stack()].tolist()

APP1= APS > 1
indices = APP1.stack().index[APP1.stack()].tolist()

APP1= ASS > 1
indices = APP1.stack().index[APP1.stack()].tolist()

APP1= ASP > 1
indices = APP1.stack().index[APP1.stack()].tolist()
#%%

def demandaCalculator(RR,RS,SR,SS,pR,pS):
    
    tamañoIdentidad=np.shape(RR)[0] + np.shape(SS)[1]
    
    identidadSuper= np.eye(tamañoIdentidad)
    
    up = np.hstack((RR,RS))
    
    down = np.hstack((SR,SS))
    
    ASuper= np.vstack((up,down))
    
    restaSuper = identidadSuper - ASuper
    
    productoSuper = np.hstack((pR,pS))
    
    demandaSuper= restaSuper @ productoSuper
    
    #######SEPARAMOS LAS DEMANDAS R,S ##########
    
    demandatamaño=np.shape(demandaSuper)[0]
    
    dR= demandaSuper[0:int(demandatamaño/2)]
    
    dS= demandaSuper[int(demandatamaño/2):demandatamaño]
    
    
    #return restaSuper, productoSuper
    return dR,dS


dP,dS=demandaCalculator(APP, APS, ASP, ASS, outputP, outputS)

#resta,producto=demandaCalculator(APP, APS, ASP, ASS, output01P, output01S)

#%%

#SHOCK SIMULATION!

def deltaD_generator(demanda, shocks):
    demandaPrima = demanda.copy()
    for i in range(len(shocks)):
        demandaPrima[shocks[i][0]-1] = demanda[shocks[i][0]-1] + demanda[shocks[i][0]-1] * shocks[i][1]
        
    deltaD = demanda - demandaPrima
    return deltaD

shocks = [[5,-0.1],[6,0.33],[7,0.33],[8,0.33]]

deltaD = deltaD_generator(dP, shocks)


#%%


