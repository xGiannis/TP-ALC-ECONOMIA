#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Eliminacion Gausianna
"""
import numpy as np


def elim_gaussianaConPermutaciones(A):
    cant_op = 0
    m=A.shape[0]
    n=A.shape[1]
    Ac = A.copy()
    
    if m!=n:
        print('Matriz no cuadrada')
        return
    
    ## desde aqui -- CODIGO A COMPLETAR
    L = np.eye(n) #inciamos L de nxn

    matricesPermutacion=[]

    print("esto es n--->",n)

    #Gaussian Elimination
    for i in range(n):
        if Ac[i,i]==0:
            # ...
            # aca supongo q tengo q hacer la permutacion 
            print(f"Pivote {i} es nulo")
            
            listaDePivotes:list[tuple] = []
            for j in range(i,n):
                listaDePivotes.append((Ac[j,i],j))   #listaDePivotes.append((elemento,indice)), tupla = (elemento,indice)

            print(listaDePivotes)

            newPivoteIndice= (0,i)
            for k in range(len(listaDePivotes)-1):
                if listaDePivotes[k][0]<= listaDePivotes[k+1][0]:
                    newPivoteIndice = listaDePivotes[k+1]

            #ahora resta crear una identidad con las filas intercambiadas por los indices newPivote[]


            identidad = np.eye(n)
            permutacion = np.eye(n)     #la transformo mas abajo, arranco con la Identidad
            permutacionFila1 = identidad[i]
            permutacionFila2 = identidad[newPivoteIndice[1]]
            permutacion[i] = permutacionFila2
            permutacion[newPivoteIndice[1]]=permutacionFila1 
            matricesPermutacion.append(permutacion) #estaran en la fila en un orden inverso al propuesto para conseguir P



            Ac= permutacion @ Ac




        
        for j in range(i+1, n):
            factor = Ac[j,i] / Ac[i,i]#aca va q multiplico entre las matrices para q de 0.
            L[j,i] = factor # guardamos el factor en la matriz L 
            Ac[j,i:] = Ac[j,i:]  - factor*Ac[i,i:]
            #cant_operaciones #... n por algo
            print(f"Matriz L despuse del paso ({j},{i})")
            print(L) 

    P=multiplicarPermutaciones(matricesPermutacion)
    ## hasta aqui
            
    #L = np.tril(Ac,-1) + np.eye(A.shape[0]) 
    print("esto queda despues de todos los pasos --->\n",Ac)
    U = np.triu(Ac)
    


    return L, U, cant_op, P

def multiplicarPermutaciones(matricesPermutaciones:list):
    P=matricesPermutaciones[len(matricesPermutaciones)-1]
    for i in range(len(matricesPermutaciones)-2,-1,-1):         #el menos 1 del segundo item es porque es -1 exclusive
        Pnmenos1 = matricesPermutaciones[i]
        P = P @ Pnmenos1

    return P


matrizEjemplo = np.array([[1,2,3,4],[1,2,5,6],[1,3,5,2],[1,3,5,4]])


def ejMatriz(matriz):
        


    print('Matriz B \n', matriz)
    
    L,U,cant_oper,P = elim_gaussianaConPermutaciones(matriz)
    
    print('Matriz L \n', L)
    print('Matriz U \n', U)
    print('Cantidad de operaciones: ', cant_oper)
    print('B=LU? ' , 'Si!' if np.allclose(np.linalg.norm(matriz - L@U, 1), 0) else 'No!')
    print('PB=LU? ' , 'Si!' if np.allclose(np.linalg.norm(P@matriz - L@U, 1), 0) else 'No!')
    print('Norma infinito de U: ', np.max(np.sum(np.abs(U), axis=1)) )

    

ejMatriz(matrizEjemplo)

