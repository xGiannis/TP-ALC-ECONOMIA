import numpy as np

def calcularLU(A):
    """Realiza la descomposicón PA = LU de una matriz A. Como parametro rescibe una matriz A y devuelve L (triangular inferior
       con 1s en la diagonal), U (triangular superior) y P (matriz de permutación)."""
    m=A.shape[0]
    n=A.shape[1]
    Ac = A.copy()
    
    if m!=n:
        print('Matriz no cuadrada')
        return
    
    L = np.eye(n) #inciamos L de nxn

    matricesPermutacion=[]

    

    #Gaussian Elimination
    for i in range(n):
        if Ac[i,i]==0:
            
            
            listaDePivotes:list[tuple] = []
            for j in range(i+1,n):
                listaDePivotes.append((Ac[j,i],j))   #listaDePivotes.append((elemento,indice)), tupla = (elemento,indice)

            

            newPivoteIndice= listaDePivotes[0]
            
            for k in range(len(listaDePivotes)-1):
                    
                if newPivoteIndice[0]<= listaDePivotes[k][0]:
                    newPivoteIndice = listaDePivotes[k]
                
            #ahora resta crear una identidad con las filas intercambiadas por los indices newPivote[]


            identidad = np.eye(n)
            permutacion = np.eye(n)     #la transformo mas abajo, arranco con la Identidad
            permutacionFila1 = identidad[i]
            permutacionFila2 = identidad[newPivoteIndice[1]]
            permutacion[i] = permutacionFila2
            permutacion[newPivoteIndice[1]]=permutacionFila1 
            matricesPermutacion.append(permutacion) #estaran en la fila en un orden inverso al propuesto para conseguir P



            Ac= permutacion @ Ac
            
            for k in range (i):
                L[0:,k] = (permutacion @ L[0:,k])
                    

        for j in range(i+1, n):
            factor = Ac[j,i] / Ac[i,i]#aca va q multiplico entre las matrices para q de 0.
            L[j,i] = factor # guardamos el factor en la matriz L 
            Ac[j,i:] = Ac[j,i:]  - factor*Ac[i,i:]
            


    if(len(matricesPermutacion)>=1):
        P=multiplicarPermutaciones(matricesPermutacion)
    else:
        P = np.eye(n)
            
    U = np.triu(Ac)
    
    return L, U, P


def multiplicarPermutaciones(matricesPermutaciones:list):
    """Como parametro recibe una lista de matrices de permutacion y devuelve el producto de ellas"""
    P=matricesPermutaciones[len(matricesPermutaciones)-1]
    for i in range(len(matricesPermutaciones)-2,-1,-1):         #el menos 1 del segundo item es porque es -1 exclusive
        Pnmenos1 = matricesPermutaciones[i]
        P = P @ Pnmenos1

    return P

def back_substitution(A_aug):
    """
    Realiza la sustitución hacia atrás para obtener la identidad en el lado izquierdo
    de la matriz aumentada y la inversa en el lado derecho.
    """
    n = A_aug.shape[0]
    
    # Desde la última fila hacia la primera
    for i in range(n-1, -1, -1):
        # Normalizar la fila de pivote
        A_aug[i] = A_aug[i] / A_aug[i, i]
        
        # Hacer ceros en las filas superiores
        for j in range(i):
            A_aug[j] -= A_aug[i] * A_aug[j, i]
    
    return A_aug[:, n:]



def inversaLU(L,U, P = None):
    """Calcula la inversa de una matriz, recibiendo como parámetros L, U y P provenientes de la descomposición PA = LU de la
       matriz A. """
    #PA = LU -> A = P^-1 LU -> A^-1 = (LU)^-1 P -> A^-1 = U^-1 L^-1 P 
    
    Id = np.eye(np.shape(U)[0])
    
    #inversa U
    U_aumentada = np.c_[U,Id]
    U_inv = back_substitution(U_aumentada)
    
    #inversa L
    L_aumentada = np.c_[L,Id]
    L_inv = escalonar_filas(L_aumentada)

    augmento = int((L_inv.shape[1])/2) 

    L_inv = L_inv[:, augmento:] #Tomo solo la parte aumentada de la matriz que estoy operando
    
    Inv = U_inv @ L_inv @ P
   
    
    return Inv



#como es triangular inferior no se hacen 0s en la diagonal entonces no hay que permutar filas.

def escalonar_filas(M):
    """Escalona las filas de una matriz M, sirve para los casos donde no se generan 0 en la diagonal."""
    n= np.shape(M)[0]
    
    for i in range(n):
        for j in range(i+1, n):
            factor = M[j,i] / M[i,i]#aca va q multiplico entre las matrices para q de 0.
            M[j,:] = M[j,:]  - factor*M[i,:]
    
    return M