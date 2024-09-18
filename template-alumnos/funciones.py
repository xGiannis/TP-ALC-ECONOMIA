import numpy as np

def calcularLU(A):
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
            print(f"Pivote {i} es nulo")
            
            listaDePivotes:list[tuple] = []
            for j in range(i+1,n):
                listaDePivotes.append((Ac[j,i],j))   #listaDePivotes.append((elemento,indice)), tupla = (elemento,indice)

            print(listaDePivotes)

            newPivoteIndice= listaDePivotes[0]
            
            for k in range(len(listaDePivotes)-1):
                    
                if newPivoteIndice[0]<= listaDePivotes[k][0]:
                    newPivoteIndice = listaDePivotes[k]
                print(newPivoteIndice)
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
                print(f'L despues de permutar \n {L}')    

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

def escalonar_filas(M):
    """ 
        Retorna la Matriz Escalonada por Filas 
    """
    A = np.copy(M)
    if (issubclass(A.dtype.type, np.integer)):
        A = A.astype(float)

    # Si A no tiene filas o columnas, ya esta escalonada
    f, c = A.shape
    if f == 0 or c == 0:
        return A

    # buscamos primer elemento no nulo de la primera columna
    i = 0
    
    while i < f and A[i,0] == 0:
        i += 1

    if i == f:
        # si todos los elementos de la primera columna son ceros
        # escalonamos filas desde la segunda columna
        B = escalonar_filas(A[:,1:])
        
        # y volvemos a agregar la primera columna de zeros
        return np.block([A[:,:1], B])


    # intercambiamos filas i <-> 0, pues el primer cero aparece en la fila i
    if i > 0:
        A[[0,i],:] = A[[i,0],:]
        
    
    # intercambiamos filas i <-> 0, pues el primer cero aparece en la fila i
    if i > 0:
        A[[0,i],:] = A[[i,0],:]

    # PASO DE TRIANGULACION GAUSSIANA:
    # a las filas subsiguientes les restamos un multiplo de la primera
    A[1:,:] -= (A[0,:] / A[0,0]) * A[1:,0:1]

    # escalonamos desde la segunda fila y segunda columna en adelante
    B = escalonar_filas(A[1:,1:])

    # reconstruimos la matriz por bloques adosando a B la primera fila 
    # y la primera columna (de ceros)
    return np.block([ [A[:1,:]], [ A[1:,:1], B] ])


def inversaLU(L,U, P = None):
    #PA = LU -> A = P^-1 LU -> A^-1 = (LU)^-1 P -> A^-1 = U^-1 L^-1 P 
    
    Id = np.eye(np.shape(U)[0])
    
    #inversa U
    U_aumentada = np.c_[U,Id]
    U_inv = back_substitution(U_aumentada)
    
    #inversa L
    L_aumentada = np.c_[L,Id]
    L_inv = escalonar_filas(L_aumentada)
    L_inv = back_substitution(L_inv)
    
    Inv = U_inv @ L_inv @ P
    print(np.shape(Inv))
    
    return Inv

def resolverSistema(A,d):
    ML = np.eye(np.shape(A)[0]) - A

    #calculamos L, U, P y la inversa de ML:
    L, U, P = calcularLU(ML)
    ML_inv = inversaLU(L,U,P)

    #calculamos p:
    p = ML_inv @ d
    
    return p