#Programa para conocer la solucion de vectores por medio de Gauss Jordan
#Algebra Lineal
#2BM1
#Perez Lucio Kevyn Alejandro
import numpy as np
#  --------------------------------------------------------------------------------------------------------------------
# FUNCIONES

def print_matrix(info, matrix):
    print (info)
    for i in range( 0, matrix.shape[0]):
        print('[', end='')
        for j in range( 0, matrix.shape[1]):
            if(j == matrix.shape[1] - 1):# Vector b
                print( '|', end=''),# Separe la matriz de coeficientes y el vector b con "|"
            print("%5.2f" %matrix[i][j], end=' ')
            if j == matrix.shape[1] - 1:# Elemento de matriz de salida m [i] [j]
                print(']', end=' ')
                print('\n')                                

def InvertirFilasColumnas(matrixa, filas):
    n = filas
    rn = filas
    # Para cada fila en AB
    for i in range(0,n-1,1):
        # columna desde diagonal i en adelante
        columna = abs(matrixa[i:,i])
        dondemax = np.argmax(columna)
        # dondemax no está en diagonal
        if (dondemax !=0):
            # intercambia filas
            temporal = np.copy(matrixa[i,:])
            matrixa[i,:] = matrixa[dondemax+i,:]
            matrixa[dondemax+i,:] = temporal 
    for i in range(n):
        suma = 0
        suma = sum(matrixa[i])
        if (suma == 0):
            matrixa = np.delete(matrixa, [i], axis=0)
            rn = n-1
    return matrixa
    return rn
    

def DiagonalAbajo(matrix, filas, columnas):
    for i in range(0,filas-1,1):
        pivote = matrix[i][i]
        adelante = i + 1
        for k in range(adelante,filas,1):
            factor = matrix[k][i]/pivote
            matrix[k][:] = matrix[k][:] - matrix[i][:]*factor
    matrixa = np.copy(matrix)
    return matrixa

def DiagonalArriba(matrix, filas, columnas):
    ultfila = filas-1
    ultcolumna = columnas-1
    for i in range(ultfila,0-1,-1):
        pivote = matrix[i][i]
        atras = i-1 
        for k in range(atras,0-1,-1):
            factor = matrix[k][i]/pivote
            matrix[k][:] = matrix[k][:] - matrix[i][:]*factor
        # diagonal a unos
        matrix[i][:] = matrix[i][:]/matrix[i][i]
    X = np.copy(matrix[:,ultcolumna])
    X = np.transpose([X])
    return X

def EliminarFilas(matrixa, filas):
    n = filas
    for i in range(n):
        suma = 0
        suma = sum(matrixa[i])
        if (suma == 0):
            matrixa = np.delete(matrixa, [i], axis=0)
    return matrixa

def FilasCeros(matrixa, filas, columnas):
    DiagonalAbajo(matrixa, filas, columnas)
    n = filas
    rn = filas
    for i in range(n):
        suma = 0
        suma = sum(matrixa[i])
        if (suma == 0):
            matrixa = np.delete(matrixa, [i], axis=0)
            rn = n-1
    return rn

def SolucionesInfinitas(matrix, literales, num_param, param):
    comienzo = len(matrix[0]) - num_param
    p_l = 0
    bucles = 0
    for indice, ecuacion in enumerate(matrix):
        print(f"{literales[indice]} = {matrix[indice][-1]:.4f}", end="")
        for p in range(comienzo, len(matrix[0])):
            if matrix[indice][p - 1] < 0:
                print("+", end="")
            print(f"{-matrix[indice][p - 1]:.4f}{param[p_l]}", end="")
            p_l = p_l + 1
        print("\n")
        p_l = 0
    p_l = 0
    for g in range(0, bucles):
        literales.pop(0)

def SolucionSolInfinitasCL(matrix, vcl, vn, rn):
    filas = len(matrix)
    print("Una posible combinación lineal suponiendo que los parametros libres valen 1 es:")
    print("")
    print("< ", end="")
    for i in range(rn):
        print("%0.2f " %(vcl[i]), end="")
    print("> = ")
    for i in range(filas):
        suma = 0
        suma = sum(matrix[i]) - 1
        print("+", end="")
        print("(", end="")
        print(suma, end="")
        print(")", end="")
        print("V%d" %(i+1), end="")
    if filas<vn:
        for i in range(vn-filas):
            indice = filas + (i+1)
            print("+(1)V" +str(indice), end="")
        
#  --------------------------------------------------------------------------------------------------------------------
# PROGRAMA PRINCIPAL

print("____PROGRAMA PARA SOLUCIONAR PROBLEMAS DE ESPACIOS VECTORIALES POR GAUSS JORDAN____")
print("")
rn = int(input("Ingresa la dimesion del vector: "))
vn = int(input("Cuantos vectores va a ingresar: "))

vcl = np.zeros(rn)
matrixa = np.zeros((rn, vn+1))
matrixc = np.zeros((rn, vn))
matrixcopia = np.zeros((rn, vn+1))
variables = []
parametros = [" λ1", " λ2", " λ3", " λ4", " λ5", " λ6", " λ7", " λ8", " λ9"]

for i in range(rn):
  vcl[i] = float(input("Ingrese el componente "+str(i+1)+" del vector u: "))

for i in range(rn):
  matrixa[i][vn] = vcl[i]

for j in range(vn):
  for i in range(rn):
    matrixa[i][j] = float(input("Ingrese el componente "+str(i+1)+" del vector "+str(j+1)+": "))

for j in range(vn):
  for i in range(rn):
    matrixc[i][j] = matrixa[i][j]
    
for j in range(vn):
  for i in range(rn):
    matrixcopia[i][j] = matrixa[i][j]
    
for i in range(rn):
    variables.append(f"x{i}")
variables = np.array(variables)

for j in range(vn):
  print("< ", end="")
  for i in range(rn):
    print("%0.2f " %(matrixa[i][j]), end="")
  print("> V%d" %(j+1))

print_matrix('La matriz de coeficientes es:', matrixc)
print("")

print_matrix('La matriz aumentada es:', matrixa)
print("")

print("")
print("COMBINACIÓN LINEAL\n")
rango_coeficientes = np.linalg.matrix_rank(matrixc)
print("R(A) = ",rango_coeficientes)
print("")
rango_aumentada = np.linalg.matrix_rank(matrixa)
print("R(A*) = ",rango_aumentada)
print("")

# COMBINACIÓN LINEAL
if rango_aumentada == rango_coeficientes and rango_aumentada == vn:
    InvertirFilasColumnas(matrixa, rn)
    matrixA2 = np.copy(matrixa)
    f_eliminadas = FilasCeros(matrixa, rn, vn)
    n = f_eliminadas
    DiagonalAbajo(matrixA2, rn, vn)
    matrixA3 = EliminarFilas(matrixa, rn)
    DiagonalArriba(matrixA3, n, vn)
    print_matrix("\nSolución única-->Es una combinación lineal", matrixA3)
    row = matrixA3.shape[0]
    col = matrixA3.shape[1]
    print("< ", end="")
    for i in range(n):
        print("%0.2f " %(vcl[i]), end="")
    print("> = ")
    for i in range(0, row):
        print("(%4.2f)V%d" %(matrixA3[i][col - 1], i+1), end="  ")
    print("") 

elif rango_aumentada == rango_coeficientes and rango_aumentada < vn:
    matrixA2 = np.copy(matrixa)
    print("\nEl sistema tiene infinidad de soluciones-->Es una combinación lineal")
    numero_parametros = vn - rango_aumentada
    f_eliminadas = FilasCeros(matrixa, rn, vn)
    print(f_eliminadas)
    n = f_eliminadas
    DiagonalAbajo(matrixA2, rn, vn)
    matrixA3 = EliminarFilas(matrixA2, rn)
    DiagonalArriba(matrixA3, n, vn)
    print_matrix('La matriz es:', matrixA3)
    print("")
    print("En términos de parámetros libres la solución es: ")
    SolucionesInfinitas(matrixA3, variables, numero_parametros, parametros)
    SolucionSolInfinitasCL(matrixA3, vcl, vn, n)
    print("")

elif rango_aumentada != rango_coeficientes:
    print("\nEl sistema no tiene solución-->No es combinación lineal")

# GENERADOR DE ESPACIO VECTORIAL
gv = 0
arrayrandom = np.random.randint(10000, size=(rn)) 
for i in range(rn):
  matrixcopia[i][vn] = arrayrandom[i]
  
print("") 
print("GENERADOR DE ESPACIO VECTORIAL\n")    
rango_coeficientes = np.linalg.matrix_rank(matrixc)
print("R(A) = ",rango_coeficientes)
print("")
rango_aumentada = np.linalg.matrix_rank(matrixcopia)
print("R(A*) = ",rango_aumentada)
print("")

if rango_aumentada == rango_coeficientes and rango_aumentada == vn:
    print("\nTiene Solución Única--> Si genera a R", (rn))
    gv = 1

elif rango_aumentada == rango_coeficientes and rango_aumentada < vn:
    print("\nEl sistema tiene infinidad de soluciones--> Si genera a R", (rn))
    gv = 1
    
elif rango_aumentada != rango_coeficientes:
    print("\nEl sistema no tiene solución--> No genera a R", (rn))
    gv = 0

# INDEPENDENCIA LINEAL
il = 0
arrayzeros = np.zeros(rn)
for i in range(rn):
  matrixcopia[i][vn] = arrayzeros[i]

print("")    
print("INDEPENDENCIA LINEAL\n") 
rango_coeficientes = np.linalg.matrix_rank(matrixc)
print("R(A) = ",rango_coeficientes)
print("")
rango_aumentada = np.linalg.matrix_rank(matrixcopia)
print("R(A*) = ",rango_aumentada)
print("")

if rango_aumentada == rango_coeficientes and rango_aumentada == vn:
    print("\nTiene Solución Única\n")
    print("{ ", end="")
    for i in range(vn):
        print("u%d " %(i+1), end="")
    print("} ", end="")
    print("--> Es un conjunto linealmente independiente")
    il = 1

elif rango_aumentada == rango_coeficientes and rango_aumentada < vn:
    print("\nEl sistema tiene infinidad de soluciones--> Es linealmente dependiente")
    il = 0
    
elif rango_aumentada != rango_coeficientes:
    print("\nEl sistema no tiene solución--> No es linealmente independiente")
    il = 0

#BASE DE UN ESPACIO VECTORIAL
print("")    
print("BASE DE UN ESPACIO VECTORIAL\n") 

if gv==1 and il==1:
    print("{ ", end="")
    for i in range(vn):
        print("V%d " %(i+1), end="")
    print("} ", end="")
    print("Es base de R", (rn))
    
else:
    print("No es base de R", (rn))
    
del(vcl)
del(matrixa)
del(matrixc)
del(matrixcopia)
del(variables)
del(parametros)

#FIN DEL PROGRAMA