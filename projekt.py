# indeks: 175948

import matplotlib.pyplot as plt
import math
import numpy as np
import time


def my_solve(matr, vec):

    var = np.zeros(len(vec))
    var[0] = vec[0]/matr.item(0, 0)
    for i in range(1, len(vec)):
        wart = 0.0
        for j in range(0, i):
            wart += matr.item(i, j) * var[j]
        var[i] = (vec[i] - wart)/matr.item(i, i)

    return var


def diagonal_matrix_from_diagonal(vec):
    matr = np.zeros((len(vec), len(vec)))
    for i in range(len(vec)):
        matr[i, i] = vec[i]
    return matr


def diagonal_matrix_from_chosen_diagonal(vec, number):
    matr_new = np.zeros(((len(vec) + abs(number)), (len(vec) + abs(number))))
    for i in range(len(vec)):
        if number < 0:
            matr_new[i+abs(number), i] = vec[i]
        if number > 0:
            matr_new[i,i+abs(number)] = vec[i]
    return matr_new


def my_tril(matr):
    matr_new = np.zeros((len(matr[0]), len(matr[0])))
    for i in range(len(matr[0])):
        for j in range(len(matr[0])):
            if i > j:
                matr_new[i, j] = matr[i, j]
    return matr_new


def my_triu(matr):
    matr_new = np.zeros((len(matr[0]), len(matr[0])))
    for i in range(len(matr[0])):
        for j in range(len(matr[0])):
            if i < j:
                matr_new[i, j] = matr[i, j]
    return matr_new


def my_norm(vec):

    s = 0
    for l in range(len(vec)):
        s += np.power(vec[l], 2)
    norm = np.sqrt(s)

    return norm


def my_back_solve(matr, vec):

    var = np.zeros(len(vec))
    for i in range(len(vec) - 1, -1, -1):
        tmp = vec[i]
        for j in range(len(vec) - 1, i, -1):
            tmp -= var[j] * matr.item(i, j)
        var[i] = tmp / matr.item(i, i)

    return var


def diagonal(matr):
    diag = []
    for i in range(len(matr[0])):
        diag.append(matr.item(i, i))
    diag = np.array(diag)
    return diag


# ZADANIE A


a = 5 + 9
c = -1
d = -1
N = 948
b = []

for i in range(0, N):
    b.append(math.sin(i*6))

b = np.array(b)
a1 = np.array([a])
a1 = np.repeat(a1, N)
a2 = np.array([c])
a2 = np.repeat(a2, N-1)
a3 = np.array([d])
a3 = np.repeat(a3, N-2)


matrix1 = diagonal_matrix_from_diagonal(a1)
matrix2 = diagonal_matrix_from_chosen_diagonal(a2, 1)
matrix3 = diagonal_matrix_from_chosen_diagonal(a2, -1)
matrix4 = diagonal_matrix_from_chosen_diagonal(a3, 2)
matrix5 = diagonal_matrix_from_chosen_diagonal(a3, -2)
matrix = np.add(matrix1, matrix2)
matrix = np.add(matrix, matrix3)
matrix = np.add(matrix, matrix4)

A = np.add(matrix, matrix5)

# ZADANIE B METODA JACOBIEO


res = []
res.append(1)
res = np.array(res)
j = 0
norm_res = []
x = np.ones(N)
counter = 0
iterations = 0
t0 = time.time()

L = my_tril(A)


U = my_triu(A)
D = diagonal_matrix_from_diagonal(diagonal(A))
right = my_solve(D, b)

while my_norm(res) > 1e-9:
    norm_res.append(my_norm(res))
    res = np.dot(A, x) - b
    variab = np.dot((L+U), x)
    x = my_solve(-D, variab) + right
    counter = counter + 1
    j = j+1

time_jacobi = time.time() - t0
iterations = counter
print('\nZADANIE B METODA JACOBIEGO')
print('Iteracje:  ', iterations)
print('Czas:  ', time_jacobi)


# ZADANIE B METODA GAUSSA SIEDLA


res = []
res.append(1)
res = np.array(res)
j = 0
norm_res = []
x = np.ones(N)
counter = 0
iterations = 0

t0 = time.time()
L = -my_tril(A)
U = -my_triu(A)
D = diagonal_matrix_from_diagonal(diagonal(A))
right = my_solve(D-L, b)

while my_norm(res) > 1e-9:
    norm_res.append(my_norm(res))
    res = np.subtract(np.dot(A, x), b)
    variab1 = np.dot(U, x)
    variab2 = (D - L)
    x = my_solve(variab2, variab1) + right
    counter = counter + 1
    j = j + 1

time_gs = time.time() - t0
iterations = counter
print('\nZADANIE B METODA GAUSSA-SEIDLA')
print('Iteracje:  ', iterations)
print('Czas:  ', time_gs)


# ZADANIE C


a = 3
c = -1
d = -1
N = 948
b = []

for i in range(0, N):
    b.append(math.sin(i*6))

b = np.array(b)
a1 = np.array([a])
a1 = np.repeat(a1, N)
a2 = np.array([c])
a2 = np.repeat(a2, N-1)
a3 = np.array([d])
a3 = np.repeat(a3, N-2)

matrix1 = diagonal_matrix_from_diagonal(a1)
matrix2 = diagonal_matrix_from_chosen_diagonal(a2, 1)
matrix3 = diagonal_matrix_from_chosen_diagonal(a2, -1)
matrix4 = diagonal_matrix_from_chosen_diagonal(a3, 2)
matrix5 = diagonal_matrix_from_chosen_diagonal(a3, -2)
matrix = np.add(matrix1, matrix2)
matrix = np.add(matrix, matrix3)
matrix = np.add(matrix, matrix4)

A = np.add(matrix, matrix5)


# C METODA JACOBIEGO


res = []
res.append(1)
res = np.array(res)
j = 0
norm_res = []
x = np.ones(N)
counter = 0
iterations = 0

t0 = time.time()
L = my_tril(A)
U = my_triu(A)
D = diagonal_matrix_from_diagonal(diagonal(A))
right = my_solve(D, b)

while my_norm(res) > 1e-9 and counter < 100:
    norm_res.append(my_norm(res))
    res = np.dot(A, x) - b
    variab = np.dot((L+U), x)
    x = my_solve(-D, variab) + right
    counter = counter + 1
    j = j+1

time_jacobi = time.time() - t0
iterations = counter
print('\nZADANIE C METODA JACOBIEGO')
print('Iteracje:  ', iterations)
print('Czas:  ', time_jacobi)


# C METODA GAUSSA SIEDLA


res = []
res.append(1)
res = np.array(res)
j = 0
norm_res = []
x = np.ones(N)
counter = 0
iterations = 0

t0 = time.time()
L = -my_tril(A)
U = -my_triu(A)
D = diagonal_matrix_from_diagonal(diagonal(A))
right = my_solve(D-L, b)

while my_norm(res) > 1e-9 and counter < 100:
    norm_res.append(my_norm(res))
    res = np.subtract(np.dot(A, x), b)
    variab1 = np.dot(U, x)
    variab2 = (D - L)
    x = my_solve(variab2, variab1) + right
    counter = counter + 1
    j = j + 1

time_gs = time.time() - t0
iterations = counter
print('\nZADANIE C METODA GAUSSA-SEIDLA')
print('Iteracje:  ', iterations)
print('Czas:  ', time_gs)


# ZADANIE D METODA BEZPOŚREDNIEGO ROZWIĄZANIA


a = 3.0
c = -1.0
d = -1.0
N = 948
b = []

for i in range(1, N+1):
    b.append(math.sin(i*6))

b = np.array(b)
a1 = np.array([a])
a1 = np.repeat(a1, N)
a2 = np.array([c])
a2 = np.repeat(a2, N-1)
a3 = np.array([d])
a3 = np.repeat(a3, N-2)

matrix1 = diagonal_matrix_from_diagonal(a1)
matrix2 = diagonal_matrix_from_chosen_diagonal(a2, 1)
matrix3 = diagonal_matrix_from_chosen_diagonal(a2, -1)
matrix4 = diagonal_matrix_from_chosen_diagonal(a3, 2)
matrix5 = diagonal_matrix_from_chosen_diagonal(a3, -2)
matrix = np.add(matrix1, matrix2)
matrix = np.add(matrix, matrix3)
matrix = np.add(matrix, matrix4)

A = np.add(matrix, matrix5)

t0 = time.time()
x = np.ones(N)
y = np.zeros(N)
U = A.copy()
L = np.eye(N)

for k in range(0, N-1):
    for j in range(k+1, N):
        L[j, k] = U[j, k]/U[k, k]
        U[j, k:N] = U[j, k:N] - L[j, k]*U[k, k:N]

y = my_solve(L, b)
x = my_back_solve(U, y)
res = np.subtract(np.dot(A, x), b)
norm_res = my_norm(res)

time_lu = time.time() - t0
print('\nZADANIE D METODA FAKTORYZACJI LU')
print('Norma residuum ', norm_res)
print('Czas:  ', time_lu)


# ZADANIE E JACOBI


N = np.array([100, 500, 1000, 2000, 3000])
time_compare_jacobi = np.zeros(5)

for q in range(0, 5):

    a = 5 + 9
    c = -1
    d = -1
    b = []

    for i in range(0, N[q]):
        b.append(math.sin(i * 6))

    b = np.array(b)
    a1 = np.array([a])
    a1 = np.repeat(a1, N[q])
    a2 = np.array([c])
    a2 = np.repeat(a2, N[q] - 1)
    a3 = np.array([d])
    a3 = np.repeat(a3, N[q] - 2)

    matrix1 = diagonal_matrix_from_diagonal(a1)
    matrix2 = diagonal_matrix_from_chosen_diagonal(a2, 1)
    matrix3 = diagonal_matrix_from_chosen_diagonal(a2, -1)
    matrix4 = diagonal_matrix_from_chosen_diagonal(a3, 2)
    matrix5 = diagonal_matrix_from_chosen_diagonal(a3, -2)
    matrix = np.add(matrix1, matrix2)
    matrix = np.add(matrix, matrix3)
    matrix = np.add(matrix, matrix4)

    A = np.add(matrix, matrix5)

    res = []
    res.append(1)
    res = np.array(res)
    j = 0
    norm_res = []
    x = np.ones(N[q])
    counter = 0
    iterations = 0

    t0 = time.time()
    L = my_tril(A)
    U = my_triu(A)
    D = diagonal_matrix_from_diagonal(diagonal(A))
    right = my_solve(D, b)

    while my_norm(res) > 1e-9:
        norm_res.append(my_norm(res))
        res = np.dot(A, x) - b
        variab = np.dot((L + U), x)
        x = my_solve(-D, variab) + right
        counter = counter + 1
        j = j + 1

    time_compare_jacobi[q] = time.time() - t0
    iterations = counter

print('\nZADANIE E METODA JACOBIEGO')
print('Czas:  ', time_compare_jacobi)
plt.plot(N, time_compare_jacobi)
plt.ylabel('czas [s]')
plt.xlabel('liczba niewiadomych')
plt.title('Wykres zależności czasu trwania algorytmu Jacobiego od liczby niewiadomych')
plt.suptitle('Metoda Jacobiego')
plt.show()


# E GAUSS SIEGLE


time_compare_gs = np.zeros(5)
N = np.array([100, 500, 1000, 2000, 3000])
for q in range(0, 5):

    a = 5 + 9
    c = -1
    d = -1
    b = []

    for i in range(0, N[q]):
        b.append(math.sin(i * 6))

    b = np.array(b)
    a1 = np.array([a])
    a1 = np.repeat(a1, N[q])
    a2 = np.array([c])
    a2 = np.repeat(a2, N[q] - 1)
    a3 = np.array([d])
    a3 = np.repeat(a3, N[q] - 2)

    matrix1 = diagonal_matrix_from_diagonal(a1)
    matrix2 = diagonal_matrix_from_chosen_diagonal(a2, 1)
    matrix3 = diagonal_matrix_from_chosen_diagonal(a2, -1)
    matrix4 = diagonal_matrix_from_chosen_diagonal(a3, 2)
    matrix5 = diagonal_matrix_from_chosen_diagonal(a3, -2)
    matrix = np.add(matrix1, matrix2)
    matrix = np.add(matrix, matrix3)
    matrix = np.add(matrix, matrix4)

    A = np.add(matrix, matrix5)

    res = []
    res.append(1)
    res = np.array(res)
    j = 0
    norm_res = []
    x = np.ones(N[q])
    counter = 0
    iterations = 0

    t0 = time.time()
    L = -my_tril(A)
    U = -my_triu(A)
    D = diagonal_matrix_from_diagonal(diagonal(A))
    right = my_solve(D - L, b)

    while my_norm(res) > 1e-9:

        norm_res.append(my_norm(res))
        res = np.subtract(np.dot(A, x), b)
        variab1 = np.dot(U, x)
        variab2 = (D - L)
        x = my_solve(variab2, variab1) + right
        counter = counter + 1
        j = j + 1

    time_compare_gs[q] = time.time() - t0
    iterations = counter

print('\nZADANIE E METODA GAUSSA-SEIDLA')
print('Czas:  ', time_compare_gs)
plt.plot(N, time_compare_gs)
plt.ylabel('czas [s]')
plt.xlabel('liczba niewiadomych')
plt.title('Wykres zależności czasu trwania algorytmu Gaussa-Seidla od liczby niewiadomych')
plt.suptitle('Metoda Gaussa Seidla')
plt.show()




# ZADANIE E FAKTORYZACJA LU


N = np.array([100, 500, 1000, 2000, 3000])
time_compare_lu = np.zeros(5)

for q in range(0, 5):

    a = 5.0 + 9.0
    c = -1.0
    d = -1.0
    b = []

    for i in range(1, N[q] + 1):
        b.append(math.sin(i * 6))

    b = np.array(b)
    a1 = np.array([a])
    a1 = np.repeat(a1, N[q])
    a2 = np.array([c])
    a2 = np.repeat(a2, N[q] - 1)
    a3 = np.array([d])
    a3 = np.repeat(a3, N[q] - 2)

    matrix1 = diagonal_matrix_from_diagonal(a1)
    matrix2 = diagonal_matrix_from_chosen_diagonal(a2, 1)
    matrix3 = diagonal_matrix_from_chosen_diagonal(a2, -1)
    matrix4 = diagonal_matrix_from_chosen_diagonal(a3, 2)
    matrix5 = diagonal_matrix_from_chosen_diagonal(a3, -2)
    matrix = np.add(matrix1, matrix2)
    matrix = np.add(matrix, matrix3)
    matrix = np.add(matrix, matrix4)

    A = np.add(matrix, matrix5)

    t0 = time.time()
    x = np.ones(N[q])
    y = np.zeros(N[q])
    U = A.copy()
    L = np.eye(N[q])

    for k in range(0, N[q] - 1):
        for j in range(k + 1, N[q]):
            L[j, k] = U[j, k] / U[k, k]
            U[j, k:N[q]] = U[j, k:N[q]] - L[j, k] * U[k, k:N[q]]

    y = my_solve(L, b)
    x = my_back_solve(U, y)
    res = np.subtract(np.dot(A, x), b)
    norm_res = my_norm(res)

    time_compare_lu[q] = time.time() - t0

print('\nZADANIE E METODA FAKTORYZACJI LU')
print('Czas:  ', time_compare_lu)
plt.plot(N, time_compare_lu)
plt.ylabel('czas [s]')
plt.xlabel('liczba niewiadomych')
plt.title('Wykres zależności czasu trwania algorytmu FAKTORYZACJI LU od liczby niewiadomych')
plt.suptitle('Metoda Faktoryzacji LU')
plt.show()

