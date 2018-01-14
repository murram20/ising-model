# ising-model# -*- coding: utf-8 -*-
"""
Created on Sun Jan 14 16:48:35 2018

@author: Murray
"""

import numpy as np
import random as rd
import matplotlib.pyplot as plt

'''
def collinear_matrix(matrix,dimension_size):
    spindown = -1
    for i in range(0,dimension_size):
        matrix.append(spin_down) #every second row is spin down
        i = i+2 #every second row is the same (collinear)
    return np.array(matrix) 
#print collinear_matrix(10) # m = 10 so its a square 10 by 10 matrix
#-------------------------------------------------------------
'''
L=10 #L = width of the lattice
D = 2 # D = lattice dimension 
def random_matrix(m):
    if D == 1:
        ran_matrix = np.zeros(m)
        for i in range(0,m):
            ran_matrix[i] = 1
    elif D == 2:
        ran_matrix = np.zeros((m,m))
        for i in range(0,m):
            for j in range(0,m):
                ran_matrix[i,j] = 1
    elif D == 3:
        ran_matrix = np.zeros((m,m,m))
        for i in range(0,m):
            for j in range(0,m):
                for z in range(0,m):
                    ran_matrix[i,j,z] = 1
    return ran_matrix
'''
def random_matrix(m): #creates a matrix of size m same as before but 
#randomly assigns a 1 or -1 to each point in the matrix 
    ran_matrix = [] #random matrix
    x = 1.0 #spin up
    y = -1.0 # spin down
    limit = 0.5 # if value is above this it gives a +1 and if value is below it gives a -1
    for i in range(m):
        R=[]
        for j in range(m):
            if rd.random()>limit: #Selects a random number between 0 and 1, when random value > 0.5 we get a +1
                    R.append(x) # x = +1
            else:
                    R.append(y) # when random value is less than 0.5 we get a -1 = y
        ran_matrix.append(R) 
    return np.array(ran_matrix)    #returns the random spin up or spin down matrix
#print random_matrix(10) #printing random matrix
   ''' 
#---------------------------------------------------------------------------------------
    
def sweeps(matrix, T):    #function that sweeps through the matrix changing spins if Hamiltonian < 0 when spin is changed
    # or if a random number is < P_flip
    Energy = 0 
    Magnetism_list = []
    Energy_list = []
    for z in range(L**D * 50):     #for z in range(no. of sweeps). This for loop repeats this 'no. of sweeps' times
        if D == 1:
            i = rd.randint(0, L-1)    #L is matrix size. L = 10 means a 10 by 10 matrix. Starts counting from 0  
            Boundary=matrix[(i+1)%L]+matrix[(i-1)%L]
            matrix1 = matrix[i]
            coordinate = i
        if D == 2:
            i = rd.randint(0, L-1)
            j = rd.randint(0, L-1)    # i is the x coordinates or columns j is the rows or y coordinates
            Boundary=matrix[(i+1)%L,j]+matrix[i,(j+1)%L]+matrix[(i-1)%L,j] + matrix[i, (j-1)%L]
            matrix1 = matrix[i,j]
            coordinate = (i,j)
        if D == 3:
            i = rd.randint(0, L-1)
            j = rd.randint(0, L-1)
            k = rd.randint(0, L-1)
            Boundary=matrix[(i+1)%L,j,k]+matrix[i,(j+1)%L, k]+matrix[(i-1)%L,j, k] + matrix[i, (j-1)%L, k] + matrix[i, j, (k+1) % L] + matrix[i, j, (k-1) % L]
            matrix1 = matrix[i, j, k]
            coordinate = (i,j,k)
        Energy_diff=2*matrix1*Boundary       #formula given in slides. Energy is the Hamiltonian
        P_flip = np.exp(-Energy_diff/T) #the probability of the spin flipping
        if Energy_diff < 0:    #If a flip is energy efficient (a flip will end up being less energy) then we flip it    
            matrix1 *= -1     #matrix gets flipped   
        elif rd.random() < P_flip:  # a random number between 0 and 1 is chosen. If it is < P_flip then we flip the spin. If not then we leave it 
            matrix1 *= -1 # The spin is flipped if random number is < P_flip
        matrix[coordinate] = matrix1
    for x in range(L**D*5):
        if D == 1:
            i = rd.randint(0, L-1)    #L is matrix size. L = 10 means a 10 by 10 matrix. Starts counting from 0  
            Boundary=matrix[(i+1)%L]+matrix[(i-1)%L]
            matrix1 = matrix[i]
            coordinate = i
        if D == 2:
            i = rd.randint(0, L-1)
            j = rd.randint(0, L-1)    # i is the x coordinates or columns j is the rows or y coordinates
            Boundary=matrix[(i+1)%L,j]+matrix[i,(j+1)%L]+matrix[(i-1)%L,j] + matrix[i, (j-1)%L]
            matrix1 = matrix[i,j]
            coordinate = (i,j)
        if D == 3: 
            i = rd.randint(0, L-1)
            j = rd.randint(0, L-1)
            k = rd.randint(0, L-1)
            Boundary=matrix[(i+1)%L,j,k]+matrix[i,(j+1)%L, k]+matrix[(i-1)%L,j, k] + matrix[i, (j-1)%L, k] + matrix[i, j, (k+1) % L] + matrix[i, j, (k-1) % L]
            coordinate = (i,j,k)
        Energy_list.append(-Boundary * matrix[coordinate])
        Energy += -Boundary * matrix[coordinate]
        Magnetism_list.append(matrix[coordinate])        
    Energy_var = np.var(Energy_list)    
    Energy /= L**D*5  # averaging the energy of the lattice
    Magnetism_variance = np.var(Magnetism_list)
    return matrix, Energy, Energy_var, Magnetism_variance
#print sweeps(collinear_matrix(L), 2000)    
#------------------------------------------------------------------------------------------

#Trying to add up the elements in the matrices to find the magnetisation of each matrix
T=0.1 #T = starting temperature
#print collinear_matrix(L)   #original collinear matrix
#print sweeps(collinear_matrix(L),T)  #updated collinear matrix
list_of_elements=np.array(sweeps(random_matrix(L),T)) # makes an array or list out of the matrix elements so they can be added up     
#magnetisation=np.sum(list_of_elements) #adding up all the elements in the array/list = the magnetisation                     
#print magnetisation
#Creating a while loop to calculate magnetisations for each matrix as T increases in steps of 0.1
magnetisation = []
temperature = []
energy_list = []
simplematrix = random_matrix(L)
energy_variance_list = []
specific_heat = []
magnetic_susceptibility = []
magnetic_variance = []
while (T < 6.0): #while loop to find magnetisation of numerous matrices as T increases and plot them
    equilibriated_lattice, energy , energy_var, mag_var = sweeps(simplematrix, T)
    list_of_elements1=np.array(equilibriated_lattice) #makes the elements into an array so they can be added up
    magnetisation1=np.sum((list_of_elements1)/L**D) #magnetisation = the sum of all the elements in a matrix. We divide by L**2 to get our results between 0 and 1
    magnetisation2=abs(magnetisation1) # We dont want any negative values of magnetisation so we use the absolute value
    print T # prints the temperature steps
    print magnetisation1 #Prints the magnetisation of different matrices of increasing temperature in steps of 0.1
    temperature.append(T)
    magnetisation.append(magnetisation2)
    energy_list.append(energy)
    energy_variance_list.append(energy_var)
    magnetic_variance.append(mag_var)
    T += 0.01 #T increases in steps of 0.1
    
for i in range(0,len(temperature)):
    #temperature[i] /= 3.4
    specific_heat.append(energy_variance_list[i] / (temperature[i]**2))
    magnetic_susceptibility.append(magnetic_variance[i] / temperature[i])
plt.figure(1)
plt.plot(temperature, magnetisation) # magnetisation Vs Temperature
plt.ylabel('Mean Magnetisation (magnetisation/no. of elements)')
plt.xlabel('Temperature (J/K_B)')
plt.title('Mean Magnetisation Vs Temperature')
plt.grid(True)
plt.figure(2)
plt.plot(temperature, energy_list)
plt.ylabel('Energy (joules)')
plt.xlabel('Temperature (J/K_B)')
plt.title('Energy Vs Temperature')
plt.grid(True)
plt.figure(3)
plt.plot(temperature, specific_heat)
plt.ylabel('Specific Heat (J/kg K)')
plt.xlabel('Temperature (J/K_B)')
plt.title('Specific Heat Vs Temperature')
plt.grid(True)
plt.figure(4)
plt.plot(temperature, magnetic_susceptibility)
plt.ylabel('Magnetic Susceptibility')
plt.xlabel('Temperature (J/K_B)')
plt.title('Magnetic Susceptibility Vs Temperature')
plt.grid(True)
plt.show()
