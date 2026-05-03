# -*- coding: utf-8 -*-
"""
Ejercicio 5

"""
import numpy as np
from funciones_ef import k_elemental,k_global
#Datos
E= 29e6 #Psi
I=200 #in4
#L=6 #m
gl=2 #dy y dphi
elemento='viga'
w=200/12 #lb/ft paso a pulgada

MC=np.array([[0,1],[1,2]])
MN=np.array([[0,0],[15*12,0],[30*12,0]])

kel,Le =k_elemental(MC,MN,gl,elemento)

for e in range(len(MC)):
    kel[e]=(E*I/Le[e]**3)*kel[e]

K=k_global(MC,MN,gl,kel)
#El resorte solo actua sobre el desplazamiento vertical del nodo 3, o sea d3y.

F=np.zeros(gl*len(MN))

#Distribcuion de cargas
for e in range(len(MC)):
    L=Le[e]
    n1=MC[e,0] #nodo 1 del elemento
    n2=MC[e,1] #nodo 2
   
    d1=n1*gl #como tengo q hacer d y phi, d0 phi0y asi..
    phi1=n1*gl+1
    d2=n2*gl
    phi2=n2*gl+1
    
    F[d1]+=-w*L/2
    F[phi1]+=-w*L**2/12
    
    F[d2]+=-w*L/2
    F[phi2]+=w*L**2/12

#vinculos
s=[0,1,4]#d0=0,phi0=0, d1=0
r=[2,3,5]#phi1, d2, phi2

u=np.zeros(len(K))

krr=K[np.ix_(r,r)]
krs=K[np.ix_(r,s)]

b=F[r]-krs@u[s]
Ku=np.linalg.solve(krr,b)

u[r]=Ku

Fint=K@u
Fint[np.abs(Fint)<1e-10]=0

R=Fint-F
R[np.abs(R)<1e-10]=0

U=u.reshape(len(MN),gl)
Fnod=Fint.reshape(len(MN),gl)
Reac=R.reshape(len(MN),gl)

print("Le=",Le)
print("K=")
print(K)
print("u=",u)
print("U nodal [dy,phi]=")
print(U)
print("F aplicada=",F)
print("Fint=",Fint)
print("Reacciones=")
print(Reac)

