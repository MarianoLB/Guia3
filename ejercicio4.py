# -*- coding: utf-8 -*-
"""
Ejercicio 4:
Determine los desplazamientos y rotaciones, y fuerzas y torques de vínculos para el sistema de la figura. Tome 
   E= 210 Gpa I= 2e-4 m4 
"""
import numpy as np
from funciones_ef import k_elemental,k_global
#Datos
E= 210e9 #Pa
I=2e-4 #m4
#L=6 #m
gl=2 #dy y dphi
elemento='viga'
kres=200e3
P=50e3

MC=np.array([[0,1],[1,2]])
MN=np.array([[0,0],[3,0],[6,0]])

kel,Le =k_elemental(MC,MN,gl,elemento)

for e in range(len(MC)):
    kel[e]=(E*I/Le[e]**3)*kel[e]

K=k_global(MC,MN,gl,kel)
#El resorte solo actua sobre el desplazamiento vertical del nodo 3, o sea d3y.
K[4,4]+=kres

F=np.zeros(gl*len(MN))
F[4]=-P
#vinculos
s=[0,1,2]#d0=0,phi0=0, d1=0
r=[3,4,5]#phi1, d2, phi2

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
print(Reac) #el del resorte no aparece pero si hago kres*u(4) y sumo todo me da 0