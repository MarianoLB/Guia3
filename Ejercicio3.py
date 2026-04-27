# -*- coding: utf-8 -*-
"""
Considere una barra con una carga axial aplicada de 
, como se muestra en la figura. Determine el desplazamiento axial y la tensión. Tome 
, E=30e6 N/cm2  A= 2cm2 y L= 60cm
. Use primero uno y luego dos elementos. Intente generalizarlo a 
 elementos y compare sus resultados con la solución teórica:
"""
import numpy as np
import matplotlib.pyplot as plt
from funciones_ef import k_elemental,k_global

#Datos
C=-10 #carga lineal con la distancia
L=60
E=30e6
A=2

#1 elemento
gl=1
MC=np.array([[0,1]])
MN=np.array([[0,0],[L,0]])

n=len(MC)
h=L/n #ej si fuera un elemento 60/1 =60 pero si ahora son dos, el L cambia y da L/30
Ft=(1/2) *C*h**2
F=np.zeros(len(MN))

#Hago la distrib de F en nodos
#Para los triangulitos
for t in range(len(F)):
    if t==0:
        F[t]+=(1/3)*Ft
    elif t==len(F)-1:
        F[t]+=(2/3)*Ft
    else:
        F[t]+=1*Ft
        
#Para los cuadraditos va cambiando el F

for c in range(len(F)):
    if c==0:
        F[c]+=0
    elif c==len(F)-1:
        Fcuadrado_izq=C*MN[c-1,0]*h
        F[c]+=(1/2)*Fcuadrado_izq
    else:
        Fcuadrado_izq=C*MN[c-1,0]*h
        Fcuadrado_der=C*MN[c,0]*h
        F[c]+=(1/2)*Fcuadrado_izq+(1/2)*Fcuadrado_der

#----
kel, Le = k_elemental(MC,MN,gl)   

for e in range(len(MC)):
    kel[e] = (E*A/Le[e]) * kel[e]
K = k_global(MC, MN, gl, kel)


#Vinculos
s=[1] #u1=0, nodo derecho fijo
r=[0] #u0 incognita

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

print("u=",u)
print("F aplicada=",F)
print("Fint=",Fint)
print("R=",R)    
#Lo probe para 2 elementos y anda, y parece q para mas elementos tambien.

#Ssolucion teorica en los nodos
T=abs(C)

# puntos para comparar dentro de la barra
xcomp=np.linspace(0,L,11)

# solucion teorica
dteo=(T/(6*A*E))*(xcomp**3-L**3)

# como hay un 1 elemento, interpolo lineal entre u[0] y u[1]
dFEM=u[0]+(u[1]-u[0])*(xcomp/L)

#----- GRAFICOS-----------------
plt.figure()

plt.plot(xcomp,dteo,"o-",label="Solucion teorica")
plt.plot(xcomp,dFEM,"o-",label="1 elemento")

plt.xlabel("x [cm]")
plt.ylabel("d(x) [cm]")
plt.grid(True)
plt.legend()
plt.show()