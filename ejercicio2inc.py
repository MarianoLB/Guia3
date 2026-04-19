
"""
Ejercicio 2 Guia 3:
  Considere el entramado mostrado en la figura, con una fuerza aplicada de 20kN
. Calcule los desplazamientos de cada uno de los nodos y las tensiones que sufre cada elemento. Todos los elementos tienen 
 E=210 GPA y una sección de  10CM2, excepgto el elemento 3 q tiene 20cm2
, excepto el elemento 3, que tiene una sección de 
. Los elementos 2 y 5 tienen una longitud de 8 metros y el elemento 3 de 4 metros.
"""
import numpy as np

#DATOS
E=210e9 #Gpa
A= 0.001#m2 o sea 10 cm2
A3=0.002 #m2
gl=2 # cada nodo tiene tiene 2 gl, x e y
nodxel=2 #nodo por elemento
#Long de los elementos
L1=8
L2=4
MN=np.array([[0,0],[L1,L2],[L1,0],[L1+L1,0]])
MC=np.array([[0,1],[2,0],[1,2],[1,3],[3,2]]) #Matriz de conectividad, los nodos de los elenentos, tengo 5 el.
#Matrices ayuda
x=np.zeros([len(MC),nodxel])
y=np.zeros([len(MC),nodxel])
L=np.zeros(len(MC))
ang=np.zeros(len(MC))
k=np.zeros(len(MC))
kel=np.zeros([nodxel*gl,nodxel*gl])
kglob=np.zeros([gl*len(MN),gl*len(MN)])


for e in range(len(MC)):
    for i in range (len(MC[e])):
       #print(i)
       x[e,i]= np.array(MN[MC[e,i],0])
       y[e,i]= np.array(MN[MC[e,i],1])
     
    ang[e]=np.arctan2((y[e,1]-y[e,0]),(x[e,1]-x[e,0]))
    L[e]=np.sqrt((x[e,1]-x[e,0])**2 + (y[e,1]-y[e,0])**2 )
    if e==2: #pa el elemento 3 q seria mi 2 tien dif area.
        k[e]=E*A3/L[e] #Estas son mis k
    else:
        k[e]=E*A/L[e]
#Agrupo las matrices- Armo la general.
for e in range(len(MC)):
    c= np.cos(ang[e])
    s= np.sin(ang[e])
    #kel esta escritco en coordendas locales del elemento, las tengo q pasar a global
    kel = k[e] * np.array([ [ c**2,  c*s, -c**2, -c*s],[ c*s,  s**2, -c*s,  -s**2], 
                            [-c**2, -c*s,  c**2,  c*s],[-c*s, -s**2,  c*s,   s**2]])
    kel[np.abs(kel)<1e-10]=0
    for i in range(len(MC[e])):
        rangoi = np.linspace(i*gl,(i+1)*gl-1,gl,dtype=int) #este me da los indices locales tipo x y de cada nodo i dentro de kel.
        rangoni = np.linspace(MC[e,i]*gl,(MC[e,i]+1)*gl-1,gl,dtype=int) # armo este para los indices globales de los gl, onda a q posicion de la matriz global corresponde el nodo local i
        for j in range(len(MC[e])):
            rangoj = np.linspace(j*gl, (j+1)*gl-1, gl,dtype=int)
            rangonj = np.linspace(MC[e,j]*gl, (MC[e,j]+1)*gl-1, gl,dtype=int)
            kglob[np.ix_(rangoni, rangonj)] = kglob[np.ix_(rangoni, rangonj)] + kel[np.ix_(rangoi, rangoj)]
        
kglob[np.abs(kglob)<1e-10]=0
#kglob=kglob/(E*A)

#s los vinculos
s=[0,1,7] #Vinculos o sea en ux0 y uy0 es 0 y en uy3y q es el 7
r=[2,3,4,5,6] #Incognitas
u=np.zeros(len(kglob)) #u0x,u0y, u1x,u1y...
F=np.zeros(len(kglob)) #F0x F0y,...
F[5]=-20e3 # este seria mi F2y

krr=kglob[np.ix_(r,r)]
krs=kglob[np.ix_(r,s)]

b=F[r]-krs@u[s]
Ku=np.linalg.solve(krr,b)
u[r]=Ku
F=kglob@u
F[np.abs(F)<1e-10]=0

#Los quiero pasar en 2 columna para ver los efectos en x e y en cada nodo.
U=u.reshape(len(MN), gl)
Fnod=F.reshape(len(MN), gl)
Umm=U*1000 # paso a mm

print("Desplazamientos de nodos en [mm] \n [  ux  ,          uy]:")
print(Umm)
print("Fuerzas nodales \n [  Fx  ,    Fy]:")
print(Fnod)

#-------------Me faltaria calcular las tensiones-------------------