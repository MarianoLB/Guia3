# -*- coding: utf-8 -*-
"""
Funciones para FEM
"""
import numpy as np

def k_elemental(MC,MN,gl):
    
    nxel = len(MC[0])
    x = np.zeros([len(MC), nxel])
    y = np.zeros([len(MC), nxel])
    L = np.zeros(len(MC))

    if gl==1:
        kel = np.zeros([len(MC), nxel*gl, nxel*gl])   # una kel por elemento
        for e in range(len(MC)):
            for i in range(len(MC[e])):
                x[e,i] = MN[MC[e,i],0]
                y[e,i] = MN[MC[e,i],1]
                
                L[e]=np.sqrt((x[e,1]-x[e,0])**2 + (y[e,1]-y[e,0])**2 )  

                kel[e] = np.array([ [1,-1],
                                    [-1,1]  ])
                
                kel[e][np.abs(kel[e]) < 1e-10] = 0    
    
    elif gl==2:
        nxel = len(MC[0])
        ang = np.zeros(len(MC))
        L = np.zeros(len(MC))
        kel = np.zeros([len(MC), nxel*gl, nxel*gl])   # una kel por elemento

        for e in range(len(MC)):
            for i in range(len(MC[e])):
                x[e,i] = MN[MC[e,i],0]
                y[e,i] = MN[MC[e,i],1]
                
                ang[e] = np.arctan2((y[e,1]-y[e,0]), (x[e,1]-x[e,0]))
                L[e]=np.sqrt((x[e,1]-x[e,0])**2 + (y[e,1]-y[e,0])**2 )  
                c = np.cos(ang[e])
                s = np.sin(ang[e])
                kel[e] = np.array([ [ c**2,   c*s,  -c**2,  -c*s],[ c*s,   s**2,  -c*s,   -s**2],
                                   [-c**2, -c*s,   c**2,    c*s],  [-c*s,  -s**2,  c*s,    s**2]])

            kel[e][np.abs(kel[e]) < 1e-10] = 0

    return kel,L

#-----------------------------------------------------------------------
def k_global(MC, MN, gl, kel):
    kglob = np.zeros([gl*len(MN), gl*len(MN)])

    for e in range(len(MC)):
        for i in range(len(MC[e])):
            rangoi = np.linspace(i*gl, (i+1)*gl-1, gl, dtype=int)
            rangoni = np.linspace(MC[e,i]*gl, (MC[e,i]+1)*gl-1, gl, dtype=int)
            for j in range(len(MC[e])):
                rangoj = np.linspace(j*gl, (j+1)*gl-1, gl, dtype=int)
                rangonj = np.linspace(MC[e,j]*gl, (MC[e,j]+1)*gl-1, gl, dtype=int)
                kglob[np.ix_(rangoni, rangonj)] = kglob[np.ix_(rangoni, rangonj)] + kel[e][np.ix_(rangoi, rangoj)]

    kglob[np.abs(kglob) < 1e-10] = 0
    return kglob