import numpy as np

from math import *
import matplotlib.pyplot as plt
G=1*10**9
IxTable=[]
IyTable=[]
Ix2Table=[]
Iy2Table=[]
ztab=[]


chord= 3.451
a= 0.45 * chord
b= 0.1011 * chord
c= 0.1011 * chord
d=0.45 * chord
z= 0
t= 0.001 #thin walled

centx =[]
centz = []

# Stringer dimension
Sh = 0.05
Sb = 0.02
St = 0.01
#AS = Sh * St + Sb * St
AS = 0.021000
n =1

t2 = 0.001
AS2 = 0.0138


while z<10.8695:
    taper = (1- (0.6*z)/10.8695)
    a= a*taper
    b= b*taper
    c=c*taper
    d=d*taper
    
    
    #Areas
    Aa = a*t
    Ab = b*t
    Ac = c*t
    Ad = d*t

    #center x
    Xa= a/2
    Xb= 0
    Xc= a
    Xd= d/2

    #center y
    Ya= 0
    Yb= b/2
    Yc= c/2
    Yd= b

    

    #---------------------------------------------------------------------------------------------
    #centroid
    x= (Aa* Xa + Ab* Xb +Ac* Xc +Ad* Xd )/(Aa+Ab+Ac+Ad)
    centx.append(x)

    y= (Aa* Ya + Ab* Yb +Ac* Yc +Ad* Yd )/(Aa+Ab+Ac+Ad)
    centz.append(y)

    # Stringers x
    SDy = b / 2
    SIx = AS * SDy ** 2
    SIx2 = AS2 * SDy **2
    

    #---------------------------------------------------------------------------------------------
    #moment of inertia x
    # moment of inertia x
    Ix = Aa * (Ya - y) ** 2 + 1 / 12 * t * c ** 3 + Ac * (Yc - y) ** 2 + 1 / 12 * t * b ** 3 + Ab * (Yb - y) ** 2 + Ad * (Yd - y) ** 2 + n * SIx
    Ix2 = Aa * (Ya - y) ** 2 + 1 / 12 * t2 * c ** 3 + Ac * (Yc - y) ** 2 + 1 / 12 * t2 * b ** 3 + Ab * (Yb - y) ** 2 + Ad * (Yd - y) ** 2 + n * SIx2
    #moment of inertia y
    Iy= 1/12 * t * a**3 + Aa* (Xa-x)**2 + Ab*(Xb-x)**2 + Ac*(Xc-x)**2 + 1/12 * t* d**3 + Ad*(Xd-x)**2
    Iy2= 1/12 * t2 * a**3 + Aa* (Xa-x)**2 + Ab*(Xb-x)**2 + Ac*(Xc-x)**2 + 1/12 * t2* d**3 + Ad*(Xd-x)**2
    

    #--------------------------------------------------------------------------------------------
    IxTable.append(Ix)
    Ix2Table.append(Ix2)
    ztab.append(z)
    IyTable.append(Iy)
    Iy2Table.append(Iy2)
    z=z+0.018116
    a=0.45 * chord
    b=0.1011 * chord
    c=0.1011 * chord
    d=0.45 * chord










centx = np.array(centx)
centz = np.array(centz)

Ixx = np.array(IxTable)
Izz = np.array(IyTable)

Ixx2 = np.array(Ix2Table)
Izz2 = np.array(Iy2Table)


# ---------------------------------------------------------------------------------------------


# Weight sihiteght
dens = 2810
Topplate = t*(1.5525+0.4*1.5525)*10.8695/2
Sideplate = t*(0.36+0.4*0.36)*10.8695/2
Stringer = AS*10.8695
Weight = dens*(Topplate*2+Sideplate*2+Stringer*n)
#print(Weight)
