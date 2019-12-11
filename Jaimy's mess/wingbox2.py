import pandas as pd
import numpy as np
import scipy as sp
from scipy import interpolate
import matplotlib.pyplot as plt
import scipy.integrate as integrate
from Inertia import Inertia_Calc
from LOADS import span, c_interpolate, N_v, N_m, T_v, T_m, IntM, N_v2, N_m2, T_v2, T_m2, IntM2

#Wing box properties
t_flange = 0.003
t_spar = 0.004
AS = np.zeros(600)
AS_lower = np.zeros(600)
b = np.linspace(0.349, 0.1396, 600)
bshort = np.linspace(0.3399, 0.1359644, 600)
a = np.linspace(1.55, 0.620, 600)
Am = (0.101*c_interpolate(span)+0.0985*c_interpolate(span))*(0.45*c_interpolate(span)-0)/2
S = 0.101*c_interpolate(span)+0.0985*c_interpolate(span)+2*np.sqrt((0.008625*c_interpolate(span))**(2)+(0.45*c_interpolate(span))**(2))

#Material properties
dens = 2810
E1 = 71*10**9
E2 = 107*10**9
G1 = 26*10**9
G2 = 42.1*10**9
vp = 0.33
yieldstress = 483*10**6

def stress_at_any_point(x, y, I_xx, I_yy, centx, centy, Nm, Tm):
    x = np.linspace(x, 0.4*x,600)
    y = np.linspace(y, 0.4 * y, 600)
    x = x - centx
    y = y - centy
    stress = ((Nm *y)/I_xx) - ((Tm *x)/I_yy)
    return stress

#Rib locations
Rib_locations = np.arange(0, 10.8695, 0.17)
Rib_locations[-1] = 10.8695
Lengths = []
spanlocations = []
for i in range(1, len(Rib_locations)):
    Length = Rib_locations[i] - Rib_locations[i-1]
    Lengths.append(Length)
    spanlocations.append(int(Rib_locations[i-1]/0.018116))
Lengths = np.array(Lengths)
spanlocations.append(600)


#Rib section start, Rib section stop, n_s, A_s, P_s, H_s, b_s, Pitch, Actual pitch
Stringer_range_upper = np.arange(17, 30 ,2)
Stringer_range_lower = np.arange(17, 30 ,2)
#Stringer_range[0] = 2

#Upper skin, column and skin buckling
def Stringers_upper(Ribs, AS, AS_lower):
    sections = np.array([0,0,0,0,0,0])
    for x in range(np.size(Rib_locations)-1):

        Length = Lengths[x]
        spanloc_start = spanlocations[x]
        spanloc_end = spanlocations[x+1]
        
        for i in Stringer_range_upper:
            A_s = 0.000568
            AS_section = A_s * i
            AS[spanloc_start:spanloc_end] = AS_section
            Areas = 2*a*t_flange+2*b*t_spar+AS+AS_lower
            kc = 4
            t_s = 0.004
            j = 1/4
            H_s = A_s/((1+j)*t_s)
            b_s = j*H_s
            centroid_s = (0.5*(H_s**2)*t_s)/A_s
            Ixx_stringers = (1/12*(H_s**3)*t_s+A_s*(centroid_s)**2)*i

            centx, centy, Ixx, Iyy = Inertia_Calc(t_flange, t_spar, AS, AS_lower)
            stress = stress_at_any_point(0.775, 0.1745, Ixx, Iyy, centx, centy, N_m, T_m)
            I = (abs(stress[spanloc_start])*Length**2)/(0.25*np.pi**2*E1)
            pitch_stringer = t_flange*((abs(stress[spanloc_start])*Areas[spanloc_start]*12*(1-vp**2))/(np.pi**2*kc*E1))**-0.5
            Apitch = a[spanloc_start]/i
            distance = a[spanloc_start] - b_s*i
            if Ixx_stringers > I:
                section = np.array([Rib_locations[x], Rib_locations[x+1], i, AS_section, distance, Apitch/pitch_stringer])
                sections = np.vstack((sections, section))
                break
    sections = np.delete(sections, 0,0)

    df = pd.DataFrame(sections)
    df.columns = ['Section start', 'Section end', 'Amount','AS_section', 'Distance', 'APitch/RPitch']
    return df, AS

#Lower skin, column and skin buckling
def Stringers_lower(Ribs, AS, AS_lower):

    sections_low = np.array([0,0,0,0,0,0])
    for x in range(np.size(Rib_locations)-1):

        Length = Lengths[x]
        spanloc_start = spanlocations[x]
        spanloc_end = spanlocations[x+1]

        #stress = stress_at_any_point(0.775, 0.1745, Ixx, Iyy, N_m, T_m)
        for m in Stringer_range_lower:
            A_s = 0.0005
            AS_lowsection = A_s * m
            AS_lower[spanloc_start:spanloc_end] = AS_lowsection
            Areas = 2*a*t_flange+2*b*t_spar+AS+AS_lower
            kc = 4
            t_s = 0.004
            j = 1/4
            H_s = A_s/((1+j)*t_s)
            b_s = j*H_s
            centroid_s = (0.5*(H_s**2)*t_s)/A_s
            Ixx_stringers = (1/12*(H_s**3)*t_s+A_s*(centroid_s)**2)*m

            centx, centy, Ixx, Iyy = Inertia_Calc(t_flange, t_spar, AS, AS_lower)
            stress_low = stress_at_any_point(0.775, -0.1745, Ixx, Iyy, centx, centy, N_m2, T_m2)
            I = (abs(stress_low[spanloc_start])*Length**2)/(0.25*np.pi**2*E1)
            pitch_stringer_low = t_flange*((abs(stress_low[spanloc_start])*Areas[spanloc_start]*12*(1-vp**2))/(np.pi**2*kc*E1))**-0.5
            Apitch = a[spanloc_start]/m
            distance = a[spanloc_start] - b_s*m
            if Ixx_stringers > I:
                section_low = np.array([Rib_locations[x], Rib_locations[x+1], m, AS_lowsection, distance, Apitch/pitch_stringer_low])
                sections_low = np.vstack((sections_low, section_low))
                break
    sections_low = np.delete(sections_low, 0,0)

    df_low = pd.DataFrame(sections_low)
    df_low.columns = ['Section start', 'Section end', 'Amount','AS_section', ' Distance', 'APitch/RPitch']
    return df_low, AS_lower

#Iteration process
df, AS = Stringers_upper(Rib_locations, AS, AS_lower)
df_low, AS_lower = Stringers_lower(Rib_locations, AS, AS_lower)

df, AS = Stringers_upper(Rib_locations, AS, AS_lower)
df_low, AS_lower = Stringers_lower(Rib_locations, AS, AS_lower)

df, AS = Stringers_upper(Rib_locations, AS, AS_lower)
df_low, AS_lower = Stringers_lower(Rib_locations, AS, AS_lower)

A_s = 0.0005
t_s = 0.004
j = 1/4
H_s = A_s/((1+j)*t_s)
b_s = j*H_s
print(H_s, b_s)

print('Upper skin')
print(df)
print('Lower skin')
print(df_low)

#Check if maximum stress does not exceed yield stress
centx, centy, Ixx, Iyy = Inertia_Calc(t_flange, t_spar, AS, AS_lower)
stressmaxtopright = stress_at_any_point(0.775, 0.1745, Ixx, Iyy, centx, centy, N_m, T_m) 
stressmaxbottomleft = stress_at_any_point(0.775, 0.1745, Ixx, Iyy, centx, centy, N_m, T_m)
stress2maxbottomright = stress_at_any_point(0.775, -0.1745, Ixx, Iyy, centx, centy, N_m, T_m) 
stress2maxtopleft = stress_at_any_point(-0.775, 0.1745, Ixx, Iyy, centx, centy, N_m, T_m)

if abs(stressmaxtopright[0]) > yieldstress:
    print("Yield stress exceeded top right load case 1")
if abs(stressmaxbottomleft[0]) > yieldstress:
    print("Yield stress exceeded bottom left load case 1")

if abs(stress2maxbottomright[0]) > yieldstress:
    print("Yield stress exceeded bottom right load case 2")
if abs(stress2maxtopleft[0]) > yieldstress:
    print("Yield stress exceeded top left load case 2")
    
#Moment of inertia top plate
'''
n_s = 25
A_s = AS/n_s 
t_s = 0.004
j = 1/4
H = A_s/((1+j)*t_s)
b = j*H

#L_s = (A_s/2)/t_s
#centroid_s = (0.5*(H**2)*t_s)/A_s
#Ixx_stringers = (1/12*(H**3)*t_s+A_s*(centroid_s)**2)*n_s
'''

#Web buckling
kv = 1.5
sf = 1.1
ks = 5.5

tau_ave_shear = N_v/(b*t_spar+bshort*t_spar)

tau_max_shear = kv*tau_ave_shear
tau = IntM/(2*Am*t_spar)

tau_max_max = tau + tau_max_shear
pitch_stiffener = t_spar*((12*(1-vp**2)*tau_max_max*sf)/(np.pi**2*ks*E1))**-0.5

#Varying pitch
pitch_locations = [0, 2.5, 5.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.8695]
pitch_Lengths = []
pitch_span = []
for i in range(1, len(pitch_locations)):
    pitch_Length = pitch_locations[i] - pitch_locations[i-1]
    pitch_Lengths.append(pitch_Length)
    pitch_span.append(int(pitch_locations[i-1]/0.018116))
pitch_span.append(600)
stiffeners = 0
for i in range(len(pitch_locations)-1):
    stiffener = pitch_Lengths[i]/pitch_stiffener[pitch_span[i]]
    stiffeners = stiffeners + stiffener
print("Amount of stiffeners =", stiffeners*2)

#Deflection
upper_change = []
lower_change = []
#Stringer reduction positions
for i in range(len(AS)-1):
    if AS[i+1] < AS[i]:
        upper_change.append(i)

for j in range(len(AS_lower)-1):
    if AS_lower[j+1] < AS_lower[j]:
        lower_change.append(j)
upper_change.append(600)
lower_change.append(600)

theta = integrate.cumtrapz(N_m, span, initial = 0)
theta2 = integrate.cumtrapz(N_m2, span, initial = 0)
v = -(1/(E1 *Ixx))*integrate.cumtrapz(theta, span, initial = 0)
v2 = -(1/(E1 *Ixx))*integrate.cumtrapz(theta2, span, initial = 0)

for i in range(len(upper_change)-1):
    k1 = upper_change[i]

    k2 = upper_change[i+1]
    v[k1+1:k2+1] =  v[k1+1:k2+1] - (v[k1+1]-v[k1+0])
    v2[k1+1:k2+1] =  v2[k1+1:k2+1] - (v2[k1+1]-v2[k1+0])
for i in range(len(lower_change)-1):
    
    k1 = lower_change[i]
    k2 = lower_change[i+1]
    v[k1+1:k2+1] =  v[k1+1:k2+1] - (v[k1+1]-v[k1+0])
    v2[k1+1:k2+1] =  v2[k1+1:k2+1] - (v2[k1+1]-v2[k1+0])
print("Wing tip deflection =", v[-1])
plt.plot(span, v)
plt.plot(span, v2)
plt.axis('equal')
plt.show()

#Angle of twist

dtheta1 = IntM/(4*Am**2*G1)*(2*a/t_flange+2*b/t_spar)
totaltheta1 = (integrate.cumtrapz(dtheta1, span, initial = 0)) *180/np.pi # Degrees
print("Twist wing tip =", totaltheta1[-1])
plt.plot(span, IntM)
plt.show()
print(IntM[506])
#Torsional Rigidity
JG = G1*(2*a/t_flange+2*b/t_spar)/Am
print(JG)

#Weight

Topplate = t_flange*(1.5525+0.4*1.5525)*10.8695/2
Sideplate = t_spar*(0.36+0.4*0.36)*10.8695/2
Stringers = integrate.trapz(AS, span)
Stringers_lower = integrate.trapz(AS_lower, span)
Stiffener_volume = 0.000568*b[0]*stiffeners
TotalArea = 0
for x in range(np.size(Rib_locations)-1):
    Area = Am[spanlocations[x]]*0.001
    TotalArea = TotalArea + Area

Weight = dens*(Topplate*2+Sideplate*2+Stringers+Stringers_lower+TotalArea+Stiffener_volume)
print("Weight =", Weight)

