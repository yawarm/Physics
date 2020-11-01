import numpy as np

#Parametere som går igjen (felles parametere)
#Alle enheter er SI-enheter
g = 9.81 #tyngdeakselerasjon
dt = 0.0001 #tidsintervall (steglengde?) - ikke noe bra navn... 
r_bane = 0.5 #bane radius
v_start = 0.0001 #Start hastighet (hold denne lav, objekter sluppet uten start fart)

#Justerbare parametere
#Alle enheter er SI-enheter
r_objekt = 0.0285 #radius på objektet som modelleres (i oppgave 1 settes den til null)
c = 2/5 #Objektets treghetsmoment konstant (vet ikke formell navn?)
µ_s = 1 #Statisk friksjonskoeffisient
µ_k = 0.6 #Kinetisk friksjonskoeffisient
theta_s = 0 #Start vinkel relativ 0
m = 1 #Objektets masse

#Analytiske uttrykk som går igjen
#Alle enheter er SI-enheter
Theta = [] #Vinkel langs banen målt relativt null - en liste
Theta.append(theta_s) #Legger til startvinkel

V = [] #Fart
V.append(v_start) #Legger til starthastighet

V_x = [] #Fartskomponent i x-retning
V_x.append(v_start*np.cos(theta_s)) #Start fart x-retning

V_y = [] #Fartkomponent i y-retning
V_y.append(-v_start*np.sin(theta_s)) #Start fart y-retning

W = [] #Vinkelhastighet om origo
W.append(v_start/(r_objekt+r_bane)) #start vinkelhastighet

X = [] #x-posisjon
X.append((r_objekt + r_bane)*np.sin(theta_s)) # start x-posisjon

Y = [] #y-posisjon
Y.append((r_objekt + r_bane)*np.cos(theta_s)) # start y-posisjon

N = [] #Normalkraft
N.append(m*g*np.cos(theta_s)) #Normalkraft ved start

A_tan_1 = [] #Tangentsiell akselerasjon
A_tan_1.append(g*np.sin(theta_s))

A_tan_2 = [] #Tangentsiell akselerasjon oppgave 2
A_tan_2.append((g*np.sin(theta_s))/(c+1)) #start tangentsiell akselerasjon oppgave 1

Tid = [] #tid
Tid.append(dt)

F_max = [] # Maksimal friksjonskraft
F_max.append(µ_s*m*g*np.cos(theta_s)) #start maksimal friksjonskraft

Alfa = [] #Vinkelakselerasjon
Alfa.append((g*np.sin(theta_s))/(r_objekt+r_bane))

F_rull = [] #Friksjonskraften under ren rulling
F_rull.append((c*m*g*np.sin(theta_s))/(c+1))

F_slur = [] #Friksjonskraft under sluring

A_tan_slur = [] #Tangtsiell akselerasjon under sluring

N_slur = [] #Normalkraft under sluring

Potensiell = [] #Potensiell energi
Potensiell.append(m*g*((r_objekt + r_bane)*np.cos(theta_s)))

K_trans = [] #Kinetisk energi grunnet translasjon
K_trans.append(m*0.5*v_start**2)

K_rot = [] #Kinetisk energi grunnet rotasjon
K_rot.append(0.5*c*m*(v_start)**2)

Mekanisk_E = [] #Total mekanisk energi i systemet
Mekanisk_E.append(Potensiell[0] + K_trans[0] + K_rot[0])

F_work = [] #Arbeid utført av friksjonskraften
F_work.append(0)

Total_E = [] #Total energi i systemet - dette er en konstant
Total_E.append(Mekanisk_E[0] + F_work[0])

V_rel = [] #Farten massesenter har relativ underlag under sluring
#Hensiktsmessig å definere denne størrelsen, siden den viser seg å være nyttig for å regne friksjonsarbeid










#Oppgave 1 (Inkludert energibetraktning)

i = 0 #Indeks for loop
while N[i] > 0:
    v = V[i] + Tid[0]*A_tan_1[i] #ny fart
    V.append(v)
    w = V[i+1]/r_bane #ny vinkelhastighet
    W.append(w)
    theta = Theta[i] + W[i+1]*Tid[0] #ny vinkel
    Theta.append(theta)
    n = m*g*(3*np.cos(Theta[i+1]) - 2*np.cos(theta_s)) #ny normalkraft
    N.append(n)
    a_tan = g*np.sin(Theta[i+1]) #ny tangetsiell akselerasjon
    A_tan_1.append(a_tan)
    v_x = V[i]*np.cos(Theta[i+1]) #ny fart x-komponent
    V_x.append(v_x)
    v_y = -V[i]*np.sin(Theta[i+1])#ny fart y-komponent
    V_y.append(v_y)
    x = (r_objekt + r_bane)*(np.sin(Theta[i+1])) #ny x-posisjon
    X.append(x)
    y = (r_objekt + r_bane)*(np.cos(Theta[i+1])) #ny-posisjon
    Y.append(y)
    potensiell = m*g*Y[i+1] #ny potensiell energi 
    Potensiell.append(potensiell)
    k_trans = 0.5*m*(V[i+1])**2 #ny kinetisk energi translasjon
    K_trans.append(k_trans)
    mekanisk_e = potensiell + k_trans #ny mekanisk energi
    Mekanisk_E.append(mekanisk_e)
    tid = dt*i
    Tid.append(tid)
    i += 1





#Oppgave 2 (inkludert energibetraktning)

if len(Theta) > 1:
    del Theta[1:], V[1:], V_x[1:], V_y[1:], W[1:], N[1:], X[1:], Y[1:], Tid[1:], A_tan_1[1::]

#Koden over kjøres for å fjerne eksisterende verdier i listene. 

i = 0 #Indeks for loop
while N[i] > 0:
    v = V[i] + Tid[0]*A_tan_2[i] #ny fart
    V.append(v)
    w = V[i+1]/(r_bane+r_objekt) #ny vinkelhastighet
    W.append(w)
    theta = Theta[i] + W[i+1]*Tid[0] #ny vinkel
    Theta.append(theta)
    n = m*g*np.cos(Theta[i+1]) + (-2*m*g*(np.cos(theta_s) - np.cos(Theta[i+1])))/(1+c) #ny normalkraft
    N.append(n)
    f_max = µ_s*N[i+1] #ny max friksjonskraft
    F_max.append(f_max)
    a_tan = (g*np.sin(Theta[i+1]))/(c+1)#ny tangetsiell akselerasjon
    A_tan_2.append(a_tan)
    alfa = A_tan_2[i+1]/(r_bane+r_objekt) #Ny vinkelakselerasjon
    Alfa.append(alfa)
    v_x = V[i]*np.cos(Theta[i+1]) #ny fart x-komponent
    V_x.append(v_x)
    v_y = -V[i]*np.sin(Theta[i+1])#ny fart y-komponent
    V_y.append(v_y)
    x = (r_objekt + r_bane)*(np.sin(Theta[i+1])) #ny x-posisjon
    X.append(x)
    y = (r_objekt + r_bane)*(np.cos(Theta[i+1])) #ny-posisjon
    Y.append(y)
    f_rull = (m*g*c*np.sin(Theta[i+1]))/(c + 1)
    F_rull.append(f_rull)
    potensiell = m*g*Y[i+1] #ny potensiell energi 
    Potensiell.append(potensiell)
    k_trans = 0.5*m*(V[i+1])**2 #ny kinetisk energi translasjon
    K_trans.append(k_trans)
    k_rot = 0.5*c*m*(V[i+1])**2 #ny kinetisk energi rotasjon
    K_rot.append(k_rot)
    mekanisk_e = potensiell + k_trans + k_rot #ny mekanisk energi
    Mekanisk_E.append(mekanisk_e)
    total_e = mekanisk_e
    Total_E.append(mekanisk_e)
    tid = dt*i #Ny tid
    Tid.append(tid)
    i += 1





#oppgave 3 (inkludert energibetraktninger)

if len(Theta) > 1:
    del Theta[1:], V[1:], V_x[1:], V_y[1:], W[1:], N[1:], X[1:], Y[1:], Tid[1:], A_tan_1[1::]

#Koden over kjøres for å fjerne eksisterende verdier i listene. 

i = 0 #Indeks for loop
while N[-1] > 0:
    v = V[i] + Tid[0]*A_tan_2[i] #ny fart
    V.append(v)
    w = V[i+1]/(r_bane+r_objekt) #ny vinkelhastighet
    W.append(w)
    theta = Theta[i] + W[i+1]*Tid[0] #ny vinkel
    Theta.append(theta)
    n = m*g*np.cos(Theta[i+1]) + (-2*m*g*(np.cos(theta_s) - np.cos(Theta[i+1])))/(1+c) #ny normalkraft
    N.append(n)
    f_max = µ_s*N[i+1] #ny max friksjonskraft
    F_max.append(f_max)
    a_tan = (g*np.sin(Theta[i+1]))/(c+1)#ny tangetsiell akselerasjon
    A_tan_2.append(a_tan)
    alfa = A_tan_2[i+1]/(r_bane+r_objekt) #Ny vinkelakselerasjon
    Alfa.append(alfa)
    v_x = V[i+1]*np.cos(Theta[i+1]) #ny fart x-komponent
    V_x.append(v_x)
    v_y = -V[i+1]*np.sin(Theta[i+1])#ny fart y-komponent
    V_y.append(v_y)
    x = (r_objekt + r_bane)*(np.sin(Theta[i+1])) #ny x-posisjon
    X.append(x)
    y = (r_objekt + r_bane)*(np.cos(Theta[i+1])) #ny-posisjon
    Y.append(y)
    f_rull = (m*g*c*np.sin(Theta[i+1]))/(c + 1)
    F_rull.append(f_rull)
    potensiell = m*g*Y[i+1] #ny potensiell energi 
    Potensiell.append(potensiell)
    k_trans = 0.5*m*(V[i+1])**2 #ny kinetisk energi translasjon
    K_trans.append(k_trans)
    k_rot = 0.5*c*m*(V[i+1])**2 #ny kinetisk energi rotasjon
    K_rot.append(k_rot)
    mekanisk_e = potensiell + k_trans + k_rot #ny mekanisk energi
    Mekanisk_E.append(mekanisk_e)
    total_e = mekanisk_e
    Total_E.append(mekanisk_e)
    if f_rull > f_max: #Betingelse for overgang fra ren rulling til sluring
        theta_kritisk = Theta[i]
        k = i #lagrer indexen vi går over fra ren rulling til sluring
        j = 0
        while N[-1] > 0:
            f_slur = µ_k*N[-1] #friksjonskraft sluring
            F_slur.append(f_slur)
            w_rull = W[i+1] + (F_slur[-1]*dt)/(c*m)
            W.append(w_rull)
            a_tan_slur = (g*np.sin(Theta[-1])) - F_slur[j] #tangentsiell akselerasjon sluring
            A_tan_slur.append(a_tan_slur)
            v = V[i] + dt*A_tan_slur[j] #ny fart
            V.append(v)
            v_rel = (V[-1]) - (W[-1]*(r_objekt))
            V_rel.append(v_rel)
            theta_slur = Theta[-1] + (v*dt)/(r_objekt+r_bane)
            Theta.append(theta_slur)
            n_slur = m*g*np.cos(Theta[-1]) - m*(v**2)/(r_bane + r_objekt)
            N.append(n_slur)
            x = (r_objekt + r_bane)*(np.sin(Theta[i+1])) #ny x-posisjon
            X.append(x)
            y = (r_objekt + r_bane)*(np.cos(Theta[i+1])) #ny-posisjon
            Y.append(y)
            f_work = F_slur[-1]*V_rel[-1]*dt #arbeid gjort av friksjonskraft
            F_work.append(f_work)
            k_trans = 0.5*m*(V[-1])**2 #Kinetisk energi translasjon
            K_trans.append(k_trans)
            k_rot = 0.5*c*m*(V[-1])**2
            K_rot.append(k_rot)
            potensiell = m*g*Y[-1]
            Potensiell.append(potensiell)
            mekanisk_e = K_trans[-1] + K_rot[-1] + Potensiell[-1]
            Mekanisk_E.append(mekanisk_e)
            total_e = Mekanisk_E[-1] + F_work[-1]
            Total_E.append(total_e)
            tid = dt*i
            Tid.append(tid)
            j += 1
            i += 1
    else:
        tid = dt*i
        Tid.append(tid)
        i += 1

