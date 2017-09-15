
# coding: utf-8

# In[1]:

import numpy as np
import matplotlib.pyplot as plt
get_ipython().magic('matplotlib inline')
from scipy import stats

plt.style.use("seaborn-whitegrid")

# lineare Regression: Fehler
def lin_sigma(x, y, m, b):
    n = len(x)
    sum1, sum2, sum3 = 0, 0, 0
    for i in range(n):
        sum1 += (y[i]-b-m*x[i])**2
        sum2 += x[i]**2
        sum3 += x[i]
    sigma_m_b = [np.sqrt(n*sum1/((n-2)*(n*sum2-sum3**2))),np.sqrt((sum2*sum1/((n-2)*(n*sum2-sum3**2))))]
    return sigma_m_b

def gew_mittel(x, s_x):    
    sum1, sum2 = 0, 0
    for i in range(len(x)):
        sum1 += x[i]/s_x[i]**2
        sum2 += 1/s_x[i]**2
    x = [sum1/sum2, np.sqrt(1/sum2)]
    return x
    

# Daten aus txt-Datei einlesen
daten = open("Magnetdaten.txt", mode="r")
inhalt = daten.readlines()
daten.close()
inhalt = [x.strip() for x in inhalt]

# Spulendurchmesser (lang, kurz, Helmholtz, Induktion) in m
D = [float(x) for x in inhalt[2].split(",")]
# Windungszahl
N = [float(x) for x in inhalt[4].split(",")]
# Länge (lang, kurz) in m
L = [float(x) for x in inhalt[6].split(",")]
# Drahtdurchmesser (lang, kurz)
DD = [float(x)*10**(-3) for x in inhalt[8].split(",")]
# Innenwiderstand des Stromintegrators in Ohm
Rs = 1000#0    # LP sagt auch 1 kOhm, ich habe es gerade geändert.
f_Rs = 30
R1 = 10000
R2 = 3.3

# Eichungszeiten in s
t = [float(x)*10**(-3) for x in inhalt[12].split(",")]
# Zahl am Integrator
zt = [float(x)*(-1) for x in inhalt[14].split(",")]
f_zt = float(inhalt[16])

# Messung per Induktionsspule
# Abstand zur Spulenmitte in m
ab1 = [float(x)*0.01 for x in inhalt[20].split(",")]
f_ab1 = float(inhalt[22])*0.01
# Zahl
z1 = [float(x)*(-1) for x in inhalt[24].split(",")]

# Messung per Hallsonde
# Abstand zur Spulenmitte in m
ab2 = [i*0.01 for i in range(int(inhalt[28])+1)]
ab2h = [i*0.01 for i in range(int(inhalt[30])+1)]
# Magnetfeld in Tesla (1 G ~ 10**(-4) T)
B = [[-1*float(x)*10**(-4) for x in inhalt[i].split(",")] for i in [32, 33, 34]]
f_B = 0.2*10**(-4)
Boff = -1*float(inhalt[38])*10**(-4)
B = [[B[i][j]-Boff for j in range(len(B[i]))]for i in range(3)]


# In[5]:

# Eichung
def k(rs, r1, r2, u, dt, z, f_rs, f_z):
    k = r2/(rs*(r2+r1)+r1*r2)*u*dt/z
    sigma = np.sqrt((f_rs*rs*(r1+r2)*u*dt*r2/(z*(rs*(r1+r2)*r1*r2)))**2 + (f_z*r2*u*dt/(z**2*(rs*(r1+r2)+r1*r2)))**2)
    #sigma = np.sqrt(sigma)
    return [k, sigma]

a = [k(Rs, R1, R2, 2, t[i], zt[i], f_Rs, f_zt) for i in range(len(t))]

print("Gewichtetes Mittel für die Propkonstante ", gew_mittel([a[i][0] for i in range(len(t))],[a[i][1] for i in range(len(t))]))

# Neuer Versuch der Eichung
# Lineare Regression von Q = k*Z, wir haben 
def Q_theorie(rs, r1, r2, u, dt):
    q = (r2+r2**2/rs)/(rs*(r2+r1)+r1*r2)*u*dt 
    return q


Q1 = [Q_theorie(Rs, R1, R2, 2, t[i]) for i in range(len(t))]
print("Erwartete Ladungen am Kondensator: ",Q1)

fig = plt.figure()
ax = plt.axes()
x = np.linspace(-900,0, 10000)
m, b, r, p, s = stats.linregress(zt,Q1)

f_m = lin_sigma(zt, Q1, m, b)[0]
f_b = lin_sigma(zt, Q1, m, b)[1]
print("Ergebnisse der Regression:")
print("m = ", m, "\pm", f_m)
print("b = ", b, "\pm", f_b)

plt.errorbar(zt, Q1, xerr = f_zt, fmt='.', color = "blue", label="Gemessene Werte")
plt.plot(x, m*x+b, color="green", label=r"Lineare Regression")

plt.xlabel(r"$Z$ in Skt.")
plt.ylabel(r"$Q$ in C")
#plt.axis([390000, 1610000,470,1000]) # ([x_Achsenstart, x-Achsenstopp, y-Achsenstart, y-Achsenstopp])
ax.legend(loc="lower left",frameon=True)
#plt.savefig("k.png")
plt.show()

#m *= 100
#b *= 100
#f_m *= 100
#f_b *= 100

def Q(m, b, z, f_m, f_b):
    q = (m*z+b)
    f_q = np.sqrt((f_m*z)**2+(f_b)**2)
    return [q, f_q]


# In[6]:

# Messung mit der Induktionsspule
q1 = [Q(m, b, z1[i], f_m, f_b) for i in range(len(z1))]

def B_mess(Rs, Q, n, d, f_Q):
    A = 1/4*np.pi*d**2
    B = [Rs*Q/(n*A)] # Wert
    B.append(f_Q*Rs/(n*A)) #Fehler
    return B

B1_lang = [B_mess(Rs, q1[i][0], N[0], D[0], q1[i][1]) for i in range(len(q1))]

def b_theorie(z, n, I, l, R):
    mu = 4*np.pi*10**(-7)
    b = mu*n*I/(2*l)*((z+l/2)/np.sqrt(R**2 + (l/2)**2)-(z-l/2)/np.sqrt(R**2+(z+l/2)**2))
    return b

#b(x, N[0], 0.5, L[0], D[0]/2)
mu = 4*np.pi*10**(-7)

fig = plt.figure()
ax = plt.axes()
x = np.linspace(0,0.3, 60)
plt.plot(x, mu*N[0]*0.5/(2*L[0])*((x+L[0]/2)/np.sqrt((D[0]/2)**2 + (x+L[0]/2)**2)-(x-L[0]/2)/np.sqrt((D[0]/2)**2+(x-L[0]/2)**2)),color="green", label=r"Theorie: B($d$)")
plt.errorbar(ab1, [B1_lang[i][0] for i in range(len(B1_lang))], xerr = f_ab1, yerr = [B1_lang[i][1] for i in range(len(B1_lang))], fmt='.', color = "blue", label="Induktionsspule")
plt.errorbar(ab2, B[0], xerr = f_ab1, yerr = f_B, fmt='.', color = "purple", label="Hallsonde")
plt.axhline(y=mu*N[0]*0.5/(L[0]), color="red",label=r"Theorie: B(0)")
plt.ylabel(r"$B$ in T")
plt.xlabel(r"Abstand $d$ in m")
#plt.axis([390000, 1610000,470,1000]) # ([x_Achsenstart, x-Achsenstopp, y-Achsenstart, y-Achsenstopp])
ax.legend(loc="lower left",frameon=True)
#plt.savefig("messung1mal100.png")
plt.show()

#Abweichung (Hall - Induktion, Hall-Theorie, Theorie-Induktion)
print(gew_mittel([B[0][2*i]-B1_lang[i][0] for i in range(9)],[f_B+B1_lang[i][1] for i in range(len(B1_lang))]))
print(gew_mittel([B[0][i]-b_theorie(i*0.01,N[0],0.5,L[0],D[0]/2) for i in range(18)],[f_B for i in range(18)]))
print(gew_mittel([b_theorie(i*0.02,N[0],0.5,L[0],D[0]/2)-B1_lang[i][0] for i in range(9)],[B1_lang[i][1] for i in range(len(B1_lang))]))


# In[4]:

fig = plt.figure()
ax = plt.axes()
x = np.linspace(0,0.18,30)
plt.errorbar(ab2, B[1], xerr = f_ab1, yerr = f_B, fmt='.', color = "red", label="Messwerte: kurze Spule")
plt.plot(x, mu*N[1]*0.5/(2*L[1])*((x+L[1]/2)/np.sqrt((D[1]/2)**2 + (x+L[1]/2)**2)-(x-L[1]/2)/np.sqrt((D[1]/2)**2+(x-L[1]/2)**2)),color="orange", label=r"Theorie: B($d$)")
plt.axhline(y=mu*N[1]*0.5/(L[1]), color = 'purple', label=r"Theorie: B(0)")
plt.ylabel(r"$B$ in T")
plt.xlabel(r"Abstand $d$ in m")
#plt.axis([390000, 1610000,470,1000]) # ([x_Achsenstart, x-Achsenstopp, y-Achsenstart, y-Achsenstopp])
ax.legend(loc="lower left",frameon=True)
#plt.savefig("kurzespule.png")
plt.show()

fig = plt.figure()
ax = plt.axes()
x = np.linspace(0, 0.09, 30)
plt.errorbar(ab2h, B[2], xerr = f_ab1, yerr = f_B, fmt='.', color = "green", label="Messwerte: Helmholtzspulenpaar")
plt.plot(x, mu*N[2]*0.25*(D[2])**2*(((x+D[2]/2)**2+(D[2])**2)**(-3/2)+((x-D[2]/2)**2+(D[2])**2)**(-3/2)),label="Theorie: B($d$)")
plt.axhline(y=(4/5)**(3/2)*N[2]*mu*0.25/(D[2]/2), color = 'purple', label=r"Theorie: $\approx$B($d$)")
plt.ylabel(r"$B$ in T")
plt.xlabel(r"Abstand $d$ in m")
#plt.axis([390000, 1610000,470,1000]) # ([x_Achsenstart, x-Achsenstopp, y-Achsenstart, y-Achsenstopp])
ax.legend(loc="lower left",frameon=True)
#plt.savefig("helmholtzspule.png")
plt.show()


# In[5]:

# Gruppenplot
fig = plt.figure()
ax = plt.axes()
#x = np.linspace(0, 0.30, 30)
plt.errorbar(ab2h, B[2], xerr = f_ab1, yerr = f_B, fmt='.', color = "blue", label="Helmholtzspulenpaar")
plt.errorbar(ab2, B[1], xerr = f_ab1, yerr = f_B, fmt='.', color = "red", label="Kurze Spule")
plt.errorbar(ab2, B[0], xerr = f_ab1, yerr = f_B, fmt='.', color = "purple", label="Lange Spule")
#plt.plot(x, mu*N[2]*0.25*(D[2])**2*(((x+D[2]/2)**2+(D[2])**2)**(-3/2)+((x-D[2]/2)**2+(D[2])**2)**(-3/2)),label="Theorie")
plt.axhline(y=(4/5)**(3/2)*N[2]*mu*0.25/(D[2]/2), color = 'green', label=r"Theorie:$\approx$ B($d$)")
plt.ylabel(r"$B$ in T")
plt.xlabel(r"Abstand $d$ in m")
plt.axis([-0.025,0.3,0,0.0025]) # ([x_Achsenstart, x-Achsenstopp, y-Achsenstart, y-Achsenstopp])
ax.legend(loc="upper right",frameon=True)
#plt.savefig("gruppenplot.png")
plt.show()
print(gew_mittel([(4/5)**(3/2)*N[2]*mu*0.25/(D[2]/2)-B[2][i] for i in range(len(B[2]))],[f_B for i in range(len(B[2]))]))
#print("sigma_B1 =", 0.05*mu*N[1]*0.5/(2*L[1])*((x+L[1]/2)/np.sqrt((D[1]/2)**2 + (x+L[1]/2)**2)-(x-L[1]/2)/np.sqrt((D[1]/2)**2+(x-L[1]/2)**2)))


# In[6]:

# Bestimmung von mu_0
def h_theorie(z, n, I, l, R, f_I,f_z):
#    mu = 4*np.pi*10**(-7)
    h = []
    h.append(n*I/(2*l)*((z+l/2)/np.sqrt(R**2 + (l/2)**2)-(z-l/2)/np.sqrt(R**2+(z+l/2)**2)))
    a1=(n*f_I/(2*l)*((z+l/2)/np.sqrt(R**2 + (l/2)**2)-(z-l/2)/np.sqrt(R**2+(z+l/2)**2)))**2
    a2= f_z*n*I/(2*l)*(1/(np.sqrt(R**2+(z+l/2)**2))-1/(np.sqrt(R**2+(z-l/2)**2))-(z+l/2)**2/(R**2+(z+l/2)**2)**(5/2)+(z-l/2)**2/(R**2+(z-l/2)**2)**(5/2))
    h.append(np.sqrt(a1+a2**2))
    return h

def h_helmholtz(z,f_z,f_i):
    h = []
    h.append(N[2]*0.25*(D[2])**2*(((z+D[2]/2)**2+(D[2])**2)**(-3/2)+((z-D[2]/2)**2+(D[2])**2)**(-3/2)))
    h.append(np.sqrt((f_i*N[2]*(D[2])**2*(((z+D[2]/2)**2+(D[2])**2)**(-3/2)+((z-D[2]/2)**2+(D[2])**2)**(-3/2)))**2
                    +(f_z*0.25*N[2]*(D[2])**2*3*((z+D[2]/2)*((z+D[2]/2)**2+D[2]**2)**(-5/2)+(z-D[2]/2)*((z-D[2]/2)**2+D[2]**2)**(-5/2))**2)))
    return h
             
def mu(b,h,f_b,f_h):
    mu = []
    mu.append(b/h)
    mu.append(np.sqrt((f_b/h)**2+(f_h*b/h**2)**2))
    return mu

mu_0 = [[mu(B1_lang[i][0], h_theorie(ab1[i],N[0],0.5,L[0],D[0]/2,0.05,f_ab1)[0],
            B1_lang[i][1], h_theorie(ab1[i],N[0],0.5,L[0],D[0]/2,0.05,f_ab1)[1]) for i in range(len(B1_lang))],
        [mu(B[0][i],h_theorie(ab2[i],N[0],0.5,L[0],D[0]/2,0.05,f_ab1)[0],f_B,h_theorie(ab2[i],N[0],0.5,L[0],D[0]/2,0.05,f_ab1)[1]) for i in range(len(B[0]))],
        [mu(B[1][i],h_theorie(ab2[i],N[1],0.5,L[1],D[1]/2,0.05,f_ab1)[0],f_B,h_theorie(ab2[i],N[1],0.5,L[1],D[1]/2,0.05,f_ab1)[0])for i in range(len(B[1]))],
       [mu(B[2][i],h_helmholtz(ab2h[i],f_ab1,0.05)[0],f_B,h_helmholtz(ab2h[i],f_ab1,0.05)[1])for i in range(len(B[2]))]]
a =[gew_mittel([mu_0[i][j][0] for j in range(len(mu_0[i]))],[mu_0[i][j][1]for j in range(len(mu_0[i]))]) for i in range(4)]
print(a)
print(gew_mittel([a[i][0] for i in range(4)],[a[i][1] for i in range(4)]))


# In[16]:

# normierter Plot
f_abn = [f_ab1/L[0],f_ab1/L[1], f_ab1/D[0]]
f_Bn = []
f_Bn.append([np.sqrt((f_B/max(B[0]))**+(f_B*B[0][i]/max(B[0])**2)**2) for i in range(len(B[0]))])
f_Bn.append([np.sqrt((f_B/max(B[1]))**+(f_B*B[1][i]/max(B[1])**2)**2) for i in range(len(B[1]))])
f_Bn.append([np.sqrt((f_B/max(B[2]))**+(f_B*B[2][i]/max(B[2])**2)**2) for i in range(len(B[2]))])
fig = plt.figure()
ax = plt.axes()
#x = np.linspace(0, 0.30, 30)
plt.errorbar([ab2h[i]/D[2] for i in range(len(ab2h))], [B[2][i]/max(B[2]) for i in range(len(B[2]))], xerr = f_abn[2], fmt='.', color = "blue", label="Helmholtzspulenpaar")
plt.errorbar([ab2[i]/L[1] for i in range(len(ab2))], [B[1][i]/max(B[1]) for i in range(len(B[1]))], xerr = f_abn[1], fmt='.', color = "red", label="Kurze Spule")
plt.errorbar([ab2[i]/L[0] for i in range(len(ab2))], [B[0][i]/max(B[0]) for i in range(len(B[0]))], xerr = f_abn[0], fmt='.', color = "purple", label="Lange Spule")
#plt.plot(x, mu*N[2]*0.25*(D[2])**2*(((x+D[2]/2)**2+(D[2])**2)**(-3/2)+((x-D[2]/2)**2+(D[2])**2)**(-3/2)),label="Theorie")
#plt.axhline(y=(4/5)**(3/2)*N[2]*mu*0.25/(D[2]/2), color = 'green', label=r"Theorie:$\approx$ B($d$)")
plt.ylabel(r"$\frac{B}{B_{\mathrm{max}}}$")
plt.xlabel(r"$\frac{d}{l}$")
plt.axis([-0.1,0.9,0.3,1.1]) # ([x_Achsenstart, x-Achsenstopp, y-Achsenstart, y-Achsenstopp])
ax.legend(loc="lower left",frameon=True)
#plt.savefig("normiert.png")
plt.show()
#print(gew_mittel([(4/5)**(3/2)*N[2]*mu*0.25/(D[2]/2)-B[2][i] for i in range(len(B[2]))],[f_B for i in range(len(B[2]))]))
#print("sigma_B1 =", 0.05*mu*N[1]*0.5/(2*L[1])*((x+L[1]/2)/np.sqrt((D[1]/2)**2 + (x+L[1]/2)**2)-(x-L[1]/2)/np.sqrt((D[1]/2)**2+(x-L[1]/2)**2)))


# In[ ]:




# In[ ]:



