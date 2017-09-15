
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
daten = open("Transformatordaten.txt", mode="r")
inhalt = daten.readlines()
daten.close()
inhalt = [x.strip() for x in inhalt]

# Fehler
f_U = 0.1
f_I = 0.05

# Primärspannung in V
U1 = [float(x) for x in inhalt[1].split(",")]
# Primärstrom in A
I1 = [float(x) for x in inhalt[3].split(",")]

# Sekundärspannung Teil 1 in V
U21 = [float(x) for x in inhalt[8].split(",")]

# Sekundärspannung Teil 2 in V
U22 = [float(x) for x in inhalt[13].split(",")]

# Strom I2 in A
I22 = [float(x) for x in inhalt[16].split(",")]
# Strom I1 in A
I21 = [float(x) for x in inhalt[18].split(",")]
# Strom I_ges in A
Iges = [float(x) for x in inhalt[20].split(",")]

# Primärstrom in A
I_prim = float(inhalt[23])
# Primärspannung in V
U_prim = float(inhalt[25])

# Längenverhältnisse
a = [float(inhalt[29+i].split(",")[0]) for i in range(6)]
x0 = [float(inhalt[29+i].split(",")[1]) for i in range(6)]
b = [float(inhalt[29+i].split(",")[2]) for i in range(6)]
y0 = [float(inhalt[29+i].split(",")[3]) for i in range(6)]

f_abl = float(inhalt[37])


# In[2]:

# U1 gegen I1 auftragen, Verlauf diskutieren (soll linear sein)

fig = plt.figure()
ax = plt.axes()
#x = np.linspace(0, 0.09, 30)
plt.errorbar(I1, U1, xerr = f_I, yerr = f_U, fmt='.', color = "green", label="Messwerte")
#plt.plot(x, mu*N[2]*0.25*(D[2])**2*(((x+D[2]/2)**2+(D[2])**2)**(-3/2)+((x-D[2]/2)**2+(D[2])**2)**(-3/2)),label="Theorie: B($d$)")
#plt.axhline(y=(4/5)**(3/2)*N[2]*mu*0.25/(D[2]/2), color = 'purple', label=r"Theorie: $\approx$B($d$)")
plt.ylabel(r"$U_1$ in V")
plt.xlabel(r"$I_1$ in A")
plt.axis([0, 2.5, 0, 16]) # ([x_Achsenstart, x-Achsenstopp, y-Achsenstart, y-Achsenstopp])
ax.legend(loc="upper left",frameon=True)
#plt.savefig("UgegenI.png")
plt.show()


# U2 gegen U1 auftragen, lineare Regression -> Steigung m = U2/U1 = 1/u, u soll ausgerechnet werden

m1, b1, r1, p1, s1 = stats.linregress(U1,U21)

f_m1 = lin_sigma(U1, U21, m1, b1)[0]
f_b1 = lin_sigma(U1, U21, m1, b1)[1]
print("Ergebnisse der Regression:")
print("m = ", m1, "\pm", f_m1)
print("b = ", b1, "\pm", f_b1)
print("u = $", 1/m1, "\pm", f_m1/m1**2, "$")

fig = plt.figure()
ax = plt.axes()
x = np.linspace(0, 16, 30)
plt.errorbar(U1, U21, xerr = f_U, yerr = f_U, fmt='.', color = "green", label="Messwerte")
plt.plot(x, m1*x+b1,color = "blue", label="Lineare Regression")
#plt.axhline(y=(4/5)**(3/2)*N[2]*mu*0.25/(D[2]/2), color = 'purple', label=r"Theorie: $\approx$B($d$)")
plt.ylabel(r"$U_2$ in V")
plt.xlabel(r"$U_1$ in V")
#plt.axis([0, 8, 0, 16]) # ([x_Achsenstart, x-Achsenstopp, y-Achsenstart, y-Achsenstopp])
ax.legend(loc="upper left",frameon=True)
#plt.savefig("U2gegenU1.png")
plt.show()


# U1 gegen U2 auftragen, lineare Regression -> Steigung m = U1/U2 = 1/u, wieder u ausrechnen

m2, b2, r2, p2, s2 = stats.linregress(U22,U1)

f_m2 = lin_sigma(U22, U1, m2, b2)[0]
f_b2 = lin_sigma(U22, U1, m2, b2)[1]
print("Ergebnisse der Regression:")
print("m = ", m2, "\pm", f_m2)
print("b = ", b2, "\pm", f_b2)
print("u = ", 1/m2, "\pm", f_m2/m2**2)

fig = plt.figure()
ax = plt.axes()
x = np.linspace(0, 16, 30)
plt.errorbar(U22, U1, xerr = f_U, yerr = f_U, fmt='.', color = "green", label="Messwerte")
plt.plot(x, m2*x+b2,color = "blue", label="Lineare Regression")
#plt.axhline(y=(4/5)**(3/2)*N[2]*mu*0.25/(D[2]/2), color = 'purple', label=r"Theorie: $\approx$B($d$)")
plt.ylabel(r"$U_1$ in V")
plt.xlabel(r"$U_2$ in V")
plt.axis([0, 16, 0, 16]) # ([x_Achsenstart, x-Achsenstopp, y-Achsenstart, y-Achsenstopp])
ax.legend(loc="upper left",frameon=True)
#plt.savefig("U1gegenU2.png")
plt.show()

# Gewichteter Mittelwert für u

# Ein Stückchen der "Auswertung" von Tim abschreiben

# Phasenverschiebung berechnen

# Phasenverschiebung berechnen

# Phasenverschiebung und Theoriekurve in einen Plot malen

# Komische Handyladegerätaufgabe rechnen


# In[3]:

# Phasenverschiebung aus phi = 2*arccos(I_ges/(2*I_1))
phi, f_phi, f_phi2 = [],[],[]
for i in range(len(I22)-1):
    phi.append(2*np.arccos(Iges[i]/(2*I21[i])))
    f_phi.append(2*f_I*np.sqrt(1/(4*I21[i]**2)+Iges[i]**2/(4*I21[i]**4))*(1/(1-Iges[i]**2/(4*I21[i]**2))))
    f_phi2.append(f_I*np.sqrt((I21[i]**2+Iges[i]**2)/(I21[i]**4-1/4*I21[i]**2*Iges[i]**2))) # aus Luis Protokoll
    
phi.append(0) # Iges/(2I1) ist ungefähr 1
f_phi2.append(f_I*np.sqrt(abs((I21[5]**2+Iges[5]**2)/(I21[5]**4-1/4*I21[5]**2*Iges[5]**2))))

print("Phasenverschiebung φ zwischen  der  Spannung  und  dem  Gesamtstrom.")

for i in range(len(I22)):
    print("$", I22[i], "$&$", phi[i], "\pm",f_phi2[i], "$\\\\\hline")


# In[4]:

# Phasenverschiebung aus Lissajousfiguren

theta = []
f_theta = []
for i in range(6):
    theta.append([np.arcsin(y0[i]/b[i]),np.arcsin(x0[i]/a[i])])
    f_theta.append([f_abl*np.sqrt(1/(b[i]**2-y0[i]**2)+y0[i]**2/(b[i]**4-y0[i]**2*b[i]**2)),f_abl*np.sqrt(1/(a[i]**2-x0[i]**2)+x0[i]**2/(a[i]**4-x0[i]**2*a[i]**2))])

print("Phasenverschiebung φ zwischen  der  Spannung  und  dem  Gesamtstrom  aus den Lissajousfiguren")
for i in range(len(I22)):
    print("$", I22[i], "$&$", theta[i][0], "\pm",f_theta[i][0],"$&$", theta[i][1], "\pm",f_theta[i][1], "$\\\\\hline")
    
print("\n\n")
    
# Gewichtete Mittelwerte
winkel_gesamt, f_winkel_gesamt = [],[]
for i in range(6):
    winkel_gesamt.append([phi[i], theta[i][0], theta[i][1]])
    f_winkel_gesamt.append([f_phi2[i], f_theta[i][0], f_theta[i][1]])
mittelphi = [gew_mittel(winkel_gesamt[i], f_winkel_gesamt[i]) for i in range(6)]

print("Gewichtete  Mittelwerte  der  Phasenverschiebung φ zwischen  der  Spannung und dem Gesamtstrom.")
for i in range(len(I22)):
    print("$", I22[i], "$&$", mittelphi[i][0], "\pm",mittelphi[i][1], "$\\\\\hline")


# In[5]:

fig = plt.figure()
ax = plt.axes()
x = np.linspace(-0.2, 1.6, 30)
plt.errorbar(I22, phi, xerr = f_I, yerr = f_phi2, fmt='.', color = "red", label=r"$\phi$ aus $I_{\mathrm{s}}$ und $I_1$")
plt.errorbar(I22, [theta[i][0] for i in range(6)], xerr = f_I, yerr =  [f_theta[i][0] for i in range(6)], fmt='.', color = "orange", label=r"$\phi$ aus $y_0$ und $b$")
plt.errorbar(I22,  [theta[i][1] for i in range(6)], xerr = f_I, yerr =  [f_theta[i][1] for i in range(6)], fmt='.', color = "purple", label=r"$\phi$ aus $x_0$ und $a$")
plt.errorbar(I22,  [mittelphi[i][0] for i in range(6)], xerr = f_I, yerr =  [mittelphi[i][1] for i in range(6)], fmt='.', color = "pink", label=r"Gew. Mittel von $\phi$")
plt.plot(x, np.arctan(0.32*np.sin(0.51)/(m2*x + 0.32*np.cos(0.51))),label="Theorie")
#plt.axhline(y=(4/5)**(3/2)*N[2]*mu*0.25/(D[2]/2), color = 'purple', label=r"Theorie: $\approx$B($d$)")
plt.ylabel(r"$\phi$ in rad.")
plt.xlabel(r"$I_2$ in A")
plt.axis([-0.2, 1.6, 0, 2]) # ([x_Achsenstart, x-Achsenstopp, y-Achsenstart, y-Achsenstopp])
ax.legend(loc="upper right",frameon=True)
#plt.savefig("phitheorie.png")
plt.show()


# In[6]:

# Wirk- und Verlustleistung bei 1.5 A und 200 V.
print(1.5*200*np.sin(mittelphi[5][0]), "\pm", mittelphi[5][1]*1.5*200*np.cos(mittelphi[5][1]))
print(1.5*200*np.cos(mittelphi[5][0]), "\pm", mittelphi[5][1]*1.5*200*np.sin(mittelphi[5][0]))
print(365*24*I21[0]*200*np.sin(mittelphi[5][0])*0.25*0.001, "\pm", 365*24*f_I*200*np.sin(mittelphi[5][0])*0.25*0.001)
print(365*24*I21[0]*200*np.cos(mittelphi[5][0])*0.25*0.001, "\pm", 365*24*f_I*200*np.cos(mittelphi[5][0])*0.25*0.001) # Nachdem Merten meinte, sin und cos wären vertauscht
#print(I21[0], "\pm", f_I, "V")
print("\n", mittelphi[5][0], "\pm", mittelphi[5][1])


# In[ ]:




# In[ ]:




# In[ ]:



