
# coding: utf-8

import numpy as np
from numpy.polynomial import polynomial as P
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
daten = open("Diaten.txt", mode="r")
inhalt = daten.readlines()
daten.close()
inhalt = [x.strip() for x in inhalt]

# Spule: Windungszahl
n = float(inhalt[1])
# Spulenstrom in A
Isp = float(inhalt[3])
# Massen der Probekörper in kg
m = [float(x)*0.001 for x in inhalt[5].split(",")]
# Magnetfeld Offset in T
Boff = -1*float(inhalt[7])
# Fehler von B in T
f_B = float(inhalt[9])
# Fehler von I in A
f_I = float(inhalt[11])
# Fehler der Waage in kg
f_W = 0.001*float(inhalt[13])
# Abstand in m
ab1 = [i*0.001 for i in range(0,int(inhalt[15])+int(inhalt[17]),int(inhalt[17]))]
# Magnetfeld in T
B1 = [-1*float(x)-Boff for x in inhalt[19].split(",")]

# Position in m ohne Magnetfeld ([b1], [b2], [b3])
poso = [[float(x)*0.001 for x in inhalt[23+i].split(",")] for i in range(3)]
# Position in m mit Magnetfeld ([b1], [b2], [b3])
posm = [[float(x)*0.001 for x in inhalt[27+i].split(",")] for i in range(3)]

# Waagenanzeige in kg ([b1], [b2], [b3]) ohne Magnetfeld
wo1 = [[float(x)*0.001 for x in inhalt[33+i].split(",")] for i in range(3)]
# Waagenanzeige in kg ([b1], [b2], [b3]) mit Magnetfeld
wm1 = [[float(x)*0.001 for x in inhalt[37+i].split(",")] for i in range(3)]
# Abstand in m
ab2 = [float(x)*0.001 for x in inhalt[43].split(",")]
# Magnetfeld in T bei (0.8, 1, 1.2, 1.4)A (Abstand in der Zeile variiert)
B2 = [[-1*float(x)-Boff for x in inhalt[45+i].split(",")] for i in range(4)]

# Strom in A
Strom = [0.8, 1, 1.2, 1.4]
rho = [16.6*10**3, 5*10**3, 9.8*10**3]

# Messung der Kraft auf den Körper B2
# Ohne Magnetfeld in g bei (0.8, 1, 1.2, 1.4) A je drei Werte
wo2 = [[float(x)*0.001 for x in inhalt[52+i].split(",")] for i in range(4)]
# Mit Magnetfeld in g bei (0.8, 1, 1.2, 1.4) A je drei Werte
wm2 = [[float(x)*0.001 for x in inhalt[57+i].split(",")] for i in range(4)]


# Plot des Ortsverlaufes von B(h)
p, cov = np.polyfit(ab1,B1, 5, cov=True)
for i in range(6):
    print("c%d = %.f \pm %.f " % (i, p[i], np.sqrt(abs(cov[i][i]))))
#p = np.polyfit(ab1, B1, 5)
for i in range(5):
    print("%d*c%d = %.f \pm %.f " % (5-i,i, (5-i)*p[i],(5-i)*np.sqrt(abs(cov[i][i]))))
p2 = [(5-i)*p[i] for i in range(5)]
f_p2 = [(5-i)*np.sqrt(abs(cov[i][i])) for i in range(5)]
f_p = [np.sqrt(abs(cov[i][i])) for i in range(6)]

def reg(x, p, f_p):
    reg = [p[0]*x**5+p[1]*x**4+p[2]*x**3+p[3]*x**2+p[4]*x+p[5]]
    reg.append(np.sqrt((f_p[0]*x**5)**2+(f_p[1]*x**4)**2+(f_p[2]*x**3)**2+(f_p[3]*x**2)**2+(f_p[4]*x)**2+f_p[5]**2))
    return reg

def produkt(x, reg, p, f_p, p2, f_p2):
    prd = [reg[0]*(p2[0]*x**4+p2[1]*x**3+p2[2]*x**2+p2[3]*x+p2[4])]
    prd.append(np.sqrt(abs((reg[1]*(p2[0]*x**4+p2[1]*x**3+p2[2]*x**2+p2[3]*x+p2[4]))**2
                      +(reg[0]*(f_p2[0]*x**4))**2
                       +(reg[0]*(f_p2[1]*x**3))**2
                      +(reg[0]*(f_p2[2]*x**2)**2
                      +(reg[0]*(f_p2[3]*x))**2
                       +(reg[0]*(f_p2[4]))**2))))
    return prd

def dbdh(x, p2):
    dbdh = (p2[0]*x**4+p2[1]*x**3+p2[2]*x**2+p2[3]*x+p2[4])
    return dbdh

fig = plt.figure()
ax = plt.axes()
x = np.linspace(0,0.07, 60)
plt.plot(x,p[0]*x**5+p[1]*x**4+p[2]*x**3+p[3]*x**2+p[4]*x+p[5],color="green", label=r"Regression 5. Ordnung")
#plt.errorbar(ab1, [B1_lang[i][0] for i in range(len(B1_lang))], xerr = f_ab1, yerr = [B1_lang[i][1] for i in range(len(B1_lang))], fmt='.', color = "blue", label="Induktionsspule")
plt.errorbar(ab1, B1, yerr = f_B, fmt='.', color = "blue", label=r"Messwerte für B(h)")
#plt.axhline(y=mu*N[0]*0.5/(L[0]), color="red",label=r"Theorie: B(0)")
plt.ylabel(r"$B$ in T")
plt.xlabel(r"Höhe $h$ in m")
plt.axis([0, 0.07,-0.6,-0.15]) # ([x_Achsenstart, x-Achsenstopp, y-Achsenstart, y-Achsenstopp])
ax.legend(loc="lower left",frameon=True)
plt.savefig("B(h).png")
plt.show()


f_bdbdh = [np.sqrt((f_B*((p2[0]*ab1[i]**4+p2[1]*ab1[i]**3+p2[2]*ab1[i]**2+p2[3]*ab1[i]+p2[4])))**2
                  +(f_p2[0]*ab1[i]**4*B1[i])**2+(f_p2[1]*ab1[i]**3*B1[i])**2
                  +(f_p2[2]*ab1[i]**2*B1[i])**2+(f_p2[3]*ab1[i]*B1[i])**2
                  +(f_p2[4]*B1[i])**2) for i in range(len(ab1))]

fig = plt.figure()
ax = plt.axes()
x = np.linspace(0,0.07, 60)
plt.plot(x,(p[0]*x**5+p[1]*x**4+p[2]*x**3+p[3]*x**2+p[4]*x+p[5])*(p2[0]*x**4+p2[1]*x**3+p2[2]*x**2+p2[3]*x+p2[4]),
         color="orange", label=r"$B$ Regression: $B\cdot\frac{\partial B}{\partial h}$")
#plt.errorbar(ab1, [B1_lang[i][0] for i in range(len(B1_lang))], xerr = f_ab1, yerr = [B1_lang[i][1] for i in range(len(B1_lang))], fmt='.', color = "blue", label="Induktionsspule")
plt.errorbar(ab1, [B1[i]*(p2[0]*ab1[i]**4+p2[1]*ab1[i]**3+p2[2]*ab1[i]**2+p2[3]*ab1[i]+p2[4]) for i in range(len(B1))], fmt='.', color = "red", label=r"$B$ Messwerte: $B\cdot\frac{\partial B}{\partial h}$")
#plt.axhline(y=mu*N[0]*0.5/(L[0]), color="red",label=r"Theorie: B(0)")
plt.ylabel(r"$B\cdot\frac{\partial B}{\partial h}$ in $\mathrm{T}^2\mathrm{m}^{-1}$")
plt.xlabel(r"Höhe $h$ in m")
plt.axis([0, 0.07,0,7]) # ([x_Achsenstart, x-Achsenstopp, y-Achsenstart, y-Achsenstopp])
ax.legend(loc="upper left",frameon=True)
plt.savefig("B(h)dBdh.png")
plt.show()


a, hoehe, prd = [],[],[]
for j in range(3):
    a.append(gew_mittel([reg(posm[j][i],p,f_p)[0] for i in range(len(posm[j]))],[reg(posm[j][i],p,f_p)[1] for i in range(len(posm[j]))]))
    hoehe.append(gew_mittel(posm[j], [0.000000001 for i in range(len(posm[j]))]))
    prd.append(gew_mittel([produkt(posm[j][i], a[j], p, f_p, p2, f_p2)[0] for i in range(len(posm[j]))],[produkt(posm[j][i], a[j], p, f_p, p2, f_p2)[1] for i in range(len(posm[j]))]))
    #print(poso[j][i]*1000, a , produkt(poso[j][i], a, p, f_p, p2, f_p2))
    
print("Gewichtete Mittelwerte für die Höhe, B und B·∂B∂h an den Messpositionen der Probekörper")
    
for i in range(3):
    print("&$%.1f$ & $%.5f \pm %.5f$ &$%.5f \pm %.5f$ \\\\\hline" % (hoehe[i][0]*1000, a[i][0], a[i][1], prd[i][0],prd[i][1]))


F1 = []
print("Gewichtete Mittelwerte für die Kraft auf die Probekörper.")
for i in range(len(wm1)):
    F1.append(gew_mittel([(wm1[i][j]-wo1[i][j])*9.81 for j in range(3)],[np.sqrt(2)*f_W*9.81 for j in range(3)]))
    print("&$%.4f \pm %.4f$\\\\\hline *10^{-4}" % (F1[i][0]*10**4,F1[i][1]*10**4))


def chi(F, rho, m, Bdbdh, f_F, f_bdbdh):
    mu = 4*np.pi*10**(-7)
    chi = [mu*F*rho/(m*Bdbdh)]
    #chi.append(np.sqrt((mu*f_F*rho/(m*Bdbdh))**2 + (f_bdbdh*mu*F*rho/(m*Bdbdh)**2)**2))
    chi.append((mu*f_F*rho/(m*Bdbdh)))
    return chi

print("Magnetische  und  spezifische  Suszeptibilität  (Fehler  der  Regression  vernachlässigt).")

for i in range(3):
    print("& $%.3f \pm %.3f$ & $%.3f \pm %.3f$\\\\\hline"%(chi(F1[i][0], rho[i], m[i], prd[i][0], F1[i][1], prd[i][1])[0]*10**4,
                                                        chi(F1[i][0], rho[i], m[i], prd[i][0], F1[i][1], prd[i][1])[1]*10**4,
                                                        chi(F1[i][0], rho[i], m[i], prd[i][0], F1[i][1], prd[i][1])[0]/rho[i]*10**6,
                                                        chi(F1[i][0], rho[i], m[i], prd[i][0], F1[i][1], prd[i][1])[1]/rho[i]*10**6))



# Für verschiedene Stromstärken den Gradienten von B bestimmen
fig = plt.figure()
ax = plt.axes()
x = np.linspace(0.04,0.07, 100)

print("Ergebnisse der Regression:")
Bm, Bcov = [],[]
for i in range(4):
    Bm1, Bcov1 = np.polyfit(ab2, B2[i], 2, cov=True)
    Bm.append(Bm1), Bcov.append(Bcov1)
    #f_Bm.append([np.sqrt(abs(Bcov[i][j][j])) for j in range(3)])
    print("$%.1f$ & $%.3f$ & $%3f$& $%3f$ \\\\\hline"
          %(Strom[i], Bm[i][0], Bm[i][1],Bm[i][2]))

plt.errorbar(ab2, B2[0], yerr = f_B, fmt='.', color = "blue", label="$B$ bei %.1f A"%(Strom[0]))
plt.plot(x, Bm[0][0]*x**2+Bm[0][1]*x+Bm[0][2], color="green", label=r"Reg. für %.1f A"%(Strom[0]))
plt.errorbar(ab2, B2[1], yerr = f_B, fmt='.', color = "orange", label="$B$ bei %.1f A"%(Strom[1]))
plt.plot(x, Bm[1][0]*x**2+Bm[1][1]*x+Bm[1][2], color="yellow", label=r"Reg. für %.1f A"%(Strom[1]))
plt.errorbar(ab2, B2[2], yerr = f_B, fmt='.', color = "red", label="$B$ bei %.1f A"%(Strom[2]))
plt.plot(x, Bm[2][0]*x**2+Bm[2][1]*x+Bm[2][2], color="pink", label=r"Reg. für %.1f A"%(Strom[2]))
plt.errorbar(ab2, B2[3], yerr = f_B, fmt='.', color = "purple", label="$B$ bei %.1f A"%(Strom[3]))
plt.plot(x, Bm[3][0]*x**2+Bm[3][1]*x+Bm[3][2], color="black", label=r"Reg. für %.1f A"%(Strom[3]))

plt.xlabel(r"Höhe $h$ in m")
plt.ylabel(r"$B$ in T")
plt.axis([0.04,0.07,-0.6, -0.15]) # ([x_Achsenstart, x-Achsenstopp, y-Achsenstart, y-Achsenstopp])
ax.legend(frameon=False, loc='lower center', ncol=2)
plt.savefig("B_strom.png")
plt.show()


fig = plt.figure()
ax = plt.axes()
x = np.linspace(0.04,0.07, 100)

plt.plot(x, 2*Bm[0][0]*x+Bm[0][1], color="blue", label=r"Gradient für %.1f A"%(Strom[0]))
plt.plot(x, 2*Bm[1][0]*x+Bm[1][1], color="green", label=r"Gradient für %.1f A"%(Strom[1]))
plt.plot(x, 2*Bm[2][0]*x+Bm[2][1], color="orange", label=r"Gradient für %.1f A"%(Strom[2]))
plt.plot(x, 2*Bm[3][0]*x+Bm[3][1], color="red", label=r"Gradient für %.1f A"%(Strom[3]))

plt.xlabel(r"Höhe $h$ in m")
plt.ylabel(r"$\frac{\partial B}{\partial h}$ in T/m")
plt.axis([0.04,0.07,1, 5.5]) # ([x_Achsenstart, x-Achsenstopp, y-Achsenstart, y-Achsenstopp])
ax.legend(frameon=False, loc='lower center', ncol=2)
plt.savefig("dB_strom.png")
plt.show()

for i in range(4):
    print("$%.1f$ & $%.3f h + %.3f$ \\\\\hline"
          %(Strom[i], 2*Bm[i][0], Bm[i][1]))


Delta_m = [gew_mittel([wm2[i][j]-wo2[i][j] for j in range(3)],[np.sqrt(2)*f_W for j in range(3)]) for i in range(4)]
FI = [[Delta_m[i][0]*9.81, Delta_m[i][1]*9.81] for i in range(4)]
print(FI)

mf, b, r, p, s = stats.linregress(Strom, [FI[i][0] for i in range(4)])
f_mf = lin_sigma(Strom, [FI[i][0] for i in range(4)], mf, b)[0]
f_b = lin_sigma(Strom, [FI[i][0] for i in range(4)], mf, b)[1]
print("""m & $%.4f \pm %.4f$ \\\\\hline
          b & $%.4f \pm %.4f$ \\\\\hline""" % (mf*10**5, f_mf*10**5, b*10**5, f_b*10**5))

fig = plt.figure()
ax = plt.axes()
x = np.linspace(0.7,1.5, 100)
plt.errorbar(Strom, [FI[i][0] for i in range(4)], yerr = [FI[i][0] for i in range(4)],
             fmt=".", color = "blue", label="Kraft $F$")
plt.plot(x, mf*x+b, color = "green", label="Lineare Regression")
plt.xlabel(r"$I$ in A")
plt.ylabel(r"$F$ in N")
plt.axis([0.7,1.5,-0.00005, 0.0003]) # ([x_Achsenstart, x-Achsenstopp, y-Achsenstart, y-Achsenstopp])
ax.legend(frameon=False, loc='lower center', ncol=2)
plt.savefig("F_strom.png")
plt.show()


dbdh1 = dbdh(ab2[2], p2)
print(dbdh1)
print(ab2[2])
chi1 = 7.116*10**(-4)

def FeldausF(dbdh, F_m, F_b, chi, I, rho, m):
    B = (F_m*I+F_b)*4*np.pi*10**(-7)*rho/(m*chi*dbdh)
    return B

fig = plt.figure()
ax = plt.axes()
x = np.linspace(0.8,1.4, 100)
#y = [FeldausF(dbdh1, mf, b, chi1, wert, rho[1], m[1]) for wert in x]
plt.errorbar(Strom, [B2[i][2] for i in range(4)], yerr = f_B,
             fmt=".", color = "blue", label="Messwerte $B$")
plt.plot(x,(mf*x+b)*4*np.pi*10**(-7)*rho[1]/(m[1]*chi1*dbdh1) , color = "green", label="Berechnung aus $F$")
plt.xlabel(r"$I$ in A")
plt.ylabel(r"$B$ in T")
#plt.axis([0.75,1.45,-0.55,-0.2]) # ([x_Achsenstart, x-Achsenstopp, y-Achsenstart, y-Achsenstopp])
ax.legend(frameon=True, loc='lower left')
plt.savefig("FeldausF.png")
plt.show()




