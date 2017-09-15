
# coding: utf-8

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
    

# Daten aus txt-Datei einlesen
daten = open("Wechseldaten.txt", mode="r")
inhalt = daten.readlines()
daten.close()
inhalt = [x.strip() for x in inhalt]

# Innenwiderstand und Fehler Amperemeter in Ohm
Ra = float(inhalt[1])
f_Ra = float(inhalt[3])
# Ohmscher Widerstand in Ohm
Ro = float(inhalt[5])
f_Ro = float(inhalt[7])
#Spulenwiderstand in Ohm
Rs = float(inhalt[9])
# Wicklungszahl n und Durchmesser d (in m) der Spule
n = int(inhalt[11])
d = float(inhalt[13])
# Kapazität in Farad
c = float(inhalt[13])*10**(-6)
f_c = float(inhalt[15])*0.01*c

#Versuchsteil 1
#Frequenzen in Hz
F1 = [float(x) for x in inhalt[21].split(",")]
# Strom in A
I1 = [float(x)*0.001 for x in inhalt[23].split(",")]
f_I1 = float(inhalt[25])*0.001
# Spannung in V
U1 = [float(x) for x in inhalt[27].split(",")]
f_U1 = float(inhalt[29])
# Phasenverschiebung
phi1 = [np.pi/float(x) for x in inhalt[31].split(",")]
f_phi1 = np.pi/float(inhalt[33])

#Versuchsteil 2
#Frequenzen in Hz
F2 = [float(x) for x in inhalt[37].split(",")]
#Strom in Ampere
I2 = [float(x)*0.001 for x in inhalt[39].split(",")]
f_I2 = float(inhalt[41])*0.001
# U in V
U2 = [float(x) for x in inhalt[43].split(",")]
f_U2 = float(inhalt[45])
# Phasenverschiebung
phi2 = []
for x in inhalt[47].split(","):
    if float(x) == 0:
        phi2.append(0)
    else:
        phi2.append(1/float(x))
f_phi2 = 1/float(inhalt[49])
# Uc in V
Uc = [float(x) for x in inhalt[51].split(",")]
f_Uc = float(inhalt[53])
# Ulc in V
Ulc = [float(x) for x in inhalt[55].split(",")]
f_Ulc = float(inhalt[57])

# Versuchsteil 3
#Frequenzen in Hz
F3 = [float(x) for x in inhalt[61].split(",")]
# Strom in A
I3 = [float(x)*0.001 for x in inhalt[63].split(",")]
f_I3 = float(inhalt[65])*0.001
# Spannung in V
U3 = [float(x) for x in inhalt[67].split(",")]
f_U3 = float(inhalt[69])

# Plot Impedanz Z^2 = (U/I)^2 gegen w^2= (2 pi f)^2
w2 = [(2*np.pi*x)**2 for x in F1]
Z2 = [(U1[i]/I1[i])**2 for i in range(len(F1))]
f_z = [np.sqrt((f_U1/I1[i])**2 + (U1[i]*f_I1/I1[i]**2)) for i in range(len(U1))]

fig = plt.figure()
ax = plt.axes()
x = np.linspace(4,16, 10000)
m, b, r, p, s = stats.linregress(w2, Z2)

print("m = %.6f \pm %.6f" % (m, lin_sigma(w2, Z2, m, b)[0]))
print("b = %.2f \pm %.2f" % (b, lin_sigma(w2, Z2, m, b)[1]))
print("r=", r)

w2_plot = [w2[i]*10**(-5) for i in range(len(w2))]
Z2_plot = [Z2[i]*10**(-5) for i in range(len(w2))]
f_z_plot = [f_z[i]*10**(-5) for i in range(len(f_z))]
m_plot = m*10**(-5)
b_plot = b*10**(-5)
print("m_plot = ", m_plot,"b_plot=", b_plot)

plt.errorbar(w2_plot, Z2_plot, yerr = f_z_plot, fmt='.', color = "blue", label="Daten")
plt.plot(x, m*x+b_plot, color="green", label="Lineare Regression")

plt.ylabel(r"$Z^2$ in $10^{5}\Omega^2$")
plt.xlabel(r"$\omega^2$ in $10^{5}\frac{1}{\mathrm{s}^2}$")
#plt.axis([390000, 1610000,470,1000]) # ([x_Achsenstart, x-Achsenstopp, y-Achsenstart, y-Achsenstopp])
ax.legend(loc="upper left",frameon=True)
plt.savefig("w2z2.png")
plt.show()

L = np.sqrt(m)
f_L = 1/(2*np.sqrt(m))*lin_sigma(w2, Z2, m, b)[0]
R = np.sqrt(b)
f_R = 1/(2*np.sqrt(b))*lin_sigma(w2, Z2, m, b)[1]
print("L &= %.4f \pm %.4f \Omega\cdot \mathrm{s}" % (L, f_L))
print("R &= %.3f \pm %.4f \Omega" % (R, f_R))

# Plot Impedanz Z = (U/I) gegen w= (2 pi f)
w1 = [(2*np.pi*x) for x in F2]
Z1 = [(U2[i]/I2[i]) for i in range(len(F2))]
f_z1 = [np.sqrt((f_U2/I2[i])**2 + (U2[i]*f_I2/I2[i]**2)) for i in range(len(U2))]
print(Z1[4:7]) # Minimum suchen...
print("Z1_min = ", Z1[5], "\pm", f_z1[5])
print(w1[5])
fig = plt.figure()
ax = plt.axes()
print("Abweichung: ", Z1[4]-Z1[6], Z1[6]-Z1[5])

plt.errorbar(w1, Z1, yerr = f_z1, fmt='.', color = "blue")

plt.ylabel(r"$Z$ in $\Omega$")
plt.xlabel(r"$\omega$ in $\frac{1}{\mathrm{s}}$")
#plt.axis([390000, 1610000,470,1000]) # ([x_Achsenstart, x-Achsenstopp, y-Achsenstart, y-Achsenstopp])
ax.legend(loc="upper left",frameon=True)
plt.savefig("wz.png")
plt.show()


# Plot Phasenverschiebung gegen w

fig = plt.figure()
ax = plt.axes()
plt.errorbar(w1, phi2, yerr = f_phi2, fmt='.', color = "blue")

plt.ylabel(r"$\phi$ in $\pi$")
plt.xlabel(r"$\omega$ in $\frac{1}{\mathrm{s}}$")
#plt.axis([390000, 1610000,470,1000]) # ([x_Achsenstart, x-Achsenstopp, y-Achsenstart, y-Achsenstopp])
#ax.legend(loc="upper left",frameon=True)
plt.savefig("phi1.png")
plt.show()
print(inhalt[55])

# gewichteter Mittelwert des Gesamtwiderstandes
def gew_mittel(x, s_x):    
    sum1, sum2 = 0, 0
    for i in range(len(x)):
        sum1 += x[i]/s_x[i]**2
        sum2 += 1/s_x[i]**2
    x = [sum1/sum2, np.sqrt(1/sum2)]
    return x

a = gew_mittel([R, Z1[5]], [f_R, 40])
b = gew_mittel([890, 890], [60, 10])
print("R_gesamt = ", a)
print("w_r = ",b)

RL = a[0] - Ro - Ra
f_RL = np.sqrt(a[1]**2 + f_Ro**2 + f_Ra**2)
print("Spulenwiderstand R_L = ", RL, "\pm", f_RL)
print("Angabe: ", Rs)

C = 1/(b[0]**2*L)
s_C = np.sqrt((f_L/(b[0]**2*L**2))**2+(2*b[1]/(L*b[0]**3))**2)
print("Kapazität der Spule C = ", C, "\pm", s_C)

# gemeinsamer Gruppenplot

fig = plt.figure()
ax = plt.axes()
plt.errorbar(w1, U2, yerr = f_U2, fmt='.', color = "blue", label=r"$U$")
plt.errorbar(w1, Uc, yerr = f_Uc, fmt='.', color = "red", label = r"$U_C$")
plt.errorbar(w1, Ulc, yerr = f_Ulc, fmt='.', color = "green", label = r"$U_{R+L}$")

plt.ylabel(r"$U$ in V")
plt.xlabel(r"$\omega$ in $\frac{1}{\mathrm{s}}$")
#plt.axis([600, 1300,0,3]) # ([x_Achsenstart, x-Achsenstopp, y-Achsenstart, y-Achsenstopp])
ax.legend(loc="upper right",frameon=True)
plt.savefig("UUU.png")
plt.show()


print(U2[5], Uc[5], Ulc[5])
print("\pm", f_U2, f_Uc, f_Ulc)
cp = U2[5]/Ulc[5]
f_cp = np.sqrt(np.sqrt((f_U2/Ulc[5])**2 + (U2[5]*f_Ulc/Ulc[5]**2)))
print("phi = ", np.arccos(cp), "\pm", np.sqrt(f_U2**2/(1-U2[5]/Ulc[5])+ f_Ulc**2/Ulc[5]**4*1/(1-U2[5]/Ulc[5])))

f_phi2 =1/(1+(b[0]*L/R)**2)*np.sqrt(b[1]**2 + f_L**2 + f_R**2/R**4)
print("phi2 = ", np.arctan(890*0.794/64.2), "\pm",f_phi2)

plt.rc("text", usetex=False)
soa = np.array([[0, 0, 1, 0], [0, 0, 0, 1], [0, 0, cp, np.sin(np.arccos(cp))]])
X, Y, U, V = zip(*soa)
plt.figure(figsize=(5,5))
ax = plt.gca()
ax.quiver(X, Y, U, V, angles='xy', scale_units='xy', scale=1)
ax.set_xlim([-0.5, 1.5])
ax.set_ylim([-0.5, 1.5])
plt.xticks([])
plt.yticks([])
plt.annotate("U",xy=(1,-0.1))
plt.annotate("U_{L+R}",xy=(0.3,1))
plt.annotate("-U_C", xy=(-0.2,1))
plt.savefig("zeiger.png")
plt.draw()
plt.show()

# Parallelkreis
w3 = [2*np.pi*F3[i] for i in range(len(F3))]
z3 = [U3[i]/I3[i] for i in range(len(U3))]
f_z3 = [np.sqrt((f_U3/I3[i])**2 + (U3[i]*f_I3/I3[i]**2)) for i in range(len(U3))]

fig = plt.figure()
ax = plt.axes()
plt.errorbar(w3, z3, yerr = f_z3, fmt='.', color = "blue")

plt.ylabel(r"$Z$ in $\Omega$")
plt.xlabel(r"$\omega$ in $\frac{1}{\mathrm{s}}$")
#plt.axis([390000, 1610000,470,1000]) # ([x_Achsenstart, x-Achsenstopp, y-Achsenstart, y-Achsenstopp])
#ax.legend(loc="upper left",frameon=True)
plt.savefig("wz_parallel.png")
plt.show()
