
# coding: utf-8

# In[26]:

# Versuch 7 Coulombsches Gesetz

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

plt.style.use("seaborn-whitegrid")

# Daten aus txt-Datei einlesen
daten = open("C_Daten.txt", mode="r")
inhalt = daten.readlines()
daten.close()
inhalt = [x.strip() for x in inhalt]

# Beladespannung in kVolt
bs = [float(x) for x in inhalt[1].split(",")]
#Abstand in m
ab = [float(x) for x in inhalt[4].split(",")]
# Spannung am Messgerät
sp = []
for i in range(10):
    sp.append([float(x) for x in inhalt[7+i].split(",")])
# Messbereich wechseln [[1->2],[2->3]]
mw = [[int(x) for x in inhalt[19].split(",")],[int(x) for x in inhalt[21].split(",")]]
#Fehler in Bereich [1, 2, 3, Abstand]
fehler = [float(inhalt[24]),float(inhalt[27])]

# Gemessene Spannung zuFeldstärke umrechnen: Messbereich 1:

for i in range(10):
    for j in range(10):
        if(j>mw[1][i]):# Messbereich 3
            sp[i][j] *= 10
        elif (j<=mw[0][i]): # Messbereich 1
            sp[i][j] = sp[i][j]*0.1

#Kugelradius
d = [float(x)*0.5 for x in inhalt[30].split(",")]
# Beladespannung in kiloVolt
bsp = [float(x) for x in inhalt[33].split(",")]
# Gemessene Spannung in Volt
gsp = [[float(x) for x in inhalt[36].split(",")],
       [float(x) for x in inhalt[37].split(",")],
       [float(x) for x in inhalt[38].split(",")]]
# Messbereich gewechselt
mw2 = [float(x) for x in inhalt[41].split(",")]

for i in range(3):
    for j in range(5):
        if (j<=mw2[i]): # Messbereich 1
            gsp[i][j] = gsp[i][j]*0.1


# In[33]:

# Versuchsteil 1
            
slope, intercept, r_value, p_value, std_err = [0 for i in range(10)], [0 for i in range(10)], [0 for i in range(10)], [0 for i in range(10)], [0 for i in range(10)]
        
for i in range(10):
    slope[i], intercept[i], r_value[i], p_value[i], std_err[i] = stats.linregress(bs, sp[i])

print("\\begin{table}[H]\n\centering\n\\begin{tabular}{|c|c|}\hline\n Abstand in m & Steigung in $\\frac{1}{m}$ \\\\\hline \n")
for i in range(10):
    print("$%.2f$ & $%.2f \pm %.2f$ \\\\\hline" % (ab[i], slope[i], std_err[i]))
print("\end{tabular}\n\caption{Steigung der linearen Regressionen}\n \label{tab:steigung}\n\end{table}")

fig = plt.figure()
ax = plt.axes()
x = np.linspace(0,10000,10000)


plt.errorbar(bs, sp[0], yerr = fehler[0], fmt='.', color = "blue", label=str(ab[0])+"m")
plt.plot(slope[0]*x+intercept[0], color = "blue")
plt.errorbar(bs, sp[1], yerr = fehler[0], fmt='.', color = "green", label=str(ab[1])+"m")
plt.plot(slope[1]*x+intercept[1], color = "green")
plt.errorbar(bs, sp[2], yerr = fehler[0], fmt='.', color = "yellow", label=str(ab[2])+"m")
plt.plot(slope[2]*x+intercept[2], color = "yellow")
plt.errorbar(bs, sp[3], yerr = fehler[0], fmt='.', color = "orange", label=str(ab[3])+"m")
plt.plot(slope[3]*x+intercept[3], color = "orange")
plt.errorbar(bs, sp[4], yerr = fehler[0], fmt='.', color = "red", label=str(ab[4])+"m")
plt.plot(slope[4]*x+intercept[4], color = "red")
plt.errorbar(bs, sp[5], yerr = fehler[0], fmt='.', color = "pink", label=str(ab[5])+"m")
plt.plot(slope[5]*x+intercept[5], color = "pink")
plt.errorbar(bs, sp[6], yerr = fehler[0], fmt='.', color = "purple", label=str(ab[6])+"m")
plt.plot(slope[6]*x+intercept[6], color = "purple")
plt.errorbar(bs, sp[7], yerr = fehler[0], fmt='.', color = "grey", label=str(ab[7])+"m")
plt.plot(slope[7]*x+intercept[7], color = "grey")
plt.errorbar(bs, sp[8], yerr = fehler[0], fmt='.', color = "black", label=str(ab[8])+"m")
plt.plot(slope[8]*x+intercept[8], color = "black")
plt.errorbar(bs, sp[9], yerr = fehler[0], fmt='.', color = "brown", label=str(ab[9])+"m")
plt.plot(slope[9]*x+intercept[9], color = "brown")

plt.ylabel("Feldstärke in kV/m")
plt.xlabel("Beladespannung in kV")
plt.axis([-1, 10.5,0,11]) # ([x_Achsenstart, x-Achsenstopp, y-Achsenstart, y-Achsenstopp])
ax.legend(loc="upper left",frameon=True)
plt.savefig("Plot_1.png")
plt.show()


# In[28]:

# Doppeltlogarithmisch Steigung gegen Abstand

fig = plt.figure()
ax = plt.axes()
x = np.linspace(0,100,100)

ln_bs = [np.log(x) for x in bs]
ln_slope = [-1*np.log(x) for x in slope]
ln_std_err = [std_err[i]/slope[i] for i in range(10)]

m, b, r, p, s = stats.linregress(ln_bs, ln_slope)

print("m = $%.2f \pm %.2f$" % (m, s))
print("b = %.2f" % b)

plt.errorbar(ln_bs, ln_slope, yerr = ln_std_err, fmt='.', color = "blue")
plt.plot(m*x+b, color = "blue")

plt.ylabel(r"ln(m) in $ln(m)$")
plt.xlabel(r"ln(r) in $ln(m)$")
plt.axis([-0.2, 2.5,-0.5,2.2]) # ([x_Achsenstart, x-Achsenstopp, y-Achsenstart, y-Achsenstopp])
ax.legend(loc="upper left",frameon=True)
plt.savefig("Plot_2.png")
plt.show()


# In[32]:

# Feldstärke gegen Beladespannung

fig = plt.figure()
ax = plt.axes()
x = np.linspace(0,100,100)

m1, b1, r1, p1, s1 = [0 for x in range(5)], [0 for x in range(5)], [0 for x in range(5)], [0 for x in range(5)], [0 for x in range(5)]

for i in range(3):
    m1[i], b1[i], r1[i], p1[i], s1[i] = stats.linregress(bsp, gsp[i])

plt.errorbar(bsp, gsp[0], yerr = fehler[0], fmt='.', color = "blue", label="R=1cm")
plt.plot(m1[0]*x+b1[0], color = "blue")
plt.errorbar(bsp, gsp[1], yerr = fehler[0], fmt='.', color = "green", label="R=2cm")
plt.plot(m1[1]*x+b1[1], color = "green")
plt.errorbar(bsp, gsp[2], yerr = fehler[0], fmt='.', color = "red", label="R=6cm")
plt.plot(m1[2]*x+b1[2], color = "red")

plt.xlabel("Beladespannung in kV")
plt.ylabel(r"Feldstärke in $\frac{kV}{m}$")
plt.axis([0.5, 5.2,0,8]) # ([x_Achsenstart, x-Achsenstopp, y-Achsenstart, y-Achsenstopp])
ax.legend(loc="upper left",frameon=True)
plt.savefig("Plot_3.png")
plt.show()

prop = []
for i in range(3):
    prop.append(m1[i]*(0.25)**2/d[i])
    
prop_fehler = []
for i in range(3):
    prop_fehler.append(1/d[i]*np.sqrt(s1[i]**2*(0.25)**2 + fehler[1]**2*4*(0.25)**2*m1[i]**2))

print("\\begin{table}[H]\n\centering\n\\begin{tabular}{|c|c|}\hline\n Kugelradius in $10^{-2}m$ & Steigung in $\\frac{1}{m}$ \\\\\hline \n")
print("1 & $%.2f \pm %.2f$ \\\\\hline \n 2 & $%.2f \pm %.2f$ \\\\\hline \n 6 & $%.2f \pm %.2f$ \\\\\hline \n" % (prop[0], prop_fehler[0], prop[1], prop_fehler[1], prop[2], prop_fehler[2]))
print("\end{tabular}\n\caption{Proportionalitätsfaktor}\n \label{tab:prop}\n\end{table}")


# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:



