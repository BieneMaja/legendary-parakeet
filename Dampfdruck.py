import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats
plt.style.use("seaborn-whitegrid")

daten = open("Dampfdaten.txt", mode="r")
inhalt = daten.readlines()
daten.close()
inhalt = [x.strip() for x in inhalt]
widerstand1 =[int(x) for x in inhalt[1].split(",")]
widerstand2 =[int(x) for x in inhalt[3].split(",")]
widerstand2.sort()

def temp(R):
    A = 3.9083*0.001
    B = -5.775*10**(-7)
    R_0=1000
    t = -A/(2*B)-np.sqrt((A/(2*B))**2+(R/R_0-1)/B) + 273.115
    return t

Temp1 = [temp(widerstand1[i]) for i in range(len(widerstand1))]
Temp2 = [temp(widerstand2[i]) for i in range(len(widerstand2))]
Temp_fehler1 = [0.3+0.005*x for x in Temp1]
Temp_fehler2 = [0.3+0.005*x for x in Temp2]
Druck = [0.000001]
for i in range(len(Temp1)-1):
    Druck.append(i+1)
Druck_plot= [np.log(x*10**5) for x in Druck] # von Bar zu Pascal
Temp_fehler_plot1 = [Temp_fehler1[i]*1/Temp1[i]**2 for i in range(len(Temp1))]
Temp_fehler_plot2 = [Temp_fehler2[i]*1/Temp2[i]**2 for i in range(len(Temp2))]
Temp_plot1 = [1/x for x in Temp1]
Temp_plot2 = [1/x for x in Temp2]

m1, b1, r_value1, p_value1, std_err1 = stats.linregress(Temp_plot1[4:],Druck_plot[4:])
print(m1, b1, r_value1, p_value1, std_err1)
m2, b2, r_value2, p_value2, std_err2 = stats.linregress(Temp_plot2,Druck_plot[5:])
print(m2, b2, r_value2, p_value2, std_err2)

fig = plt.figure()
ax = plt.axes()
x = np.linspace(0,100,100)
plt.plot(m1*x+b1+0.1, label="Regression")
plt.errorbar(Temp_plot1,Druck_plot,xerr = Temp_fehler_plot1,fmt='.',label="Daten (aufwärmend)")
plt.xlabel(r"$\frac{1}{T} / \frac{1}{\mathrm{K}}$")
plt.ylabel(r"$ln(p) /$ Pa")
plt.axis([0.0019,0.0028,11,15.3]) # ([x_Achsenstart, x-Achsenstopp, y-Achsenstart, y-Achsenstopp])
ax.legend(loc="upper right",frameon=True)
#plt.savefig("P_T_1.png")
plt.show()
fig = plt.figure()
ax = plt.axes()
plt.plot(m2*x+b2+0.1, label="Regression")
plt.errorbar(Temp_plot2,Druck_plot[5:],xerr = Temp_fehler_plot2,fmt='.',label="Daten (abfallend)")
plt.xlabel(r"$\frac{1}{T} / \frac{1}{ \mathrm{K}}$")
plt.ylabel(r"$ln(p) /$ Pa")
plt.axis([0.0019,0.0024,13,15.3]) # ([x_Achsenstart, x-Achsenstopp, y-Achsenstart, y-Achsenstopp])
ax.legend(loc="upper right",frameon=True)
plt.savefig("P_T_2.png")
plt.show()

def lambdav(m):
    R = 8.314
    lambdav = -m*R
    return lambdav
def s_lambdav(lambdav,sm):
    R = 8.314
    sigma = abs(sm*R)
    return sigma
def t_s(m,p,b):
    t = m/(np.log(p)-b)
    return t
def s_t_s(sm,p,b):
    sigma = abs(sm/(np.log(p)-b))
    return sigma
def P_s(m,T,b):
    ps = np.exp(m/T +b)
    return ps
def s_P_s(m,sm,T,b):
    sigma = abs(sm* 1/T*np.exp(m/T+b))
    return sigma

Lambdav = [lambdav(m1),lambdav(m2)]
s_Lambdav = [s_lambdav(Lambdav[0],std_err1),s_lambdav(Lambdav[1],std_err2)]
T_s = [t_s(m1,101300,b1),t_s(m2,101300,b2)]
s_T_s = [s_t_s(std_err1,101300,b1),s_t_s(std_err2,101300,b2)]
p_s = [P_s(m1,273.115,b1),P_s(m2,273.115,b2)]
s_p_s = [s_P_s(m1,std_err1,271.115,b1),s_P_s(m2,std_err2,273.115,b2)]
tabelle = "Anstieg&$"+"{:.1f}".format(Lambdav[0])+"\\pm"+"{:.1f}".format(abs(s_Lambdav[0]))+"$&$"+"{:.1f}".format(T_s[0])+"\\pm"+"{:.1f}".format(abs(s_T_s[0]))+"$&$"+"{:.1f}".format(p_s[0])+"\\pm"+"{:.1f}".format(s_p_s[0])+"$\\"+"\\"+"\\hline "+"Abfall&$"+"{:.1f}".format(Lambdav[1])+"\\pm"+"{:.1f}".format(s_Lambdav[1])+"$&$"+"{:.1f}".format(T_s[1])+"\\pm"+"{:.1f}".format(abs(s_T_s[1]))+"$&$"+"{:.1f}".format(p_s[1])+"\\pm"+"{:.1f}".format(s_p_s[1])+"$\\"+"\\"+"\\hline"
print(tabelle)

def siede(T_0,rho,h,p_0,R,V):
    g= 9.81
    T = 1/(1/T_0+(rho* g* h)/(p_0)*R/V)
    return T
print(siede(373.15,1.293,2962,101300,8.13,40805)-273.15)
