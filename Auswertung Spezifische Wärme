from uncertainties import unumpy
from uncertainties import ufloat
import matplotlib.pyplot as plt
import numpy as np
plt.style.use("seaborn-whitegrid")

daten = open("Messdaten_spezi.txt", mode="r")
inhalt = daten.readlines()
daten.close()
inhalt = [x.strip() for x in inhalt] 
spannung =[int(x) for x in inhalt[1].split(",")]
ddruck = [[float(x) for x in inhalt[i].split(",")] for i in range(3,13)]
dp = unumpy.uarray(ddruck,[[2,2,2],[2,2,2],[2,2,2],[2,2,2],[2,2,2],[2,2,2],[2,2,2],[5,5,5],[5,5,5],[5,5,5]])
gas = [inhalt[22],inhalt[36],inhalt[47]]
schwingungen = [int(x) for x in inhalt[24].split(",")]
messung = [float(x) for x in inhalt[20].split(",")]
schdauer = [[[int(x) for x in inhalt[i].split(",")]for i in range(26,31)],
           [[int(x) for x in inhalt[i].split(",")]for i in range(37,42)],
           [[int(x) for x in inhalt[i].split(",")]for i in range(48,53)]]
amplitude = [float(x) for x in [inhalt[32],inhalt[43],inhalt[54]]]
s_amplitude = [float(x) for x in[inhalt[34],inhalt[45],inhalt[56]]]

barsd = [np.mean([0.001*schdauer[i][j][k]/schwingungen[j] for j in range(5)
                  for k in range(len(schdauer[i][j]))]) for i in range(3)]
schdauer_stat = [np.sqrt(1/(22*21)*sum([(0.001*schdauer[i][j][k]/schwingungen[j]-barsd[i])**2 for j in range(5)
                                               for k in range(len(schdauer[i][j]))])) for i in range(3)]
h1 = [[int(x) for x in inhalt[i].split(",")] for i in [63,68,73]]
s_h = 5
h2 = [[float(x) for x in inhalt[i].split(",")] for i in [65,70,75]]
vol = 10
b = float(inhalt[79])

def Druck(b,m,A):
    g=9.81
    p=b+m*g/A
    return p
def meff(m,A,p,a):
    meff = m+A*p*2*a
    return meff
def s_meff(A,p,a,sa):
    sigma = sa*(A*p*2)
    return sigma
# Freiheitsgrade der Luft
def frei(q,p,dp,V,dv):
    f= 2*(q-p*dv)/(V*dp+p*dv)
    return f
# Molwärme der Luft
def molw(f):
    R=8.314
    cv = f*R*0.5
    return cv
# cp/cv nach Rüchard (Adiabatenexponent)
def cpcv1(m,V,T,p,d):
    cpcv = (64*m*V)/(T**2*p*d**4)
    return cpcv
#Adiabatenexponent nach Clement-Desormes
def cpcv2(h1,h2):
    cpcv = h1 / (h1-h2)
    return cpcv
def sigma_cpcv2(h1,h2,s_h):
    sigma = np.sqrt(s_h**2*(((h1-h2-1)*h1/(h1-h2)**2)**2+(-h1/(h1-h2)**2)**2))
    return sigma
def s_cpcv1(m,V,T,p,d,s_t,s_m):
    sigma = np.sqrt(s_t**2*(-(2*T*p*d**4)*(64*m*V)/(T**2*p*d**4)**2)**2+s_m**2*((64*V)/(T**2*p*d**4))**2)
    return sigma
def freicpcv(cpcv):
    f=2/(cpcv-1)
    return f
def sigma_freicpcv(cpcv,s_cpcv):
    sigma = s_cpcv*(-2)/(cpcv-1)**2
    return sigma
def gew_mittel(x,s):
    mittel = (sum([x[i]/s[i]**2 for i in range(len(x))]))/(sum([1/s[i]**2 for i in range(len(x))]))
    sigma = np.sqrt(1/sum([1/s[i]**2 for i in range(len(x))]))
    return mittel, sigma

#Versuch A
R,h = float(input("Radius des Zylinders (A) in m: ")),float(input("Höhe des Zylinders (A) in m: "))
c= 10**(-5)
r1,r2=0.002,0.092
VZ = np.pi*R**2*h
dV = np.pi*r1**2*dp/1000.
dq = [1/2.*c*U**2 for U in spannung]
DP = dp*9.81*(1+r1**2/r2**2)

fig = plt.figure()
ax = plt.axes()
for i in range(3):
    plt.errorbar(dq,
                 [DP[j][i].n for j in range(10)],
                 yerr=[DP[j][i].s for j in range(10)], fmt='.k')
plt.xlabel(r"$\Delta Q /$ J")
plt.ylabel(r"$\Delta p /$ Pa")
plt.savefig("dQdP.png")
plt.show()
f=[[frei(dq[j],103150,DP[j][i],VZ,dV[j][i]) for i in range(3)]for j in range(10)]
f_mittel=np.mean([np.mean(f[i]) for i in range(10)])
molwa = [[molw(x) for x in f[i]] for i in range(10)]
molwa_mittel = np.mean(molwa[i])
longstring = ["$"+str("{:.4f}".format(dq[j]))+"$&$"+str("{:.4f}".format(10**7*dV[j][i].n))+"\\pm"+str("{:.4f}".format(10**7*dV[j][i].s))+"$&$"+str("{:.4f}".format(DP[j][i].n))+"\\pm"+str("{:.4f}".format(DP[j][i].s))+"$&$"+str("{:.4f}".format(f[j][i].n))+"\\pm"+str("{:.4f}".format(f[j][i].s))+"$&$"+str("{:.4f}".format(molwa[j][i].n))+"\\pm"+str("{:.4f}".format(molwa[j][i].s))+"$\\"+"\\" for j in range(10)for i in range(3) ]
verylongstring = """\\begin{table}[h!]
    \\begin{tabular}{|c|c|c|c|c|}\hline
        $\Delta Q$ / J & $\Delta V$ / $10^{-7} m^3$ & $\Delta p$/Pa & $f$ & Molwärme /$\frac{\mathrm{J}}{\mathrm{molK}}$ \\\hline"""+"".join(longstring)+"""\\hline
    \end{tabular}
    \caption{Zwischenwerte zur Berechnung der spezifischen Wärme}
    \label{zwischen1}
\end{table}"""
print("Volumen des Zylinders (A) in m^3 = ",VZ,"\n Zwischenwerte als Latex-Tabelle : ")
print(verylongstring)
print("Mittelwert der Freiheitsgrade = ",f_mittel, "\n Mittelwert der Molwärme (J/(mol*K)) = ", molwa_mittel)

#Versuch B
p=Druck(messung[3]*100,messung[0]*0.001,np.pi*(messung[1]*0.001)**2/4.)
Dichte = [1.293,1.98,1.784]
meff = [meff(messung[0]*0.001,np.pi*(messung[1]*0.001)**2/4.,
             Dichte[i],0.001*amplitude[i])
        for i in range(3)]
s_meff = [s_meff(np.pi*(messung[1]*0.001)**2/4.,Dichte[i],0.001*amplitude[i],0.0005) for i in range(3)]

cpcv=[[[cpcv1(meff[i],messung[2]*0.000001,schdauer[i][j][k]/(1000*schwingungen[j]),p,messung[1]*0.001)
       for k in range(len(schdauer[i][j]))] for j in range(5)] for i in range(3)]
s_cpcv = [[[s_cpcv1(meff[i],messung[2]*0.000001,schdauer[i][j][k]/(1000*schwingungen[j]),
                    p,messung[1]*0.001,schdauer_stat[i],s_meff[i])
            for k in range(len(schdauer[i][j]))] for j in range(5)] for i in range(3)]
a=[[np.mean(cpcv[j][i]) for i in range(5)]for j in range(3)]
mittel_cpcv = [np.mean(a[0]),np.mean(a[1]),np.mean(a[2])]
mittel_s_cpcv=[np.mean([np.mean(s_cpcv[i][j]) for j in range(5)]) for i in range(3)]
longstring = ["$"+str("{:.3f}".format(cpcv[0][j][k]))+"\\pm"+str("{:.3f}".format(s_cpcv[0][j][k]))+"$"+"&"+"$"+str("{:.3f}".format(cpcv[1][j][k]))+"\\pm"+str("{:.3f}".format(s_cpcv[1][j][k]))+"$"+"&"+"$"+str("{:.3f}".format(cpcv[2][j][k]))+"\\pm"+str("{:.3f}".format(s_cpcv[2][j][k]))+"$"+"\\"+"\\" for j in range(5) for k in range(len(cpcv[0][j]))]
verylongstring ="""\\begin{table}
\\begin{tabular}{|c|c|c|}\hline
    $\kappa$ Luft & $\kappa \, CO_2$  & $\kappa$ Argon  \\\hline"""+"".join(longstring)+"""\\hline
    \end{tabular}
    \caption{Zwischenwerte der Berechnung der Adiabatenexponenten}
    \label{zwischen2}
\end{table}"""
print("Effektive Massen (Luft, CO_2,Argon) in kg= ", meff,"\\pm",s_meff)
print("Mittel cp/cv (Luft, CO_2,Argon)=",mittel_cpcv,"\\pm",mittel_s_cpcv)
print("Mittel Freiheitsgrade = ", [freicpcv(mittel_cpcv[i]) for i in range(3)],"\\pm",[sigma_freicpcv(mittel_cpcv[i],mittel_s_cpcv[i]) for i in range(3)])
print("Zwischenwerte als Latex-Tabelle : ",verylongstring)

#Versuch C

cpcv2=[[cpcv2(h1[i][j]*0.001,h2[i][j]*0.001) for j in range(len(h1[i]))] for i in range(len(h1))]
s_cpcv2=[[sigma_cpcv2(h1[i][j]*0.001,h2[i][j]*0.001,s_h*0.001) for j in range(len(h1[i]))]for i in range(len(h1))]
mittel_cpcv2 = [np.mean(cpcv2[i]) for i in range(len(h1))]
mittel_s_cpcv2 = [np.mean(s_cpcv2[i]) for i in range(len(h1))]

print("C: Adiabatenexponent cp/cv = "+ "".join([str(cpcv2[i][j]) + "\\pm"+ str(s_cpcv2[i][j]) + "\n" for i in range(3) for j in range(5)]))
print("C: Mittel Adiabatenexponent cp/cv = ", np.mean(mittel_cpcv2),"\\pm",np.mean(mittel_s_cpcv2))


DP1 = [[(h1[i][j])*9.81*2 for j in range(len(h1[i]))]for i in range(3)]
DP2 = [[(h2[i][j])*9.81*2 for j in range(len(h1[i]))]for i in range(3)]

dV = [[np.pi*0.005**2*0.5*(h1[i][j])/1000. for j in range(len(h1[i]))]for i in range(3)]
dV2 = [[np.pi*0.005**2*0.5*(h2[i][j])/1000. for j in range(len(h1[i]))]for i in range(3)]

V = [[0.001*vol+a for a in dV[i]]for i in range(3)]
P1 = [[100*b+a for a in DP1[i]]for i in range(3)]
P2 = [[100*b+a for a in DP2[i]]for i in range(3)]

fig = plt.figure()
ax = plt.axes()
plt.scatter(0.001*vol,100*b, color="blue",label=r"$b$")
plt.scatter(V,P1,color="red",label = r"$b+\Delta p_1$")
plt.scatter(V,P2,color="green",label=r"$b+\Delta p_2$")
ax.legend(loc='upper left', frameon=True)
plt.xlabel("Volumen / $m^3$")
plt.ylabel("Druck / Pa")
plt.show()
fig.savefig("P_V.png")

x =[a for a in mittel_cpcv]
for a in mittel_cpcv2:
    x.append(a)
s = [a for a in mittel_s_cpcv]
for a in mittel_s_cpcv2:
    s.append(a)

print("Gewichtetes Mittel des Adiabatenexponenten aus B und C (Mittel,Fehler): ",gew_mittel(x,s))
