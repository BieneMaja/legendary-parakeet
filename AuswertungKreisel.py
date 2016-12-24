from uncertainties import unumpy
import numpy as np
from uncertainties import ufloat
import matplotlib.pyplot as plt
print("""Bitte gib die Werte in SI-Einheiten (m,kg,s) ein, sonst stimmen die Ergebnisse nicht.
Falls die Eingabe im Plural ist, gib die Werte wie folgt mit einem Komma getrennt ein: 0.1,1.34,2,09 .
Falls Fehler auftreten, kopier bitte die Fehlermeldung und schicke sie mir - oder versuche, es selber zu lösen :)""")
mes1 = [float(input("Raddurchmesser: "))*0.5,float(input("Radgewicht: ")),float(input("Masse des Zusatzgewichtes: ")),float(input("Raddicke: "))]
Sd = unumpy.uarray([[float(x) for x in input("Schwingungsdauer links: ").split(",")],
                    [float(x) for x in input("Schwingungsdauer rechts: ").split(",")]],
                   [float(input("Fehler: "))])
def thm(M,r,fac):
    J = fac*M*r**2
    return J
def thms(sd,m,z):
    T = sd*0.1
    g = 9.81
    J = (T**2*g*z*m)/(4*np.pi**2)-m*z**2
    return J
Thms=[[],[]]
for i in range(2):
    for j in range(len(Sd[0])):
        J = thms(Sd[i][j],mes1[2],mes1[0])
        Thms[i].append(J)
print(thm(mes1[1],mes1[0],0.5),r" = Thm horizontal ($kgm^2 s^{-1}$) = 0.5*M*r**2")
print(thm(mes1[1],mes1[0],0.25)+thm(mes1[1],0.03,1/12.),r" = Thm vertikal ($kgm^2 s^{-1}$) = 0.25*M*r**2")
print(Thms," = Thm aus Schwingungsdauer = (T**2*g*z*m)/(4*np.pi**2)-m*z**2")

m = [float(x) for x in input("Massen: ").split(",")]
mes2 = unumpy.uarray([float(input("Masse Ausgleichgewicht: ")),
                      float(input("Dicke Ausgleichgewicht: ")),
                     float(input("Durchmesser Ausgleichgewicht: ")),
                    float(input("Abstand Ausgleichgewicht zu Unterstützpunkt: ")),
                      float(input("Abstand Kerbe - Drehachse: ")),
                      float(input("Stablänge: "))],[float(x) for x in input("Fehler: ").split(",")])
a=[]
c=[]
for i in range(len(m)):
    a.append([float(x) for x in input("Rotationsperioden " + str(m[i])+" kg:").split(",")])
b=[float(input("Fehler der Lichtschranke: "))]
rotp1 = unumpy.uarray(a,b) 
for i in range(len(m)):
    c.append([float(x) for x in input("Halbe Präzessionsperioden "+ str(m[i])+" kg:").split(",")])
d=[float(input("Fehler der Stoppuhr + Reaktionszeit: "))]
prazf = unumpy.uarray(c,d)
wr = 1/(2*np.pi)*(rotp1)
wrp = [[a[i][j]/(2*np.pi) for j in range(len(a[i]))] for i in range(len(m))]
wp = [np.pi/(prazf[i][j]) for i in range(len(m)) for j in range(len(prazf[i]))]
wpp = [[np.pi/c[i][j] for j in range(len(c[i]))] for i in range(len(m))]
wppf = [[prazf[i][j].s for j in range(len(prazf[i]))] for i in range(len(m))]
wrpf = [[rotp1[i][j].std_dev for j in range(len(rotp1[i]))] for i in range(len(m))]

f = plt.figure()
ax = plt.axes()
for i in range(len(m)):
    plt.errorbar(wpp[i],wrp[i], yerr=wppf[i],xerr=wrpf[i], fmt='.',label=str(m[i])+" kg");
#plt.axis([0,3.5,0,0.05]) # ([x_Achsenstart, x-Achsenstopp, y-Achsenstart, y-Achsenstopp])
plt.ylabel(r"$ \frac{1}{\omega_R} /s$")
plt.xlabel(r"$\omega_P / s^{-1}$")
ax.legend(loc="upper left",frameon=True)
plt.savefig("wpwr.png")
plt.show()

def thm3(r,m,wr,wp):
    g = 9.81
    J = (r*m*g)/(wr*wp)
    return J
for i in range(len(m)):
    print(thm3(mes2[4],m[i],np.mean(wr[i]),np.mean(wp[i]))," = Thm für "+str(m[i])+" = (r*m*g)/(wr*wp)")

rotpp = [[float(x) for x in input("Rotationsperioden "+str(m[i])+" kg: ").split(",")] for i in range(len(m))]
nutpp = [[float(x) for x in input("10*Nutationsperioden "+str(m[i])+" kg: ").split(",")] for i in range(len(m))]
for i in range(len(m)):
    nutpp[i] = [x*0.1 for x in nutpp[i]]
f = plt.figure()
ax = plt.axes()
for i in range(len(m)):
    plt.errorbar(rotpp[i],nutpp[i], yerr=0.001,xerr=0.001, fmt='.',label=str(m[i])+" kg");
plt.ylabel(r"$\omega_N /s^{-1}$")
plt.xlabel(r"$\omega_R / s^{-1}$")
ax.legend(loc="upper left",frameon=True)
plt.savefig("wnwr.png")
plt.show()
print("""Hier habe ich als Fehler (für x und y) 0.001 eingestellt. Das war geschätzt...
Die Plots wurden unter den Namen 'wpwr.png' und 'wnwr.png' gespeichert. Falls du das alles mehrmals gemacht hast,
ist wahrscheinlich nur die neueste Version gespeichert und der Rest überschrieben.
Viel Spaß beim Suchen und ich hoffe, es war hilfreich :)""")
