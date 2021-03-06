import matplotlib.pyplot as plt
plt.style.use("seaborn-whitegrid")
import numpy as np
from uncertainties import unumpy

#Aufnahme der Messdaten
Radius_kr = [float(x) for x in [input("Radius der Kugel: "),
                                            input("Radius des Zylinders: "),
                                            input("Radius der Scheibe: "),
                                            input("Radius des Hohlzylinders (innen): "),
                                            input("Radius des Hohlzylinders (außen): "),
                                            input("Abstand der Hanteln: "),
                                            input("Kantenlänge des Würfels: "),
                                            input("Länge des Stabes: "),]]
Radius_k = unumpy.uarray(Radius_kr,
                         [float(x) for x in input("Fehler: ").split(",")])
Masse_kr = [float(x) for x in [input("Masse der Kugel: "),
                                            input("Masse des Zylinders: "),
                                            input("Masse der Scheibe: "),
                                            input("Masse des Hohlzylinders: "),
                                            input("Masse des Hohlzylinders: "),
                                            input("Masse der Hanteln: "),
                                            input("Masse des Würfels: "),
                                            input("Masse des Stabes: "),]]
Masse_k = unumpy.uarray(Masse_kr,
                         [float(x) for x in input("Fehler: ").split(",")])
Jfaktor = [2./5., 1./2., 1./2., 1./2. ,1./2., 1./3. , 1./6., 1./12. ]

Masse_l = [float(x) for x in input("Massen links: ").split(",")]
winkel_lp = [float(x) for x in input("Winkelausschlag (°) links: ")]
winkel_l = unumpy.uarray(winkel_lp,[float(input("Fehler: "))])
Masse_r = [float(x) for x in input("Massen rechts: ").split(",")]
winkel_rp = [float(x) for x in input("Winkelausschlag (°) rechts: ")]
winkel_r = unumpy.uarray(winkel_lp,[float(input("Fehler: "))])

winkel_lp = [winkel_lp[i]/180 *np.pi for i in range(len(winkel_lp))]
winkel_rp = [winkel_rp[i]/180 *np.pi for i in range(len(winkel_rp))]

def Drehmoment(m,r):
    d = 9.81*m*r
    return d
Dl=[Drehmoment(Masse_l[i],Radius_k[2]) for i in range(len(Masse_l))]
Dlp = [x.nominal_value for x in Dl]
Dr = [Drehmoment(Masse_r[i],Radius_k[2]) for i in range(len(Masse_r))]
Drp = [x.nominal_value for x in Dr]
Fehler = [[x.s for x in Dl],[x.s for x in Dr]]
x= np.linspace(0,3.5, 1000)
m,b = np.polyfit(winkel_lp,Dlp,1)
m2,b2 = np.polyfit(winkel_rp,Drp,1)
fig = plt.figure()
ax = plt.axes()
plt.errorbar(winkel_lp,Dlp, xerr=0.01745,yerr=Fehler[0] ,fmt='.',label="links");
plt.errorbar(winkel_rp,Drp, xerr=0.01745,yerr=Fehler[1] ,fmt='.',label="rechts");
plt.plot(x, m*x + b, linestyle='solid',label='Regression links');
plt.plot(x, m2*x + b2, linestyle='solid',label='Regression rechts');
#plt.axis([0,3.5,0,0.05]) # ([x_Achsenstart, x-Achsenstopp, y-Achsenstart, y-Achsenstopp])
#plt.title("Winkelausschlag mit verschiedenem angreifenden Drehmoment")
plt.ylabel("Drehmoment ($kg^2\cdot m^2\cdot s^{-2}$)")
plt.xlabel("Winkelausschlag/$\pi$")
ax.legend(loc="upper left",frameon=True)
#plt.savefig("winkelreg.png")
plt.show()
print("Rückstellmoment links durch Regression: "+str(m) + " Einheit?")
print("Rückstellmoment rechts durch Regression: "+str(m2) + " Einheit?")

Sd_w = unumpy.uarray([[0.1*float(x) for x in input("Schwingungsdauer der Kugel: ").split(",")],
                      [0.1*float(x) for x in input("Schwingungsdauer des Zylinders: ").split(",")],
                      [0.1*float(x) for x in input("Schwingungsdauer der Scheibe: ").split(",")],
                      [0.1*float(x) for x in input("Schwingungsdauer des Hohlzylinders: ").split(",")],
                      [0.1*float(x) for x in input("Schwingungsdauer der Hanteln: ").split(",")],
                      [0.1*float(x) for x in input("Schwingungsdauer des Würfels (Mitte): ").split(",")],
                      [0.1*float(x) for x in input("Schwingungsdauer des Würfels (Ecke): ").split(",")],
                      [0.1*float(x) for x in input("Schwingungsdauer des Stabes im Schwerpunkt: ").split(",")],
                      [0.1*float(x) for x in input("Schwingungsdauer des Stabes parrallel dazu: ").split(",")]],
                     [0.1*float(input("Fehler: "))])
#Berechnetes Thm
def Teta(Masse,Lange,Faktor):
    J = Faktor * Masse * (Lange)**2
    return J
Theta = []
for i in range(8):
    if i == 3:
        Theta.append(Teta(Masse_kr[i+1],Radius_kr[i+1],Jfaktor[i+1])+Teta(Masse_kr[i],Radius_kr[i],Jfaktor[i]))
    elif i == 4:
        Theta.append(":)")
    elif i ==  5:
        Theta.append(Teta(Masse_kr[i]/3.,Radius_kr[i],Jfaktor[i])+2*Masse_kr[i]/3.*(Radius_kr[i]/2.)**2)
    else:
        Theta.append(Teta(Masse_kr[i],Radius_kr[i],Jfaktor[i]))
print("""Thm ausgerechnet, Reihenfolge: Kugel, Zylinder, Scheibe, Hohlzylinder,:), Hanteln, Würfel, Stab.
      Der Würfel hat mit beiden Drehachsen das gleiche Trägheitsmoment und 
      mit dem Steinerschen Satz muss das Trägheitsmoment des Stabes mit der verschobenen Drehachse ausgerechnet werden.""")
print(Theta , "kg * m**2")
#Trägheitsmoment aus der Schwingdauer, man rechnet mit dem Mittelwert des Rückstellmoments W weiter.
Wl = [m,m2]
W = np.mean(Wl)
def Teta_schwingung(W,T):
    J = W*(T**2/((2*np.pi)**2))
    return(J)

Thm_schw = [Teta_schwingung(W,np.mean(Sd_w[i])) for i in range(len(Sd_w))]
print("""Trägheitsmoment aus der Schwingung berechnet (Reihenfolge wie in der Eingabe) in kg * m**2:\n""",Thm_schw)

Sd_Tisch = unumpy.uarray([[float(x) for x in input("1) 10*Schwingungsdauer des Tischchens (13 Werte, Abstand 15°): ").split(",")],
                          [float(x) for x in input("2) 10*Schwingungsdauer des Tischchens (13 Werte, Abstand 15°): ").split(",")]],
                         [float(input("Fehler: "))])
Grad_Tisch = [0,15,30,45,60,75,90,105,120,135,150,165,180]
Thm_tisch1 =[Teta_schwingung(W,Sd_Tisch[0][i]) for i in range(len(Grad_Tisch))]
Thm_tisch2 =[Teta_schwingung(W,Sd_Tisch[1][i]) for i in range(len(Grad_Tisch))]
Thm_Tisch =[[Thm_tisch1[i],Thm_tisch2[i]] for i in range(len(Grad_Tisch))]
Thm_Tischm =[np.mean(Thm_Tisch[i]) for i in range(len(Grad_Tisch))]

#reziproke Quadratwurzeln
Grad_Tischp=[]
rq=[]
for i in range(len(Grad_Tisch)):
    rq.append(1/(np.sqrt(Thm_Tischm[i].nominal_value)))
    Grad_Tischp.append(Grad_Tisch[i]/180. *np.pi)
for i in range(len(Grad_Tisch)):
    rq.append(1/(np.sqrt(Thm_Tischm[i].nominal_value)))
    Grad_Tischp.append(np.pi+Grad_Tisch[i]/180. *np.pi)
print("Reziproke Quadratwurzeln: ",rq)
print("Maximum: ",max(rq),"Minimum: ",min(rq))

fig = plt.figure()
ax = plt.subplot(111, projection='polar')
ax.plot(Grad_Tischp, rq, color='r', linewidth=1)
ax.set_rmax(70)
ax.grid(True)
#plt.savefig('Trägheitsmoment-Ellipse.png')
plt.show()

print("Teil B")
#Reihenfolge 0 Radius Felge, 1 Radius Rad(Faden), 2 Masse Gewicht, 3 Abstand Gewicht-Drehachse
Messungen = unumpy.uarray([float(x) for x in [input("Radius Felge: "),
                                             input("Radius kleines Rad: "),
                                             input("Masse des Gewichtes (physikalisches Pendel): "),
                                             input("Abstand Gewicht zur Drehachse: ")]],
                          [float(input("Fehler: "))])

Massen = [float(x) for x in input("Angehängte Massen: ").split(",")]
abstand = [[float(x) for x in input("Abstand der Sekundenmarken " +str(Massen[i])+" kg: ").split(",")]for i in range(len(Massen))]
abstandfehler = [float(x) for x in input("Fehler (4 Werte erwartet): ").split(",")]
#abstand_mit_fehler = unumpy.uarray(abstand,abstandfehler)
ax = plt.axes()
fig = plt.figure()
beschl = []
for i in range(4):
    x = np.linspace(0,len(abstand[i])-1,len(abstand[i]))
    m,b = np.polyfit(x,abstand[i],1)
    beschl.append(m)
    plt.errorbar(x, abstand[i], yerr=abstandfehler[i], fmt='.', label=str(i));
    plt.plot(x, m*x + b, linestyle='solid',label='Regression '+str(i));
#    plt.axis([-0.5,len(abstand_p[i]),0,numpy.max(abstand_p[i])+1])
    plt.xlabel("Zeit (s)")
    plt.ylabel("Abstand der Markierungen "+str(i+1))
    ax.legend(loc="upper left",frameon=True)
    #plt.savefig("abstandreg"+str(i)+".png")
    plt.show()
print("Beschleunigung des großen Rades durch Regression: ",beschl,"m*s**-1")
#a_klein = [beschl[i]*Messungen[1]/Messungen[0] for i in range(len(Massen))]
#print("Beschleunigung des kleinen Rades: ",a,"m*s**-1")
J=[(Massen[i]*9.81*Messungen[1]*Messungen[0])/beschl[i] - Massen[i]*Messungen[1]**2]
print("Thm aus der Beschleunigung berechnet: ",J,"m*s**-1")
print("Mittelwert des Thm: ",np.mean(J))

# Aus der Schwingungsdauer Thm berechnen: 
Sd_v3 = unumpy.uarray([float(x) for x in input("Schwingungsdauer links (2x), rechts (2x): ").split(",")],[float(input("Fehler: "))])
abstand_end1 = []
def Thm(Sd,m,z):
    T=Sd*0.2
    J = (T**2*9.81*z*m)/(4*np.pi**2) - m*z**2
    return J
ThetaB2 = [Thm(Sd_v3[i],Messungen[2],Messungen[3]) for i in range(4)]
u=np.mean(ThetaB2[0:2])
v=np.mean(ThetaB2[2:])
w=np.mean(ThetaB2[4:])
print("Thm links: ",ThetaB2[:2],"Thm rechts: ",ThetaB2[2:4],
      "Mittelwert links: ",u, "Mittelwert rechts: ",v,
      "Mittelwert gesamt: ",w)
