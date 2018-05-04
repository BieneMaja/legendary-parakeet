import urllib.request #um Zeug aus dem Internet zu fischen
from bs4 import BeautifulSoup # um aus html einen schönen str zu machen

daten = open("Gibtsneues.txt", mode="r+") # "r+" heißt lesen und schreiben
inhalt = daten.readlines() # Liste aus str, die jeweils in einer Zeile stehen
# In der Zeile steht jeweils, wie oft "Blatt", "Skript" oder ähnliches das letzte Mal auf der Seite stand.

seiten = ["https://www.physik.uni-muenchen.de/lehre/vorlesungen/sose_18/T3_-Elektrodynamik/uebungen/index.html",
         "https://www.physik.uni-muenchen.de/lehre/vorlesungen/sose_18/E4-E4p-Atom-und-Molekuelphysik/uebungen/index.html",
         "https://www.physik.uni-muenchen.de/lehre/vorlesungen/sose_18/Quantum-Field-Theory/tutorial/index.html",
         "https://www.mathematik.uni-muenchen.de/~leidl/NuMaSoSe2018/"]

neues = [0 for i in range(len(seiten))]
suchen = [ "Übungsblatt", "Übungsblatt", "Exercise sheet","Blatt"]
wort = ["Blatt", "Blatt", "Blatt", "Blatt"]
fach = ["T3","E4","QFT","M4"]

for i in range(len(seiten)):
    html_doc1 = urllib.request.urlopen(seiten[i]) #holt sich html Dokument aus seiten
    soup1 = BeautifulSoup(html_doc1, 'html.parser')
    html_doc1.close()
    
    a = soup1.get_text() # macht einen schönen str daraus

    neues[i] = a.count(suchen[i])
        
    if neues[i] > int(inhalt[i]): # also wenn es öfter auf der Seite steht als das letzte Mal
        print("Es gibt ein neues {} in {}.\n siehe {}".format(wort[i], fach[i], seiten[i]))
        inhalt[i] = "".join([str(neues[i]), "\n"])  # Zahl der vorkommenden Wörter wird erneuert
    else:
        print("Es gibt kein neues %s in %s." % (wort[i], fach[i]))
    
daten.seek(0) # damit der Pointer auf den Anfang der Datei zeigt und die Zeilen überschrieben werden
daten.writelines(''.join(inhalt)) # Neue Zahlen werden in die txt Datei geschrieben

daten.close()
