Um die Fluideigenschaften von CO2 (Dichte, Viskosität, Enthalpie) als 
konstitutive Beziehungen in MUFTE_UG zu erhalten, müssen die folgenden
Schritte durchgeführt werden:

1. Ermittlung des für das gestellte Problem notwenigen Druck- und
   Temperaturbereiches und der nötigen Auflösung der Stützstellen für
   Druck und Temperatur.
2. Anwendung der Beziehung von Span & Wagner (in diesem Verzeichnis bzw. in
   pml_CO2): 
   das Programm "extractproperties" muss entsprechend des p- und T-Bereiches und der 
   Auflösung angepasst werden. Die Ausführung von "extractproperties" erzeugt den 
   Datensatz (datensatz.dat) mit Dichte, Enthalpie, und Viskosität.  
3. Für die Erzeugung der c-Funktionen aus dem Datensatz werden PERL-Scripts
   verwendet. Das Skript data.prl nimmt den Datensatz datensatz.dat und 
   erzeugt die Ausgabe. Im Pearl-Skript muss lediglich der Parameter $which
   angepasst werden (1=density, 2=enthalpy, 3=solubility) um die entsprechende 
   Funktion zu schreiben.
4. Überprüfen ob Datenarrays richtig geschriebn wurden (manchmal werden 2 Zahlen 
   zusammen geschrieben, obwohl ein Komma+Lehrzeichen dazwischen gehört ??) 
5. Ein Test der neuen Funktion kann im Verzeichnis ./TEST durchgeführt werden.
   Hierzu muss die erzeugte Funktion density.c an die Datei testconstrel.c 
   angehängt werden (cat density.c >> TEST/testconstrel.c). Vorsicht: ist 
   die alte Funktion von density.c noch in testconstrel.c enthalten?
6. Die Ausgabe density.c kann nun an unser constrel-File angehängt werden, z.B.:
	cat density.c >> pml/constrel_CO2.c
7. Die alte Funktion im constrel muss dann natürlich gelöscht werden. Falls
   es sich um eine zusätzliche Funktion handelt, muss diese im constrel_CO2.h
   deklariert werden.


Dieses Vorgehen ermöglicht es, den für das gestellte Problem notwendigen 
Druck- und Temperaturbereich vorzugeben, schnell in das Programm zu
implementieren und somit ein aufwendiges "mitschleppen" zusätzlicher Werte,
die nicht benötigt werden, zu vermeiden. Dasselbe gilt auch für die Wahl der
Stützstellen.


Hinweise
--------

- Die Datei datensatz.dat muss das folgende Format haben, da andernfalls
  Probleme bei der Umformatierung im Pearl-Script auftreten können:
	
	Druck 		 Temperatur	 Dichte 	  Enthalpie	Viskosität	  

	1.000000e+05     2.830000e+02    1.881498e+00    -1.370482e+01  1.0E-5
	1.000000e+05     2.840000e+02    1.874743e+00    -1.286816e+01  2.0E-5

  Die Pearl-Scripts geben keine Warnungen oder Fehlermeldungen aus, wenn
  bei der Sortierung etwas schiefgehen sollte, deswegen sollte das Ergebnis
  auch nochmal überprüft werden.


Historie und Fehler
-------------------

- Ursprünglich wollte ich die Auflösung der p-, T-Stützstellen für den
  kritischen Bereich feiner wählen als in den Gebieten, in denen nicht 
  so viel passiert. Da bei der Interpolation zwischen zwei Drücken und
  zwei Temperaturen gemittelt wird, kam es an dem Übergang zwischen den
  unterschiedlichen Funktionen zu Unstetigkeiten. Deshalb entschied ich
  mich erstmal für eine gleichbleibende Auflösung über den gesamten
  p-, T-Bereich.

- Ursprünglich habe ich die Felder folgendermaßen deklariert und definiert: 
  
	Deklaration:

	DOUBLE T[41];
	DOUBLE p[110];
	DOUBLE rho[41][110];

	Definition:

	p[0] = 1.000000e+05;
	T[0] = 2.830000e+02;
	rho[0][0] = 1.881498e+00;

  Aus programmiertechnischen Gründen habe ich es dann so versucht:
	
	DOUBLE T[41]={273.15,283.15,...};
	DOUBLE p[110]={1.0E5,2.0E5,...};
	DOUBLE rho[41][110]={1.88e+00, 3.78e+00, ..}{...};

