![HSWT Logo](https://upload.wikimedia.org/wikipedia/commons/thumb/8/8b/HSWT_Logo_gruen.png/960px-HSWT_Logo_gruen.png)
## Modul: Programmierung für Datenanalyse, Bildverarbeitung und Simulation 
## Projekt 8: Modellierung von Fermentationsdaten - Rühmann M., Stehno C.
## Inhaltsverzeichnis
- [Projektbeschreibung](#projektbeschreibung)
- [Programmablauf und -architektur](#programmablauf-und--architektur)
- [Installation und Verwendung](#installation-und-verwendung)
- [Modellschwächen und Limitationen](#modellschwächen-und-limitationen)
- [Wissenschaftliche Quellen](#wissenschaftliche-quellen)
- [Nutzung von Sprachmodellen](#nutzung-von-sprachmodellen)
- [Autoren](#autoren)

## Projektbeschreibung

Dieser Chinese Hamster Ovary (CHO) Fermentations-Simulator bietet eine webbasierte Anwendung zur **Simulation**, **Analyse** und **Visualisierung** eines **Batch-Fermentationsprozesses** in einem **Bioreaktor**. Die zugrunde liegende **numerische Modellierung** bildet den Fermentationsverlauf über **12 Tage (288 Stunden)** realitätsnah ab und ermöglicht die systematische Untersuchung von Prozessparametern auf die Antikörperproduktion.

Der Simulator basiert auf **wissenschaftlich fundierten kinetischen Modellen** und **empirischen Daten**, die im Rahmen eines **Zellkultur-Praktikums im Wintersemester** gewonnen wurden. Dabei wurden ausgewählte **Parameter direkt aus experimentellen Ergebnissen** abgeleitet. Ziel des zugrundeliegenden Praktikums war die Produktion eines **monoklonalen Antikörpers** gegen den Insulin-like Growth Factor 1-Rezeptor (IGF-1R) mithilfe der **CHO-Zelllinie DG44** [1].

## Programmablauf und -architektur

### Kernklasse: CHOFermentationSimulator
- **Initialisierung:** Grundbedingungen für Simulation werden eingestellt
- **Simulation:** Euler-Integration für stündliche Berechnung von Zellwachstum, Substratverbrauch und Antikörperproduktion
- **Stress-Modellierung:** Gauß-Funktionen für Umweltstress (Temperatur, pH, Sauerstoff, Glukose) [2,3]
- **Substratkinetik:** Monod-Gleichungen mit optionaler Haldane-Inhibition [4,5]


### Interaktives Dashboard
- **Sidebar-Parametereingabe**
   - Programm startet mit optimalen Prozessparametern [2,3]
- **Drei-Tab-Layout:** Zeitverläufe, Daten, Ergebnisse
   - **Zeitverläufe (Tab 1):**
      - Zellwachstum (viable/tote Zellen) und Viabilität
      - Substratverbrauch vs. Zelldichte
      - Antikörperproduktion über Zeit
      - Automatische KPI-Berechnung
   - **Daten (Tab 2):**
      - Darstellung der Rohdaten
   - **Ergebnisse (Tab 3):**
      - Tabellatische Darstellung der KPIs und der Simulations-Parameter
      - Zusammenhänge werden mittels Korrelationsmatrix dargestellt 

## Installation und Verwendung

### Voraussetzungen
```bash
pip install streamlit numpy pandas matplotlib seaborn
```

**Bibliotheken-Übersicht:**
- [Streamlit](https://docs.streamlit.io/) – Web-App-Framework für interaktive Anwendungen
- [NumPy](https://numpy.org/doc/) – Numerische Berechnungen und Array-Operationen
- [Pandas](https://pandas.pydata.org/docs/) – Datenanalyse und Tabellenverwaltung
- [Matplotlib](https://matplotlib.org/stable/contents.html) – Wissenschaftliche Visualisierung
- [Seaborn](https://seaborn.pydata.org/) – Statistische Datenvisualisierung

### Anwendung starten
```bash
streamlit run CHOFermentation.py
```

### Bedienungsanleitung

1. **Parameter einstellen** (Sidebar)
2. **Simulation starten**
3. **Ergebnisse analysieren**
   - **Tab "Zeitverläufe"**: Grafische Darstellung mit KPI-Übersicht
   - **Tab "Daten"**: Tabellarische Rohdaten (6-Stunden-Intervalle)
   - **Tab "Ergebnisse"**: Multi-Simulationsvergleich mit Korrelationsanalyse

## Modellschwächen und Limitationen

- **Vereinfachte Stoffwechselwege**: Laktatproduktion und komplexe Metabolitinteraktionen nicht berücksichtigt [6]
- **Konstante kinetische Parameter**: Temperatur- und pH-Abhängigkeit der Wachstumsraten nicht dynamisch modelliert
- **Fehlende Nährstofflimitierungen**: Keine Berücksichtigung von Aminosäure-, Vitamin- oder Spurenelementmangel
- **Statische Stress-Parameter**: Zelladaptation und Toleranzveränderungen über Zeit nicht berücksichtigt
- **Osmotischer Stress vereinfacht**: Glukose-Stress-Modell erfasst komplexe osmotische Regulationsmechanismen nur näherungsweise


## Wissenschaftliche Quellen

[1] Augustin, I. (WS 2024/2025). Praktikum Zellkulturen in der Biotechnologie. Skript zur Produktion eines humanen monoklonalen Antikörpers. Hochschule Weihenstephan-Triesdorf.

[2] Ritacco, F.V., Wu, Y., Khetan, A. (2018). Cell culture media for recombinant protein expression in Chinese hamster ovary (CHO) cells. *Biotechnology and Bioengineering*, 115(6), 1493-1507.

[3] Trummer, E., et al. (2006). Process parameter shifting: Part I. Effect of DOT, pH, and temperature on the performance of Epo-Fc expressing CHO cells cultivated in controlled batch bioreactors. *Biotechnology and Bioengineering*, 94(6), 1073-1088.

[4] Monod, J. (1949). The growth of bacterial cultures. *Annual Review of Microbiology*, 3, 371-394.

[5] Haldane, J.B.S. (1930). Enzymes. Longmans, Green and Co., London.

[6] Lao, M.S., Toth, D. (1997). Effects of ammonium and lactate on growth and metabolism of a recombinant Chinese hamster ovary cell culture. *Biotechnology Progress*, 13(5), 688-691.


## Nutzung von Sprachmodellen

Bei der Entwicklung der Anwendung wurden KI-gestützte Sprachmodelle zur Ideenfindung, Codeentwicklung und Dokumentationsunterstützung verwendet.


[Claude.ai](https://claude.ai/)

abgerufen am 03.06.2025:
1. Erstelle ein Python-Skript zur Simulation einer CHO-Zellfermentation im Batch-Verfahren.
2. Welche kinetischen und numerischen Modelle sind am besten geeignet, um das Wachstum von CHO-Zellen in einem Batch-Fermentationsprozess zu beschreiben?
3. Wie lassen sich Umwelteinflüsse wie Temperatur, pH und Sauerstoff mathematisch modellieren, um realistische Bedingungen für die CHO-Zellfermentation zu simulieren?
4. Wie können optimale Betriebsbedingungen für CHO-Zellkulturen durch mathematische Funktionen beschrieben werden?

abgerufen am 06.06.2025:

5. Wie lässt sich die Visualisierung der Simulationsergebnisse mit Streamlit umsetzen? Implementiere diese Funktion in unser Skript.
6. Welche Ansätze gibt es, um mehrere Messwerte mit unterschiedlichen Skalierungen gleichzeitig in einem Diagramm darzustellen?
7. Wie kann eine Sekundärachse in einem Diagramm implementiert werden, um unterschiedliche Messwerte darzustellen? Implementiere diese Funktion in unser Skript.

abgerufen am 12.06.2025:

8. Wie können Key Performance Indicators (KPIs) neben den Diagrammen in der Visualisierung angezeigt werden? Implementiere diese Funktion in unser Skript.
9. Wie lassen sich Parameter-Sensitivitätsanalysen in Simulationsstudien implementieren und visualisieren?
10. Wie kann eine Korrelations-Heatmap aus KPIs und Prozessparametern erstellt werden? Implementiere diese Funktion in unser Skript.
11. Entspricht der erstellte Programmcode dem PEP 8-Formatierungsstil?


## Autoren

**Hochschule Weihenstephan-Triesdorf**
- Milena Rühmann (1389270)
- Christian Stehno (1386026)

*Entwickelt im Rahmen des Masterstudiengangs Biotechnologie an der Hochschule Weihenstephan-Triesdorf*
