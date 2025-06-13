##### Hochschule Weihenstephan-Triesdorf
##### Programmierung für Datenanalyse, Bildverarbeitung und Simulation 
##### Betreuerin: Prof. Dr. Kristina Eisen

#### Projekt 8: Modellierung von Fermentationsdaten 

### Ersteller/-in: Milena Rühmann, Christian Stehno
### Datum: 13.06.2025

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import streamlit as st


class CHOFermentationSimulator:
    """Simuliert CHO-Zellfermentation mit Monod-Kinetik."""
    
    def __init__(self, duration=288, time_step=1.0):
        """Initialisiert den Fermentationssimulator.
        
        Args:
            duration: Gesamtdauer der Simulation in Stunden.
            time_step: Zeitschritt in Stunden.
        """
        # Zeitparameter für die Simulation setzen (288h = 12 Tage Standard)
        self.duration = duration
        self.time_step = time_step
        self.time_points = np.arange(0, duration + time_step, time_step)
        
        # Kinetische Parameter für das Monod-Modell (Wachstum und Substrat)
        self.max_growth_rate = 0.035
        self.substrate_affinity = 2.0
        self.yield_coefficient = 0.4  # Biomasse/Substrat-Verhältnis
        self.maintenance_coefficient = 0.02  # Erhaltungsstoffwechsel [1/h]
        self.death_rate = 0.005  # Grundsterberate [1/h]
        
        # Produktbildungsparameter (Antikörper und Inhibierung)
        self.antibody_productivity = 0.8
        self.inhibition_constant = 50.0
        
        # Optimale Kultivierungsbedingungen für CHO-Zellen
        self.optimal = {'temp': 37.0, 'ph': 7.2, 'oxygen': 50.0,
                        'glucose': 20}
        
        # Toleranzbereiche (Standardabweichungen) für jeden Parameter
        self.sigma = {'temp': 2.0, 'ph': 0.4, 'oxygen': 20.0,
                      'glucose': 40.0}
    
    def stress_gaussian(self, deviation, sigma):
        """Berechnet die Normalverteilung für die Umgebungsbedingungen."""
        return np.exp(-0.5 * (deviation / sigma)**2)
    
    def evaluate_environment(self, temperature, ph, oxygen_saturation,
                             glucose):
        """Bewertet Umgebungsbedingungen und berechnet Effekte.
        
        Berechnet sowohl positive Aktivitätseffekte als auch
        negative Stresseffekte auf die Zellen.
        """
        # Abweichungen von optimalen Bedingungen berechnen
        deviations = {
            'temp': temperature - self.optimal['temp'],
            'ph': ph - self.optimal['ph'],
            'oxygen': oxygen_saturation - self.optimal['oxygen'],
            'glucose': glucose - self.optimal['glucose']
        }
        
        # Berechnet multiplikativen Aktivitätseffekt auf Wachstum
        activity_factors = [self.stress_gaussian(dev, self.sigma[key])
                           for key, dev in deviations.items()]
        activity_effect = np.prod(activity_factors)
        
        # Berechnet Stresseffekt für erhöhte Sterberate
        stress_factors = [1.0 / self.stress_gaussian(dev, self.sigma[key])
                         for key, dev in deviations.items()]
        combined_stress = np.prod(stress_factors)
        
        return activity_effect, self.death_rate * combined_stress
    
    def monod_kinetics(self, substrate, k_s, k_i=None):
        """Implementiert Monod- oder Haldane-Kinetik.
        
        Monod: Sättigung bei hohen Substratkonzentrationen
        Haldane: Zusätzliche Inhibierung bei sehr hohen Konzentrationen
        """
        # Schutz vor negativen Substratkonzentrationen
        if substrate <= 0:
            return 0
            
        # Einfache Monod-Kinetik ohne Substratinhibierung
        if k_i is None:
            return substrate / (k_s + substrate)
            
        # Haldane-Kinetik mit Substratinhibierung bei hohen Konzentrationen
        return substrate / (k_s + substrate + (substrate**2 / k_i))
    
    def simulate(self, initial_glucose, initial_vcd, temperature,
                 ph_constant, oxygen_saturation):
        """Führt die Hauptfermentationssimulation durch.
        
        Berechnet zeitliche Entwicklung aller Kultivierungsparameter
        mittels numerischer Integration (Euler-Verfahren).
        """
        n = len(self.time_points)
        
        # Initialisiert Datenarrays für alle Messgrößen
        keys = ['glucose', 'vcd', 'tcd', 'dead_cells', 'viability',
                'antibody_titer']
        data = {k: np.zeros(n) for k in keys}
        
        # Setzt Anfangswerte für t=0 (zu Beginn alle Zellen lebendig)
        data['glucose'][0] = initial_glucose
        data['vcd'][0] = data['tcd'][0] = initial_vcd
        data['viability'][0] = 100.0  # 100% Viabilität am Start
        
        # Hauptsimulationsschleife mit numerischer Integration (Euler)
        for i in range(1, n):
            dt = self.time_step
            prev_glucose = data['glucose'][i-1]
            prev_vcd = data['vcd'][i-1]
            
            # Bewertet aktuelle Umgebungsbedingungen
            activity_effect, death_rate = self.evaluate_environment(
                temperature, ph_constant, oxygen_saturation, prev_glucose)
            
            # Berechnet Substratlimitierung nach Monod/Haldane
            substrate_factor = self.monod_kinetics(
                prev_glucose, self.substrate_affinity,
                self.inhibition_constant)
            
            # Berechnet aktuelle spezifische Wachstumsrate
            growth_rate = (self.max_growth_rate * substrate_factor *
                          activity_effect if prev_vcd > 0 else 0)
            
            # Aktualisiert Zellkonzentrationen (Euler-Integration)
            vcd_growth = growth_rate * prev_vcd * dt
            vcd_death = death_rate * prev_vcd * dt
            
            # Neue lebende Zellzahl (nicht unter 0)
            vcd_current = max(0, prev_vcd + vcd_growth - vcd_death)
            data['vcd'][i] = vcd_current
            
            # Akkumuliert tote Zellen (sterben nicht ab)
            data['dead_cells'][i] = data['dead_cells'][i-1] + vcd_death
            data['tcd'][i] = vcd_current + data['dead_cells'][i]
            
            # Berechnet Viabilität als Anteil lebender Zellen
            data['viability'][i] = (vcd_current / data['tcd'][i] * 100 
                                   if data['tcd'][i] > 0 else 0)
            
            # Aktualisiert Glukoseverbrauch für Wachstum + Erhaltung
            glucose_consumption = (vcd_growth / self.yield_coefficient +
                                 self.maintenance_coefficient * prev_vcd * dt)
            data['glucose'][i] = max(0, prev_glucose - glucose_consumption)
            
            # Berechnet Antikörperproduktion durch lebende Zellen
            if vcd_current > 0:
                antibody_production = (self.antibody_productivity *
                                     vcd_current * dt * activity_effect)
                data['antibody_titer'][i] = (
                    data['antibody_titer'][i-1] + antibody_production)
            else:
                # Keine Produktion ohne lebende Zellen
                data['antibody_titer'][i] = data['antibody_titer'][i-1]
        
        # Erstellt strukturierten DataFrame mit gerundeten Werten
        return pd.DataFrame({
            'Zeit (h)': np.round(self.time_points, 2),
            'Glukose (g/L)': np.round(data['glucose'], 2),
            'VCD (10^6 Zellen/mL)': np.round(data['vcd'], 2),
            'TCD (10^6 Zellen/mL)': np.round(data['tcd'], 2),
            'Viabilität (%)': np.round(data['viability'], 2),
            'Antikörper-Titer (µg/mL)': np.round(data['antibody_titer'], 2),
            'Temperatur (°C)': np.full(n, temperature),
            'pH': np.full(n, ph_constant),
            'Sauerstoff (%)': np.full(n, oxygen_saturation),
        })


def combine_legends(ax1, ax2, loc='upper right'):
    """Hilfsfunktion zum Kombinieren von Legenden zweier Achsen.
    
    Vereint Legenden von primärer und sekundärer y-Achse
    für übersichtliche Darstellung in einem Plot.
    """
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc=loc)


def create_plots(data):
    """Erstellt wissenschaftliche Visualisierungen der Fermentationsdaten.
    
    Generiert drei subplot-Grafiken für verschiedene Aspekte:
    1. Zellwachstum und Viabilität, 2. Substratverbrauch, 3. Produktbildung
    """
    # Konvertiert Stunden in Tage für bessere Lesbarkeit der Achsen
    time_days = data['Zeit (h)'] / 24.0
    fig, axes = plt.subplots(3, 1, figsize=(12, 10))

    # Plot 1: Zellwachstum und Viabilität über Zeit
    axes[0].plot(time_days, data['VCD (10^6 Zellen/mL)'], 'g-',
                 label='VCD', linewidth=2)
    axes[0].plot(time_days, data['TCD (10^6 Zellen/mL)'], 'orange',
                 label='TCD', linewidth=2)
    
    # Sekundäre y-Achse für Viabilität (unterschiedliche Einheit)
    ax0_twin = axes[0].twinx()
    ax0_twin.plot(time_days, data['Viabilität (%)'], 'purple',
                  label='Viabilität (%)', linewidth=2)
    
    # Achsenbeschriftungen und Formatierung
    axes[0].set_ylabel('Zellkonzentration (10⁶ Zellen/mL)')
    ax0_twin.set_ylabel('Viabilität (%)')
    axes[0].set_title('Zellwachstum und Viabilität')
    axes[0].set_xlabel('Zeit (Tage)')
    axes[0].set_xlim(0, 12)  # 12 Tage Darstellung
    axes[0].grid(True, alpha=0.3)
    combine_legends(axes[0], ax0_twin)

    # Plot 2: Glukoseverbrauch gekoppelt mit Zellwachstum
    axes[1].plot(time_days, data['Glukose (g/L)'], 'b-',
                 label='Glukose', linewidth=2)
    ax1_twin = axes[1].twinx()
    ax1_twin.plot(time_days, data['VCD (10^6 Zellen/mL)'], 'g--',
                  label='VCD', linewidth=1.5, alpha=0.7)
    
    # Zweite y-Achse für VCD zur Korrelationsdarstellung
    axes[1].set_ylabel('Glukose (g/L)')
    ax1_twin.set_ylabel('VCD (10⁶ Zellen/mL)')
    axes[1].set_title('Substratverbrauch')
    axes[1].set_xlabel('Zeit (Tage)')
    axes[1].set_xlim(0, 12)
    axes[1].grid(True, alpha=0.3)
    combine_legends(axes[1], ax1_twin)

    # Plot 3: Antikörperproduktion als Hauptzielgröße
    axes[2].plot(time_days, data['Antikörper-Titer (µg/mL)'], 'r-',
                 label='Antikörper-Titer', linewidth=2)
    axes[2].set_ylabel('Antikörper-Titer (µg/mL)')
    axes[2].set_title('Antikörperproduktion')
    axes[2].set_xlabel('Zeit (Tage)')
    axes[2].set_xlim(0, 12)
    axes[2].legend()
    axes[2].grid(True, alpha=0.3)

    # Optimiert Layout um Überlappungen zu vermeiden
    plt.tight_layout()
    return fig


def nummer_eingabe(label, value, fmt="%.2f", step=0.1, min_value=0.0):
    """Hilfsfunktion für einheitliche Streamlit number_input Erstellung.
    
    Standardisiert Parametereingabe-Widgets mit einheitlicher Formatierung
    und verhindert negative Werte durch min_value-Beschränkung.
    """
    return st.number_input(label, value=value, format=fmt, step=step, 
                          min_value=min_value)


def calculate_kpis(data):
    """Berechnet Key Performance Indicators aus Simulationsdaten.
    
    Extrahiert wichtigste Kennzahlen für Prozessbewertung:
    Endtiter, maximale Zelldichte, Viabilitätsstatistiken.
    """
    return {
        'final_titer': data['Antikörper-Titer (µg/mL)'].iloc[-1],
        'max_vcd': data['VCD (10^6 Zellen/mL)'].max(),
        'avg_viability': data['Viabilität (%)'].mean(),
        'min_viability': data['Viabilität (%)'].min()
    }


def main():
    """Hauptfunktion für die Streamlit-Anwendung.
    
    Steuert die gesamte Web-App: Layout, Parametereingabe,
    Simulation, Ergebnisdarstellung und Vergleichsanalysen.
    """
    # Konfiguriert Streamlit-Seite mit breitem Layout
    st.set_page_config(
        page_title="CHO-Zell Fermentations-Simulator", layout="wide")
    
    # Header-Bereich mit HSWT-Logo
    col1, col2 = st.columns([5, 3])
    with col1:
        st.title("CHO-Zell Fermentations-Simulator")
    with col2:
        logo_url = ("https://upload.wikimedia.org/wikipedia/commons/thumb/"
                   "8/8b/HSWT_Logo_gruen.png/960px-HSWT_Logo_gruen.png")
        st.image(logo_url, use_container_width=True)
    
    # Initialisiert Session State für persistente Datenhaltung
    if 'data' not in st.session_state:
        st.session_state.data = None
    if 'results_df' not in st.session_state:
        st.session_state.results_df = pd.DataFrame()
    
    # Erstellt Simulator-Instanz mit Standardparametern
    simulator = CHOFermentationSimulator()
    
    # Sidebar für strukturierte Parametereingabe
    with st.sidebar:
        st.header("Fermentationsparameter")
        
        # Startbedingungen der Kultivierung
        st.subheader("Startbedingungen")
        initial_glucose = nummer_eingabe("Anfangsglukose (g/L):", 25.00,
                                        step=1.0)
        initial_vcd = nummer_eingabe("Anfangs-VCD (10^6 Zellen/mL):", 0.50)
        
        # Umgebungsbedingungen während der Kultivierung
        st.subheader("Prozessbedingungen")
        temperature = nummer_eingabe("Temperatur (°C):", 37.0, "%.1f")
        ph_constant = nummer_eingabe("pH-Wert (konstant):", 7.20,
                                    step=0.01)
        oxygen_saturation = nummer_eingabe("Sauerstoffsättigung (%):", 50.0,
                                          "%.1f", 1.0)
        
        # "Simulation starten"-Button mit Fortschrittsanzeige
        if st.button("Simulation starten", type="primary"):
            with st.spinner("Simulation läuft..."):
                # Führt Hauptsimulation mit eingegebenen Parametern durch
                st.session_state.data = simulator.simulate(
                    initial_glucose, initial_vcd, temperature,
                    ph_constant, oxygen_saturation)
                
                # Berechnet KPIs und speichert Ergebnisse für Vergleich
                kpis = calculate_kpis(st.session_state.data)
                
                # Erstellt neuen Ergebnisdatensatz für Vergleichstabelle
                new_result = pd.DataFrame({
                    'Anfangs-Glukose (g/L)': [initial_glucose],
                    'Anfangs-VCD (10^6 Zellen/mL)': [initial_vcd],
                    'Temperatur (°C)': [temperature],
                    'pH': [ph_constant],
                    'Sauerstoff (%)': [oxygen_saturation],
                    'Antikörper-Gesamtausbeute (mg)': [kpis['final_titer']],
                    'Max VCD (10^6 Zellen/mL)': [kpis['max_vcd']],
                    'Durchschn. Viabilität (%)': [kpis['avg_viability']]
                })
                
                # Fügt neues Ergebnis zur Vergleichstabelle hinzu
                st.session_state.results_df = pd.concat(
                    [st.session_state.results_df, new_result],
                    ignore_index=True)
            
            st.success("Simulation abgeschlossen!")
        
        st.divider()
    
    # Hauptinhaltsbereich mit Tab-Navigation nach Simulation
    if st.session_state.data is not None:
        data = st.session_state.data
        kpis = calculate_kpis(data)
        
        # Strukturiert Ergebnisse in übersichtlichen Tabs
        tab1, tab2, tab3 = st.tabs(["Zeitverläufe", "Daten", "Ergebnisse"])
        
        # Tab 1: Grafische Zeitverlaufs-Darstellung
        with tab1:
            col1, col2 = st.columns([3, 1])
            with col1:
                # Erstellt und zeigt wissenschaftliche Plots
                fig = create_plots(data)
                st.pyplot(fig)
            with col2:
                # KPI-Dashboard neben den Plots für schnelle Bewertung
                st.markdown("### 📊 **Key Performance Indicators**")
                st.metric("🧬 **Max. VCD**", f"{kpis['max_vcd']:.2f}")
                st.metric("💚 **Ø Viabilität**",
                         f"{kpis['avg_viability']:.1f}%")
                st.metric("💚 **Min. Viabilität**",
                         f"{kpis['min_viability']:.1f}%")
                st.metric("🎯 **Antikörper-Gesamtausbeute**", 
                         f"{kpis['final_titer']:.2f} mg")
                    
        # Tab 2: Tabellarische Rohdaten-Darstellung
        with tab2:
            # Zeigt reduzierte Daten (alle 6h) für bessere Übersicht
            display_data = data.iloc[::6]
            st.dataframe(display_data, use_container_width=True,
                        hide_index=True)
        
        # Tab 3: Vergleichsanalyse mehrerer Simulationsläufe
        with tab3:
            if not st.session_state.results_df.empty:
                st.subheader("Simulationsergebnisse")
                
                # Fügt Laufnummer für bessere Orientierung hinzu
                display_results = st.session_state.results_df.copy()
                display_results.insert(0, 'Nr.', range(1,
                                       len(display_results) + 1))
                st.dataframe(display_results, use_container_width=True,
                            hide_index=True)
                
                # Zeigt Korrelationsanalyse bei mehreren Simulationsläufen
                if len(st.session_state.results_df) > 1:
                    st.subheader("Parameter-KPI Korrelationen")
                    
                    # Definiert Input-Parameter und Output-KPIs
                    parameters = ['Anfangs-Glukose (g/L)',
                                 'Anfangs-VCD (10^6 Zellen/mL)',
                                 'Temperatur (°C)', 'pH', 'Sauerstoff (%)']
                    kpis_cols = ['Antikörper-Gesamtausbeute (mg)',
                                'Max VCD (10^6 Zellen/mL)',
                                'Durchschn. Viabilität (%)']
                    
                    # Berechnet Korrelationsmatrix zwischen Parametern und KPIs
                    corr_data = st.session_state.results_df[parameters +
                                                           kpis_cols]
                    param_kpi_corr = corr_data.corr()
                    param_to_kpi = param_kpi_corr.loc[parameters, kpis_cols]
                    
                    # Erstellt Heatmap für Korrelationsvisualisierung
                    fig1, ax1 = plt.subplots(figsize=(10, 6))
                    sns.heatmap(param_to_kpi, annot=True, cmap='RdBu_r',
                               center=0, fmt='.2f', ax=ax1,
                               cbar_kws={'label': 'Korrelation'})
                    plt.xticks(rotation=45, ha='right')
                    plt.yticks(rotation=0)
                    plt.tight_layout()
                    st.pyplot(fig1)
            else:
                st.info("Führen Sie mehrere Simulationen durch, um "
                       "Ergebnisse zu vergleichen.")
    else:
        # Anweisungen für Erstnutzer anzeigen
        st.info("Bitte starten Sie die Simulation, um die Ergebnisse "
               "anzuzeigen.")


if __name__ == "__main__":
    main()
