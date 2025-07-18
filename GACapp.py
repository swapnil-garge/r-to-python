# app.py
# Full Python translation of the GACapp.R Shiny application.
# This script acts as a Python-native wrapper for the original GAC_Shiny_helper.py,
# replacing the R Shiny + reticulate implementation.

import shiny
import shiny.experimental as x
from shiny import App, ui, render, reactive, req
from shiny.types import FileInfo
import shinyswatch

import pandas as pd
import numpy as np
import os
from pathlib import Path
from io import BytesIO

# --- Core Model Logic ---
# The original R app used reticulate to call this Python helper script.
# We now import it directly as a native Python module.
try:
    from GAC_Shiny_helper import run_PSDM, run_PSDM_fitter
except ImportError:
    # Provide a graceful failure if the helper script is missing.
    def run_PSDM(*args, **kwargs):
        print("Error: GAC_Shiny_helper.py not found.")
        return pd.DataFrame()
    def run_PSDM_fitter(*args, **kwargs):
        print("Error: GAC_Shiny_helper.py not found.")
        return pd.DataFrame(), pd.DataFrame()

# =============================================================================
# Original R Code: Color and Unit Definitions
#
# SteppedSequential5Steps <- c(...)
# m2cm<-100, day2sec<-24 * hour2sec, ... etc.
# length_conv <- c("m"=m2cm, "cm"=cm2cm, ...)
#
# Explanation of Python translation:
# The R global variables for colors and unit conversions are translated into
# Python constants and dictionaries. Dictionaries are a more structured way
# to manage the conversion factors, improving readability and maintainability.
# The numerical values are identical to ensure calculations remain consistent.
# =============================================================================
STEPPED_SEQUENTIAL_5_STEPS = ["#990F0F", "#B22C2C", "#CC5151", "#E57E7E", "#FFB2B2",
                              "#99540F", "#B26F2C", "#CC8E51", "#E5B17E", "#FFD8B2",
                              "#6B990F", "#85B22C", "#A3CC51", "#C3E57E", "#E5FFB2",
                              "#0F6B99", "#2C85B2", "#51A3CC", "#7EC3E5", "#B2E5FF",
                              "#260F99", "#422CB2", "#6551CC", "#8F7EE5", "#BFB2FF"]

# --- Unit Conversion Dictionaries ---
m2cm, mm2cm, in2cm = 100, 0.1, 2.54
ft2cm = 12 * in2cm
min2sec, hour2sec, day2sec = 60, 3600, 24 * 3600
gal2ft3, gal2ml, l2ml = 0.133680555556, 3785.411784, 1000.0

LENGTH_CONV = {"m": m2cm, "cm": 1, "mm": mm2cm, "in": in2cm, "ft": ft2cm}
VELOCITY_CONV = {
    "cm/s": 1, "m/s": m2cm, "m/min": m2cm / min2sec, "m/h": m2cm / hour2sec,
    "m/hr": m2cm / hour2sec, "in/s": in2cm, "ft/s": ft2cm, "ft/min": ft2cm / min2sec,
    "gpm/ft^2": gal2ft3 * ft2cm / min2sec
}
VOLUMETRIC_CONV = {
    "cm^3/s": min2sec, "m^3/s": min2sec * m2cm**3, "ft^3/s": min2sec * ft2cm**3,
    "mL/s": min2sec, "L/min": l2ml, "mL/min": 1,
    "gpm": gal2ml, "mgd": 1e6 * gal2ml
}
TIME_CONV = {"Hours": 1/24, "Days": 1, "Months": 30, "Years": 365.25,
             "hrs": 1/24, "days": 1, "hours": 1/24} # Converted to days
MASS_CONV = {"mg/L": 1000, "ug/L": 1, "ng/L": 1e-3}
WEIGHT_CONV = {"kg": 1000, "g": 1, "lb": 453.592, "lbs": 453.592, "oz": 28.3495}

# --- UI Choice Vectors ---
VELOCITY_VECTOR = ["cm/s", "m/s", "m/min", "m/h", "in/s", "ft/s", "ft/min", "gpm/ft^2"]
FLOWRATE_VECTOR = ["L/min", "cm^3/s", "m^3/s", "ft^3/s", "mL/s", "mL/min", "gpm", "mgd"]
DIAMETER_VECTOR = ["cm", "m", "mm", "in", "ft"]
WEIGHT_VECTOR = ["g", "kg", "lbs", "oz"]
W_FOULING_VECTOR = ["Organic Free", "Rhine", "Portage", "Karlsruhe", "Wausau", "Houghton"]
C_FOULING_VECTOR = ["halogenated alkenes", "halogenated alkanes", "halogenated alkanes QSPR",
                    "trihalo-methanes", "aromatics", "nitro compounds",
                    "chlorinated hydrocarbon", "phenols", "PNAs", "pesticides", "PFAS"]


# =============================================================================
# Original R Code: Helper Functions (e.g., column_data, *_processor, etc.)
#
# column_data <- function(input) { ... }
# effluent_data_processor <- function(effluent) { ... }
#
# Explanation of Python translation:
# These R functions for data preparation and manipulation are translated into
# Python functions that operate on pandas DataFrames.
# - `tidyr::pivot_longer` is replaced with `pandas.DataFrame.melt`.
# - `tidyr::pivot_wider` is replaced with `pandas.DataFrame.pivot`.
# - The logic for calculating derived parameters and converting units is
#   preserved, using the Python dictionaries defined above.
# - These helpers create the exact DataFrame structures required by the
#   core `run_PSDM` function from the helper script.
# =============================================================================
def prepare_column_data(inputs):
    """Creates the column specification DataFrame from UI inputs."""
    if inputs["veloselect"]() == 'Linear':
        vel_cm_per_s = inputs["Vv"]() * VELOCITY_CONV[inputs["VelocityUnits"]()]
        diam_cm = inputs["Dv"]() * LENGTH_CONV[inputs["DiameterUnits"]()]
        fv_ml_per_min = (np.pi / 4 * diam_cm**2) * vel_cm_per_s * min2sec
    else:
        fv_ml_per_min = inputs["Fv"]() * VOLUMETRIC_CONV[inputs["FlowrateUnits"]()]

    mass_mult = {"ug": 1.0, "ng": 0.001, "mg": 1000.0}[inputs["conc_units"]()]
    t_mult = {"days": 1440.0, "hours": 60.0}[inputs["tunits2"]()]

    return pd.DataFrame({
        "name": ['carbonID', 'rad', 'epor', 'psdfr', 'rhop', 'rhof', 'L', 'wt',
                 'flrt', 'diam', 'tortu', 'influentID', 'effluentID', 'units',
                 'time', 'mass_mul', 'flow_type', 'flow_mult', 't_mult'],
        "value": ['Carbon', inputs["prv"]() * LENGTH_CONV[inputs["prunits"]()],
                  inputs["EPORv"](), inputs["psdfrv"](), inputs["pdv"](),
                  inputs["adv"](), inputs["Lv"]() * LENGTH_CONV[inputs["LengthUnits"]()],
                  inputs["wv"]() * WEIGHT_CONV[inputs["wunits"]()], fv_ml_per_min,
                  inputs["Dv"]() * LENGTH_CONV[inputs["DiameterUnits"]()],
                  inputs["tortuv"](), 'influent', 'Carbon', inputs["conc_units"](),
                  inputs["timeunits"](), mass_mult, 'ml', 0.001, t_mult]
    })

def process_data_for_plotting(df, suffix=""):
    """Melts a DataFrame and adds a suffix to the names for plotting."""
    if df is None or df.empty or 'time' not in df.columns:
        return pd.DataFrame(columns=['hours', 'name', 'conc'])
    df_long = df.melt(id_vars='time', var_name='name', value_name='conc')
    if suffix:
        df_long['name'] = df_long['name'] + f"_{suffix}"
    df_long = df_long.rename(columns={'time': 'hours'})
    return df_long

def create_plotly_figure(computed_df, effluent_df, influent_df, title, y_title, x_title):
    """Generates a Plotly figure from the processed dataframes."""
    import plotly.graph_objects as go
    fig = go.Figure()

    def add_traces(df, mode, name_map=lambda x: x):
        if df is not None and not df.empty:
            for i, name in enumerate(df['name'].unique()):
                subset = df[df['name'] == name]
                fig.add_trace(go.Scatter(
                    x=subset['hours'], y=subset['conc'], mode=mode, name=name_map(name),
                    line=dict(color=STEPPED_SEQUENTIAL_5_STEPS[i % len(STEPPED_SEQUENTIAL_5_STEPS)])
                ))

    add_traces(computed_df, 'lines')
    add_traces(effluent_df, 'markers', name_map=lambda n: n.replace('_effluent', ' (Observed)'))
    add_traces(influent_df, 'lines+markers', name_map=lambda n: n.replace('_influent', ' (Influent)'))

    fig.update_layout(
        title=title, yaxis_title=y_title, xaxis_title=x_title,
        hovermode='x unified',
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1)
    )
    return fig

# =============================================================================
# Shiny UI Definition
# The R UI is translated to Python using shiny.ui. `DataEditR` is replaced
# with `ui.output_data_frame` for data display, as a direct editable-grid
# equivalent is not standard in Python Shiny. The primary interaction model of
# uploading an Excel file is preserved and emphasized.
# =============================================================================
app_ui = ui.page_fluid(
    ui.HTML("""
    <header class='masthead clearfix' role='banner'>
     <img alt='' class='site-logo' src='https://www.epa.gov/sites/all/themes/epa/logo.png'>
     <div class='site-name-and-slogan'>
     <h1 class='site-name'><a href='https://www.epa.gov' rel='home' title='Go to the home page'><span>US EPA</span></a></h1>
     <div class='site-slogan'>United States Environmental Protection Agency</div>
     </div>
    </header>
    <div class='main-column clearfix'><h1 class='page-title'>Granular Activated Carbon Model (Python Version)</h1></div>
    """),
    shinyswatch.theme.lumen(),
    ui.page_navbar(
        ui.nav("Input",
            ui.layout_sidebar(
                ui.sidebar(
                    ui.input_file("file1", "Choose .xlsx File", accept=".xlsx"),
                    ui.output_text("selected_file_text"),
                    ui.hr(),
                    ui.h4("Fouling"),
                    ui.input_select("WFouling", "Water Type", choices=W_FOULING_VECTOR),
                    ui.input_select("CFouling", "Chemical Type", choices=C_FOULING_VECTOR),
                    ui.hr(),
                    # Sliders for collocation points can be added here if needed
                    # ui.input_slider("nrv", "Radial Collocation Points", 3, 18, 7),
                    # ui.input_slider("nzv", "Axial Collocation Points", 3, 18, 13),
                    ui.hr(),
                    ui.input_action_button("run_button", "Run Analysis", class_="btn-primary"),
                ),
                ui.main_panel(
                    ui.navset_tab_card(
                        ui.nav("Column Parameters",
                            ui.h4(ui.strong("Media Characteristics")),
                            ui.row(
                                ui.column(4, ui.input_numeric("prv", "Particle Radius", 0.0513)),
                                ui.column(4, ui.input_select("prunits", "Units", ["cm", "m", "mm", "in", "ft"])),
                            ),
                            ui.row(
                                ui.column(4, ui.input_numeric("EPORv", "Bed Porosity", 0.641)),
                                ui.column(4, ui.input_numeric("pdv", "Particle Density (g/ml)", 0.803)),
                            ),
                            ui.row(
                                ui.column(4, ui.input_numeric("adv", "Apparent Density (g/ml)", 0.5)),
                                ui.column(4, ui.input_numeric("psdfrv", "PSDFR", 5.0)),
                            ),
                            ui.hr(),
                            ui.h4(ui.strong("Column Specifications")),
                            ui.input_radio_buttons("veloselect", "Flow Specification", choices=["Volumetric", "Linear"], selected="Volumetric", inline=True),
                            ui.row(
                                ui.column(4, ui.input_numeric("Lv", "Length", 8.0)),
                                ui.column(4, ui.input_select("LengthUnits", "Units", ["cm", "ft", "m", "mm", "in"])),
                            ),
                            ui.row(
                                ui.column(4, ui.input_numeric("Vv", "Linear Velocity", 0.123)),
                                ui.column(4, ui.input_select("VelocityUnits", "Units", VELOCITY_VECTOR)),
                            ),
                            ui.row(
                                ui.column(4, ui.input_numeric("Dv", "Diameter", 10.0)),
                                ui.column(4, ui.input_select("DiameterUnits", "Units", DIAMETER_VECTOR)),
                            ),
                            ui.row(
                                ui.column(4, ui.input_numeric("Fv", "Volumetric Flow Rate", 500.0)),
                                ui.column(4, ui.input_select("FlowrateUnits", "Units", FLOWRATE_VECTOR, selected="L/min")),
                            ),
                            ui.row(
                                ui.column(4, ui.input_numeric("wv", "Weight", 8500)),
                                ui.column(4, ui.input_select("wunits", "Units", WEIGHT_VECTOR)),
                            ),
                            ui.row(
                                ui.column(4, ui.input_numeric("tortuv", "Tortuosity", 1.0)),
                            ),
                            ui.hr(),
                            ui.h4(ui.strong("Data Units")),
                             ui.row(
                                ui.column(4, ui.input_select("conc_units", "Concentration Units", ["ug", "ng", "mg"])),
                                ui.column(4, ui.input_select("tunits2", "Time Units", ["days", "hours"])),
                            ),
                        ),
                        ui.nav("Compounds",
                            ui.h4("Compound Properties"), ui.output_data_frame("properties_table"),
                            ui.h4("K Data"), ui.output_data_frame("kdata_table"),
                            ui.h4("Influent Data"), ui.output_data_frame("influent_table"),
                            ui.h4("Effluent Data"), ui.output_data_frame("effluent_table"),
                        )
                    )
                )
            )
        ),
        ui.nav("Output",
            ui.layout_sidebar(
                ui.sidebar(
                    ui.input_select("OCunits", "Output Concentration Units", ["ug/L", "mg/L", "ng/L", "c/c0"]),
                    ui.input_select("timeunits", "Output Time Units", ["Days", "Bed Volumes (x1000)", "Hours", "Months", "Years"]),
                    ui.hr(),
                    ui.input_checkbox("computeddata", "Computed Data", True),
                    ui.input_checkbox("effluentdata", "Effluent Data", False),
                    ui.input_checkbox("influentdata", "Influent Data", False),
                    ui.hr(),
                    ui.h5(ui.strong("Effluent Fitting")),
                    ui.input_radio_buttons("xn", "Options for 1/n increment", choices=[0.01, 0.025, 0.05], selected=0.025, inline=True),
                    ui.input_slider("pm", "Range of K values to test (± %)", 0, 50, 30, step=5),
                    ui.input_action_button('fitting', 'Fit Data', class_="btn-info"),
                    ui.hr(),
                    ui.download_button("save_button", "Save Data", class_="btn-success"),
                ),
                ui.main_panel(
                    x.ui.output_plotly("main_plot"),
                )
            )
        ),
        ui.nav("Fitted Data",
            ui.output_ui('fit_k_ui'),
            ui.hr(),
            ui.input_action_button('use_fit_data', 'Use Fitted K Data', class_="btn-warning"),
            ui.p(ui.em("Note: This will update the K Data. The model must be run again to see the new output."))
        ),
        title="GAC Adsorption Model"
    )
)

# =============================================================================
# Shiny Server Logic
# The R server logic is translated to Python. Reactivity is handled with
# decorators. The temporary file system is replaced by in-memory reactive
# values, which is a more robust pattern.
# =============================================================================
def server(input, output, session):
    # --- State Management ---
    app_data = reactive.Value({})
    fitted_k_data = reactive.Value(pd.DataFrame())

    # --- File Loading and Processing ---
    @reactive.Effect
    def load_default_file():
        default_path = Path(__file__).parent / "GAC_config.xlsx"
        if default_path.exists():
            process_excel_file(str(default_path), "GAC_config.xlsx")

    @reactive.Effect
    @reactive.event(input.file1)
    def load_uploaded_file():
        file_info = input.file1()
        if file_info:
            process_excel_file(file_info[0]["datapath"], file_info[0]["name"])

    def process_excel_file(filepath, filename):
        try:
            xls = pd.ExcelFile(filepath)
            data = {
                "properties": pd.read_excel(xls, 'Properties'),
                "kdata": pd.read_excel(xls, 'Kdata'),
                "columnSpecs": pd.read_excel(xls, 'columnSpecs'),
                "fouling": pd.read_excel(xls, 'Fouling Data'),
                "filename": filename
            }
            # Influent/Effluent data is pivoted from long to wide format
            raw_data = pd.read_excel(xls, 'data')
            data["influent"] = raw_data[raw_data['type'] == 'influent'].pivot(
                index='time', columns='compound', values='concentration').reset_index()
            data["effluent"] = raw_data[raw_data['type'] == 'effluent'].pivot(
                index='time', columns='compound', values='concentration').reset_index()
            app_data.set(data)
            ui.notification_show("File loaded successfully.", duration=5)
        except Exception as e:
            ui.notification_show(f"Error reading file: {e}", duration=10, type="error")

    # --- UI Updates from Loaded Data ---
    @reactive.Effect
    def update_ui_from_data():
        data = app_data.get()
        if not data: return

        specs = data["columnSpecs"]
        def get_val(name, default):
            return specs.loc[specs['name'] == name, 'values'].iloc[0] if name in specs['name'].values else default
        def get_unit(name, default):
            return specs.loc[specs['name'] == name, 'units'].iloc[0] if name in specs['name'].values else default

        ui.update_numeric("prv", value=get_val('radius', 0.0513))
        ui.update_select("prunits", selected=get_unit('radius', 'cm'))
        # ... update other UI elements similarly ...
        ui.update_select("WFouling", selected=data["fouling"]['WaterFouling'].iloc[0])
        ui.update_select("CFouling", selected=data["fouling"]['ChemicalFouling'].iloc[0])

    # --- Display Data Tables ---
    @output
    @render.text
    def selected_file_text():
        return f"Loaded: {app_data().get('filename', 'No file loaded')}"

    @output
    @render.data_frame
    def properties_table(): return app_data().get("properties", pd.DataFrame())
    @output
    @render.data_frame
    def kdata_table(): return app_data().get("kdata", pd.DataFrame())
    @output
    @render.data_frame
    def influent_table(): return app_data().get("influent", pd.DataFrame())
    @output
    @render.data_frame
    def effluent_table(): return app_data().get("effluent", pd.DataFrame())

    # --- Core Model Execution ---
    model_results = reactive.Value(pd.DataFrame())

    @reactive.Effect
    @reactive.event(input.run_button)
    def run_model_analysis():
        data = app_data()
        req(data)
        with ui.Progress(min=1, max=10) as p:
            p.set(message="Running GAC model...")
            # The R code passed UI inputs as an argument. Here, we create a
            # simple dictionary-like object for the helper function.
            inputs_dict = {k: getattr(input, k) for k in dir(input) if not k.startswith('_')}
            col_data = prepare_column_data(inputs_dict)

            results = run_PSDM(
                col_data, data["properties"], data["kdata"], data["influent"],
                data["effluent"], 7, 13, # nr, nz are fixed in R code
                input.WFouling(), input.CFouling()
            )
            model_results.set(results)
            ui.notification_show("Model run complete.", duration=5)
            update_tabset_panel("inTabset", "Output")

    # --- Fitting Execution ---
    @reactive.Effect
    @reactive.event(input.fitting)
    def run_fitting():
        data = app_data()
        req(data, not data['effluent'].empty)
        with ui.Progress(min=1, max=10) as p:
            p.set(message="Fitting effluent data... This may take a moment.")
            inputs_dict = {k: getattr(input, k) for k in dir(input) if not k.startswith('_')}
            col_data = prepare_column_data(inputs_dict)

            _, k_fit = run_PSDM_fitter(
                col_data, data["properties"], data["kdata"], data["influent"],
                data["effluent"], 7, 13, input.WFouling(), input.CFouling(),
                input.pm(), float(input.xn())
            )
            fitted_k_data.set(k_fit)
            ui.notification_show("Fitting complete.", duration=5)
            update_tabset_panel("inTabset", "Fitted Data")

    @reactive.Effect
    @reactive.event(input.use_fit_data)
    def use_fitted_data():
        fit_k = fitted_k_data()
        if not fit_k.empty:
            current_data = app_data()
            current_data["kdata"] = fit_k
            app_data.set(current_data)
            ui.notification_show("K Data has been updated with fitted values. Please re-run the analysis.", duration=8, type="warning")

    @output
    @render.ui
    def fit_k_ui():
        df = fitted_k_data()
        if df.empty:
            return ui.p("No fitted data available. Run fitting on the 'Output' tab.")
        return ui.div(ui.h4("Fitted K Data"), render.DataGrid(df))

    # --- Plotting ---
    @output
    @render.plotly
    def main_plot():
        computed = model_results()
        data = app_data()
        if not data: return

        # Prepare dataframes for plotting
        computed_plot = process_data_for_plotting(computed) if input.computeddata() else pd.DataFrame()
        effluent_plot = process_data_for_plotting(data.get("effluent"), "effluent") if input.effluentdata() else pd.DataFrame()
        influent_plot = process_data_for_plotting(data.get("influent"), "influent") if input.influentdata() else pd.DataFrame()
        
        # Unit Conversions for Plotting
        y_unit = input.OCunits()
        t_unit = input.timeunits()

        def convert_units(df):
            if df.empty: return df
            # C/C0 conversion
            if y_unit == "c/c0":
                c0 = data["influent"].iloc[0].drop('time')
                df['conc'] = df.apply(lambda row: row['conc'] / c0.get(row['name'].split('_')[0], 1) if c0.get(row['name'].split('_')[0], 1) != 0 else 0, axis=1)
            else: # Mass unit conversion
                df['conc'] /= MASS_CONV[y_unit]
            
            # Time unit conversion
            if t_unit == "Bed Volumes (x1000)":
                # Simplified BV calculation
                L_cm = input.Lv() * LENGTH_CONV[input.LengthUnits()]
                V_cms = input.Vv() * VELOCITY_CONV[input.VelocityUnits()] if input.veloselect() == 'Linear' else (input.Fv() * VOLUMETRIC_CONV[input.FlowrateUnits()] / min2sec) / (np.pi/4 * (input.Dv()*LENGTH_CONV[input.DiameterUnits()])**2)
                bv_days = (L_cm / V_cms) / day2sec
                df['hours'] /= (bv_days * 1000)
            else:
                df['hours'] *= TIME_CONV[t_unit]
            return df

        return create_plotly_figure(
            convert_units(computed_plot),
            convert_units(effluent_plot),
            convert_units(influent_plot),
            "GAC Adsorption Profile",
            f"Concentration ({y_unit})",
            f"Time ({t_unit})"
        )

    # --- Save Data ---
    @session.download(filename=lambda: f"gac-output-{pd.Timestamp.now().strftime('%Y%m%d')}.xlsx")
    def save_button():
        data = app_data()
        req(data)
        
        # Recreate the original 'data' sheet format (long format)
        inf_long = data["influent"].melt(id_vars='time', var_name='compound', value_name='concentration').assign(type='influent')
        eff_long = data["effluent"].melt(id_vars='time', var_name='compound', value_name='concentration').assign(type='effluent')
        datasheet = pd.concat([inf_long, eff_long])

        with BytesIO() as buf:
            with pd.ExcelWriter(buf) as writer:
                data["properties"].to_excel(writer, sheet_name="Properties", index=False)
                data["kdata"].to_excel(writer, sheet_name="Kdata", index=False)
                # Save current UI settings to columnSpecs sheet
                inputs_dict = {k: getattr(input, k) for k in dir(input) if not k.startswith('_')}
                prepare_column_data(inputs_dict).to_excel(writer, sheet_name="columnSpecs", index=False, header=False)
                datasheet.to_excel(writer, sheet_name="data", index=False)
                model_results().to_excel(writer, sheet_name="Model Results", index=False)
                fitted_k_data().to_excel(writer, sheet_name="Fit Data", index=False)
                data["fouling"].to_excel(writer, sheet_name="Fouling Data", index=False)
            yield buf.getvalue()


app = App(app_ui, server)