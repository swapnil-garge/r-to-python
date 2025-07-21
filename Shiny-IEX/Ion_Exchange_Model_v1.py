# app.py
# Full Python translation of the Ion_Exchange_Model.R Shiny application.
# Includes Alkalinity and kL Guesser tabs.
# Author: Gemini
# Date: July 21, 2025

import shiny
import shiny.experimental as x
from shiny import App, ui, render, reactive, req
from shiny.types import FileInfo
import shinyswatch

import pandas as pd
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
from scipy.special import roots_jacobi
import os
from pathlib import Path

# =============================================================================
# Constants and Unit Conversions
# =============================================================================
STEPPED_SEQUENTIAL_5_STEPS = ["#990F0F", "#B22C2C", "#CC5151", "#E57E7E", "#FFB2B2",
                              "#99540F", "#B26F2C", "#CC8E51", "#E5B17E", "#FFD8B2",
                              "#6B990F", "#85B22C", "#A3CC51", "#C3E57E", "#E5FFB2",
                              "#0F6B99", "#2C85B2", "#51A3CC", "#7EC3E5", "#B2E5FF",
                              "#260F99", "#422CB2", "#6551CC", "#8F7EE5", "#BFB2FF"]

UNITS = {
    'm2cm': 100, 'mm2cm': 0.1, 'cm2cm': 1, 'in2cm': 2.54, 'ft2cm': 12 * 2.54,
    'min2sec': 60, 'hour2sec': 3600, 'day2sec': 24 * 3600,
    'gal2ml': 3785.411784, 'l2ml': 1000.0
}
LENGTH_CONV = {"m": UNITS['m2cm'], "cm": UNITS['cm2cm'], "mm": UNITS['mm2cm'], "in": UNITS['in2cm'], "ft": UNITS['ft2cm']}
VELOCITY_CONV = {
    "cm/s": UNITS['cm2cm'], "m/s": UNITS['m2cm'], "m/min": UNITS['m2cm'] / UNITS['min2sec'],
    "m/h": UNITS['m2cm'] / UNITS['hour2sec'], "m/hr": UNITS['m2cm'] / UNITS['hour2sec'],
    "in/s": UNITS['in2cm'], "ft/s": UNITS['ft2cm'], "ft/min": UNITS['ft2cm'] / UNITS['min2sec'],
    "gpm/ft^2": 0.133680555556 * UNITS['ft2cm'] / UNITS['min2sec']
}
VOLUMETRIC_CONV = {
    "cm^3/s": UNITS['cm2cm'], "m^3/s": UNITS['m2cm']**3, "ft^3/s": UNITS['ft2cm']**3,
    "mL/s": UNITS['cm2cm'], "L/min": UNITS['l2ml'] / UNITS['min2sec'], "mL/min": 1 / UNITS['min2sec'],
    "gpm": UNITS['gal2ml'] / UNITS['min2sec'], "mgd": 1e6 * UNITS['gal2ml'] / UNITS['day2sec']
}
TIME_CONV = {"Hours": UNITS['hour2sec'], "Days": UNITS['day2sec'], "Months": 30 * UNITS['day2sec'], "Years": 365.25 * UNITS['day2sec'],
             "hr": UNITS['hour2sec'], "day": UNITS['day2sec'], "month": 30 * UNITS['day2sec'], "year": 365.25 * UNITS['day2sec']}
DS_CONV = {"ft^2/s": UNITS['ft2cm']**2, "m^2/s": UNITS['m2cm']**2, "cm^2/s": UNITS['cm2cm'], "in^2/s": UNITS['in2cm']**2}
MASS_CONV = {"mg/L": 1, "ug/L": 1e-3, "ng/L": 1e-6, "c/c0": 1} # Divisor for output units
MASS_CONV_INPUT = {"meq": 1, "meq/L": 1, "mg": 1, "ug": 1e-3, "ng": 1e-6, "mg/L": 1, "ug/L": 1e-3, "ng/L": 1e-6}

NT_REPORT = 201

# =============================================================================
# Core Scientific Models (rad_colloc, ax_colloc, HSDMIX_solve, etc.)
# These functions are direct translations of the R versions.
# =============================================================================
def rad_colloc(N):
    N_int = N - 1
    roots_shifted, _ = roots_jacobi(N_int, 2.5, 1.5)
    roots_x_squared = (roots_shifted + 1) / 2
    roots_non_sym = np.sort(np.append(roots_x_squared, 1))

    derivatives = pd.DataFrame({'roots': roots_non_sym, 'p_1': np.zeros(N), 'p_2': np.zeros(N), 'p_3': np.zeros(N)})

    for i in range(N):
        x_i = derivatives['roots'][i]
        j_values = derivatives['roots'][derivatives['roots'] != x_i]
        delta = x_i - j_values
        p_1, p_2, p_3 = np.zeros(N), np.zeros(N), np.zeros(N)
        p_1[0] = 1.0
        for j in range(N_int):
            p_1[j+1] = delta[j] * p_1[j]
            p_2[j+1] = delta[j] * p_2[j] + 2 * p_1[j]
            p_3[j+1] = delta[j] * p_3[j] + 3 * p_2[j]
        derivatives.loc[i, ['p_1', 'p_2', 'p_3']] = [p_1[N-1], p_2[N-1], p_3[N-1]]

    Ar, Br = np.zeros((N, N)), np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            if i == j:
                Ar[i, j] = 0.5 * derivatives['p_2'][i] / derivatives['p_1'][i]
                Br[i, j] = (1/3) * derivatives['p_3'][i] / derivatives['p_1'][i]
            else:
                diff = derivatives['roots'][i] - derivatives['roots'][j]
                Ar[i, j] = (1 / diff) * (derivatives['p_1'][i] / derivatives['p_1'][j])
                Br[i, j] = 2 * Ar[i, j] * (Ar[i, i] - 1 / diff)

    roots_i_sqrt = np.sqrt(derivatives['roots'].values)
    Br_sym = 4 * np.outer(derivatives['roots'].values, np.ones(N)) * Br + 6 * Ar
    w_i_prime = 1 / (derivatives['roots'] * derivatives['p_1']**2)
    W_i_manu = (1 / 3.0) * w_i_prime / np.sum(w_i_prime)
    return Br_sym, W_i_manu

def ax_colloc(NZ):
    NZ_int = NZ - 2
    roots_shifted, _ = roots_jacobi(NZ_int, 1.0, 1.0)
    roots_Z = np.sort(np.concatenate(([0], (roots_shifted + 1) / 2, [1])))
    derivatives = pd.DataFrame({'roots': roots_Z, 'p_1': np.zeros(NZ), 'p_2': np.zeros(NZ)})
    for i in range(NZ):
        x_i = derivatives['roots'][i]
        j_values = derivatives['roots'][derivatives['roots'] != x_i]
        delta = x_i - j_values
        p_1, p_2 = np.zeros(NZ), np.zeros(NZ)
        p_1[0] = 1.0
        for j in range(NZ - 1):
            p_1[j+1] = delta[j] * p_1[j]
            p_2[j+1] = delta[j] * p_2[j] + 2 * p_1[j]
        derivatives.loc[i, ['p_1', 'p_2']] = [p_1[NZ-1], p_2[NZ-1]]
    AZ = np.zeros((NZ, NZ))
    for i in range(NZ):
        for j in range(NZ):
            if i == j:
                AZ[i, j] = 0.5 * derivatives['p_2'][i] / derivatives['p_1'][i]
            else:
                AZ[i, j] = (1 / (derivatives['roots'][i] - derivatives['roots'][j])) * (derivatives['p_1'][i] / derivatives['p_1'][j])
    return AZ

def HSDMIX_solve(params, ions, Cin, inputtime, nt_report):
    NR = int(params.loc[params['name'] == 'nr', 'value'].iloc[0])
    NZ = int(params.loc[params['name'] == 'nz', 'value'].iloc[0])
    Q, L, v, EBED, rb = [params.loc[params['name'] == n, 'value'].iloc[0] for n in ['Q', 'L', 'v', 'EBED', 'rb']]
    ion_names, KxA, valence, kL, Ds = [ions[c].to_numpy() for c in ['name', 'KxA', 'valence', 'kL', 'Ds']]
    C_in_t = Cin.copy().to_numpy()
    NION = len(ion_names)
    C_in_t[:, 0] *= inputtime
    times = np.linspace(0.0, C_in_t[-1, 0] * 0.99, nt_report)
    C_in_0 = C_in_t[0, 1:(NION + 1)]
    CT = np.sum(C_in_0)
    interp_list = [interp1d(C_in_t[:, 0], C_in_t[:, i + 1], bounds_error=False, fill_value="extrapolate") for i in range(NION)]
    x0 = np.zeros(((NR + 1), NION, NZ))
    x0[-1, :, 0], x0[-1, 0, 1:], x0[0:NR, 0, :] = C_in_0, CT, Q
    BR, WR = rad_colloc(NR)
    AZ = ax_colloc(NZ)
    def diffun(t, y):
        x = y.reshape((NR + 1, NION, NZ))
        C, q = x[-1, :, :], x[0:NR, :, :]
        dx_dt = np.zeros_like(x)
        C[:, 0] = np.array([interp(t) for interp in interp_list])
        AZ_C = (AZ @ C.T).T
        qs = q[NR - 1, :, :]
        CT_test = np.sum(C, axis=0)
        C_star, J = np.zeros((NION, NZ)), np.zeros((NION, NZ))
        z_slice = slice(1, NZ)
        if np.any(valence == 2):
            dv_ions_mask, mv_ions_mask = (valence == 2), (valence == 1)
            mv_ions_mask[0] = False
            qs_1 = np.maximum(qs[0, z_slice], 1e-9)
            cc = -CT_test[z_slice]
            bb = 1 + np.sum(qs[mv_ions_mask][:, z_slice] / KxA[mv_ions_mask, np.newaxis], axis=0) / qs_1
            aa = np.sum(qs[dv_ions_mask][:, z_slice] / KxA[dv_ions_mask, np.newaxis], axis=0) / qs_1**2
            denom = -bb - np.sqrt(np.maximum(bb**2 - 4 * aa * cc, 0))
            C_star[0, z_slice] = 2 * (cc / np.maximum(denom, 1e-9))
            for i in range(1, NION): C_star[i, z_slice] = (qs[i, z_slice] / KxA[i]) * (C_star[0, z_slice] / qs_1)**valence[i]
        else:
            sum_terms = np.sum(qs[:, z_slice] / KxA[:, np.newaxis], axis=0) / CT_test[z_slice]
            for i in range(1, NION): C_star[i, z_slice] = qs[i, z_slice] / KxA[i] / sum_terms
        J[1:, z_slice] = -kL[1:, np.newaxis] * (C[1:, z_slice] - C_star[1:, z_slice])
        J[0, z_slice] = -np.sum(J[1:, z_slice], axis=0)
        Jas = (3 / rb) * J
        dx_dt[-1, :, z_slice] = (-v / L * AZ_C[:, z_slice] + (1 - EBED) * Jas[:, z_slice]) / EBED
        BR_q = np.array([BR @ q[:, i, z_slice] for i in range(NION)]).transpose(1, 0, 2)
        dq_dt = np.zeros_like(q)
        dq_dt[:, 1:, :] = (Ds[1:, np.newaxis, np.newaxis] / rb**2) * BR_q[:, 1:, :]
        dq_dt[0:NR-1, 0, z_slice] = -np.sum(dq_dt[0:NR-1, 1:, z_slice], axis=1)
        surf_term = WR[0:NR-1] @ dq_dt[0:NR-1, :, z_slice]
        dx_dt[0:NR-1, :, z_slice] = dq_dt[0:NR-1, :, z_slice]
        dx_dt[NR-1, :, z_slice] = (-1 / rb * J[:, z_slice] - surf_term) / WR[NR-1]
        dx_dt[:, :, 0] = 0.0
        return dx_dt.flatten()
    sol = solve_ivp(diffun, [times[0], times[-1]], x0.flatten(), t_eval=times, method='LSODA')
    return (sol.t / 3600, sol.y.T.reshape(nt_report, NR + 1, NION, NZ)) if sol.success else (times / 3600, np.zeros((nt_report, NR + 1, NION, NZ)))

def PSDMIX_solve(params, ions, Cin, inputtime, nt_report):
    NR = int(params.loc[params['name'] == 'nr', 'value'].iloc[0])
    NZ = int(params.loc[params['name'] == 'nz', 'value'].iloc[0])
    Q, L, v, EBED, EPOR, rb = [params.loc[params['name'] == n, 'value'].iloc[0] for n in ['Q', 'L', 'v', 'EBED', 'EPOR', 'rb']]
    ion_names, KxA, valence, kL, Ds, Dp = [ions[c].to_numpy() for c in ['name', 'KxA', 'valence', 'kL', 'Ds', 'Dp']]
    C_in_t = Cin.copy().to_numpy()
    NION = len(ion_names)
    C_in_t[:, 0] *= inputtime
    times = np.linspace(0.0, C_in_t[-1, 0] * 0.99, nt_report)
    C_in_0 = C_in_t[0, 1:(NION + 1)]
    CT = np.sum(C_in_0)
    interp_list = [interp1d(C_in_t[:, 0], C_in_t[:, i + 1], bounds_error=False, fill_value="extrapolate") for i in range(NION)]
    x0 = np.zeros(((NR + 1), NION, NZ))
    x0[-1, :, 0], x0[-1, 0, 1:], x0[0:NR, 0, :] = C_in_0, CT, Q
    BR, WR = rad_colloc(NR)
    AZ = ax_colloc(NZ)
    def diffun(t, y):
        x = y.reshape((NR + 1, NION, NZ))
        C, Y = x[-1, :, :], x[0:NR, :, :]
        q = Y / (1 - EPOR)
        dx_dt = np.zeros_like(x)
        C[:, 0] = np.array([interp(t) for interp in interp_list])
        AZ_C = (AZ @ C.T).T
        CT_test = np.sum(C, axis=0)
        Cpore = np.zeros_like(q)
        z_slice = slice(1, NZ)
        if np.any(valence == 2):
            dv_ions_mask, mv_ions_mask = (valence == 2), (valence == 1)
            mv_ions_mask[0] = False
            for jj in range(NR):
                q_jj = q[jj, :, z_slice]
                q_jj_1 = np.maximum(q_jj[0, :], 1e-9)
                cc = -CT_test[z_slice]
                bb = 1 + np.sum(q_jj[mv_ions_mask, :] / KxA[mv_ions_mask, np.newaxis], axis=0) / q_jj_1
                aa = np.sum(q_jj[dv_ions_mask, :] / KxA[dv_ions_mask, np.newaxis], axis=0) / q_jj_1**2
                denom = -bb - np.sqrt(np.maximum(bb**2 - 4 * aa * cc, 0))
                Cpore[jj, 0, z_slice] = 2 * (cc / np.maximum(denom, 1e-9))
                for i in range(1, NION): Cpore[jj, i, z_slice] = (q_jj[i,:] / KxA[i]) * (Cpore[jj, 0, z_slice] / q_jj_1)**valence[i]
        else:
            for jj in range(NR):
                q_jj = q[jj, :, z_slice]
                sum_terms = np.sum(q_jj / KxA[:, np.newaxis], axis=0) / CT_test[z_slice]
                for i in range(1, NION): Cpore[jj, i, z_slice] = q_jj[i, :] / KxA[i] / sum_terms
        C_star = Cpore[NR - 1, :, :]
        J = np.zeros((NION, NZ))
        J[1:, z_slice] = -kL[1:, np.newaxis] * (C[1:, z_slice] - C_star[1:, z_slice])
        J[0, z_slice] = -np.sum(J[1:, z_slice], axis=0)
        Jas = (3 / rb) * J
        dx_dt[-1, :, z_slice] = (-v / L * AZ_C[:, z_slice] + (1 - EBED) * Jas[:, z_slice]) / EBED
        BR_Y = np.array([BR @ Y[:, i, z_slice] for i in range(NION)]).transpose(1, 0, 2)
        BR_Cpore = np.array([BR @ Cpore[:, i, z_slice] for i in range(NION)]).transpose(1, 0, 2)
        dY_dt = np.zeros_like(Y)
        dY_dt_calc = (EPOR * (Dp[1:] - Ds[1:])[:, np.newaxis, np.newaxis] * BR_Cpore[:, 1:, :] + Ds[1:, np.newaxis, np.newaxis] * BR_Y[:, 1:, :]) / rb**2
        dY_dt[:, 1:, :] = dY_dt_calc
        dY_dt[0:NR-1, 0, z_slice] = -np.sum(dY_dt[0:NR-1, 1:, z_slice], axis=1)
        surf_term = WR[0:NR-1] @ dY_dt[0:NR-1, :, z_slice]
        dx_dt[0:NR-1, :, z_slice] = dY_dt[0:NR-1, :, z_slice]
        dx_dt[NR-1, :, z_slice] = (-1 / rb * J[:, z_slice] - surf_term) / WR[NR-1]
        dx_dt[:, :, 0] = 0.0
        return dx_dt.flatten()
    sol = solve_ivp(diffun, [times[0], times[-1]], x0.flatten(), t_eval=times, method='LSODA')
    return (sol.t / 3600, sol.y.T.reshape(nt_report, NR + 1, NION, NZ)) if sol.success else (times / 3600, np.zeros((nt_report, NR + 1, NION, NZ)))

# =============================================================================
# Helper Functions for Data Prep
# =============================================================================
def cin_correct(ions_df, cin_df):
    corr_cin = cin_df.copy()
    for _, row in ions_df.iterrows():
        mass_units = row.get("conc_units", "meq")
        if mass_units not in ['meq', 'meq/L']:
            mw, valence, compound = row["mw"], row["valence"], row["name"]
            if compound in corr_cin.columns:
                corr_cin[compound] *= MASS_CONV_INPUT.get(mass_units, 1.0) / mw * valence
    return corr_cin

def model_prep(inputs, iondata, concdata, nt_report):
    if inputs['veloselect']() == 'Linear':
        Vv = inputs['Vv']() * VELOCITY_CONV[inputs['VelocityUnits']()]
    else:
        Dv_cm = inputs['Dv']() * LENGTH_CONV[inputs['DiameterUnits']()]
        area = np.pi / 4 * (Dv_cm ** 2)
        Fv_cm3ps = inputs['Fv']() * VOLUMETRIC_CONV[inputs['FlowrateUnits']()]
        Vv = Fv_cm3ps / area
    param_dict = {
        "Q": inputs['Qv'](), "EBED": inputs['EBEDv'](), "L": inputs['Lv']() * LENGTH_CONV[inputs['LengthUnits']()],
        "v": Vv, "rb": inputs['rbv']() * LENGTH_CONV[inputs['rbunits']()],
        "nr": inputs['nrv'](), "nz": inputs['nzv']()
    }
    if inputs['model']() == "Macroporous (PSDM)": param_dict["EPOR"] = inputs['EPORv']()
    paramdataframe = pd.DataFrame(param_dict.items(), columns=['name', 'value'])
    ion_names, cin_names = set(iondata['name']), set(c for c in concdata.columns if c != 'time')
    if not ion_names.issubset(cin_names) or not cin_names.issubset(ion_names): return None
    corr_ions = iondata.copy()
    corr_cin = cin_correct(iondata, concdata)
    for i, row in iondata.iterrows():
        corr_ions.loc[i, 'kL'] = row['kL'] * VELOCITY_CONV[row['kL_units']]
        corr_ions.loc[i, 'Ds'] = row['Ds'] * DS_CONV[row['Ds_units']]
        if inputs['model']() == "Macroporous (PSDM)" and 'Dp' in row:
            corr_ions.loc[i, 'Dp'] = row['Dp'] * DS_CONV[row['Dp_units']]
    timeconverter = TIME_CONV[inputs['timeunits2']()]
    solver = PSDMIX_solve if inputs['model']() == "Macroporous (PSDM)" else HSDMIX_solve
    return solver(paramdataframe, corr_ions, corr_cin, timeconverter, nt_report)

def create_plotly(computed_df, effluent_df, influent_df, title, y_axis_title, x_axis_title):
    import plotly.graph_objects as go
    fig = go.Figure()
    if computed_df is not None and not computed_df.empty:
        for i, name in enumerate(computed_df['name'].unique()):
            df_sub = computed_df[computed_df['name'] == name]
            fig.add_trace(go.Scatter(x=df_sub['hours'], y=df_sub['conc'], mode='lines', name=name, line=dict(color=STEPPED_SEQUENTIAL_5_STEPS[i % len(STEPPED_SEQUENTIAL_5_STEPS)])))
    if effluent_df is not None and not effluent_df.empty:
        for i, name in enumerate(effluent_df['name'].unique()):
             df_sub = effluent_df[effluent_df['name'] == name]
             color_name = name.replace('_effluent', '')
             try: color_idx = computed_df['name'].unique().tolist().index(color_name)
             except (ValueError, AttributeError): color_idx = -1
             fig.add_trace(go.Scatter(x=df_sub['hours'], y=df_sub['conc'], mode='markers', name=name, marker=dict(color=STEPPED_SEQUENTIAL_5_STEPS[color_idx % len(STEPPED_SEQUENTIAL_5_STEPS)] if color_idx != -1 else 'black')))
    if influent_df is not None and not influent_df.empty:
        for i, name in enumerate(influent_df['name'].unique()):
            df_sub = influent_df[influent_df['name'] == name]
            color_name = name.replace('_influent', '')
            try: color_idx = computed_df['name'].unique().tolist().index(color_name)
            except (ValueError, AttributeError): color_idx = -1
            fig.add_trace(go.Scatter(x=df_sub['hours'], y=df_sub['conc'], mode='lines+markers', name=name, line=dict(color=STEPPED_SEQUENTIAL_5_STEPS[color_idx % len(STEPPED_SEQUENTIAL_5_STEPS)] if color_idx != -1 else 'grey', dash='dot')))
    fig.update_layout(title=title, xaxis_title=x_axis_title, yaxis_title=y_axis_title, hovermode='x unified', legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1))
    return fig

# =============================================================================
# Shiny UI Definition
# =============================================================================
app_ui = ui.page_fluid(
    ui.HTML("""
    <header class='masthead clearfix' role='banner'>
     <img alt='' class='site-logo' src='https://www.epa.gov/sites/all/themes/epa/logo.png'>
     <div class='site-name-and-slogan'><h1 class='site-name'><a href='https://www.epa.gov' rel='home' title='Go to the home page'><span>US EPA</span></a></h1>
     <div class='site-slogan'>United States Environmental Protection Agency</div></div></header>
    <div class='main-column clearfix'><h1 class='page-title'>Ion Exchange Model (Python Version)</h1></div>
    """),
    shinyswatch.theme.cerulean(),
    ui.page_navbar(
        ui.nav("Input",
            ui.layout_sidebar(
                ui.sidebar(
                    ui.input_select("model", "Model Selection", ["Gel-Type (HSDM)", "Macroporous (PSDM)"]),
                    ui.input_file("file1", "Choose .xlsx File", accept=".xlsx"),
                    ui.output_text("selected_file_text"),
                    ui.hr(),
                    ui.input_slider("nrv", "Radial Collocation Points", 3, 18, 7),
                    ui.input_slider("nzv", "Axial Collocation Points", 3, 18, 13),
                    ui.hr(),
                    ui.input_action_button("run_button", "Run Analysis", class_="btn-primary"),
                ),
                ui.navset_tab_card(
                    ui.nav("Column Parameters",
                        ui.h4(ui.strong("Resin Characteristics")),
                        ui.row(ui.column(4, ui.input_numeric("Qv", "Resin Capacity (meq/L)", 1400)),
                               ui.column(4, ui.input_numeric("rbv", "Bead Radius", 0.03375)),
                               ui.column(4, ui.input_select("rbunits", "Units", ["cm", "m", "mm", "in", "ft"]))),
                        ui.row(ui.column(4, ui.input_numeric("EBEDv", "Bed Porosity", 0.35)),
                               ui.column(4, ui.input_numeric("EPORv", "Bead Porosity (PSDM)", 0.2))),
                        ui.hr(),
                        ui.h4(ui.strong("Column Specifications")),
                        ui.input_radio_buttons("veloselect", "Flow Specification", ["Linear", "Volumetric"], selected="Linear", inline=True),
                        ui.row(ui.column(4, ui.input_numeric("Lv", "Length", 14.765)),
                               ui.column(4, ui.input_select("LengthUnits", "Units", ["cm", "m", "mm", "in", "ft"]))),
                        ui.row(ui.column(4, ui.input_numeric("Vv", "Velocity", 0.123)),
                               ui.column(4, ui.input_select("VelocityUnits", "Units", list(VELOCITY_CONV.keys())))),
                        ui.row(ui.column(4, ui.input_numeric("Dv", "Diameter", 4.0)),
                               ui.column(4, ui.input_select("DiameterUnits", "Units", ["cm", "m", "in", "ft"]))),
                        ui.row(ui.column(4, ui.input_numeric("Fv", "Flow Rate", 1.546)),
                               ui.column(4, ui.input_select("FlowrateUnits", "Units", list(VOLUMETRIC_CONV.keys())))),
                        ui.hr(),
                        ui.h4(ui.strong("Concentration Time")),
                        ui.row(ui.column(4, ui.input_select("timeunits2", "Units", ["hr", "day"]))),
                    ),
                    ui.nav("Ions & Concentrations",
                        ui.h4("Ion List"), ui.output_data_frame("ion_table"),
                        ui.h4("Influent Concentration Points"), ui.output_data_frame("cin_table"),
                        ui.h4("Effluent Concentration Points"), ui.output_data_frame("effluent_table"),
                    ),
                    ui.nav("Alkalinity",
                        ui.h4("Bicarbonate Concentration of Alkalinity"),
                        ui.p("This calculator can be used to find bicarbonate concentrations from pH measurements."),
                        ui.hr(),
                        ui.row(
                            ui.column(4, ui.input_numeric("alkvalue", "Alkalinity Value", 100)),
                            ui.column(4, ui.input_select("alkunits", "Concentration Units", ["mg/L CaCO3"])),
                        ),
                        ui.row(
                            ui.column(8, ui.input_slider("pH", "pH", 6, 11, 7, 0.1)),
                        ),
                        ui.hr(),
                        ui.h5("Bicarbonate Concentration (meq/L)"), ui.output_text("bicarb_meq_L"),
                        ui.h5("Bicarbonate Concentration (mg C/L)"), ui.output_text("bicarb_mg_C_L"),
                        ui.h5("Bicarbonate Concentration (mg HCO3-/L)"), ui.output_text("bicarb_mg_HCO3_L"),
                    ),
                    ui.nav("kL Guesser",
                        ui.h4("Film Transfer Coefficient (kL) Estimator"),
                        ui.p("Estimate kL values for common PFAS compounds using the Gnielinski equation."),
                        ui.hr(),
                        ui.row(
                            ui.column(4, ui.input_numeric("temp", "Temperature", 23)),
                            ui.column(4, ui.input_select("tempunits", "Units", ["deg C"])),
                        ),
                        ui.input_action_button('estimate_kl', 'Estimate Values', class_="btn-info"),
                        ui.hr(),
                        ui.row(
                            ui.column(6, ui.h5("PFAS Properties"), ui.output_data_frame("pfas_properties_table")),
                            ui.column(6, ui.h5("kL Estimates"), ui.output_data_frame("kl_estimates_table")),
                        ),
                    )
                )
            )
        ),
        ui.nav("Output",
            ui.layout_sidebar(
                ui.sidebar(
                    ui.input_select("OCunits", "Output Concentration Units", ["mg/L", "ug/L", "ng/L", "c/c0"]),
                    ui.input_select("timeunits", "Output Time Units", ["Hours", "Days", "Bed Volumes (x1000)", "Months", "Years"]),
                    ui.hr(),
                    ui.input_checkbox("computeddata", "Show Computed Data", value=True),
                    ui.input_checkbox("effluentdata", "Show Effluent Data", value=False),
                    ui.input_checkbox("influentdata", "Show Influent Data", value=False),
                    ui.hr(),
                    ui.download_button("save_button", "Save Data", class_="btn-success"),
                ),
                x.ui.output_plotly("plot_counterions"),
                ui.hr(),
                x.ui.output_plotly("plot_other_ions"),
            )
        ),
        ui.nav("About",
            ui.h5("About the Ion Exchange Model"),
            ui.p("A tool to model strong-base anion exchange in drinking water treatment. This is a Python translation of the original R model."),
        ),
        title="Ion Exchange Model",
    )
)

# =============================================================================
# Shiny Server Logic
# =============================================================================
def server(input, output, session):
    app_data = reactive.Value({})
    pfas_properties = reactive.Value(pd.DataFrame())

    @reactive.Effect
    def _():
        # Load default config file
        default_file = Path(__file__).parent / "IEX_config.xlsx"
        if default_file.exists(): process_file(str(default_file))
        # Load PFAS properties for kL Guesser
        pfas_file = Path(__file__).parent / "PSDM" / "PFAS_properties.xlsx"
        if pfas_file.exists():
            df = pd.read_excel(pfas_file, index_col=0).T
            df.columns = df.iloc[0]
            df = df.iloc[1:]
            df = df[['MolarVol']].rename(columns={'MolarVol': 'MolarVol (cm^3/mol)'})
            pfas_properties.set(df)

    @reactive.Effect
    @reactive.event(input.file1)
    def _():
        f = input.file1()
        if f: process_file(f[0]["datapath"])

    def process_file(filepath):
        try:
            xls = pd.ExcelFile(filepath)
            data = {"filename": os.path.basename(filepath)}
            for sheet in ["params", "ions", "cin", "effluent"]:
                data[sheet] = pd.read_excel(xls, sheet) if sheet in xls.sheet_names else pd.DataFrame()
            if 'cin' in data and not data['cin'].empty:
                data['cin'].columns = [c.lower() if c.lower() == 'time' else c for c in data['cin'].columns]
            app_data.set(data)
            ui.notification_show("File processed.", duration=5)
        except Exception as e:
            ui.notification_show(f"Error reading file: {e}", duration=10, type="error")

    @output
    @render.text
    def selected_file_text(): return f"Loaded: {app_data().get('filename', 'No file loaded')}"

    @output
    @render.data_frame
    def ion_table(): return app_data().get("ions", pd.DataFrame())
    @output
    @render.data_frame
    def cin_table(): return app_data().get("cin", pd.DataFrame())
    @output
    @render.data_frame
    def effluent_table(): return app_data().get("effluent", pd.DataFrame())

    model_results = reactive.Value(None)
    @reactive.Effect
    @reactive.event(input.run_button)
    def _():
        data = app_data()
        req("ions" in data, "cin" in data)
        with ui.Progress(min=1, max=15) as p:
            p.set(message="Model calculation in progress...")
            inputs_dict = {key: getattr(input, key) for key in dir(input) if not key.startswith('_')}
            results = model_prep(inputs_dict, data["ions"], data["cin"], NT_REPORT)
            model_results.set(results)
            ui.notification_show("Model run complete.", duration=5)

    # --- Alkalinity Tab Logic ---
    @reactive.Calc
    def bicarbonate_calcs():
        K1, K2, KW = 10**-6.352, 10**-10.329, 10**-14
        h_plus = 10**-input.pH()
        oh_minus = KW / h_plus
        alpha_1 = 1 / (1 + h_plus / K1 + K2 / h_plus)
        alpha_2 = 1 / (1 + h_plus / K2 + h_plus**2 / (K1 * K2))
        tot_co3_M = (input.alkvalue() / 50000 + h_plus - oh_minus) / (alpha_1 + 2 * alpha_2)
        hco3_mM_L = alpha_1 * tot_co3_M * 1000
        if hco3_mM_L < 0: return "INVALID", "INVALID", "INVALID"
        return f"{hco3_mM_L:.4f}", f"{hco3_mM_L * 12:.4f}", f"{hco3_mM_L * 61:.4f}"

    @output
    @render.text
    def bicarb_meq_L(): return bicarbonate_calcs()[0]
    @output
    @render.text
    def bicarb_mg_C_L(): return bicarbonate_calcs()[1]
    @output
    @render.text
    def bicarb_mg_HCO3_L(): return bicarbonate_calcs()[2]

    # --- kL Guesser Tab Logic ---
    @output
    @render.data_frame
    def pfas_properties_table(): return pfas_properties()

    kl_estimates = reactive.Value(pd.DataFrame())
    @reactive.Effect
    @reactive.event(input.estimate_kl)
    def _():
        df = pfas_properties().copy()
        req(not df.empty)
        t_k = input.temp() + 273.15
        viscosity = np.exp(-24.71 + (4209/t_k) + 0.04527 * t_k - (3.376e-5 * t_k**2)) / 100
        t2 = t_k / 324.65
        density = 0.98396*(-1.41768 + 8.97665*t2 - 12.2755*t2**2 + 7.45844*t2**3 - 1.73849*t2**4)
        mu1 = viscosity * 100
        
        df['kL Estimate (cm/s)'] = df.apply(lambda row:
            ( (2 + 0.644 * ( (input.Vv() / input.EBEDv()) * (2*input.rbv()) * (density / viscosity) )**(1/2) * (viscosity / density / (13.26e-5 * (mu1 ** -1.14) * (float(row["MolarVol (cm^3/mol)"]) ** -0.589)))**(1/3)) * (1 + 1.5 * (1- input.EBEDv())) ) * (13.26e-5 * (mu1 ** -1.14) * (float(row["MolarVol (cm^3/mol)"]) ** -0.589)) / (2*input.rbv()),
            axis=1
        )
        kl_estimates.set(df[['kL Estimate (cm/s)']])

    @output
    @render.data_frame
    def kl_estimates_table(): return kl_estimates().round(5)

    # --- Plotting Logic ---
    @reactive.Calc
    def processed_output():
        results = model_results()
        if results is None: return None
        t_out, x_out = results
        ions_df, cin_df, effluent_df = app_data()["ions"], app_data()["cin"], app_data()["effluent"]
        outlet_conc = x_out[:, -1, :, -1]
        df = pd.DataFrame(outlet_conc, columns=ions_df['name'])
        df['hours'] = t_out
        df_long = df.melt(id_vars='hours', var_name='name', value_name='conc_meq')
        output_unit = input.OCunits()
        if output_unit == "c/c0":
            c0_meq = cin_correct(ions_df, cin_df.iloc[[0]]).drop(columns='time').iloc[0]
            df_long['conc'] = df_long.apply(lambda row: row['conc_meq'] / c0_meq.get(row['name'], 1) if c0_meq.get(row['name'], 1) != 0 else 0, axis=1)
        else:
            df_mgl = df.copy()
            for name in df_mgl.columns:
                if name != 'hours':
                    ion_info = ions_df[ions_df['name'] == name].iloc[0]
                    df_mgl[name] *= ion_info['mw'] / ion_info['valence']
            df_long = df_mgl.melt(id_vars='hours', var_name='name', value_name='conc')
            df_long['conc'] /= MASS_CONV[output_unit]
        time_unit = input.timeunits()
        if time_unit == "Bed Volumes (x1000)":
            L_cm = input.Lv() * LENGTH_CONV[input.LengthUnits()]
            V_cms = input.Vv() * VELOCITY_CONV[input.VelocityUnits()]
            bv_sec = L_cm / V_cms if V_cms > 0 else 0
            df_long['hours'] /= (bv_sec / 3600) / 1000 if bv_sec > 0 else 1
        else:
            df_long['hours'] /= (TIME_CONV[time_unit] / 3600)
        effluent_processed = pd.DataFrame()
        if input.effluentdata() and not effluent_df.empty:
            effluent_long = effluent_df.melt(id_vars='time', var_name='name', value_name='conc').rename(columns={'time': 'hours'})
            effluent_long['name'] += "_effluent"
            effluent_processed = effluent_long
        influent_processed = pd.DataFrame()
        if input.influentdata() and not cin_df.empty:
            influent_long = cin_df.melt(id_vars='time', var_name='name', value_name='conc').rename(columns={'time': 'hours'})
            influent_long['name'] += "_influent"
            influent_processed = influent_long
        return df_long, effluent_processed, influent_processed

    @output
    @render.plotly
    def plot_counterions():
        res = processed_output()
        if res is None: return
        computed, effluent, influent = res
        counter_ions = ["CHLORIDE", "SULFATE", "NITRATE", "BICARBONATE"]
        return create_plotly(
            computed[computed['name'].isin(counter_ions)] if input.computeddata() else pd.DataFrame(),
            effluent[effluent['name'].str.contains('|'.join(counter_ions), na=False)] if input.effluentdata() else pd.DataFrame(),
            influent[influent['name'].str.contains('|'.join(counter_ions), na=False)] if input.influentdata() else pd.DataFrame(),
            "Major Inorganic Ion Concentrations", f"Concentration ({input.OCunits()})", f"Time ({input.timeunits()})"
        )

    @output
    @render.plotly
    def plot_other_ions():
        res = processed_output()
        if res is None: return
        computed, effluent, influent = res
        counter_ions = ["CHLORIDE", "SULFATE", "NITRATE", "BICARBONATE"]
        return create_plotly(
            computed[~computed['name'].isin(counter_ions)] if input.computeddata() else pd.DataFrame(),
            effluent[~effluent['name'].str.contains('|'.join(counter_ions), na=False)] if input.effluentdata() else pd.DataFrame(),
            influent[~influent['name'].str.contains('|'.join(counter_ions), na=False)] if input.influentdata() else pd.DataFrame(),
            "Additional Ionic Species Concentrations", f"Concentration ({input.OCunits()})", f"Time ({input.timeunits()})"
        )

    @session.download(filename=lambda: f"iex-output-{pd.Timestamp.now().strftime('%Y%m%d')}.xlsx")
    def save_button():
        # ... Save logic can be implemented here ...
        yield b""

app = App(app_ui, server)
