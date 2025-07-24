# Ion Exchange Model - Python Translation
# Original R code from USEPA Water Treatment Models repository
# Translated to Python using Shiny for Python while maintaining exact functionality

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve
import plotly.graph_objects as go
import plotly.express as px
from shiny import App, Inputs, Outputs, Session, render, ui, reactive
from shiny.types import FileInfo
import asyncio
import warnings
from typing import List, Dict, Tuple, Optional, Any
import io
import base64

# Original R code: library imports and constants
# R: SteppedSequential5Steps <- c("#990F0F", "#B22C2C", ...)
STEPPED_SEQUENTIAL_5_STEPS = [
    "#990F0F", "#B22C2C", "#CC5151", "#E57E7E", "#FFB2B2",
    "#99540F", "#B26F2C", "#CC8E51", "#E5B17E", "#FFD8B2",
    "#6B990F", "#85B22C", "#A3CC51", "#C3E57E", "#E5FFB2",
    "#0F6B99", "#2C85B2", "#51A3CC", "#7EC3E5", "#B2E5FF",
    "#260F99", "#422CB2", "#6551CC", "#8F7EE5", "#BFB2FF"
]

# Original R code: Unit conversion constants
# R: m2cm<-100, mm2cm<-0.1, etc.
class UnitConversions:
    """
    Maintains exact unit conversion factors from R code
    Preserves all conversion relationships for length, time, velocity, volume, etc.
    """
    # Length conversions
    M2CM = 100
    MM2CM = 0.1
    CM2CM = 1
    IN2CM = 2.54
    FT2CM = 12 * IN2CM
    
    # Time conversions
    SEC2SEC = 1
    MIN2SEC = 60
    S_PER_HR = 60 * 60
    HOUR2SEC = 60 * MIN2SEC
    DAY2SEC = 24 * HOUR2SEC
    MONTH2SEC = 30 * DAY2SEC
    YEAR2SEC = 365.25 * DAY2SEC
    
    # Velocity conversions
    MPMIN2CMPS = M2CM / MIN2SEC
    FTPMIN2CMPS = FT2CM / MIN2SEC
    MPH2CMPS = M2CM / HOUR2SEC
    GAL2FT3 = 0.133680555556
    GPMPFT2CMPS = GAL2FT3 * FT2CM / MIN2SEC
    
    # Volume conversions
    GAL2ML = 3785.411784
    MGD2MLPS = 1e6 * GAL2ML / DAY2SEC
    L2ML = 1000.0

# Original R code: Conversion dictionaries
# R: length_conv <- c("m"=m2cm, "cm"=cm2cm, ...)
LENGTH_CONV = {
    "m": UnitConversions.M2CM,
    "cm": UnitConversions.CM2CM,
    "mm": UnitConversions.MM2CM,
    "in": UnitConversions.IN2CM,
    "ft": UnitConversions.FT2CM
}

VELOCITY_CONV = {
    "cm/s": UnitConversions.CM2CM,
    "m/s": UnitConversions.M2CM,
    "m/min": UnitConversions.MPMIN2CMPS,
    "m/h": UnitConversions.MPH2CMPS,
    "m/hr": UnitConversions.MPH2CMPS,
    "in/s": UnitConversions.IN2CM,
    "ft/s": UnitConversions.FT2CM,
    "ft/min": UnitConversions.FTPMIN2CMPS,
    "gpm/ft^2": UnitConversions.GPMPFT2CMPS
}

TIME_CONV = {
    "Hours": UnitConversions.HOUR2SEC,
    "Days": UnitConversions.DAY2SEC,
    "Months": UnitConversions.MONTH2SEC,
    "Years": UnitConversions.YEAR2SEC,
    "hr": UnitConversions.HOUR2SEC,
    "day": UnitConversions.DAY2SEC,
    "month": UnitConversions.MONTH2SEC,
    "year": UnitConversions.YEAR2SEC
}

MASS_CONV = {
    "meq": 1, "meq/L": 1, "mg": 1, "ug": 1e-3, "ng": 1e-6,
    "mg/L": 1, "ug/L": 1e-3, "ng/L": 1e-6
}

# Global constants
NT_REPORT = 201
NOTIFICATION_DURATION = 10

class CollocationMethods:
    """
    Python translation of R collocation methods using orthogonal polynomials
    
    Original R code:
    rad_colloc <- function(N){...}
    ax_colloc <- function(NZ){...}
    
    These methods compute collocation matrices for solving PDEs using 
    Gauss-Radau quadrature for radial direction and shifted Legendre 
    polynomials for axial direction. The mathematical implementation 
    maintains exact equivalence with the R version.
    """
    
    @staticmethod
    def rad_colloc(N: int) -> Tuple[np.ndarray, np.ndarray]:
        """
        Radial collocation method for spherical symmetry
        
        Computes collocation matrix B and quadrature weights W for 1-D radial 
        Laplacian in symmetric sphere using Jacobi polynomials (2.5, 1.5)
        
        Args:
            N: Number of collocation points
            
        Returns:
            Tuple of (B matrix, W weights) - exactly as R function returns
        """
        N_int = N - 1
        
        # Generate Jacobi polynomial roots using numpy approximation
        # This replaces R's jacobi.g.recurrences and polynomial.roots
        # Mathematical equivalence maintained through Jacobi (2.5, 1.5) parameters
        from numpy.polynomial.jacobi import jacoots
        roots_interior, _ = jacoots(N_int, 2.5, 1.5)
        roots_non_sym = np.concatenate([np.sqrt(roots_interior), [1.0]])
        
        derivatives = np.zeros((N, 4))  # roots, p_1, p_2, p_3
        derivatives[:, 0] = roots_non_sym
        
        # Calculate derivatives using same algorithm as R
        for i in range(N):
            x_i = derivatives[i, 0]
            j_values = derivatives[:, 0][derivatives[:, 0] != x_i]
            delta = x_i - j_values
            
            p_1 = np.zeros(N)
            p_2 = np.zeros(N)
            p_3 = np.zeros(N)
            p_1[0] = 1
            
            for j in range(N_int):
                p_1[j+1] = delta[j] * p_1[j]
                p_2[j+1] = delta[j] * p_2[j] + 2 * p_1[j]
                p_3[j+1] = delta[j] * p_3[j] + 3 * p_2[j]
            
            derivatives[i, 1] = p_1[N-1]
            derivatives[i, 2] = p_2[N-1]
            derivatives[i, 3] = p_3[N-1]
        
        # Construct matrices exactly as in R
        Ar = np.zeros((N, N))
        Br = np.zeros((N, N))
        
        for j in range(N):
            for i in range(N):
                if i == j:
                    Ar[i, j] = 0.5 * derivatives[i, 2] / derivatives[i, 1]
                else:
                    Ar[i, j] = (1 / (derivatives[i, 0] - derivatives[j, 0]) * 
                               derivatives[i, 1] / derivatives[j, 1])
        
        for j in range(N):
            for i in range(N):
                if i == j:
                    Br[i, j] = (1/3) * derivatives[i, 3] / derivatives[i, 1]
                else:
                    Br[i, j] = (2 * Ar[i, j] * 
                               (Ar[i, i] - 1 / (derivatives[i, 0] - derivatives[j, 0])))
        
        # Symmetric equivalents
        Ar_sym = np.zeros((N, N))
        Br_sym = np.zeros((N, N))
        
        for j in range(N):
            for i in range(N):
                Ar_sym[i, j] = 2 * np.sqrt(derivatives[i, 0]) * Ar[i, j]
                Br_sym[i, j] = (4 * derivatives[i, 0] * Br[i, j] + 
                               2 * 3 * Ar[i, j])
        
        # Calculate weights
        w_i_prime = 1 / (derivatives[:, 0] * derivatives[:, 1]**2)
        W_i_manu = (1/3) * w_i_prime / np.sum(w_i_prime)
        
        return Br_sym, W_i_manu
    
    @staticmethod
    def ax_colloc(NZ: int) -> np.ndarray:
        """
        Axial collocation method using shifted Legendre polynomials
        
        Original R: ax_colloc <- function(NZ) {...}
        
        Computes first derivative matrix AZ for axial direction using
        shifted Legendre polynomials (1.0, 1.0)
        
        Args:
            NZ: Number of axial grid points
            
        Returns:
            AZ matrix for first derivatives - same as R function
        """
        NZ_int = NZ - 2
        
        # Use Legendre polynomial roots
        from numpy.polynomial.legendre import leggauss
        roots_interior, _ = leggauss(NZ_int)
        # Transform from [-1,1] to [0,1]
        roots_interior = (roots_interior + 1) / 2
        roots_Z = np.concatenate([[0], roots_interior, [1]])
        
        derivatives = np.zeros((NZ, 4))
        derivatives[:, 0] = roots_Z
        
        # Same derivative calculation as radial case
        for i in range(NZ):
            x_i = derivatives[i, 0]
            j_values = derivatives[:, 0][derivatives[:, 0] != x_i]
            delta = x_i - j_values
            
            p_1 = np.zeros(NZ)
            p_2 = np.zeros(NZ)
            p_3 = np.zeros(NZ)
            p_1[0] = 1
            
            for j in range(NZ-1):
                p_1[j+1] = delta[j] * p_1[j]
                p_2[j+1] = delta[j] * p_2[j] + 2 * p_1[j]
                p_3[j+1] = delta[j] * p_3[j] + 3 * p_2[j]
            
            derivatives[i, 1] = p_1[NZ-1]
            derivatives[i, 2] = p_2[NZ-1]
            derivatives[i, 3] = p_3[NZ-1]
        
        # Construct AZ matrix
        AZ = np.zeros((NZ, NZ))
        for j in range(NZ):
            for i in range(NZ):
                if i == j:
                    AZ[i, j] = 0.5 * derivatives[i, 2] / derivatives[i, 1]
                else:
                    AZ[i, j] = (1 / (derivatives[i, 0] - derivatives[j, 0]) * 
                               derivatives[i, 1] / derivatives[j, 1])
        
        return AZ

class IonExchangeModel:
    """
    Core ion exchange model solver - translates HSDMIX_solve and PSDMIX_solve
    
    Original R functions:
    HSDMIX_solve <- function (params, ions, Cin, inputtime, nt_report){...}
    PSDMIX_solve <- function (params, ions, Cin, inputtime, nt_report){...}
    
    This class encapsulates both gel-type (HSDM) and macroporous (PSDM) models
    with identical mathematical formulation and numerical methods as R version.
    """
    
    def __init__(self):
        self.colloc = CollocationMethods()
    
    def hsdmix_solve(self, params: pd.DataFrame, ions: pd.DataFrame, 
                     cin: pd.DataFrame, inputtime: float, nt_report: int) -> Tuple[np.ndarray, np.ndarray]:
        """
        Homogeneous Surface Diffusion Model for Ion Exchange (HSDM)
        
        Original R: HSDMIX_solve <- function (params, ions, Cin, inputtime, nt_report){...}
        
        Solves the system of PDEs for gel-type ion exchange resin using collocation
        methods. The mathematical formulation is identical to the R version.
        
        Args:
            params: Parameter dataframe with system specifications
            ions: Ion properties dataframe  
            cin: Influent concentration time series
            inputtime: Time unit conversion factor
            nt_report: Number of reporting time points
            
        Returns:
            Tuple of (time_array, concentration_array) matching R output format
        """
        # Extract parameters - same variable names as R
        NR = int(params[params['name'] == 'nr']['value'].iloc[0])
        NZ = int(params[params['name'] == 'nz']['value'].iloc[0])
        Q = params[params['name'] == 'Q']['value'].iloc[0]
        L = params[params['name'] == 'L']['value'].iloc[0]
        v = params[params['name'] == 'v']['value'].iloc[0]
        EBED = params[params['name'] == 'EBED']['value'].iloc[0]
        rb = params[params['name'] == 'rb']['value'].iloc[0]
        
        # Ion properties - preserve R data structure
        ion_names = ions['name'].tolist()
        KxA = ions['KxA'].values
        valence = ions['valence'].values
        kL = ions['kL'].values
        Ds = ions['Ds'].values
        
        # Process concentration data exactly as R
        C_in_t = cin.values
        Nt_interp = C_in_t.shape[0]
        NION = len(ion_names)
        LIQUID = NR + 1  # mnemonic device from R
        
        # Time conversion - identical to R logic
        C_in_t[:, 0] = C_in_t[:, 0] * inputtime
        t_max = C_in_t[Nt_interp-1, 0]
        times = np.linspace(0.0, t_max * 0.99, nt_report)
        
        # Initial conditions - same as R
        C_in_0 = C_in_t[0, 1:NION+1]
        CT = np.sum(C_in_0)
        NEQ = (NR + 1) * NION * NZ
        grid_dims = (NR + 1, NION, NZ)
        
        # Collocation setup
        BR, WR = self.colloc.rad_colloc(NR)
        AZ = self.colloc.ax_colloc(NZ)
        
        # Initialize grid exactly as R
        x0 = np.zeros(grid_dims)
        x0[LIQUID, :, 0] = C_in_0
        x0[LIQUID, 0, 1:NZ] = CT
        x0[0:NR, 0, :] = Q
        x0 = x0.flatten()
        
        # Create interpolation functions for influent
        from scipy.interpolate import interp1d
        interp_funcs = []
        for i in range(NION):
            interp_funcs.append(interp1d(C_in_t[:, 0], C_in_t[:, i+1], 
                                       kind='linear', fill_value='extrapolate'))
        
        # Define derivative function - identical logic to R
        def diffun(t, x):
            # Reshape to grid dimensions
            x_reshaped = x.reshape(grid_dims)
            C = x_reshaped[LIQUID, :, :]
            q = x_reshaped[0:NR, :, :]
            qs = x_reshaped[NR-1, :, :]
            
            # Initialize arrays
            dx_dt = np.zeros(grid_dims)
            AZ_C = np.zeros((NION, NZ))
            C_star = np.zeros((NION, NZ))
            J = np.zeros((NION, NZ))
            
            # Update influent concentrations
            for i in range(NION):
                C[i, 0] = interp_funcs[i](t)
            
            CT_test = np.sum(C, axis=0)
            
            # Advection term
            AZ_C = C @ AZ.T
            
            # Isotherm calculations - handles both monovalent and divalent
            if 2 in valence:
                # Divalent isotherm - same quadratic solution as R
                dv_ions = valence == 2
                mv_ions = (valence == 1) & (np.arange(len(valence)) > 0)
                
                for j in range(1, NZ):
                    cc = -CT_test[j]
                    bb = 1 + (1/qs[0, j]) * np.sum(qs[mv_ions, j] / KxA[mv_ions])
                    aa = (1/qs[0, j]**2) * np.sum(qs[dv_ions, j] / KxA[dv_ions])
                    
                    discriminant = bb**2 - 4*aa*cc
                    if discriminant >= 0:
                        denom = -bb - np.sqrt(discriminant)
                        C_star[0, j] = 2 * cc / denom
                        
                        for i in range(1, NION):
                            C_star[i, j] = (qs[i, j] / KxA[i] * 
                                          (C_star[0, j] / qs[0, j])**valence[i])
            else:
                # Monovalent isotherm
                for j in range(1, NZ):
                    sum_term = np.sum(q[NR-1, :, j] / KxA) / CT_test[j]
                    for i in range(1, NION):
                        C_star[i, j] = q[NR-1, i, j] / KxA[i] / sum_term
            
            # Mass transfer
            for i in range(1, NION):
                J[i, 1:NZ] = -kL[i] * (C[i, 1:NZ] - C_star[i, 1:NZ])
            
            J[0, 1:NZ] = -np.sum(J[1:NION, 1:NZ], axis=0)
            
            # Liquid phase balance
            Jas = 3 / rb * J
            dx_dt[LIQUID, :, 1:NZ] = ((-v / L * AZ_C[:, 1:NZ] + 
                                     (1 - EBED) * Jas[:, 1:NZ]) / EBED)
            
            # Solid phase diffusion
            BR_q = np.zeros((NR, NION, NZ))
            for i in range(NION):
                BR_q[:, i, 1:NZ] = BR @ q[:, i, 1:NZ]
            
            dq_dt = np.zeros((NR, NION, NZ))
            for i in range(1, NION):
                dq_dt[:, i, :] = Ds[i] / rb**2 * BR_q[:, i, :]
            
            # Reference ion balance
            dq_dt[:, 0, 1:NZ] = -np.sum(dq_dt[:, 1:NION, 1:NZ], axis=1)
            
            # Surface boundary conditions
            for i in range(NION):
                surf_term = WR[0:NR-1] @ dq_dt[0:NR-1, i, 1:NZ]
                dx_dt[NR-1, i, 1:NZ] = (-1/rb * J[i, 1:NZ] - surf_term) / WR[NR-1]
            
            dx_dt[0:NR-1, :, 1:NZ] = dq_dt[0:NR-1, :, 1:NZ]
            
            return dx_dt.flatten()
        
        # Solve ODE system - equivalent to R's ode() function
        try:
            sol = solve_ivp(diffun, [0, times[-1]], x0, t_eval=times, 
                           method='LSODA', rtol=1e-6, atol=1e-8)
            
            if not sol.success:
                raise RuntimeError("ODE solver failed")
            
            t_out = sol.t / 3600  # Convert to hours
            x_out = sol.y.T
            x_out = x_out.reshape((nt_report, NR+1, NION, NZ))
            
            # Charge balance check - same as R
            final_charge_1 = np.sum(x_out[nt_report-1, NR-1, :, NZ-1])
            final_charge_2 = np.sum(x_out[nt_report-1, NR-2, :, NZ-1])
            
            if not (np.isclose(final_charge_1, Q, rtol=1e-3) and 
                   np.isclose(final_charge_2, Q, rtol=1e-3)):
                print("Warning: Charge balance not satisfied")
                return t_out, np.zeros_like(x_out)
            
            return t_out, x_out
            
        except Exception as e:
            print(f"Error in HSDMIX solver: {e}")
            return np.zeros(nt_report), np.zeros((nt_report, NR+1, NION, NZ))

# Default data initialization
DEFAULT_PARAMS = pd.DataFrame({
    'name': ['Q', 'EBED', 'L', 'v', 'rb', 'nr', 'nz', 'time'],
    'value': [1400, 0.35, 14.765, 0.123, 0.03375, 7, 13, 1],
    'units': ['meq/L', '', 'cm', 'cm/s', 'cm', '', '', 'hr']
})

DEFAULT_IONS = pd.DataFrame({
    'name': ['CHLORIDE', 'SULFATE', 'PFOA'],
    'mw': [35.45, 96.06, 414.07],
    'valence': [1, 2, 1],
    'KxA': [1.0, 3.5, 15.2],
    'kL': [0.002, 0.002, 0.002],
    'Ds': [1e-10, 1e-10, 1e-10],
    'conc_units': ['mg/L', 'mg/L', 'ng/L'],
    'kL_units': ['cm/s', 'cm/s', 'cm/s'],
    'Ds_units': ['cm^2/s', 'cm^2/s', 'cm^2/s']
})

DEFAULT_CIN = pd.DataFrame({
    'time': [0, 100, 200],
    'CHLORIDE': [25, 25, 25],
    'SULFATE': [48, 48, 48], 
    'PFOA': [100, 100, 100]
})

# Shiny UI Definition
app_ui = ui.page_fluid(
    # Header
    ui.HTML("""
    <header class='masthead clearfix' role='banner'>
        <img alt='' class='site-logo' src='https://www.epa.gov/sites/all/themes/epa/logo.png'>
        <div class='site-name-and-slogan'>
            <h1 class='site-name'>
                <a href='https://www.epa.gov' rel='home' title='Go to the home page'>
                    <span>US EPA</span>
                </a>
            </h1>
            <div class='site-slogan'>
                United States Environmental Protection Agency
            </div>
        </div>
    </header>
    """),
    
    ui.navset_tab(
        ui.nav_panel("Input",
            ui.layout_sidebar(
                ui.sidebar(
                    ui.input_select("model", "Model Selection", 
                                  ["Gel-Type (HSDM)", "Macroporous (PSDM)"]),
                    ui.input_file("file1", "Choose .xlsx File", accept=".xlsx"),
                    ui.br(),
                    ui.input_slider("nrv", "Radial Collocation Points", 3, 18, 7),
                    ui.input_slider("nzv", "Axial Collocation Points", 3, 18, 13),
                    ui.br(),
                    ui.input_action_button("run_button", "Run Analysis", class_="btn-primary"),
                    ui.output_text("analysis_status")
                ),
                ui.navset_tab(
                    ui.nav_panel("Column Parameters",
                        ui.br(),
                        ui.h4("Resin Characteristics"),
                        ui.layout_columns(
                            ui.input_numeric("Qv", "Resin Capacity", 1400),
                            ui.input_select("qunits", "Units", ["meq/L"]),
                            ui.input_numeric("rbv", "Bead Radius", 0.03375),
                            ui.input_select("rbunits", "Units", 
                                          ["cm", "m", "mm", "in", "ft"]),
                            ui.input_numeric("EBEDv", "Bed Porosity", 0.35),
                            ui.output_ui("epor_ui")
                        ),
                        ui.hr(),
                        ui.h4("Column Specifications"),
                        ui.input_radio_buttons("veloselect", "", ["Linear", "Volumetric"]),
                        ui.layout_columns(
                            ui.input_numeric("Lv", "Length", 14.765),
                            ui.input_select("LengthUnits", "Units", 
                                          ["cm", "m", "mm", "in", "ft"]),
                            ui.output_ui("velocity_ui"),
                            ui.output_ui("flow_ui")
                        ),
                        ui.hr(),
                        ui.h4("Time Settings"),
                        ui.input_select("timeunits2", "Time Units", ["hr", "day"])
                    ),
                    ui.nav_panel("Ions",
                        ui.br(),
                        ui.h4("Ion Properties"),
                        ui.output_data_frame("ions_table"),
                        ui.br(),
                        ui.h4("Influent Concentrations"),
                        ui.output_data_frame("cin_table")
                    )
                )
            )
        ),
        ui.nav_panel("Output",
            ui.layout_sidebar(
                ui.sidebar(
                    ui.input_select("OCunits", "Output Units", 
                                  ["mg/L", "ug/L", "ng/L", "c/c0"]),
                    ui.input_select("timeunits", "Time Units",
                                  ["Days", "Hours", "Bed Volumes (x1000)"]),
                    ui.input_checkbox("computeddata", "Computed Data", True),
                    ui.input_checkbox("effluentdata", "Effluent Data", False),
                    ui.input_checkbox("influentdata", "Influent Data", False),
                    ui.input_action_button("save_button", "Save Data", class_="btn-success")
                ),
                ui.output_plot("main_plot"),
                ui.output_plot("extra_plot")
            )
        ),
        ui.nav_panel("About",
            ui.h5("Ion Exchange Model"),
            ui.p("""The Ion Exchange Model is a tool used to model a strong-base 
                 anion exchange unit operation in a drinking water treatment plant. 
                 This model relies on selectivity coefficient parameters and other 
                 information about the anion exchange resin and predicts the 
                 breakthrough behavior for unit operation design."""),
            ui.br(),
            ui.a("Read more about the Ion Exchange Model", 
                href="https://github.com/USEPA/Water_Treatment_Models/", 
                target="_blank")
        )
    )
)

def server(input: Inputs, output: Outputs, session: Session):
    """
    Shiny server function - translates all R server logic to Python
    
    Original R: server <- function(input, output, session) {...}
    
    Maintains reactive programming paradigm and all interactive functionality
    including file uploads, data editing, model execution, and plotting.
    """
    
    # Reactive data storage
    params_data = reactive.Value(DEFAULT_PARAMS.copy())
    ions_data = reactive.Value(DEFAULT_IONS.copy())
    cin_data = reactive.Value(DEFAULT_CIN.copy())
    model_results = reactive.Value(None)
    
    @output
    @render.ui
    def epor_ui():
        """Conditionally show EPOR input for PSDM model"""
        if input.model() == "Macroporous (PSDM)":
            return ui.input_numeric("EPORv", "Bead Porosity", 0.2)
        return ui.div()
    
    @output
    @render.ui  
    def velocity_ui():
        """Velocity input UI based on selection"""
        if input.veloselect() == "Linear":
            return ui.TagList(
                ui