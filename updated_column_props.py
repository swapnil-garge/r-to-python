ui.nav("Column Parameters",
    ui.h4(ui.strong("Media Characteristics")),
    ui.row(
        ui.column(4, ui.input_numeric("prv", "Particle Radius", 0.0513)),
        ui.column(4, ui.input_select("prunits", "Units", ["cm", "m", "mm", "in", "ft"])),
    ),
    ui.row(
        ui.column(4, ui.input_numeric("EPORv", "Bed Porosity", 0.641)),
    ),
    ui.row(
        ui.column(4, ui.input_numeric("pdv", "Particle Density", 0.803)),
        ui.column(4, ui.input_select("pdunits", "Units", ["g/ml"])),
    ),
    ui.row(
        ui.column(4, ui.input_numeric("adv", "Apparent Density", 0.5)),
        ui.column(4, ui.input_select("adunits", "Units", ["g/ml"])),
    ),
    ui.row(
        ui.column(4, ui.input_numeric("psdfrv", "PSDFR", 5.0)),
    ),
    ui.hr(),
    ui.h4(ui.strong("Column Specifications")),
    ui.row(
        ui.column(4, ui.input_numeric("Lv", "Length", 8.0)),
        ui.column(4, ui.input_select("LengthUnits", "Units", ["cm", "ft", "m", "mm", "in"])),
    ),
    ui.row(
        ui.column(8, ui.input_radio_buttons("veloselect", "Flow Specification", choices=["Volumetric", "Linear"], selected="Volumetric", inline=True)),
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
    ),
    ui.row(
        ui.column(4, ui.input_select("tunits2", "Time Units", ["days", "hours"])),
    ),
)