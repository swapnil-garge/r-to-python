@reactive.Effect
    def update_ui_from_data():
        """
        Updates all relevant UI input controls with values from the loaded
        Excel file's 'columnSpecs' and 'Fouling Data' sheets.
        This function is a direct translation of the `observe` block in the
        original R server logic.
        """
        data = app_data.get()
        if not data or "columnSpecs" not in data:
            return

        specs = data["columnSpecs"]
        fouling = data.get("fouling", pd.DataFrame())

        # Helper to safely get a value or unit from the specs DataFrame
        def get_spec(name, column='values', default=None):
            if name in specs['name'].values:
                return specs.loc[specs['name'] == name, column].iloc[0]
            return default

        # Update numeric inputs
        ui.update_numeric("prv", value=get_spec('radius', default=0.0513))
        ui.update_numeric("EPORv", value=get_spec('porosity', default=0.641))
        ui.update_numeric("pdv", value=get_spec('particleDensity', default=0.803))
        ui.update_numeric("adv", value=get_spec('apparentDensity', default=0.5))
        ui.update_numeric("psdfrv", value=get_spec('psdfr', default=5.0))
        ui.update_numeric("Lv", value=get_spec('length', default=8.0))
        ui.update_numeric("wv", value=get_spec('weight', default=8500))
        ui.update_numeric("Dv", value=get_spec('diameter', default=10.0))
        ui.update_numeric("tortuv", value=get_spec('tortuosity', default=1.0))
        
        # Update linear vs. volumetric flow inputs
        if 'v' in specs['name'].values:
            ui.update_radio_buttons("veloselect", selected="Linear")
            ui.update_numeric("Vv", value=get_spec('v', default=0.123))
        elif 'flowrate' in specs['name'].values:
            ui.update_radio_buttons("veloselect", selected="Volumetric")
            ui.update_numeric("Fv", value=get_spec('flowrate', default=500.0))

        # Update select inputs with unique choices, setting the loaded value as selected
        ui.update_select("prunits", choices=list(pd.unique([get_spec('radius', 'units', 'cm'), "cm", "m", "mm", "in", "ft"])))
        ui.update_select("LengthUnits", choices=list(pd.unique([get_spec('length', 'units', 'cm'), "cm", "ft", "m", "mm", "in"])))
        ui.update_select("wunits", choices=list(pd.unique([get_spec('weight', 'units', 'g'), "g", "kg", "lbs", "oz"])))
        ui.update_select("conc_units", choices=list(pd.unique([get_spec('units', 'values', 'ug'), "ug", "ng", "mg"])))
        ui.update_select("tunits2", choices=list(pd.unique([get_spec('time', 'values', 'days'), "days", "hours"])))
        ui.update_select("DiameterUnits", choices=list(pd.unique([get_spec('diameter', 'units', 'cm'), "cm", "ft", "mm", "m", "in"])))
        ui.update_select("VelocityUnits", choices=list(pd.unique([get_spec('v', 'units', 'cm/s')] + VELOCITY_VECTOR)))
        ui.update_select("FlowrateUnits", choices=list(pd.unique([get_spec('flowrate', 'units', 'L/min')] + FLOWRATE_VECTOR)))

        # Update fouling data
        if not fouling.empty:
            ui.update_select("WFouling", selected=fouling['WaterFouling'].iloc[0])
            ui.update_select("CFouling", selected=fouling['ChemicalFouling'].iloc[0])