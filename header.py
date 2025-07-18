# This is the corrected UI definition for app.py

# The essential HTML for the EPA banner
epa_banner_html = """
<header class='masthead clearfix' role='banner'>
     <img alt='' class='site-logo' src='https://www.epa.gov/sites/all/themes/epa/logo.png'>
     <div class='site-name-and-slogan'>
         <h1 class='site-name'><a href='https://www.epa.gov' rel='home' title='Go to the home page'><span>US EPA</span></a></h1>
         <div class='site-slogan'>United States Environmental Protection Agency</div>
     </div>
</header>
"""

app_ui = ui.page_fluid(
    # Add the CSS file from the 'www' directory
    ui.tags.head(
        ui.tags.link(rel="stylesheet", type="text/css", href="style.css")
    ),

    # 1. Display the EPA banner
    ui.HTML(epa_banner_html),

    # 2. Add the main page title, styled correctly
    ui.div(
        ui.h1("Granular Activated Carbon Model (Python Version)", class_="page-title"),
        class_="main-column clearfix"
    ),
    
    # 3. Add the rest of your Shiny UI (theme and navbar)
    shinyswatch.theme.lumen(),
    ui.page_navbar(
        ui.nav("Input",
            # ... content of the Input tab
        ),
        ui.nav("Output",
            # ... content of the Output tab
        ),
        ui.nav("Fitted Data",
            # ... content of the Fitted Data tab
        ),
        title="GAC Adsorption Model"
    )
)
