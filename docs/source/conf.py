import os
import sys

sys.path.insert(0, os.path.abspath('..'))

from abrat import __version__

# -- Project information -----------------------------------------------------
project = 'AbRAT'
copyright = '2025, Christoph Kreer'
author = 'Christoph Kreer'
# Dynamische Versionsnummer aus dem Package
release = __version__
version = release

# -- General configuration ---------------------------------------------------
extensions = [
    "myst_parser",
    "sphinx_design",
    "sphinx.ext.mathjax",
]

myst_enable_extensions = [
    "colon_fence",
    "substitution",
    "dollarmath",
    "amsmath",
]

myst_substitutions = {
    "AbRAT": '<span class="abrat">AbRAT</span>',
}

templates_path = ['_templates']
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------
html_theme = 'sphinx_rtd_theme'
html_theme_options = {
    'display_version': True,
}

html_static_path = ["_static",
                    "../../abrat/gui/assets/workflow.png"]

html_css_files = [
    "custom.css",
]

myst_disable_html = False
myst_allow_html = True