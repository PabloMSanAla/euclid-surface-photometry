# Configuration file for the Sphinx documentation builder.
#
# Euclid Surface Photometry Pipeline – Sphinx configuration
# Sánchez-Alarcón et al. (2026), A&A (Euclid Q1 paper)

import os
import sys

# -- Path setup --------------------------------------------------------------
sys.path.insert(0, os.path.abspath('..'))

# -- Project information -----------------------------------------------------
project = 'Euclid Surface Photometry Pipeline'
author  = 'P. M. Sánchez-Alarcón et al. (Euclid Collaboration)'
copyright = '2026, Euclid Collaboration'
release  = '0.1.2'

# -- General configuration ---------------------------------------------------
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx.ext.mathjax',
    'sphinx.ext.intersphinx',
    'sphinx.ext.autosectionlabel',
]

templates_path   = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# -- Options for HTML output -------------------------------------------------
html_theme        = 'sphinx_rtd_theme'
html_static_path  = ['_static']
html_logo         = '_static/euclid_logo.png'  # place logo here if available

# -- Napoleon (NumPy / Google docstrings) ------------------------------------
napoleon_google_docstring = False
napoleon_numpy_docstring  = True
napoleon_use_param        = True
napoleon_use_rtype        = True

# -- AutoDoc -----------------------------------------------------------------
autodoc_member_order    = 'bysource'
autodoc_default_options = {
    'members': True,
    'undoc-members': False,
    'show-inheritance': True,
}

# -- Intersphinx mapping -----------------------------------------------------
intersphinx_mapping = {
    'python':  ('https://docs.python.org/3', None),
    'numpy':   ('https://numpy.org/doc/stable/', None),
    'astropy': ('https://docs.astropy.org/en/stable/', None),
    'scipy':   ('https://docs.scipy.org/doc/scipy/', None),
}

# -- MathJax -----------------------------------------------------------------
mathjax3_config = {
    'tex': {
        'macros': {
            'IE': r'I_\mathrm{E}',
            'YE': r'Y_\mathrm{E}',
            'JE': r'J_\mathrm{E}',
            'HE': r'H_\mathrm{E}',
        }
    }
}
