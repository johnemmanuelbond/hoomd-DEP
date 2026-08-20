# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os,sys
from unittest.mock import MagicMock
sys.path.insert(0,os.path.abspath('..'))
sys.path.insert(1,'./src')
sys.path.insert(2,'../src')
sys.path.insert(3,'../../src')

project = 'hoomd-DEP'
copyright = '2026, John E. Bond'
author = 'John E. Bond'
release = '1.0.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.intersphinx',
    'sphinx_rtd_theme',
    'sphinx.ext.napoleon',
    'nbsphinx',
]

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'hoomd': ('https://hoomd-blue.readthedocs.io/en/stable/', None),
}

templates_path = ['_templates']
exclude_patterns = []

autodoc_mock_imports = ["hoomd", "src._dep"]#, "hoomd.dep","hoomd._hover","_dep"]
# sys.modules['src._dep'] = MagicMock()


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'furo'
html_static_path = ['_static']
