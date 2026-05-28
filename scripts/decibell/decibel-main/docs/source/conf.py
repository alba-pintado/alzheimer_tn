# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(__file__ ,"../../..")))
import sphinx_rtd_theme
import recommonmark
from recommonmark.parser import CommonMarkParser

source_parsers = {".md": CommonMarkParser}

# -- Project information -----------------------------------------------------

project = 'decibel'
copyright = '2021, Olga Ibañez-Solé and Alex M. Ascensión'
author = 'Olga Ibañez-Solé and Alex M. Ascensión'

# The full version, including alpha/beta/rc tags
release = '1.0'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['sphinx.ext.autodoc',
              "recommonmark",
              "sphinx.ext.napoleon",
              "sphinx_rtd_theme",
              'sphinx.ext.autosectionlabel',]

# Napoleon settings
autosummary_generate = True
autodoc_member_order = 'bysource'
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_use_rtype = True  # having a separate entry generally helps readability
napoleon_use_param = True
napoleon_custom_sections = [('Params', 'Parameters')]
todo_include_todos = False
autosectionlabel_prefix_document = True


# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []
source_suffix = [".rst", ".md"]

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_rtd_theme"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# Chose index.rst as the documentation file, and not contents.rst
master_doc = 'index'
