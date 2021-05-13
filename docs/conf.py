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
# import os
# import sys
import os, sys
sys.path.insert(0, os.path.abspath('.'))


# -- Project information -----------------------------------------------------

project = 'aplc-optimization'
copyright = '2020, Bryony Nickson'
author = 'Bryony Nickson'

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['sphinx_rtd_theme', 'sphinx_automodapi.automodapi',
]

numpydoc_show_class_members = False

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

master_doc = 'index'

source_suffix = {
    '.rst': 'restructuredtext',
}

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
import sphinx
import stsci_rtd_theme


def setup(app):
    app.add_stylesheet("stsci.css")

# the below is not strictly necessary but helps with extensions you may use across versions
from distutils.version import LooseVersion

on_rtd = os.environ.get('READTHEDOCS', None) == 'True'
if on_rtd:
    extensions.append('sphinx.ext.mathjax')
elif LooseVersion(sphinx.__version__) < LooseVersion('1.4'):
    extensions.append('sphinx.ext.pngmath')
else:
    extensions.append('sphinx.ext.imgmath')

html_theme = 'stsci_rtd_theme'
html_theme_path = [stsci_rtd_theme.get_html_theme_path()]

html_favicon = 'images/logo.jpg'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']