# -*- coding: utf-8 -*-
#
# Configuration file for the Sphinx documentation builder.
#
# This file does only contain a selection of the most common options. For a
# full list see the documentation:
# http://www.sphinx-doc.org/en/master/config

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import shutil
import subprocess
import sys
from pathlib import Path

DOCS_DIR = Path(__file__).resolve().parents[1]
REPO_ROOT = DOCS_DIR.parents[1]
SOURCE_ROOT = REPO_ROOT / 'src'
BUILD_EXTENSION_ROOT = REPO_ROOT / '_build' / 'src'

if BUILD_EXTENSION_ROOT.is_dir():
    # The build tree supplies liblibra_core and its legacy split extensions
    # (libutil, liblinalg, ...), which many Python modules import directly.
    sys.path.insert(0, str(BUILD_EXTENSION_ROOT))
sys.path.insert(1 if BUILD_EXTENSION_ROOT.is_dir() else 0, str(SOURCE_ROOT))

# Load the package using the complete local build environment, then prefer the
# live source directory for all libra_py submodules so docs never lag edits.
import libra_py
libra_py.__path__.insert(0, str(SOURCE_ROOT / 'libra_py'))


# -- Project information -----------------------------------------------------

project = 'Libra'
copyright = '2019–2026, Libra development team'
author = 'Libra development team'

# The short X.Y version
version = '1.0'
# The full version, including alpha/beta/rc tags
release = '1.0.0'


# -- General configuration ---------------------------------------------------

# If your documentation needs a minimal Sphinx version, state it here.
#
# needs_sphinx = '1.0'

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.intersphinx',
    'sphinx.ext.mathjax',
    'sphinx.ext.ifconfig',
    'sphinx.ext.viewcode',
    'sphinx.ext.githubpages',
    'sphinx.ext.napoleon',
    'sphinx.ext.autosummary',
    'sphinx.ext.extlinks',
    'myst_parser'
]


def _generate_python_api_reference():
    """Regenerate recursive ``libra_py`` stubs before reading sources."""
    destination = Path(__file__).parent / 'reference' / 'generated'
    if destination.exists():
        shutil.rmtree(destination)
    destination.mkdir(parents=True)
    command = [
        sys.executable, '-m', 'sphinx.ext.apidoc',
        '-q', '--force', '--module-first', '--separate',
        '--maxdepth', '5', '--output-dir', str(destination),
        str(SOURCE_ROOT / 'libra_py'),
        str(SOURCE_ROOT / 'libra_py' / 'dynamics' / 'tsh' / 'recipes'),
        str(SOURCE_ROOT / 'libra_py' / 'workflows' / 'librax' / 'md.py'),
    ]
    subprocess.run(command, check=True)

    recipe_root = SOURCE_ROOT / 'libra_py' / 'dynamics' / 'tsh' / 'recipes'
    recipe_modules = sorted(
        path.stem for path in recipe_root.glob('*.py')
        if path.name != '__init__.py'
    )
    catalog = [
        'Surface-hopping recipe catalog',
        '==============================',
        '',
        f'This catalog contains all {len(recipe_modules)} generated configuration ',
        'modules in ``libra_py.dynamics.tsh.recipes``. Each module exposes the ',
        'same public function:',
        '',
        '``load(dyn_general)``',
        '   Update and return the supplied dynamics-parameter dictionary with ',
        '   the method combination encoded by the module name.',
        '',
        '.. hlist::',
        '   :columns: 3',
        '',
    ]
    catalog.extend(
        f'   * ``libra_py.dynamics.tsh.recipes.{name}``'
        for name in recipe_modules
    )
    (destination / 'recipe_catalog.rst').write_text(
        '\n'.join(catalog) + '\n', encoding='utf-8'
    )


_generate_python_api_reference()

language = 'en'

myst_enable_extensions = [
    "colon_fence",
    "deflist",
]


# Napoleon settings
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = False
napoleon_use_admonition_for_notes = False
napoleon_use_admonition_for_references = False
#napoleon_use_ivar = False
napoleon_use_param = True
napoleon_use_rtype = True
autosummary_generate = True
napoleon_use_ivar = True
autodoc_default_options = {
    'members': True,
    'undoc-members': True,
    'show-inheritance': True,
    'member-order': 'bysource',
}
autodoc_typehints = 'description'
autodoc_mock_imports = [
    'dftbplus', 'hippynn', 'lammps', 'mpi4py', 'openbabel', 'psi4', 'pyscf',
    'qchem', 'tensorflow',
]


# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
#
# source_suffix = ['.rst', '.md']
source_suffix = ['.rst', '.md']

# The master toctree document.
master_doc = 'index'

# The language for content autogenerated by Sphinx. Refer to documentation
# for a list of supported languages.
#
# This is also used if you do content translation via gettext catalogs.
# Usually you set "language" from the command line for these cases.
# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = [
    'reference/generated/**/__pycache__',
    'reference/libra_py.rst',
    'reference/libra_py/**',
]

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'default'



# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
#html_theme = 'alabaster'
html_theme = 'sphinx_rtd_theme'

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
#
html_theme_options = {
    'collapse_navigation': True,
    'navigation_depth': 2,
    'sticky_navigation': True,
    'titles_only': True,
}
html_title = 'Libra Documentation'
html_short_title = 'Libra'
html_show_sourcelink = True

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']
html_css_files = ['custom.css']

# Custom sidebar templates, must be a dictionary that maps document names
# to template names.
#
# The default sidebars (for documents that don't match any pattern) are
# defined by theme itself.  Builtin themes are using these templates by
# default: ``['localtoc.html', 'relations.html', 'sourcelink.html',
# 'searchbox.html']``.
#
# html_sidebars = {}


# -- Options for HTMLHelp output ---------------------------------------------

# Output file base name for HTML help builder.
htmlhelp_basename = 'Libra-Documentationdoc'


# -- Options for LaTeX output ------------------------------------------------

latex_elements = {
    # The paper size ('letterpaper' or 'a4paper').
    #
    # 'papersize': 'letterpaper',

    # The font size ('10pt', '11pt' or '12pt').
    #
    # 'pointsize': '10pt',

    # Additional stuff for the LaTeX preamble.
    #
    # 'preamble': '',

    # Latex figure (float) alignment
    #
    # 'figure_align': 'htbp',
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
latex_documents = [
    (master_doc, 'Libra.tex', u'Libra Documentation',
     u'Alexey V. Akimov', 'manual'),
]


# -- Options for manual page output ------------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [
    (master_doc, 'libra', u'Libra Documentation',
     [author], 1)
]


# -- Options for Texinfo output ----------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
    (master_doc, 'Libra', u'Libra Documentation',
     author, 'Libra', 'One line description of project.',
     'Miscellaneous'),
]


# -- Options for Epub output -------------------------------------------------

# Bibliographic Dublin Core info.
epub_title = project

# The unique identifier of the text. This can be a ISBN number
# or the project homepage.
#
# epub_identifier = ''

# A unique identification for the text.
#
# epub_uid = ''

# A list of files that should not be packed into the epub file.
epub_exclude_files = ['search.html']


# -- Extension configuration -------------------------------------------------

# -- Options for intersphinx extension ---------------------------------------

# Example configuration for intersphinx: refer to the Python standard library.
intersphinx_mapping = {}
