# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "hairpin2"
copyright = "2025, Alex Byrne"
author = "Alex Byrne"
release = "3.0.0"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
]

source_suffix = {".rst": "restructuredtext", ".md": "markdown"}

templates_path = ["_templates"]
exclude_patterns = ["hairpin2/infrastructure/*", "hairpin2/VERSION.py"]

add_module_names = False
# shorter type hints (e.g. List[int] not typing.List[int])
python_use_unqualified_type_names = True

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]

# RTD theme options for consistent sidebar
html_theme_options = {
    "navigation_depth": 4,
    "includehidden": False,
    "titles_only": False,
    "collapse_navigation": False,
    "sticky_navigation": True,
    "navigation_with_keys": False,
    "style_external_links": True,
}

autodoc_default_options = {
    "members": True,
    "undoc-members": False,
    "show-inheritance": False,
    "private-members": False,
}
add_module_names = False
autodoc_member_order = "bysource"

# Ensure consistent sidebar across all pages
html_sidebars = {
    "**": [
        "globaltoc.html",
        "relations.html",
        "sourcelink.html",
        "searchbox.html",
    ],
}
