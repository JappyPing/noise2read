# -*- coding: utf-8 -*-
# @Author: Pengyao Ping
# @Date:   2023-02-16 11:02:09
# @Last Modified by:   Pengyao Ping
# @Last Modified time: 2023-04-27 17:38:55
# Configuration file for the Sphinx documentation builder.
import os
import sys
from typing import Any, Dict
# -- Project information

project = 'noise2read'
copyright = '2023, Pengyao Ping'
author = 'Pengyao Ping'

release = '0.1'
version = '0.1.0'

# -- General configuration ---------------------------------------------------
#

extensions = [
    # Sphinx's own extensions
    "sphinx.ext.autodoc",
    "sphinx.ext.extlinks",
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
    "sphinx.ext.todo",
    "sphinx.ext.viewcode",
    # Our custom extension, only meant for Furo's own documentation.
    "furo.sphinxext",
    # External stuff
    "myst_parser",
    "sphinx_copybutton",
    "sphinx_design",
    "sphinx_inline_tabs",
]
templates_path = ["_templates"]

# -- Options for extlinks ----------------------------------------------------
#

extlinks = {
    "pypi": ("https://pypi.org/project/%s/", "%s"),
}

# -- Options for intersphinx -------------------------------------------------
#

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "sphinx": ("https://www.sphinx-doc.org/en/master", None),
}

# -- Options for TODOs -------------------------------------------------------
#

todo_include_todos = True

# -- Options for Markdown files ----------------------------------------------
#

myst_enable_extensions = [
    "colon_fence",
    "deflist",
]
myst_heading_anchors = 3

# -- Options for HTML output -------------------------------------------------
#

html_theme = "furo"
html_title = "ShortReadsCorrection"
language = "en"
html_logo = "../../logo/logo.svg"

html_static_path = ["_static"]
# html_css_files = ["pied-piper-admonition.css"]

html_theme_options = {
    "sidebar_hide_name": False,
    # "light_logo": "logo_transparent.svg",
    # "dark_logo": "logo_dark.svg",
    "light_css_variables": {
        "color-brand-primary": "#336790",  # "blue"
        "color-brand-content": "#336790",
    },
    "dark_css_variables": {
        "color-brand-primary": "#E5B62F",  # "yellow"
        "color-brand-content": "#E5B62F",
    },
}

html_theme_options: Dict[str, Any] = {
    "footer_icons": [
        {
            "name": "GitHub",
            "url": "https://github.com/Jappy0/noise2read",
            "html": """
                <svg stroke="currentColor" fill="currentColor" stroke-width="0" viewBox="0 0 20 20">
                    <path fill-rule="evenodd" d="M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38 0-.19-.01-.82-.01-1.49-2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13-.28-.15-.68-.52-.01-.53.63-.01 1.08.58 1.23.82.72 1.21 1.87.87 2.33.66.07-.52.28-.87.51-1.07-1.78-.2-3.64-.89-3.64-3.95 0-.87.31-1.59.82-2.15-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82.64-.18 1.32-.27 2-.27.68 0 1.36.09 2 .27 1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92.08 2.12.51.56.82 1.27.82 2.15 0 3.07-1.87 3.75-3.65 3.95.29.25.54.73.54 1.48 0 1.07-.01 1.93-.01 2.2 0 .21.15.46.55.38A8.013 8.013 0 0 0 16 8c0-4.42-3.58-8-8-8z"></path>
                </svg>
            """,
            "class": "",
        },
    ],
    "source_repository": "https://github.com/Jappy0/noise2read",
    "source_branch": "master",
    "source_directory": "docs/",
}

# if "READTHEDOCS" in os.environ:
#     html_theme_options["announcement"] = (
#         "This documentation is hosted on Read the Docs only for testing. Please use "
#         "<a href='https://pradyunsg.me/furo/'>the main documentation</a> instead."
#     )

# -- Options for theme development -------------------------------------------
# Make sure these are all set to the default values.

html_js_files = []
html_context: Dict[str, Any] = {}
# html_show_sphinx = False
# html_show_copyright = False
# html_last_updated_fmt = ""

# RTD_TESTING = False
# if RTD_TESTING or "FURO_RTD_TESTING" in os.environ:
#     del html_theme_options["footer_icons"]

#     html_css_files += [
#         "https://assets.readthedocs.org/static/css/readthedocs-doc-embed.css",
#         "https://assets.readthedocs.org/static/css/badge_only.css",
#     ]
#     html_js_files += [
#         "readthedocs-dummy.js",
#         "https://assets.readthedocs.org/static/javascript/readthedocs-doc-embed.js",
#     ]
#     html_context["READTHEDOCS"] = True
#     html_context["current_version"] = "latest"
#     html_context["conf_py_path"] = "/docs/"
#     html_context["display_github"] = True
#     html_context["github_user"] = "pradyunsg"
#     html_context["github_repo"] = "furo"
#     html_context["github_version"] = "main"
#     html_context["slug"] = "furo"

# FONT_AWESOME_TESTING = False
# if FONT_AWESOME_TESTING:
#     html_css_files += [
#         "https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.0.0/css/fontawesome.min.css",
#         "https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.0.0/css/solid.min.css",
#         "https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.0.0/css/brands.min.css",
#     ]

#     html_theme_options["footer_icons"] = [
#         {
#             "name": "GitHub",
#             "url": "https://github.com/pradyunsg/furo",
#             "html": "",
#             "class": "fa-solid fa-github fa-2x",
#         },
#     ]