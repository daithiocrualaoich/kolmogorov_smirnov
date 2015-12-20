project = u'Kolmogorov-Smirnov'
copyright = u'2015, Daithi O Crualaoich'
author = u'Daithi O Crualaoich'

version = '1.0.1'
release = '1.0.1'

extensions = [
  'sphinx.ext.pngmath',
  'sphinxcontrib.googleanalytics'
]

templates_path = ['_templates']
exclude_patterns = ['_build']

source_suffix = '.rst'
master_doc = 'index'

language = None
pygments_style = 'sphinx'

googleanalytics_id = 'UA-71626319-1'

# -- Options for HTML output ----------------------------------------------

html_theme = 'alabaster'
html_title = project

html_static_path = ['_static']

html_sidebars = {
   '**': ['localtoc.html'],
}

html_use_index = False
html_show_sourcelink = False
html_show_sphinx = False
html_show_copyright = False
