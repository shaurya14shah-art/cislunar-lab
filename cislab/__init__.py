"""cislab - a teaching-grade physics core for the restricted three-body problem,
general N-body dynamics and first-order relativity.

The package is deliberately split so that *nothing in here knows about rendering*.
The browser front end (``web/``) only ever talks to ``cislab.server``; every number
it draws comes from the functions in this package.
"""
__version__ = "1.0.0"
