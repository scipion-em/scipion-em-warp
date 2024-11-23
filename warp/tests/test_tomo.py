from pyworkflow.utils import weakImport

with (weakImport("tomo")):
    from .protocols_tomo_tests import *
