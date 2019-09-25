from amuse.community import *
from amuse.test.amusetest import TestWithMPI

from .interface import vaderInterface
from .interface import vader

class vaderInterfaceTests(TestWithMPI):
    
    def test1(self):
        instance = vader()
        instance.initialize_code()
        instance.initialize_keplerian_grid(128, True, 0.1|units.AU, 10.|units.AU, 1.|units.MSun)
        instance.stop()
    
