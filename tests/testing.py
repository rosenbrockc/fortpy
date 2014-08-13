import os
import sys

sys.path.append("/Users/trunks/codes/fortpy-dist")
import fortpy
from fortpy.testing.tester import UnitTester

t = UnitTester("./tests/staging", False, "./tests/unittests/templates", "./fortpy/templates",
               True)
t.parser.reparse(os.path.abspath("./tests/unittests/derivative_structure_generator.f90"))
t.writeall("./tests/unittests")
