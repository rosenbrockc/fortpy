import unittest as ut

class TestParameterDefinitions(ut.TestCase):
    """Unit tests the parsing functionality for a representative set of parameter
    and local variable definitions in fortran.
    """
    def setUp(self):
        """Initialize the parser instance for use by all the test cases."""
        from ancle.parsers.variable import VariableParser
        self.parser = VariableParser()
        self.cases  = [{"src": "integer, pointer    :: cList(:,:) ! INTENT(OUT): List of",
                        "res": [{"type": "integer",
                                 "kind": None,
                                 "modifiers": "pointer",
                                 "name": "cList",
                                 "default": None,
                                 "dimension": ":,:"}]},
                       {"src": "integer, dimension(size(volTable,1)) :: digCnt, digit",
                        "res": [{"type": "integer",
                                 "kind": None,
                                 "modifiers": None,
                                 "name": "digCnt",
                                 "default": None,
                                 "dimension": "size(volTable,1)"},
                                {"type": "integer",
                                 "kind": None,
                                 "modifiers": None,
                                 "name": "digit",
                                 "default": None,
                                 "dimension": "size(volTable,1)"}]},
                       {"src": "integer, allocatable :: label(:,:), a(:)",
                        "res": [{"type": "integer",
                                 "kind": None,
                                 "modifiers": "allocatable",
                                 "name": "label",
                                 "default": None,
                                 "dimension": ":,:"},
                                {"type": "integer",
                                 "kind": None,
                                 "modifiers": "allocatable",
                                 "name": "a",
                                 "default": None,
                                 "dimension": ":"}]},
                       {"src": "integer j",
                        "res": [{"type": "integer",
                                 "kind": None,
                                 "modifiers": None,
                                 "name": "j",
                                 "default": None,
                                 "dimension": ":"}]},                       
                       {"src": "type(RotPermList) :: rpl(:) ! The reduced (unique)",
                        "res": [{"type": "type",
                                 "kind": "RotPermList",
                                 "modifiers": None,
                                 "name": "rpl",
                                 "default": None,
                                 "dimension": ":"}]},
                       {"src": "real(dp) :: rd(size(pd,1),size(pd,2))",
                        "res": [{"type": "real",
                                 "kind": "dp",
                                 "modifiers": None,
                                 "name": "rd",
                                 "default": None,
                                 "dimension": "size(pd,1),size(pd,2)"}]},
                       {"src": "logical, pointer :: g(:,:) => null()",
                        "res": [{"type": "logical",
                                 "kind": None,
                                 "modifiers": "pointer",
                                 "name": "g",
                                 "default": "> null()",
                                 "dimension": ":,:"}]},
                       {"src": "character(30), public :: clusters_file = 'clustertrack.in'",
                        "res": [{}]}]

    def test_parsing(self):
        


    
