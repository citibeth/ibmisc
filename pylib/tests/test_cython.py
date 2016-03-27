# IBMisc: Misc. Routines for IceBin (and other code)
# Copyright (c) 2013-2016 by Elizabeth Fischer
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import unittest
import numpy as np
import ibmisc
import copy

class TestCython(unittest.TestCase):
    """Tests basic passing of Numpy array data through Cython to Blitz++"""

    def setUp(self):
        pass

    def test_blitz_1d(self):
        aa = np.array([float(i) for i in range(0,5)])
        bb = copy.copy(aa)
        ibmisc.example_double_blitz(aa)
        for a,b in zip(aa,bb):
            self.assertEqual(a,b*2)

    def test_blitz_2d(self):
        aa = np.random.rand(2,3)
        # Should only work on 1-D arrays
        with self.assertRaises(Exception):
            ibmisc.example_double_blitz(aa)

    def test_blitz_discontinuous(self):
        aa2 = np.random.rand(2,3)
        aa = aa2[:,0]   # Non-unit strides
        bb = copy.copy(aa)
        self.assertNotEqual(bb.strides[0], aa.strides[0])
        ibmisc.example_double_blitz(aa)
        for a,b in zip(aa,bb):
            self.assertAlmostEqual(a,b*2)


if __name__ == '__main__':
    unittest.main()
