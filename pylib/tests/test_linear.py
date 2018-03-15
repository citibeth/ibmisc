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
from numpy.testing import *
import ibmisc
import copy



answers = \
    {False : np.array([   # force_conservation=False
        25583., -17000., -17000., -17000.,
        -17000., -17000., -17000., -17000.,
        -17000., -17000., 27183., -17000.,
        -17000., -17000., -17000., -17000.,
        -17000., -17000., -17000., -17000.,
        31983., -17000., -17000., -17000.,
        -17000., -17000., -17000., -17000.,
        -17000., -17000., -17000., -17000.,
        -17000., -17000., -17000., -17000.,
        -17000., -17000., -17000.,-17000.
    ]),
    True: np.array([
           26430.75905728, -17000., -17000., -17000., -17000.,
          -17000., -17000., -17000., -17000., -17000.,
           28083.77920705, -17000., -17000., -17000., -17000.,
          -17000., -17000., -17000., -17000., -17000.,
           33042.83965637, -17000., -17000., -17000., -17000.,
          -17000., -17000., -17000., -17000., -17000.,
          -17000., -17000., -17000., -17000., -17000.,
          -17000., -17000., -17000., -17000., -17000.       
    ])}

class TestLinear(unittest.TestCase):
    """Tests ibmisc/linear from the Python interface."""

    def test_linear(self):
        for force_conservation in (False, True):
            for linear_type in ('EIGEN', 'COMPRESSED'):
                print('-------------------------', force_conservation, linear_type)
                BvA1 = ibmisc.example_linear_weighted(linear_type)
                shape = BvA1.shape
                self.assertEqual(shape, (40,50))

                ncio = ibmisc.NcIO('__pylinear_{}.nc'.format(linear_type), 'w')
                BvA1.ncio(ncio, 'BvA')
                ncio.close()

                ncio = ibmisc.NcIO('__pylinear_{}.nc'.format(linear_type), 'r')
                BvA1 = ibmisc.nc_read_weighted(ncio, 'BvA')
                ncio.close()

                aa = np.zeros(shape[1])
                for i in range(0,aa.shape[0]):
                    aa[i] = 32.*i*i - 17.

                bb1 = BvA1.apply_M(aa, fill=-17000., force_conservation=force_conservation)
                #print(bb1)
                assert_allclose(answers[force_conservation], bb1.reshape(-1))
#                self.assertTrue(np.all(answers[force_conservation] == bb1))



if __name__ == '__main__':
    unittest.main()
