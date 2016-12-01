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

# This Python2 script is meant to be run by OTHER libraries that need
# to interface between Fortran 90 and Blitz/C++ arrays.
# It generates the required Fortran interface code.

from __future__ import print_function
import string
import sys

# ---------------------------------------------------
arr_spec_tpl = string.Template('''
    type, bind(c) :: arr_spec_$rank
        type(c_ptr) :: base
        type(c_ptr) :: deltas($rank)
        integer(c_int) :: lbounds($rank)
        integer(c_int) :: ubounds($rank)
    end type''')

# ---------------------------------------------------
c_loc_tpl = string.Template('''
    function c_loc_$ctype(x)
        use, intrinsic :: iso_c_binding
        $ftype, target :: x
        type(c_ptr) :: c_loc_$ctype

        c_loc_$ctype = c_loc(x)
    end function
''')

# ---------------------------------------------------
get_spec_per_dim_tpl = string.Template('''

        ! ------- Dimension $dim
        spec%lbounds($dim) = low$dim
        spec%ubounds($dim) = low$dim + ubound(arr,$dim) - 1
        if (spec%lbounds($dim) < spec%ubounds($dim)) then
            spec%deltas($dim) = c_loc_$ctype(arr($deltas))
        else
            spec%deltas($dim) = spec%base
        end if''')

def get_spec_per_dim_body(ctype, rank) :
    out = []
    for dim in range(1,rank+1) :
        deltas = ['lbound(arr,%d)%s' % (d, '+1' if d == dim else '')
            for d in range(1,rank+1)]
        deltas = string.join(deltas, ',')
        out.append(get_spec_per_dim_tpl.substitute(ctype=ctype,
            dim=str(dim), deltas=deltas))
    return string.join(out, '')
# ---------------------------------------------------
get_spec_tpl = string.Template('''
    subroutine get_spec_${ctype}_$rank(arr, $lows, spec)
    implicit none
    $ftype, dimension($colons), target :: arr
    integer :: $lows
    type(arr_spec_$rank) :: spec

        spec%base = c_loc_$ctype( arr($lbounds) )
        $body
    end subroutine''')


def get_spec(ftype,ctype,rank) :
    d = {'ftype' : ftype, 'ctype' : ctype, 'rank' : str(rank)}
    d['lbounds'] = string.join(['lbound(arr,%d)' % i for i in range(1,rank+1)], ',')
    d['lows'] = string.join(['low%d' % i for i in range(1,rank+1)], ',')
    d['body'] = get_spec_per_dim_body(ctype, rank)
    d['colons'] = string.join([':' for i in range(1,rank+1)], ',')
    return get_spec_tpl.substitute(d)

# ---------------------------------------------------
def get_module(module_name, fctypes, ranks) :
    out = []
    out.append('''! ======== DO NOT EDIT!!!!  Machine Generated!!!!
module %s\n''' % module_name)
    out.append('''

use, intrinsic :: iso_c_binding
implicit none
    ''')

    for rank in ranks :
        out.append(arr_spec_tpl.substitute(rank=str(rank)))

    out.append('\ncontains\n')


    for (ftype,ctype) in fctypes :
        out.append(c_loc_tpl.substitute(ctype=ctype,ftype=ftype))
        for rank in ranks :
            out.append('\t! ================ Type (%s, %s), Rank %d' % (ftype, ctype, rank))
            out.append(get_spec(ftype,ctype,rank))
            out.append('\n\n')

    out.append('\nend module %s\n' % module_name)

    return string.join(out, '')
# ---------------------------------------------------
# ======================================================
get_spec_macro_tpl = string.Template('#define GET_SPEC_${uctype}_$rank(arr, spec) get_spec_${ctype}_$rank(arr, $lbounds, spec)')

def get_macros(fctypes, ranks) :
    out = []

    for (ftype,ctype) in fctypes :
        for rank in ranks :

            lbounds = string.join(['lbound(arr,%d)' % i for i in range(1,rank+1)], ',')
            out.append(get_spec_macro_tpl.substitute(ctype=ctype,uctype=string.upper(ctype),rank=str(rank),lbounds=lbounds))
    return string.join(out,'\n')
# ---------------------------------------------------

if len(sys.argv) < 4 :
    print('Usage: %s <module_name> <output.f90> <output.h>')
    sys.exit(-1)

# Definitions will be generate for the cross product of fctypes and ranks.

# The list of Fortran/C tyes to generate definitions.
# Format is (fortran-type, c-type)
fctypes = [
    ('real*8', 'double'),
    ('integer', 'int')]
ranks = [1,2,3]

module = get_module(sys.argv[1], fctypes, ranks)
macros = get_macros(fctypes, ranks)

#f90blitz_f_f90 = open('giss/f90blitz_f.f90', 'w')
f90blitz_f_f90 = open(sys.argv[2], 'w')

f90blitz_f_f90.write(module)
f90blitz_f_f90.close()

#f90blitz_f_h = open('../include/giss/f90blitz_f.h', 'w')
f90blitz_f_h = open(sys.argv[3], 'w')
f90blitz_f_h.write(macros)
f90blitz_f_h.close()
