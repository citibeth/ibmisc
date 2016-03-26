/*
 * IBMisc: Misc. Routines for IceBin (and other code)
 * Copyright (c) 2013-2016 by Elizabeth Fischer
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

        int num_match = 0;
        do {
restart_loop: ;
            // Scan forward first iterator
            for (;;++i1) {
                if (i1.eof()) {
                    _eof = true;
                    return;
                }
                if (*i1 == next_match) break;
                if (*i1 > next_match) {
                    next_match = *i1;
                    break;
                }
            }

            // NOW: *i1 == next_match

#if JOIN_RANK >= 2
            for (;;++i2) {
                if (i2.eof()) {
                    _eof = true;
                    return;
                }
                if (*i2 == next_match) break;
                if (*i2 > next_match) {
                    next_match = *i2;
                    ++i1;
                    goto restart_loop;
                }
            }

            // NOW: *i1 == *i2 == next_match

#if JOIN_RANK >= 3
            for (;;++i3) {
                if (i3.eof()) {
                    _eof = true;
                    return;
                }
                if (*i3 == next_match) break;
                if (*i3 > next_match) {
                    next_match = *i3;
                    ++i1;
                    ++i2;
                    goto restart_loop;
                }
            }

            // NOW: *i1 == *i2 == *i3 == next_match
#endif
#endif
        } while (false);
