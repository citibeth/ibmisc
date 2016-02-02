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
