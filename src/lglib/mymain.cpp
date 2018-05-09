/*
 * This file is part of FreeFem++.
 *
 * FreeFem++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * FreeFem++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
 */

/// \file

int mainff (int argc, char **argv);
namespace ffapi {void init ();
}
extern void init_ptr_parallelepmi ();

/// called by platform-dependent main() in [[file:../Graphics/sansrgraph.cpp::calling_mymain]]
int mymain (int argc, char **argv) {
	ffapi::init();	// [[file:~/ff/src/fflib/ffapi.cpp::init]]

	// Calls either [[file:~/ff/src/mpi/parallelempi.cpp::init_ptr_parallelepmi]] or
	// [[file:~/ff/src/mpi/parallelempi-empty.cpp::init_ptr_parallelepmi]]
	init_ptr_parallelepmi();

	return mainff(argc, argv);	// [[file:lg.ypp::mainff]]
}

