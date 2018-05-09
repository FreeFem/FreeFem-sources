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

#ifndef LABEL_HPP
#define LABEL_HPP

// <<Label>>
class Label {	// reference number for the physics
	friend inline ostream &operator << (ostream &f, const Label &r)
	{f << r.lab; return f;}

	friend inline istream &operator >> (istream &f, Label &r)
	{f >> r.lab; return f;}

	public:
		int lab;
		Label (int r = 0): lab(r) {}

		bool onGamma () const {return lab;}

		int operator ! () const {return !lab;}

		bool operator < (const Label &r) const {return lab < r.lab;}

		bool operator == (const Label &r) const {return lab == r.lab;}

		bool operator > (const Label &r) const {return lab > r.lab;}
};

#endif

