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

// -*- Mode : c++ -*-
//
// SUMMARY  :
// USAGE    :
// ORG      :
// AUTHOR   : Frederic Hecht
// E-MAIL   : hecht@ann.jussieu.fr
//

#include "config-wrapper.h"
#include <string>
#include <list>
#include <map>
typedef std::list<std::string> OneEnvironmentData;
typedef std::map<std::string, OneEnvironmentData> EnvironmentData;

extern EnvironmentData ffenvironment;
extern long verbosity;

bool EnvironmentInsert (std::string key, std::string item, std::string before);

void GetEnvironment ();
void EnvironmentLoad ();

