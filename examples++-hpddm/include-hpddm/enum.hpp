/*
   This file is part of HPDDM.

   Author(s): Pierre Jolivet <pierre.jolivet@enseeiht.fr>
        Date: 2012-10-04

   Copyright (C) 2011-2014 Université de Grenoble
                 2015      Eidgenössische Technische Hochschule Zürich
                 2016-     Centre National de la Recherche Scientifique

   HPDDM is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published
   by the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   HPDDM is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with HPDDM.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _HPDDM_ENUM_
#define _HPDDM_ENUM_

namespace HPDDM {
/* Enum: FetiPrcndtnr
 *
 *  Defines the FETI preconditioner used in the projection.
 *
 * NONE         - No preconditioner.
 * SUPERLUMPED  - Approximation of the local Schur complement by the diagonal of <Schur::bb>.
 * LUMPED       - Approximation of the local Schur complement by <Schur::bb>.
 * DIRICHLET    - Local Schur complement.
 *
 * See also: <Feti>. */
enum class FetiPrcndtnr : char {
    NONE, SUPERLUMPED, LUMPED, DIRICHLET
};
} // HPDDM
#endif // _HPDDM_ENUM_
