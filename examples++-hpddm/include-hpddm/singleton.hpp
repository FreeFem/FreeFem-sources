 /*
   This file is part of HPDDM.

   Author(s): Pierre Jolivet <pierre.jolivet@enseeiht.fr>
        Date: 2016-01-11

   Copyright (C) 2016-     Centre National de la Recherche Scientifique

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

#ifndef _HPDDM_SINGLETON_
#define _HPDDM_SINGLETON_

namespace HPDDM {
/* Class: Singleton
 *  A base class for creating singletons. */
class Singleton {
    protected:
        template<int N>
        class construct_key { };
        Singleton() { }
        Singleton(const Singleton&) = delete;
        template<class T, int N, class... Args>
        static std::shared_ptr<T> get(const Args&... arguments) {
            static std::shared_ptr<T> instance = std::make_shared<T>(construct_key<N>(), arguments...);
            return instance;
        }
};
} // HPDDM
#endif // _HPDDM_SINGLETON_
