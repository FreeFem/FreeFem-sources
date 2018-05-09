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

#include <direct.h>
#include <string.h>
#include <cstdlib>
#include <string>
#include <cstring>
#include <iostream>
#include <tchar.h>
#include <mmc.h>
using namespace std;

const string CB = "\"", CE = "\"";
const char SLACH = '/';
const char BACKSLACH = '\\';
const char dirsep = BACKSLACH, dirnsep = SLACH;
string escapepath (string f) {
	string r;

	for (int i = 0; i < f.length(); ++i) {
		if (f[i] == ' ' || f[i] == '(' || f[i] == ')' || f[i] == '=')
			r.append("^");

		r.append(f.begin() + i, f.begin() + i + 1);
	}

	return r;
}

string MyGetModuleFileName () {
	// const size_t lx=2048;
	// char b[lx+1];
	char fullPath[1024];
	string r, f;
	int lp = 0;
	int l = GetModuleFileName(NULL, fullPath, 1024);

	if (l > 0) {
		f = fullPath;
		r = escapepath(f);
	}

	// cout << l <<" MyGetModuleFileName " << fullPath << "\n"
	// << r << "\n" << f << "*****"<< endl;
	return r;
}

string DirName (const string &fn) {
	// cout << " dir ff++ " << fn<< " ::" << dirsep << endl;

	size_t s = fn.rfind(dirsep);

	if (s == string::npos) return string();

	// cout<< " **" << string(fn,0,s) << " s " << s << endl;
	return string(fn, 0, s) + dirsep;
}

int main (int argc, const char **argv) {
	string SC = " ";
	string addargs = SC + "-wait" + SC + "-log";
	int debug = 0;
	string pp;
	string cmd;

	cmd += DirName(MyGetModuleFileName());
	cmd += "FreeFem++.exe";

	if (debug) cout << "cmd:: " << cmd << endl;

	if (argc <= 1) {
		cerr << " Sorry no file name " << endl;
		cerr << " Drag and Drop the file icon on the application  icon or double clip on script file"
		     << endl;
		cmd += addargs;
		int ret = system(cmd.c_str());
		return 0;
	}

	for (int i = 1; i < argc; ++i) {
		if (strcmp("++d", argv[i]) == 0)
			debug = 1;
		else {
			cmd += SC;
			cmd += CB;
			cmd += argv[i];
			cmd += CE;
			if (!pp.length() && strlen(argv[i]) > 2)
				if (argv[i][1] == ':')
					pp = argv[i];

			;
			if (debug) cout << "  ffl: arg " << i << argv[i] << endl;
		}
	}

	if (pp.length()) {
		if (debug) cout << "  ffl: file:" << pp << endl;

		string dir = DirName(pp);
		if (debug)
			cout << "  ffl:  chdir to " << dir << endl;

		_chdir(dir.c_str());
	}

	cout << flush;
	cerr << "\n*" << cmd << "*\n";

	cmd += addargs;
	cerr << "\n**" << cmd << "**\n";
	if (debug)
		cout << "exec " << cmd << endl;

	int ret = system(cmd.c_str());
	return ret;
}

