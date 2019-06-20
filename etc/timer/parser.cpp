// PARSER - v0.3.1

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <filesystem>

using namespace std;
namespace fs = filesystem;
#define DEBUG
#define IND (string)"\t"

#define functionDefinition "__FUNCTION_DEFINITION__ "
#define typeDefinition "__TYPE_DEFINITION__ "
#define timeElapsed "__TIME_ELAPSED__ "
#define parameterDefinition "__PARAMETER_DEFINITION__ "

void findFlagAndExtract (const string &line, vector<string> *out);
bool processLog(const string &filePathm, vector<string> *output);
void addLogToHtml(vector<string> *input, ofstream &htmlFile);
void addTimeDiv (ofstream &file, vector<string>::iterator it, int indent);
string getFlag (const string &line);
string getValue(const string &line);
string getClass (const string &line);
string getInd (int n);

/*!
 * \brief Main function
 * \param argc Number of arguments
 * \param argv Argument(s) value(s)
 * \return Success / Failure
 */
int main (int argc, char *argv[]) {
  // Parameters
  string unitDirectory = "../FreeFem-sources/unit/";
  string flagValOutputFile = "flagVal";
  string htmlOutputFile = "logs.html";

  if (argc == 4) {
    unitDirectory = argv[1];
    flagValOutputFile = argv[2];
    htmlOutputFile = argv[3];
  }
  else if (argc != 1) {
    cout << "Correct use : ./parser.exe [unitdirectory] [flagvaloutputfile] [htmlfile]" << endl;
    return EXIT_FAILURE;
  }

  // Create html file
  ofstream flagValFile(flagValOutputFile);
  ofstream htmlFile(htmlOutputFile);

  // Retrieve and process all .log files
  for(auto &logdir : fs::directory_iterator(unitDirectory)) {
    if (!logdir.is_directory())
      continue;

    for(auto &log : fs::directory_iterator(logdir.path())) {
      if (log.path().extension() != ".log")
        continue;

      #ifdef DEBUG
      	cout << endl << log.path().filename().string() << endl;
      #endif

      // find the flags and their associated values
      vector<string> output;
      if (!processLog(log.path().string(), &output))
        return EXIT_FAILURE;

      // add flags and their associated values to flagVal file
      for (vector<string>::iterator it = output.begin(); it != output.end(); ++it)
    		flagValFile << *it << endl;

      // add formated flag and value to html file
      addLogToHtml(&output, htmlFile);
    }
  }

  flagValFile.close();
  htmlFile.close();
	return EXIT_SUCCESS;
}

/*!
 * \brief Find and extract a flag in a line and add it to output vector
 * \param line Line
 * \param type Type
 * \param out Out
 */
void findFlagAndExtract (const string &line, vector<string> *output) {
	if (!line.length())
		return;
	size_t pos = line.find(functionDefinition);
  if (pos == string::npos) pos = line.find(typeDefinition);
  if (pos == string::npos) pos = line.find(timeElapsed);
  if (pos == string::npos) pos = line.find(parameterDefinition);
	if (line[0] == '_' && pos != string::npos)
		output->push_back(line);
}

/*!
 * \brief Process log file
 * \param filePath Log file path
 * \param output String vector containing flag anv their associated values
 * \return Successfully processed or not
 */
bool processLog(const string &filePath, vector<string> *output) {
  vector<string> fileContent;
	ifstream file(filePath, ios::in); // open file, read only mode
	if (file) {
		string temp;
		while (getline(file, temp))	// read all lines, keep ones containing flags
      findFlagAndExtract(temp, output);
	} else {	// return error
		cerr << "Unable to open " << filePath << endl;
    return false;
	}
  file.close();

  #ifdef DEBUG
    cout << output->size() << (output->size() < 2 ? " flag" : " flags") << " found" << endl;
  	// Display vector content
  	for (vector<string>::iterator it = output->begin(); it != output->end(); ++it)
  		cout << *it << endl;
  #endif

  return true;
}

/*!
  * \brief Add formated log content to html file
  * \param input Log formatted content
  * \param htmlFile HTML output file
  */
void addLogToHtml(vector<string> *input, ofstream &htmlFile) {
  if (!htmlFile || !input->size())
    return;

  vector<string>::iterator it = input->begin();

  if (input->size() == 1) {
    htmlFile << "<div class=\"" << getClass(*it) << "\">" << endl;
    htmlFile << getInd(1) << "<p>" << getValue(*it) << "</p>" << endl;
    htmlFile << "</div>" << endl;
    return;
  }

  while((it != input->end()) && (getFlag(*it) == functionDefinition)) {
    htmlFile << "<div class=\"" << getClass(*it) << "\">" << endl;
    htmlFile << getInd(1) << "<p>" << getValue(*it) << "</p>" << endl;
    it++;

    flagsleft: {
      if ((getFlag(*it) != typeDefinition) && (getFlag(*it) != parameterDefinition))
        goto timeonly;

      while ((it != input->end()) && ((getFlag(*it) == typeDefinition) || (getFlag(*it) == parameterDefinition))) {
        htmlFile << getInd(1) << "<div class=\"" << getClass(*it) << "\">" << endl;
        htmlFile << getInd(2) << "<p>" << getValue(*it) << "</p>" << endl;
        it++;
        while ((it < input->end()) && (getFlag(*it) == timeElapsed)) {
          addTimeDiv(htmlFile, it, 2);
          it++;
        }
        htmlFile << getInd(1) << "</div>" << endl;
      }

      timeonly: while ((it != input->end()) && (getFlag(*it) == timeElapsed)) {
        addTimeDiv(htmlFile, it, 1);
        it++;
      }
    }

    if ((it != input->end()) && (getFlag(*it) != functionDefinition))
      goto flagsleft;

    htmlFile << "</div>" << endl;
  }
}

/*!
  * \brief Add div element containing time elapsed to html file
  * \param it String iterator
  * \param indent Div indentation
  */
void addTimeDiv (ofstream &htmlFile, vector<string>::iterator it, int indent) {
  htmlFile << getInd(indent) << "<div class=\"" << getClass(*it) << "\">" << endl;
  htmlFile << getInd(indent + 1) << "<p>" << getValue(*it) << "</p>" << endl;
  htmlFile << getInd(indent) <<"</div>" << endl;
}

/*!
  * \brief Extract flag on a line
  * \param line Line to process
  * \return Flag on line
  */
string getFlag (const string &line) {
  return line.substr(0, line.find(" ") + 1);
}

/*!
  * \brief Extract value after a flag
  * \param line Line containing a flag
  * \return Value on line
  */
string getValue (const string &line) {
  return line.substr(line.find(" ") + 1);
}

/*!
  * \brief Returns class for html format
  * \param line Line containing the flag
  * \return Html class
  */
string getClass (const string &line) {
  string flag = getFlag(line);
  if (flag == timeElapsed)
    return "time";
  if (flag == typeDefinition)
    return "type";
  if (flag == parameterDefinition)
    return "parameter";
  if (flag == functionDefinition)
    return "function";
  return "";
}

/*!
  * \brief Returns a specific number of indentations
  * \param n Number of indentations
  * \return Indentations
  */
string getInd (int n) {
  string out = "";
  for (int i = 0; i < n; i++)
    out += IND;
  return out;
}
