#include <iostream>
#include <string>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <stdlib.h>
#include "qgrs.h"
using namespace std;


string format_qgrs(G4 & g, char c) {
	string sequence = g.sequence;
	string t = "";
	for ( int i = 0; i < g.tetrads; i++) t += c;
	sequence.replace(0, g.tetrads, t);
	sequence.replace(g.tetrad2-g.start, g.tetrads, t);
	sequence.replace(g.tetrad3-g.start, g.tetrads, t);
	sequence.replace(g.tetrad4-g.start, g.tetrads, t);
	return sequence;
}

int pos_length = 5;
void print(ostream & out, string id, G4  & qgrs, bool raw, string sep, bool format, char c) {
	out << left;
	out << setw(raw ? 10 : 0) << id << sep;
	out << right;
	out << setw(raw ? pos_length : 0) << qgrs.start << sep;
	out << setw(raw ? pos_length : 0) << qgrs.tetrad2 << sep;
	out << setw(raw ? pos_length : 0) << qgrs.tetrad3 << sep;
	out << setw(raw ? pos_length : 0) << qgrs.tetrad4 << sep;

	out << setw(raw ? 3 : 0) << qgrs.tetrads << sep;
	out << setw(raw ? 4 : 0) << qgrs.gscore << sep;

	out << "  " << left;
	if (format) {
    	out << format_qgrs(qgrs, c) << sep;
	}
	else {
		out << qgrs.sequence << sep;
	}
	out << endl;
}
void print_heading(ostream & out, bool raw, string sep) {
	out << left;
	out << setw(raw ? 10 : 0) << "ID" << sep;
	out << right;
	out << setw(raw ? pos_length : 0) << "T1" << sep;
	out << setw(raw ? pos_length : 0) << "T2" << sep;
	out << setw(raw ? pos_length : 0) << "T3" << sep;
	out << setw(raw ? pos_length : 0) << "T4" << sep;

	out << setw(raw ? 3 : 0) << "TS" << sep;
	out << setw(raw ? 4 : 0) << "GS" << sep;

	out << "  " << left;
	out << "SEQ" << sep;
	out << endl;
	if ( raw ) {
		out << "--------------------------------------------------------------------------------------------" << endl;
	}
}

void help() {
	cout << "-----------------------------------------------------------------------------------------------------------" << endl;
	cout << " Command line options" << endl;
	cout << "-----------------------------------------------------------------------------------------------------------" << endl;
	cout << "-h                    help (this)" << endl;
	cout << "-csv                  output in csv mode" << endl;
	cout << "-json                 output in JSON mode" << endl;
	cout << "-i [input filename]   read input from a file as specified (defaults to stdin)" << endl;
	cout << "-o [output filename]  write output to file as specified (defaults to stdout)" << endl;
	cout << "-t [n]                filter output to QGRS with number tetrads equal to or greater than n (defaults to 2)" << endl;
	cout << "-s [n]                filter output to QGRS with g-score equal to or greater than n (defaults to 17)" << endl;
	cout << "-g [c]                replace all G characters within tetrads with given character." << endl;
	cout << "-v                    include overlapping QGRS in output (very large outputs may be generated)" << endl;
	cout << "-notitle              omit column titles in output (does not apply to JSON)" << endl;

	cout << "\n-----------------------------------------------------------------------------------------------------------" << endl;
	cout << " Output (default or csv)" << endl;;
	cout << "-----------------------------------------------------------------------------------------------------------" << endl;
	cout << "Column 1 - ID (order found in sequence).  x.y where x is the primary id, and y is number assigned overlaps.  " << endl;
	cout << "           For example, all QGRS listed as 2.y overlap QGRS listed with ID 2 - where 2 is the QGRS resulting" << endl;
	cout << "           in the highest G-Score in the group." << endl;
	cout << "Column 2 - Position of the start of the first tetrad (relative to beginning of input sequence)" << endl;
	cout << "Column 3 - Position of the start of the second tetrad (relative to beginning of input sequence)" << endl;
	cout << "Column 4 - Position of the start of the third tetrad (relative to beginning of input sequence)" << endl;
	cout << "Column 5 - Position of the start of the fourth tetrad (relative to beginning of input sequence)" << endl;
	cout << "Column 6 - Number of tetrads" << endl;
	cout << "Column 7 - G-Score" << endl;
	cout << "Column 8 - Sequence" << endl;
	cout << "\n-----------------------------------------------------------------------------------------------------------" << endl;
	cout << " Manual sequence entry)" << endl;;
	cout << "-----------------------------------------------------------------------------------------------------------" << endl;
	cout << "When you run this program without the -i option, the sequence must be entered via stdin.  Unless you are piping " << endl;
	cout << "in file contents, you will need to actually copy/paste or type the sequence data.  The program expects an EOF" << endl;
	cout << "signal at the end of the sequence - which will signal that it has received all input.  On the keyboard, type " << endl;
	cout << "Ctrl+D to enter the EOF signal!" << endl;
	cout << "-----------------------------------------------------------------------------------------------------------" << endl;
}


int main(int argc, char * argv[]) {
	string sequence = "";
	bool header = true;
	bool overlaps = false;
	int tetrads = 2;
	int gscore = 17;
	bool raw = true;
	bool json = false;
	string sep = "";
	ifstream in_fs;
	ofstream out_fs;
	bool infile = false;
	bool outfile = false;
	char sub;
	bool format = false;

	for ( int i = 1; i < argc; i++ ) {
		string arg = std::string(argv[i]);
		if ( arg == "-h" || arg == "-help") {
			help();
			return 0;
		}
		if ( arg == "-csv") {
			raw = false;
			sep = ",";
		}
		if ( arg == "-json") {
			json = true;
		}
		if ( arg == "-i") {
			infile = true;
			in_fs.open(argv[i+1]);
			if (in_fs.fail() ) {
				cout << "Failed to open " << argv[i+1] << endl;
				return 0;
			}
		}
		if ( arg == "-o") {
			outfile = true;
			out_fs.open(argv[i+1]);
			if (out_fs.fail() ) {
				cout << "Failed to open " << argv[i+1] << endl;
				return 0;
			}

		}
		if ( arg == "-t") {
			tetrads = atoi(argv[i+1]);
		}
		if ( arg == "-s") {
			gscore = atoi(argv[i+1]);
		}
		if ( arg == "-v") {
			overlaps = true;
		}
		if (arg == "-g") {
			sub = argv[i+1][0];
			format = true;
		}
		if (arg == "-notitle") {
			header = false;
		}
	}

	std::istream & in = (infile ? in_fs : std::cin);
	std::ostream & out = (outfile ? out_fs : std::cout);

	string buff;
	while (std::getline(in, buff))
	{
    	sequence += buff;
	}

	pos_length = log10(sequence.length()) + 2;

	if ( json) {
		out << findJSON(sequence, overlaps, tetrads, gscore);
	}
	else {
		out << endl;
		
		if ( header ) {
			print_heading(out, raw, sep);
		}
		
		vector<G4> results = find(sequence, overlaps, tetrads, gscore);
		int id = 1;
		if ( results.size() == 0 ) {
			cerr << "No QGRS found" << endl;
		}
		for (auto &qgrs : results) {
			print(out, std::to_string(id), qgrs, raw, sep, format, sub);
			int overlap_id = 1;
			for (auto &overlap : qgrs.overlaps) {
				string overlap_id_str = std::to_string(id) + "." + std::to_string(overlap_id++);
				print(cout, overlap_id_str, overlap, raw, sep, format, sub);

			}
			id++;
		}
	}
	if ( outfile ) out_fs.close();
}
