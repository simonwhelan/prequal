/*
 * Options.cpp
 * **
 * Sets the options from the command line
 *
 *  Created on: Oct 4, 2017
 *      Author: simon
 */

#include "SeqFilter.h"
#include "Sequence.h"

using namespace::std;

COptions::COptions(int argc, char *argv[]) {

	// Check basic input is matched
	if (argc < 2) {
		cout << "Incorrect command line. Usage: \n\t" << PROGRAM_NAME << " [options] input_file\n\t-h [all] for [full] options\n";
		exit(-1);
	}
	// Get the input file
	if (strcmp(argv[1], "-h") != 0) {
		if (argv[argc - 1][0] == '-') {
			cout << "\nInput file must be specified last. you tried input of <" << argv[argc - 1] << ">\n";
			exit(-1);
		}
		_inputFile = argv[argc-1];
		argc--;
	}
	for (int i = 1; i < argc; i++) {
		if (strcmp(argv[i], "-h") == 0) {
			cout << "----------------------------------------------------\n\t"
					<< PROGRAM_NAME << " v." << VERSION_NUMBER
					<< "\n\t\tWritten by Simon Whelan"
					<< "\n----------------------------------------------------";
			if(i + 1 < argc) {
				if(strcmp(argv[i+1],"all") == 0) {
					cout << "\n\nOptions affecting the core region and filtering:";
					cout << "\n\t-corerun X       \t: X number of high posterior sites at beginning and end before \n\t\t\t\t\ta core region is defined [DEFAULT "
						 << _runBeforeInside << "]";
					cout << "\n\t-nocore          \t: No core region will be defined [DEFAULT uses core region]";
					cout << "\n\t-removeall       \t: Remove all sites rather than those outside the core region \n\t\t\t\t\t[DEFAULT sites in core filtered with X]";
					cout << "\n\t-corefilter X    \t: The character outputted in the core for the  filtered alignment [DEFAULT is X]";
					cout << "\n\nOptions affecting output formats:";
					cout << "\n\t-outsuffix X     \t: Output file will be the original name with X as a suffix [DEFAULT .filtered]";
					cout << "\n\t-dosummary X     \t: Output summary statistics to file suffixed with X [if no X then DEFAULT .summary]";
					cout << "\n\t-dodetail X      \t: Output detailed statistics to file suffixed with X [if no X then DEFAULT .detail]";
					cout << "\n\t-noPP            \t: Stop outputting the posterior probability matrix";
					cout << "\n\nOptions affecting posterior probabilities and filtering:";
					cout << "\n\t-pptype X [Y]       \t: Specify the algorithm used to calculate posterior probabilities\n\t\t\t\t\tX = all : for all against all sequence comparisons";
					cout << "\n\t\t\t\t\tX = closest : for Y closest relatives [DEFAULT; Y = 10]\n\t\t\t\t\tX = longest : for comparing the Y longest sequences [Y = 10]";
					cout << "\n\t-filterprop X       \t: Filter the sequences so that in total X proportion (range 0.0 - 1.0) \n\t\t\t\t\tof the sequences are maintained. (DEFAULT: "
						 << _keepProportion << ")";
					cout << "\n\t-filterthresh X     \t: Filter the sequences to the posterior probabilities threshold X \n\t\t\t\t\t(range 0.0 - 1.0). Any amount of sequence can be removed. Not recommended.";
					cout << "\n\t-filterjoin X       \t: Extend filtering over regions of unfiltered sequence less than X [DEFAULT X = " << FilterRange() << "]";
					cout << "\n\t-nofilterlist X     \t: Specify a file X that contains a list of taxa names that will \n\t\t\t\t\tnot be filtered. In X one name per line.";
					cout << "\n\t-nofilterword X     \t: Specify a file X that contains a list of words and sequence names that contain \n\t\t\t\t\tthose words will not be filtered. In X one word per line.";
					cout << "\n\nUsage:\n\t" << PROGRAM_NAME << " [options] input_file\n\n";
					exit(-1);
				}
			}
			cout << "\n\nSimple options (-h all for full options)";
			cout << "\n\t-corerun X       \t: X number of high posterior sites at beginning and end before \n\t\t\t\t\ta core region is defined [DEFAULT "
				 << _runBeforeInside << "]";
			cout << "\n\t-pptype X [Y]       \t: Specify the algorithm used to calculate posterior probabilities\n\t\t\t\t\tX = all : for all against all sequence comparisons";
			cout << "\n\t-filterjoin X       \t: Extend filtering over regions of unfiltered sequence less than X [DEFAULT X = " << FilterRange() << "]";
			cout << "\n\t-nofilterlist X     \t: Specify a file X that contains a list of taxa names that will \n\t\t\t\t\tnot be filtered. In X one name per line.";
			cout << "\n\nUsage:\n\t" << PROGRAM_NAME << " [options] input_file\n\n";
			exit(-1);
		}
		/////////// Core stuff
		else if (strcmp(argv[i], "-corerun") == 0) {
			if (i + 1 == argc) {
				cout << "\nError when definining -corerun. Specify number afterwards";
				exit(-1);
			}
			_runBeforeInside = atoi(argv[++i]);
			if (!InRange(_runBeforeInside, 0, 100)) {
				cout << "\n-corerun " << _runBeforeInside
						<< " outside good limits (0,100)";
			}
		} else if (strcmp(argv[i], "-nocore") == 0) {
			_remove2Core = false;
		} else if (strcmp(argv[i], "-removeall") == 0) {
			_removeAll = true;
		} else if (strcmp(argv[i], "-corefilter") == 0) {
			if (i + 1 == argc) {
				cout << "\nError when definining -corefilter. Specify a character afterwards";
				exit(-1);
			}
			_coreFilter = argv[++i][0];
		}
		//////////// File stuff
		else if (strcmp(argv[i], "-outsuffix") == 0) {
			if (i + 1 == argc) {
				cout << "\nError when definining -outsuffix. Specify a suffix afterwards";
				exit(-1);
			}
			_outputSuffix = argv[++i];
			if (_outputSuffix.size() < 1) {
				cout << "\nMust specify a meaningful outsuffix";
			}
		} else if (strcmp(argv[i], "-dosummary") == 0) {
			_doSummary = true;
			if (argv[i + 1][0] != '-') {
				continue;
			}
			_summarySuffix = argv[++i];
			if (_summarySuffix.size() < 1) {
				cout << "\nMust specify a meaningful summary suffix";
			}
		} else if (strcmp(argv[i], "-dodetail") == 0) {
			_doDetail = true;
			if (argv[i + 1][0] != '-') {
				continue;
			}
			_detailSuffix = argv[++i];
			if (_detailSuffix.size() < 1) {
				cout << "\nMust specify a meaningful detail suffix";
			}
		} else if (strcmp(argv[i], "-noPP") == 0) {
			_doPPs = false;
		}
		///////////// Filtering and PP stuff
		else if (strcmp(argv[i], "-pptype") == 0) {
			if (i + 1 == argc ) {
				cout << "\nError when definining -pptype. Specify type afterwards (all; closest; longest)";
				exit(-1);
			}
			string pp_option = argv[++i];
			if(pp_option == "all") {
				_ppCalcs = -1;
			} else if(pp_option == "closest") {
				_ppCalcs = 0;
			} else if(pp_option == "longest") {
				_ppCalcs = 1;
			} else {
				cout << "\nUnknown option when defining -pptype : " << pp_option; exit(-1);
			}
			if(i + 1 < argc && _ppCalcs != -1) {
				if(argv[i+1][0] != '-') {
					_ppCalcNumber = atoi(argv[++i]);
					if(!InRange(_ppCalcNumber,1,1000)) {
						cout << "\n-pptype X Y : Y is not in expected range: " << _ppCalcNumber << endl; exit(-1);
					}
				}
			}
		} else if (strcmp(argv[i], "-nofilterlist") == 0) {
			if (i + 1 == argc) {
				cout << "\nError when definining -nofilterlist. Specify a file afterwards";
				exit(-1);
			}
			_noFilterList = GetStringList(argv[++i]);
			cout << "\nFilter list consists of " << _noFilterList.size()
					<< " entries";
		} else if (strcmp(argv[i], "-nofilterregex") == 0) {
			if (i + 1 == argc) {
				cout << "\nError when definining -nofilterregex. Specify a file afterwards\n";
				exit(-1);
			}
			_noFilterWord = GetStringList(argv[++i]);
			cout << "\nFilter word consists of " << _noFilterWord.size()
					<< " entries";
		} else if (strcmp(argv[i], "-filterthresh") == 0) {
			if (i + 1 < argc) {
				_filterThreshold = atof(argv[++i]);
			} else {
				_filterThreshold = -1;
			}
			if (!InRange(_filterThreshold, 0.0, 1.0)) {
				cout << "\nError when defining -filterthresh. Specify a valid number (0,1) afterwards\n";
				exit(-1);
			}
		} else if (strcmp(argv[i], "-filterprop") == 0) {
			if (_filterThreshold > 0) {
				cout << "\nCannot set -filterprop and -filterthresh at the same time.\n";
				exit(-1);
			}
			if (i + 1 < argc) {
				_keepProportion = atof(argv[++i]);
			} else {
				_keepProportion = -1;
			}
			if (!InRange(_keepProportion, 0.0, 1.0)) {
				cout << "\nError when defining -filterprop. Specify a valid number (0,1) afterwards\n";
				exit(-1);
			}
		} else if (strcmp(argv[i], "-filterjoin") == 0) {
			if (i + 1 < argc) {
				_joinFilterRange = atof(argv[++i]);
			} else {
				_joinFilterRange = -1;
			}
			if(!InRange(_joinFilterRange,0,1000)) {
				cout << "\nNeed to define a value after -filterjoin of 0 or greater \n";exit(-1);
			}
			if (!InRange(_keepProportion, 0.0, 1.0)) {
				cout << "\nError when defining -filterprop. Specify a valid number (0,1) afterwards\n";
				exit(-1);
			}
		}
		///////////// Default option
		else {
			cout << "Incorrect command line. Usage: \n\t" << PROGRAM_NAME
					<< " [options] input_file\n\t-h for full options\n";
			exit(-1);
		}

	}
}

vector<string> COptions::GetStringList(string FilterListFile) {
	ifstream input(FilterListFile.c_str());
	if (!input.good()) {
		std::cerr << "Error opening '" << FilterListFile
				<< "'. Provide a valid file." << std::endl;
		exit(-1);
	}
	string line, name;
	vector < string > FilterList;
	while (getline(input, line)) {
		if (line.size() < 1) {
			continue;
		}
		line = RemoveWhiteSpace(line);
		FilterList.push_back(line);
	}
	return FilterList;
}
