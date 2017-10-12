/*
 * SeqFilter.h
 *
 *  Created on: 12 Sep 2017
 *      Author: simon
 */

#ifndef SEQFILTER_H_
#define SEQFILTER_H_

#include <algorithm>
#include <sys/stat.h>
#include <unistd.h>
#include <unordered_map>
#include <iostream>
#include <vector>
#include <numeric>
#include <regex>

extern "C" {
#define EXTERN extern
#include "utils.h"
#include "hmm.h"
#undef EXTERN
}

/////////////// General naming stuff
const std::string PROGRAM_NAME = "seqfilter";
const std::string VERSION_NUMBER = "0.02a";

// Class holding the options for the program. Implemented in Options.cpp
class COptions {
public:
	COptions(int argc, char * argv[]);			// Creates the options from reading the command line. Only done once

	// Access functions
	// 1. File handling
	std::string Infile() { return _inputFile; }
	std::string OutSuffix() { return _outputSuffix; }
	bool DoSummary() { return _doSummary; }
	std::string SummarySuffix() { return _summarySuffix; }
	bool DoDetail() { return _doDetail; }
	std::string DetailSuffix() { return _detailSuffix; }
	bool DoPPs() { return _doPPs; }
	std::string PPSuffix() { return _ppSuffix; }
	bool Overwrite_PP() { return _overwritePP; }
	// 2. Core regions
	int RunBeforeInside() { return _runBeforeInside; }
	bool Remove2Core() { return _remove2Core; }
	bool RemoveAll() { return _removeAll; }
	char CoreFilter() { return _coreFilter; }
	// 3. Filtering options
	bool IgnoreSequence(std::string seq_name) {
		if(_noFilterList.empty() && _noFilterWord.empty()) { return false; }
		if(find(_noFilterList.begin(), _noFilterList.end(),seq_name) != _noFilterList.end()) { std::cout << "\nSequence filtered: " << seq_name; return true; }
		for(int i = 0; i < _noFilterWord.size(); i++) {
			if(seq_name.find(_noFilterWord[i]) != std::string::npos) { std::cout << "\nWord filtered " << seq_name; return true; }
		}
		return false;
	}
	bool DoKeepProportion() { if(_filterThreshold < 0) { return true; } else { return false; } }
	double KeepProportion() { return _keepProportion; }
	double KeepThreshold() { return _filterThreshold; }
	int FilterRange() { return _joinFilterRange; }
	// 4. PP calculation options
	int PPcalcNumber() { return _ppCalcs; }
	bool AllPP() { if(_ppCalcs == -1) { return true; } return false; }
	bool LongPP() { if(_ppCalcs == 1) { return true; } return false; }
	bool ClosePP() { if(_ppCalcs == 0) { return true; } return false; }
	int PPnumber() { return _ppCalcNumber; }
	// 5. Approximations for pairHMM calculations
	bool DoApprox() { return _approxPP; }
	int DefaultBound() { return 25; }

private:
	// Behaviour defining core and how the different parts are treated/filtered
	int _runBeforeInside = 3;			// Number of positively defined homologies before a sequence is consider part of the core
	bool _remove2Core = true;			// Whether parts outside the core are removed
	bool _removeAll = false;			// Whether all fails are removed
	char _coreFilter = 'X';			// The character used for filtering in the core
	// File handling
	std::string _inputFile;
	std::string _outputSuffix = ".filtered";
	bool _doSummary = true;				// General summary
	std::string _summarySuffix = ".summary";
	bool _doDetail = true;				// Detailed output per site
	std::string _detailSuffix = ".detail";
	bool _doPPs = true;					// Grid of PPs for fast downstream computation
	std::string _ppSuffix = ".PP";
	// Filtering options
	std::vector <std::string> _noFilterList;		// List of taxa names that will not be filtered at all
	std::vector <std::string> _noFilterWord;		// List of taxa regular expressions that will not be filtered at all
	double _keepProportion = 0.85;		// The proportion of sites to be kept
	double _filterThreshold = -1.0;		// Hard set the filter threshold
	int _joinFilterRange = 25;			// The maximum gap between filtered characters before the whole segment is filtered
	// Posterior probability calculation options
	int _ppCalcs = 0;					// Options are 0: closest _ppCalcNumber; 1: longest _ppCalcNumber; -1: all
	int _ppCalcNumber = 10;
	// Options that are currently disabled
	bool _overwritePP = false;
	bool _approxPP = true;				// Whether to use full or approximate PP computations according to my hack
	// Internal functions
	std::vector <std::string> GetStringList(std::string FilterListFile);
};

double mean(std::vector <double> vec);
double stdev(std::vector <double> vec);
void RunHMM(std::string outFile, bool forceOverwrite = false);		// C++ interaction function with Zorro
double TargetCutoff(double prop2Keep);							// Computes the target cutoff given a specific proportion to keep
void DoFiltering(double Threshold);								// Applies the filtering method

// Kmer stuff
void MakeKmers();
void extractKmers(std::string& seq, std::unordered_map<std::string,short>* umap);
std::vector <double> MakeKmerDistanceRow(int CurSeq);
unsigned int commonKmerCount(unsigned int i, unsigned int j);
std::vector <int> GetClosest(int Sequence, int NumberClosest);

inline bool file_exist (const std::string& name) {
  struct stat buffer;
  return (stat (name.c_str(), &buffer) == 0);
}

inline std::string read_line(std::istream &in) {
	std::string tmp;
	getline(in,tmp);
	if(!in.good()) { std::cout << "\nError reading file..."; exit(-1); }
	return tmp;
}

template <typename T>
std::vector<int> ordered(std::vector<T> const& values) {
    std::vector<int> indices(values.size());
    std::iota(begin(indices), end(indices), static_cast<int>(0));

    std::sort(
        begin(indices), end(indices),
        [&](int a, int b) { return values[a] < values[b]; }
    );
    return indices;
}


#endif /* SEQFILTER_H_ */
