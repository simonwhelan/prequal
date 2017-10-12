/*
 * Sequence.h
 *
 *  Created on: Oct 13, 2017
 *      Author: simon
 */

#ifndef SEQUENCE_H_
#define SEQUENCE_H_
#include <algorithm>
#include <sys/stat.h>
#include <unistd.h>
#include <unordered_map>
#include <vector>
#include <numeric>
#include <regex>
#include <cassert>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>


// Basic class for sequences
class CSequence {
public:
	CSequence() { };								// Blank constructor
	CSequence(std::string name, std::string seq);	// Standard constructor
//	std::string name;
//	std::string seq;
	std::vector <bool> Inside;				// Whether the character is on the inside or outside
	std::vector <bool> Remove;				// Whether to remove in the filter (true = remove)
	double PropInside;						// The proportion of the sequence labeled inside
	double PropRemoved;						// The proportion of hte sequence labeled to be removed
	void AddSequence(std::string seq);
	void AddName(std::string name);
	static void SetFilter(char filterOut) { _filterOut = filterOut; };
	int length() { return _seq.size(); }
	static int MaxLength() { return _maxLength; }
	std::string RealSeq(int pos = -1);				// Outputs the unfiltered seq
	std::string Seq(int pos = -1, bool filter = true, bool showOutside = false);		// Output the sequence (or pos i of sequence)
	std::string Name() { return _name; }
	bool Filter(int pos);					// Whether pos should be filtered/removed in any way

	std::string out() { return _name + " " + _seq; }
	void CalculateSummary() {
		int in = 0, rem = 0;
		for(int i = 0; i < length(); i++) {
			if(Inside[i]) { in++; }
			if(Remove[i]) { rem++; }
		}
		PropInside = (double) in / (double) length();
		PropRemoved = (double) rem / (double) length();
	}
private:
	static int _maxLength;		// Maximum length of the sequences examined
	std::string _name;			// The sequence
	std::string _seq;			// The name
	static char _filterOut;		// The string output on filtering

	void InitialiseFlags() {
		assert(Inside.empty() && Remove.empty());
		Inside.assign(_seq.size(),true);
		Remove.assign(_seq.size(),false);
	}

};

// File reader
std::vector <CSequence> *FASTAReader(std::string SeqFile);

// Other minor tools
template <class TRange> bool InRange(TRange Val, TRange LowerBound, TRange UpperBound) { return ( !(Val < LowerBound) && ( Val < UpperBound) ); }
#define my_min(a,b) ((a)<(b)?(a):(b))
#define my_max(a,b) ((a)>(b)?(a):(b))
std::string RemoveWhiteSpace(std::string s);
std::vector <std::string> Tokenise(std::string line);	// Tokenise a string
std::vector <std::string> Tokenise(std::string line, std::string Delim);		// Tokenise a string according to delimiter Delim
inline void ProgressSpinner(int suffix = -1) {
        static int count = 0;
        static char progress_spinner [] = "/-\\|";
        printf("\r%c",progress_spinner[count++]);
        if(suffix >= 0) { printf(" %d",suffix); }
        fflush(stdout);
        if(count == 4) { count = 0; }
};


#endif /* SEQUENCE_H_ */
