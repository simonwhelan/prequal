/*
 * Sequence.cpp
 *
 *  Created on: Oct 13, 2017
 *      Author: simon
 */

#include "Sequence.h"

int CSequence::_maxLength = 0;
char CSequence::_filterOut = 'X';

using namespace::std;

//////////////// CSequence
CSequence::CSequence(std::string name, std::string seq) {
	AddName(name);
	AddSequence(seq);
	InitialiseFlags();
}
void CSequence::AddName(std::string name) {
	if(!_name.empty()) {
		std::cout << "\nCSequence ERROR: Trying to added name to non-empty named sequence\n";
		exit(-1);
	}
	_name = name;
}
void CSequence::AddSequence(std::string seq) {
	if (!_seq.empty()) {
		std::cout << "\nCSequence ERROR: Trying to added sequence to non-empty named sequence\n";
		exit(-1);
	}
	_seq = seq;
	if(_seq.size() > _maxLength) { _maxLength = _seq.size(); }
}
std::string CSequence::Seq(int pos, bool filter, bool showOutside) {
	std::stringstream ss;
	if(pos != -1) {
		if(!showOutside && !Inside[pos]) { ss << "0"; }
		else if(filter && Remove[pos]) { if(_filterOut != '\0') { ss <<_filterOut; } }
		else { ss << _seq[pos]; }

	} else {
		for( int i = 0 ; i < _seq.size(); i++) {
			if(!showOutside && !Inside[i]) { continue; }
			if(filter && Remove[i]) { if(_filterOut != '\0') { ss <<_filterOut; } }
			else { ss << _seq[i]; }
		}
	}
	return ss.str();
}
std::string CSequence::RealSeq(int pos) {
	if(pos != -1) {
		return _seq.substr(pos,1);
	}
	return _seq;
}
bool CSequence::Filter(int pos) {
	if(Remove[pos] || !Inside[pos]) { return true; }
	return false;
}


/////////////// Minor functions
// File reader
std::vector <CSequence> *FASTAReader(std::string SeqFile) {
	std::vector <CSequence> *RetSeq = new std::vector<CSequence>();

	std::ifstream input(SeqFile.c_str(), std::ifstream::in);
    if(!input.good()){
        std::cerr << "Error opening '"<< SeqFile <<"'. Provide a valid file." << std::endl;
        exit(-1);
    }

    std::string line, name, content;

    while(getline( input, line ) ){
    	line = RemoveWhiteSpace(line);
    	if(line.empty()) { continue; }		// Skip blank lines
        if(line[0] == '>' ){ // Identifier marker
            if( !line.empty() ){
            	if(!name.empty()) { // Push the sequence back
            		RetSeq->push_back(CSequence(RemoveWhiteSpace(name),RemoveWhiteSpace(content)));
            		name.clear(); content.clear();
            	}
                name = line.substr(1);
            }
            content.clear();
        } else if( !name.empty() ){
            content += RemoveWhiteSpace(line);
        }
    }
    // Add the final sequence
    RetSeq->push_back(CSequence(RemoveWhiteSpace(name),RemoveWhiteSpace(content)));
    return RetSeq;
}

std::string RemoveWhiteSpace(std::string s) {
	s.erase( std::remove_if( s.begin(), s.end(), ::isspace ), s.end() );
	return s;
}

std::vector <std::string> Tokenise(std::string line) {
	std::string buf;
	std::stringstream in(line);
	std::vector <std::string> Toks;
	Toks.~vector();
	while(in >> buf) { Toks.push_back(buf); }
	return Toks;
}
std::vector <std::string> Tokenise(std::string line, std::string Delim)	{
	size_t i = 0, j,j1;
	std::vector <std::string> Toks;
	while(i != (int)line.size())	{
		j = line.find(Delim,i+1);
		if(j == std::string::npos) { j = j1 = (int)line.size(); } else { j1 = j+1; }
		Toks.push_back(line.substr(i,j-i));
		i = j1;
	}
	return Toks;
}

