/*
 * SeqFilter.cpp
 *
 *  Created on: 12 Sep 2017
 *      Author: simon
 *
 *      ---
 *
 *      Main functions for performing filtering
 *
 *      //
 *
 *      Current idea list:
 *       * Run together removed regions, for example XXXARGNDEXXX -> XXXXXXXXXXXX so frameshifts are better caught
 *       * Beginning and end regions are too generous -- user should define run_inside
 *       * Need a full set of options and the beginning of documentation
 *       * Store bounds for particular sequences in the pairHMM so custom bounds are internalised and multiple calls are not needed.
 */


#include <algorithm>
#include <iomanip>
#include "SeqFilter.h"
#include "Sequence.h"

using namespace::std;

// Global variables. Ugly, but easiest quick solution
COptions *options = NULL;
vector <CSequence> *data = NULL;
double **PP = NULL;
vector<unordered_map<string,short>*>* kmers;
unsigned int kmerSize = 4;

int main(int argc, char * argv[]) {

	// Collect options
	options = new COptions(argc, argv);
	CSequence::SetFilter(options->CoreFilter());
	double threshold;
	vector<double> values;

	// Read data and sort initialisation
	data = FASTAReader(options->Infile()); // Reads the sequences
	Nseq = data->size(); // Allocates the number of sequences in utils.h
	sequence = new char*[data->size()]; // Allocates sequence array from utits.h
	lens = new int[data->size()]; // Gets the memory for lengths from utils.h
	PP = new double*[data->size()]; // Initialise Posterior probabilities
	for (int i = 0; i < data->size(); i++) {
		sequence[i] = new char[data->at(i).length()];
		PP[i] = new double[data->at(i).length()];
		lens[i] = data->at(i).length();
		for (int j = 0; j < data->at(i).length(); j++) {
			if (data->at(i).Seq(j).find("-") != std::string::npos) {
				cout << "\nERROR: sequence contains gaps and will break when pushed through a pairHMM...\n";
				exit(-1);
			}
			sequence[i][j] = pep2num(data->at(i).Seq(j)[0]); // Do the weird character to int to character conversion that Zorro demands...
			PP[i][j] = 0.0;
		}
	}
	cout << "\nThere are " << data->size() << " sequences of max length " << CSequence::MaxLength();
	//	for(int i = 0; i < data->size(); i++) { cout << "\ni=" << i << "\t" << data->at(i).out(); }

	// Run the HMM if needed
	RunHMM(options->Infile() + options->PPSuffix(), options->Overwrite_PP());

	// Define the threshold
	if (options->DoKeepProportion()) {
		cout << "\n\nExamining posterior probabilities to determine appropriate thresholds to retain " << options->KeepProportion() * 100 << "% of sequence" << flush;
		threshold = TargetCutoff(options->KeepProportion());
	} else {
		threshold = options->KeepThreshold();
		cout << "\n\nThreshold set to input value";
	}
	assert(InRange(threshold, 0.0, 1.0));

	// Do the filtering
	DoFiltering(threshold);

	////////////////////////////////////////////////////////////////
	cout << "\n\nOutputting results: ";
	if (options->DoDetail()) {
		cout << "\n\tDoing detailed output to " << options->Infile() << options->DetailSuffix() << flush;
		ofstream detail_out(options->Infile() + options->DetailSuffix());
		detail_out << std::fixed;
		detail_out << std::setprecision(4);
		detail_out << "# [seq_pos]seq_character\tmaxPP\tToRemove\tInside\n";
		for (int i = 0; i < data->size(); i++) {
			values.clear();
			double max = 0.0;
			for(int j = 0; j < data->at(i).length(); j++) {
				if(PP[i][j] > max) { max = PP[i][j]; }
				values.push_back(PP[i][j]);
			}
			detail_out << data->at(i).Name();
			detail_out << "\nmean= " << mean(values) << " : stdev= " << stdev(values) << " : poscut= " << mean(values) - (4*stdev(values)) << " : max= " << max;
			for(int j = 0; j < data->at(i).length(); j++) {
				detail_out << "\n["<< j<< "]" << data->at(i).RealSeq(j) << "\t" << PP[i][j] << "\t" << data->at(i).Remove[j] << "\t" << data->at(i).Inside[j];
			}
			detail_out << endl;
		}
		detail_out.close();
		cout << " ... done" << flush;
	}
	if(options->DoSummary()) {
		cout << "\n\tDoing summary output to " << options->Infile() << options->SummarySuffix() << flush;
		ofstream summary_out(options->Infile() + options->SummarySuffix());
		summary_out << std::fixed;
		summary_out << std::setprecision(4);
		// Calculate statistics
		double rem_mean = 0, rem_max = 0, in_mean = 0, in_min = 1.0;
		int rem_index = -1, in_index = -1;
		for(int i = 0; i < data->size(); i++) {
			data->at(i).CalculateSummary();
			// Removed
			rem_mean += data->at(i).PropRemoved;
			if(data->at(i).PropRemoved > rem_max) { rem_max = data->at(i).PropRemoved; rem_index = i; }
			// Inside
			in_mean += data->at(i).PropInside;
			if(data->at(i).PropInside < in_min) { in_min = data->at(i).PropInside; in_index = i; }
		}
		rem_mean /= (double) data->size();
		in_mean /= (double) data->size();
		// Output
		summary_out << "\nThere are " << data->size() << " sequences";
		summary_out << "\nRemoval:\n\tOn average " << rem_mean * 100 << "% of sequence retained";
		if(rem_index >= 0) { summary_out << "\n\tSequence with most removed (" << rem_max * 100 << "%) is [" << rem_index << "] = "<< data->at(rem_index).Name(); }
		summary_out << "\nCore regions:\n\tOne average " << in_mean * 100 << "% of sequence is in the core region";
		if(in_index >= 0) { summary_out << "\n\tSequence with least in core (" << in_min * 100 << "%) is [" << in_index << "] = "<< data->at(in_index).Name(); }
		summary_out << "\n##";
		for(int i = 0; i < data->size(); i++) {
			summary_out << "\n["<<i<<"] " << data->at(i).Name() << " has " << data->at(i).PropRemoved * 100 << "% removed and " << data->at(i).PropInside* 100 << "% in the core";
		}
		summary_out.close();
	}
	cout << "\n\tOutputting filtered sequences to " << options->Infile() << options->OutSuffix();
	int total_char = 0;
	int output_char = 0;
	ofstream sequence_out(options->Infile() + options->OutSuffix());
	for(int i = 0; i < data->size(); i++) {
		total_char += data->at(i).length();
		sequence_out << ">" << data->at(i).Name() << endl;
		// The whole sequence when filtered
		if(options->IgnoreSequence(data->at(i).Name())) {
			output_char += data->at(i).length();
			sequence_out << data->at(i).Seq();
			continue;
		}
		// The filtered style sequence determined by the CSequence class with some added counting stuff
		string output = data->at(i).Seq();
		for(int j = 0; j < output.size(); j++) { if(output[j] != options->CoreFilter()) { output_char ++; } }
		sequence_out << output << endl;
	}
	cout << "\nIn total " << output_char << " / " << total_char << "(" << 100* ((double)output_char/(double)total_char) << "%) of sequences retained";
	sequence_out.close();

	// Clean up memory
	for(int i = 0; i < data->size(); i++) { delete [] PP[i]; delete [] sequence[i];  } delete [] sequence; delete [] PP;
	delete data;

	cout << " ... program complete!\n\n";
}

#define DEBUG_PP 0
// Run the HMM stuff if required
void RunHMM(string outFile, bool forceOverwrite) {
	assert(data != NULL);
	// If file exists read the PPs
	if(file_exist(outFile) && !forceOverwrite) {
		cout << "\nPosterior probability file <" << outFile << "> exists. \n\tReading file";
		vector <string> Toks;
		string tmp;
		ifstream in(outFile);
		// Check sequence number
		tmp = read_line(in);
		Toks = Tokenise(tmp);
		if(Toks.size() != 1) { cout << "\nNumber of sequences on first line contains multiple tokens in PP file\n"; exit(-1); }
		if(atoi(Toks[0].c_str()) != data->size()) { cout << "\nNumber of sequences in first line of PP file does not match data\n"; exit(-1); }
		cout << " ... checking names" << flush;
		// Check sequence names
		for(int i = 0 ;i < data->size(); i++) {
			tmp = read_line(in);
			Toks = Tokenise(tmp);
			if(Toks.size() != 1) { cout << "\nName of sequence " << i << " has white spaces. This is not allowed\n"; exit(-1); }
			if(Toks[0] != data->at(i).Name()) { cout << "\nName of sequence " << i << " does not match that in data. This is not allowed\n"; exit(-1); }
		}
		cout << " ... checking PPs" << flush;
		// Input PPs
		double value;
		for(int i = 0 ;i < data->size(); i++) {
			tmp = read_line(in);
			Toks = Tokenise(tmp);
			if(Toks.size() != data->at(i).length()) { cout << "\nLength of sequence " << i << " does not match PPs. This is not allowed\n"; exit(-1); }
			for(int j = 0; j < data->at(i).length(); j ++) {
				value = (double) atof(Toks[j].c_str());
				if(value < 0.0 || value > 1.0) { cout << "\nThe PPs for sequence " << i << " contains a non probability (" << value<< ")"; exit(-1); }
				PP[i][j] = value;
			}
		}
		in.close();
		cout << " ... success" << flush;
	} else {
		// Calculate max length (repeat, but that's okay...)
		vector <int> lengthDist;
		for(int i = 0; i < data->size(); i++) {
			lengthDist.push_back(data->at(i).length());
		}
		cout << "\nPrepping pairHMM" << flush;
		initHMM(CSequence::MaxLength());
		///////////////////////////////// Get the PPs according to options ////////////////////////////////////"
		if(data->size() - options->PPnumber() < 5 || options->AllPP()) { // Get all the posteriors
			cout << "\nCollecting all posterior probabilities. This may take some time for large data sets:\n";
			SimonGetPosteriors(CSequence::MaxLength(), PP, options->DoApprox(), options->DefaultBound());
		} else if(options->ClosePP()) { // The closest
			cout << "\nCollecting subset of posterior probabilities based on closest " << options->PPnumber() << " sequences determined by Kmers\n\tThis may take time for larger data sets:\n" << flush;
			int rLsize = -1;
			MakeKmers();
			cout << "\nCreating collection sets\n";
			int **runList = new int*[data->size()];
			for(int i = 0; i < data->size(); i++) {
				vector <int> list = GetClosest(i, options->PPnumber());
				ProgressSpinner(i);
				runList[i] = new int[list.size()+1];
				if(rLsize < 0) { rLsize = list.size(); }
				else if(rLsize != list.size()) { cout << "\nProblem with different sized lists of running groups..."; }
				for(int j = 0; j < list.size(); j++) { runList[i][j] = list[j]; }
			}
			assert(rLsize > 0);
			cout << " ... done\nGetting posterior probabilities:\n";
			SimonGetSubsetPosteriors(CSequence::MaxLength(), PP, runList, rLsize, options->DoApprox(), options->DefaultBound());
			for(int i = 0; i < data->size(); i++) { delete [] runList[i]; } delete [] runList;
		} else if(options->LongPP()) { // The longest
			cout << "\nCollecting subset of posterior probabilities based on longest " << options->PPnumber() << " sequences:\n" << flush;
			int *coreList = NULL;
			vector <int> newL = ordered(lengthDist);
			coreList = new int[options->PPnumber()];
			for(int i = 0; i < options->PPnumber(); i++) {
				coreList[i] = newL[data->size() - 1 - i];
			}
			int **runList = new int*[data->size()];
			for(int i =0 ; i < data->size(); i++) {
				runList[i] = coreList;
			}
			SimonGetSubsetPosteriors(CSequence::MaxLength(), PP, runList, options->PPnumber(), options->DoApprox(), options->DefaultBound());
			for(int i =0 ; i < data->size(); i++) { runList[i] = NULL; }
			delete [] runList;
			delete [] coreList;
		} else {
			cout << "\nUnknown option for calculating posteriors...\n"; exit(-1);
		}
		// Output PPs
		if (options->DoPPs()) {
			cout << "\nOutputting posterior probabilities to " << outFile
					<< flush;
			ofstream out(outFile.c_str());
			out << data->size();
			for (int i = 0; i < data->size(); i++) {
				out << "\n" << data->at(i).Name();
			}
			for (int i = 0; i < data->size(); i++) {
				out << "\n";
				for (int j = 0; j < data->at(i).length(); j++) {
#if DEBUG_PP == 1
					out <<"[" << i << ","<< j<<"]";
#endif
					out << PP[i][j] << "\t";
				}
			}
			out << "\n";
			out.close();
			cout << "... done" << flush;
		}
	}
}

// Returns the cutoff based on the empirical set of PPs in PP[][]
double TargetCutoff(double prop2Keep) {
	cout << "\n\n<<<<<<<<<<<<<<<<< TODO: Make better cut-off for filtered sequences! <<<<<<<<<<<<<<<\n";
	cout << "\nDetermining empirical cut-off to retain " << prop2Keep * 100<< "% of sequence data";
	int total_length = 0;
	vector <double> tmp_PP;
	for(int i = 0; i < data->size(); i++)  {
		total_length += data->at(i).length();
		for(int j = 0; j < data->at(i).length(); j++)  {
			tmp_PP.push_back(PP[i][j]);
		}
	}
	std::sort(tmp_PP.begin(),tmp_PP.end());
	int count_stop;
	cout << "\nHelpful cut-offs ([PropRetained] Cutoffs):";
	cout << std::fixed;
	cout << std::setprecision(4);

	int spacer = 5;
	for(double x = 1.0; x >= 0.75; x-= 0.01,spacer ++) {
		if(spacer >= 5) { cout << "\n"; spacer = 0; }
		count_stop = (int)((1.0 - x) * (double) total_length);
		cout << "\t[" << x << "] " << tmp_PP[count_stop];
	}
	count_stop = (int)((1.0 - prop2Keep) * (double) total_length);
	return tmp_PP[count_stop];
}

void DoFiltering(double threshold) {
	cout << "\n\nPerforming filtering";
	cout << "\n\tApplying standard threshold " << threshold;
	int thresholdCount = 0;
	// Apply the threshold in a simple way
	for (int i = 0; i < data->size(); i++) {
		// 2. Find the sites to be filtered
		//		cout << " filtered sites..." << flush;
		for (int j = 0; j < data->at(i).length(); j++) {
			if (PP[i][j] < threshold) {
				thresholdCount++;
				data->at(i).Remove[j] = true;
			}
		}
	}
	cout << " resulting in " << thresholdCount << " sites removed" << flush;
	// Do the joining of filtered/outside regions if options require so
	if(options->FilterRange() > 0) {
		cout << "\n\tExtending filtered regions with width of " << options->FilterRange() << " ";
		int filter_count = 0;
		for(int i = 0 ; i < data->size(); i++) {
			int lastFilter = 0;
			for(int j = 0; j < data->at(i).length(); j++) {
				if(data->at(i).Filter(j)) {
					if(j - lastFilter < options->FilterRange() && j - lastFilter > 1) {
//						cout << "\nRemoving seq["<<i<<"]["<<j<<"] " << data->at(i).Name() << "== >"<< data->at(i).Seq(j) << "< : due to value " << j - lastFilter << " range (" << lastFilter << "," << j << ")";
//						cout << "\n\t" << data->at(i).RealSeq().substr(lastFilter,j-lastFilter) << " : length " << j - lastFilter;
						for(int k = j; k > lastFilter; k--) {
							data->at(i).Remove[k] = true;
						}
						filter_count++;
					}
					lastFilter = j;
				}
			}
		}
		cout << " ... " << filter_count << " additional regions removed" << flush;
	}
	// Tidy the front and back
//	if(false) {
	if(options->RunBeforeInside() > 0) {
		cout << "\n\tApplying front/back trimming for runs of " << options->RunBeforeInside();
		int seqTrimmed = 0;
		for (int i = 0; i < data->size(); i++) {
			// 1. Get the front ...
			bool DoOutside = false;
			for(int j = my_min(data->at(i).length(),options->RunBeforeInside())-1; j > 0; j--) {
				if(data->at(i).Filter(j)) { if(!DoOutside) { seqTrimmed++;} DoOutside = true;}
				if(DoOutside) { data->at(i).Inside[j] = false; data->at(i).Remove[j] = true; }
			}
			for(int j = 0; j < data->at(i).length(); j++) {
				if(!data->at(i).Filter(j)) { break; }
				data->at(i).Inside[j] = false;
			}
			// 2. And the back ...
			DoOutside = false;
			for(int j = my_max(0,data->at(i).length() - options->RunBeforeInside()); j < data->at(i).length(); j++)  {
				if(data->at(i).Filter(j)) { if(!DoOutside) { seqTrimmed++;} DoOutside = true; }
				if(DoOutside) { data->at(i).Inside[j] = false; data->at(i).Remove[j] = true; }
			}
			for(int j = data->at(i).length() - 1; j >=0; j--) {
				if(!data->at(i).Filter(j)) { break; }
				data->at(i).Inside[j] = false;
			}
		}
		cout << " resulting in " <<seqTrimmed << " sections removed" << flush;
	}
	cout << "\n\t... done" << flush;
}

//////////////////////////////////////// Kmer stuff from PaHMM tree //////////////////////////////////////////

void MakeKmers() {
	assert(data != NULL);
	// Get the memory
	kmers = new vector<unordered_map<string,short>*>(data->size());

	for(int i = 0; i < data->size(); i++) {
		(*kmers)[i] = new unordered_map<string,short>();
		string seq = data->at(i).Seq();
		extractKmers(seq, (*kmers)[i]);
	}
}
void extractKmers(string& seq, unordered_map<string, short>* umap)
{
	string kmer;

	if(seq.size() < kmerSize)
		return;

	for(unsigned int i = 0; i< (seq.size() - kmerSize); i++)
	{
		kmer = seq.substr(i, kmerSize);
		++((*umap)[kmer]);
	}
}

vector <int> GetClosest(int Sequence, int NumberClosest) {
	vector <double> dist = MakeKmerDistanceRow(Sequence);
	vector <int> indices = ordered(dist);
	for(int i = 0; i < indices.size(); i++) {
		if(indices[i] == Sequence) {
			indices.erase(indices.begin() + i);
			break;
		}
	}
	if(indices.size() > NumberClosest) {
		indices.erase(indices.begin() + NumberClosest, indices.end());
	}
	return indices;
}

vector <double> MakeKmerDistanceRow(int CurSeq) {
	vector <double> RetDist;
	for(int i = 0; i < data->size(); i++) {
		if(i == CurSeq) { RetDist.push_back(1.0); continue; }	// For the current sequence return error
		double identity = 1.0 - commonKmerCount(i,CurSeq)/(double)my_min(data->at(i).length(),data->at(CurSeq).length());
		identity = (pow(100, identity-1.04)+0.01)/0.6;	// Marcin's AA function
		RetDist.push_back(identity);
	}
	return RetDist;
}

unsigned int commonKmerCount(unsigned int i, unsigned int j)
{
	unordered_map<string, short>* m1 = (*kmers)[i];
	unordered_map<string, short>* m2 = (*kmers)[j];
	unsigned int commonCount = 0;

	for(auto it = m1->begin(); it != m1->end(); it++)
	{
		commonCount += std::min((short)((*m2)[it->first]), (short)(it->second));
	}
	//DEBUG("Common k-mer count between seq. " << i << " and " << j << " is " << commonCount);
	return commonCount;


}


double mean(vector <double> vec) {
	double ret = 0.0;
	for(int i = 0; i < vec.size(); i++) {
		ret += vec[i];
	}
	return ret / (double)vec.size();
}

double stdev(vector <double> vec) {
	double ret = 0, ave = mean(vec);
	for(int i = 0; i < vec.size(); i++) {
		ret += (vec[i] - ave) * (vec[i] - ave);
	}
	ret /= (double) vec.size();
	return sqrt(ret);
}
