#include <map>
#include <vector>
#include <string>
#include <iomanip>
#include <fstream>
#include <iostream>

#include "TString.h"
#include "TSystem.h"
#include "TStopwatch.h"
#include "TDatime.h"

#include "AnaMaker/AnaMaker.h"

using std::map; 
using std::cout; 
using std::endl; 
using std::string; 
using std::vector; 
using std::ifstream; 

int main(int argc, char** argv) {
	/*
		Arguments: 2
		:file list: str, path to the file list
		:output file name: str, root file name of output, subfix ".root" is not included
	*/

	const char* lg = "[LOG] - Basic QA: ";

	cout << lg << "Number of arguments received: " << argc << endl;
	if (argc < 2) {
		cout << lg << "EC201: Number of arguments doesn't match!" << endl;
		return 101;
	}
	string fileList = argv[1];
	string outputName = string(argv[2]) + ".root";
	cout << lg << "File list: " << fileList << endl;
	cout << lg << "Output name: " << outputName << endl;

	ifstream list(fileList);
	if (!list.is_open()) {
    	cout << lg << "EC202: File list cannot open!" << fileList << endl;
    	return 102;
	}
	vector<string> files;
	string fName;
	while (getline(list, fName)) {
		if (fName.empty()) { continue; }
		files.push_back(fName);
	}
	list.close();

	TStopwatch timer;
  	timer.Start();
	const int nFiles = files.size();
	if (nFiles == 0) {
		cout << lg << "EC203: No file in the list!" << endl;
		return 103;
	}
	cout << lg << "Number of input files: " << nFiles << endl;
	
	cout << lg << "Initializing analysis maker" << endl;
	AnaMaker ana(files, outputName);
	ana.Init();
	const Int_t nEv = ana.GetAnaEntries();
    cout << lg << "Number of events: " << nEv << endl;
	for (Int_t iEv=0; iEv<nEv; iEv++) {
		ana.Make(iEv);
	}
	ana.Finish();
	timer.Stop();
	Double_t elapsedTime = timer.RealTime();
	int hours = elapsedTime / 3600;
	int minutes = (elapsedTime / 60) - (hours * 60);
	int seconds = elapsedTime - (hours * 3600) - (minutes * 60);

	std::cout << lg << "Run time: " << std::setfill('0') << std::setw(2) << hours << ":"
			<< std::setfill('0') << std::setw(2) << minutes << ":"
			<< std::setfill('0') << std::setw(2) << seconds << std::endl;
	return 0;

}
