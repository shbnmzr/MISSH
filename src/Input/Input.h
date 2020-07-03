#ifndef SRC_INPUTCLASS_H_
#define SRC_INPUTCLASS_H_

#include <memory>
#include <unordered_map>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include "../Utilities/Utilities.h"

enum FileType {Fasta, Fastq, UnknownFile};
enum PairType {SingleEnd, PairedEnd, UnknownPair};

using namespace std;

struct SingleEndFile
{
	typedef shared_ptr<SingleEndFile> Ptr;
	SingleEndFile();
	bool init(const string& path);
	bool isCorrect() const;
	//Get
	const string& getPath() const;
	const vector<string>& getPathParse() const;
	const string getDirectory() const;
	const string& getFilename() const;
	const string getFilenameWithoutExt() const;
	const string getExt() const;
	FileType getFileType() const;
	const char& getSequenceDelimiter() const;
	void reset();

private:
	string path;
	vector<string> path_parse;
	FileType file_type;
	char read_delimiter;
	bool open_correct;
};

struct PairFiles : public pair<SingleEndFile, SingleEndFile>
{
	PairFiles();
	bool init(const string& path_end1, const string& path_end2);
	bool isCorrect(); //return if the file or files opened isCorrect
	PairType getPairType() const;
	FileType getFileType() const;
	const string& getIdentify() const;
	void reset();

private:
	PairType pair_type;
	string identify;
	bool open_correct;
};

#endif /* SRC_INPUTCLASS_H_ */
