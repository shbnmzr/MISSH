/*
 * InputClass.cpp
 *
 *  Created on: 10/apr/2015
 *      Author: Samuele Girotto
 */

#include "Input.h"

void PairFiles::reset() {
	this->pair_type = UnknownPair;
	this->identify = "";
	this->open_correct = false;
}

PairFiles::PairFiles() {
	this->reset();
}

bool PairFiles::init(const string& path_end1, const string& path_end2) {
	this->reset();

	bool init1 = this->first.init(path_end1); //Prova prima path
	if(init1)
	{
		this->pair_type = SingleEnd;
		bool init2 = this->second.init(path_end2);
		if(init2)
			this->pair_type = PairedEnd;
	}
	else
	{
		init1 = this->first.init(path_end2); //Prova seconda path
		if(init1)
			this->pair_type = SingleEnd;
		else
			this->pair_type = UnknownPair;
	}

	if(init1)
	{
		if(this->pair_type == SingleEnd)
			this->identify = this->first.getFilename();
		else if (this->pair_type == PairedEnd)
			this->identify = this->first.getFilenameWithoutExt();
	}
	//Check correct
	this->open_correct = init1;//almeno il primo
	if(!path_end1.empty() && !path_end2.empty())//dovrebbe essere Paired end
	{
		if(this->pair_type == PairedEnd)
		{
			if(this->first.getFileType() != this->second.getFileType()) //Stesso tipo di file
			{
				cerr << "Different file type Fasta and Fastq" << endl << flush;
				this->open_correct &= false;
			}
			this->open_correct &= true;
		}
	}
	return this->open_correct;
}

bool PairFiles::isCorrect() {
	return this->open_correct;
}

PairType PairFiles::getPairType() const {
	return pair_type;
}

FileType PairFiles::getFileType() const {
	return this->first.getFileType();
}

void SingleEndFile::reset() {
	this->path = "";
	this->path_parse.clear(); this->path_parse.shrink_to_fit();
	this->file_type = UnknownFile;
	this->read_delimiter = '>';
	this->open_correct = false;
}

SingleEndFile::SingleEndFile() {
	this->reset();
}

bool SingleEndFile::init(const string& path)
{
	this->reset();

	string line; //store the read data
	ifstream file(path);
	if(file.is_open())
	{
		//Prendi prima riga e guarda primo elemento per capire la tipologia di file
		getline(file,line);
		if(!line.empty())
		{
			this->path = path;
			parseLine(this->path, this->path_parse, {"/"});
			this->open_correct = true;

			if(line[0] == '>')
			{
				this->file_type = FileType::Fasta;
				this->read_delimiter = '>';
			}
			if(line[0] == '@')
			{
				this->file_type = FileType::Fastq;
				this->read_delimiter = '@';
			}
			file.close();
			return true;
		}
		else
		{
			file.close();
			return false;
		}
	}
	else
	{
		if(path != "")
			cerr << "Fail to open: " << path << endl << flush;
		return false;
	}
}

bool SingleEndFile::isCorrect() const {
	return this->open_correct && this->file_type != UnknownFile;
}

const string& SingleEndFile::getPath() const {
	return path;
}

const vector<string>& SingleEndFile::getPathParse() const {
	return path_parse;
}

const string SingleEndFile::getDirectory() const
{
	string dir = "";
	for(size_t i = 0; i < this->path_parse.size() - 1; ++i)
		dir += this->path_parse[i] + "/";
	return dir;
}

const string& SingleEndFile::getFilename() const {
	return this->path_parse.back();
}

const string SingleEndFile::getFilenameWithoutExt() const {
	string filename_without_ext = this->path_parse.back();
	string delimiter = ".";
	size_t pos = filename_without_ext.rfind(delimiter);
	filename_without_ext.erase(pos, filename_without_ext.length());
	return filename_without_ext;
}

const string SingleEndFile::getExt() const {
	string filename_without_ext = this->path_parse.back();
	string delimiter = ".";
	size_t pos = filename_without_ext.rfind(delimiter);
	return filename_without_ext.substr(pos+1, filename_without_ext.length());
}

FileType SingleEndFile::getFileType() const {
	return file_type;
}

const string& PairFiles::getIdentify() const {
	return identify;
}

const char& SingleEndFile::getSequenceDelimiter() const {
	return read_delimiter;
}
