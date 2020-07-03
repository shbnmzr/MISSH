/*
 * FileScan.cpp
 *
 *  Created on: 09/feb/2016
 *      Author: samuele
 */

#include "FileScan.h"

FileScan::FileScan() {
}

FileScan::~FileScan() {
}

// copy assignment
FileScan& FileScan::operator=(const FileScan& other) {
	this->init(other.file.getPath());
    return *this;
}

bool FileScan::init(const string& path) {
	bool ret = this->file.init(path);
	if(ret)
		ret &= this->getStreamPosRead(this->file.getPath());
	return ret;
}

size_t FileScan::getSequenceNumber() const {
	return this->pos_read.size();
}

bool FileScan::getSequenceWithIndex(size_t index, Sequence& read, Parser parserHeader)
{
	if(index >= this->pos_read.size())
		return false;
	streampos actual = this->stream.tellg();
	if(this->stream.eof() || this->stream.fail() || actual != this->pos_read[index].start)
	{
		this->stream.clear();
		this->stream.seekg(this->pos_read[index].start);
	}
	bool parseHeader = parserHeader != NULL;
	string line = "";
	bool isLineQuality = false;
	while(actual != this->pos_read[index].end && getline(this->stream,line))
	{
		if(!line.empty())
		{
			if(line[0] == this->file.getSequenceDelimiter())
			{
				//due casi:
				//o delimitatore non fa parte di qualità se isLineQuality == false, quindi line contiene la descrizione del read
				//o delimitatore fa parte di qualità se isLineQuality == true e line contiene la riga qualità
				if(!isLineQuality)
				{
					read.setIndexFile(index);
					if(parseHeader)
						read.appendHeader(line, parserHeader);
					else
						read.appendHeader(line);
				}
				else //altrimenti line contiene la sequenza qualità che inizia per il delimitatore
				{
					read.appendQuality(line);
					isLineQuality = false;
				}
			}
			else
			{
				//caso line non inizia per il delimitatore tre casi:
				//o line contiene il + (header qualità) quindi prossima line sarà qualità, si imposta isLineQuality == true
				//o line contiene la line di qualità quindi isLineQuality == true
				//o line contiene la line di sequenza quindi isLineQuality == false
				if(this->file.getFileType() == Fastq && line[0] == '+')
				{
					read.appendHeaderQuality(line);
					isLineQuality = true; //prossima line sarà linea qualità
				}
				else
				{
					if(!isLineQuality) //Non è line qualità quindi ha la sequenza da memorizzare
						read.appendSequence(line);
					else
					{
						read.appendQuality(line);
						isLineQuality = false;
					}
				}
			}
		}
		actual = this->stream.tellg();
	}
	return true;
}

bool FileScan::getStreamPosRead(const string& path)
{
	this->stream.open(path);
	string empty = "";
	string line; //store the read data
	if(this->stream.is_open())
	{
		bool isLineQuality = false;
		streampos beforeCallGetLine = this->stream.tellg();
		while(getline(this->stream,line))
		{
			if(!line.empty())
			{
				if(line[0] == this->file.getSequenceDelimiter())
				{
					//due casi:
					//o delimitatore non fa parte di qualità se isLineQuality == false, quindi line contiene la descrizione del read
					//o delimitatore fa parte di qualità se isLineQuality == true e line contiene la riga qualità
					if(!isLineQuality)
					{
						//////////////////SAVE///////////////////
						if(!this->pos_read.empty()) //End previus is start actual
							this->pos_read[this->pos_read.size()-1].end = beforeCallGetLine;
						PositionRead pos;
						pos.start = beforeCallGetLine;
						this->pos_read.push_back(pos);
						//////////////////END_SAVE///////////////
					}
					else //altrimenti line contiene la sequenza qualità che inizia per il delimitatore
					{
						isLineQuality = false;
					}
				}
				else
				{
					//caso line non inizia per il delimitatore tre casi:
					//o line contiene il + (header qualità) quindi prossima line sarà qualità, si imposta isLineQuality == true
					//o line contiene la line di qualità quindi isLineQuality == true
					//o line contiene la line di sequenza quindi isLineQuality == false
					if(this->file.getFileType() == Fastq && line[0] == '+')
						isLineQuality = true; //prossima line sarà linea qualità
					else
						isLineQuality = false;
				}
			}
			beforeCallGetLine = this->stream.tellg(); //Sono già su nuova riga
		}
		//////////////////SAVE///////////////////
		//Save last end
		if(!this->pos_read.empty()) //End previus is start actual
			this->pos_read[this->pos_read.size()-1].end = beforeCallGetLine;
		//////////////////END_SAVE///////////////
		this->stream.close();
		this->stream.open(path, ifstream::binary);
		return true;
	}
	else
	{
		if(!this->file.getPath().empty())
			cout<<"Fail to open "<< this->file.getPath() << endl << flush;
		return false;
	}
}

const SingleEndFile& FileScan::getFile() const {
	return file;
}

bool FileScan::isCorrect() const {
	return this->file.isCorrect();
}

void FileScan::reset() {
	this->file.reset();
	ifstream().swap(this->stream);
	this->pos_read.clear(); this->pos_read.shrink_to_fit();
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////7

FilesScan::FilesScan() {
	this->reset();
}

FilesScan& FilesScan::operator=(const FilesScan& other) {
	this->init(other.first.getFile().getPath(), other.second.getFile().getPath());
	return *this;
}

void FilesScan::reset() {
	this->pair_type = UnknownPair;
	this->identify = "";
	this->open_correct = false;
	this->first.reset();
	this->second.reset();
}

bool FilesScan::init(const string& path_end1, const string& path_end2) {
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
			this->identify = this->first.getFile().getFilename();
		else if (this->pair_type == PairedEnd)
			this->identify = this->first.getFile().getFilenameWithoutExt();
	}

	//Check correct
	this->open_correct = init1;//almeno il primo
	if(!path_end1.empty() && !path_end2.empty())//dovrebbe essere Paired end
	{
		if(this->pair_type == PairedEnd)
		{
			if(this->first.getSequenceNumber() != this->second.getSequenceNumber())
			{
				this->open_correct &= false;
				cerr << "Different number of read in file " << this->first.getFile().getPath() << " and " << this->second.getFile().getPath() << endl << flush;
			}
			if(this->first.getFile().getFileType() != this->second.getFile().getFileType())
			{
				this->open_correct &= false;
				cerr << "Different file type Fasta and Fastq" << endl << flush;
			}
		}
	}
	return this->open_correct;
}

FilesScan::~FilesScan() {
	// TODO Auto-generated destructor stub
}

size_t FilesScan::getPairedReadsNumber() const {
	return this->first.getSequenceNumber();
}

size_t FilesScan::getSequencesNumber() const {
	return this->first.getSequenceNumber() + this->second.getSequenceNumber();
}

const string& FilesScan::getIdentify() const {
	return identify;
}

PairType FilesScan::getPairType() const {
	return pair_type;
}

FileType FilesScan::getFileType() const {
	return this->first.getFile().getFileType();
}

bool FilesScan::isCorrect() const {
	return this->open_correct;
}
