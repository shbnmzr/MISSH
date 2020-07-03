/*
 * Sequence.h
 *
 *  Created on: 09/feb/2016
 *      Author: samuele
 */

#ifndef INPUT_SEQUENCE_H_
#define INPUT_SEQUENCE_H_

#include <string>
using namespace std;

class Sequence {
public:
	typedef void (Sequence::*Parser)();
	Sequence();
	virtual ~Sequence();

	inline size_t getIndexFile() const {
		return indexFile;
	}

	inline void setIndexFile(size_t indexFile) {
		this->indexFile = indexFile;
	}

	inline const string& getHeader() const {
		return header;
	}

	inline void appendHeader(const string& header) {
		this->header.append(header);
	}

	inline void appendHeader(const string& header, void (Sequence::*parser)()) {
		this->header.append(header);
		(this->*parser)();
	}

	inline const string& getSequence() const {
		return sequence;
	}

	inline void appendSequence(const string& sequence) {
		this->sequence.append(sequence);
	}

	inline const string& getHeaderQuality() const {
		return headerQuality;
	}

	inline void appendHeaderQuality(const string& headerQuality) {
		this->headerQuality.append(headerQuality);
	}

	inline const string& getQuality() const {
		return quality;
	}

	inline void appendQuality(const string& quality) {
		this->quality.append(quality);
	}

	inline const string& getId() const {
		return ID;
	}

	inline void setId(const string& id) {
		ID = id;
	}

	inline const string& getFlagEnd() const {
		return flagEnd;
	}

	inline void setFlagEnd(const string& flagEnd) {
		this->flagEnd = flagEnd;
	}

	inline bool isSequenceAllN() const
	{
		for(size_t i = 0; i < this->sequence.length(); ++i)
			if(this->sequence[i] != 'N')
				return false;
		return true;
	}

	inline bool haveSequencePercent_N(double perc) const
	{
		size_t contN = 0;
		for(size_t i = 0; i < this->sequence.length(); ++i)
			if(this->sequence[i] == 'N')
				++contN;
		double percN = (double)contN/(double)this->sequence.length();
		return percN >= perc;
	}

	inline void parser1() {
		string copy_header = this->header;
		//Find next ID
		string delimiter = ">r";
		size_t pos = copy_header.find(delimiter);
		copy_header.erase(0, pos + delimiter.length());
		delimiter = ".";
		pos = copy_header.find(delimiter);

		//Save next ID
		this->setId(copy_header.substr(0, pos));

		//Erase until IDFile
		copy_header.erase(0, pos + delimiter.length());

		//The sequence ID and the SOURCE are the same for both end(if paired End)
		this->setFlagEnd(copy_header.substr(0, 1));
	}

	inline void parser2()	{
		string copy_header = this->header;
		//Find next ID
		string delimiter = "@";
		size_t pos = copy_header.find(delimiter);
		copy_header.erase(0, pos + delimiter.length());
		delimiter = "_";
		pos = copy_header.find(delimiter);

		//Save next ID
		this->setId(copy_header.substr(0, pos));

		//Erase until IDFile
		copy_header.erase(0, pos + delimiter.length());

		//The sequence ID and the SOURCE are the same for both end(if paired End)
		this->setFlagEnd("" + copy_header.length() - 1);
	}

private:
	size_t indexFile;
	string header;
	string sequence;
	string headerQuality;
	string quality;
	
	string ID;
	string flagEnd;
};

#endif /* INPUT_SEQUENCE_H_ */
