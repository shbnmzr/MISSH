/*---------------------------------------//
Authors: Le Van Vinh
Faculty of Information technology
Hochiminh City University of Technology and Education
and
Faculty of Computer Science and Engineering
Hochiminh City Univeristy of Technology
//---------------------------------------*/
#ifndef	__UTILITIES_CPP__
#define __UTILITIES_CPP__

#include "Utilities.h"

bool file_exist(string& path)
{
    ifstream infile(path);
    bool ret = infile.good();
    infile.close();
    return ret;
}

bool getLines(string& path, vector<string>& lines) {
	ifstream infile(path);
	if(infile.is_open())
	{
		string line = "";
		while(getline(infile, line))
			lines.push_back(line);
		infile.close();
		return true;
	}
	else
		return false;
}

void createDirAndSubDir(string path) {
	string dir = "";
	string delimiter = "/";
	size_t pos = path.find(delimiter);
	while(pos != path.npos)
	{
		dir += path.substr(0, pos + delimiter.length());
		path.erase(0, pos + delimiter.length());
		mkdir(dir.c_str(), S_IRWXU);
		pos = path.find(delimiter);
	}
}

void parseLine(string line, vector<string>& lineV, vector<string> delimiter)
{
	//Find min
	size_t index_min = numeric_limits<size_t>::max();
	size_t min_pos = numeric_limits<size_t>::max();
	for(size_t i = 0; i < delimiter.size(); ++i)
	{
		size_t tmp_pos = line.find(delimiter[i]);
		if(min_pos > tmp_pos)
		{
			index_min = i;
			min_pos = tmp_pos;
		}
	}

	size_t pos = line.find(delimiter[index_min]);
	while(pos != line.npos)
	{
		lineV.push_back(line.substr(0,pos));
		line.erase(0, pos + delimiter[index_min].length());
		pos = line.find(delimiter[index_min]);
	}
	lineV.push_back(line);
}

string LCSubstr(const string& x, const string& y){
    int m = x.length(), n=y.length();

    uint LCSuff[m][n];

    for(int j=0; j<=n; j++)
        LCSuff[0][j] = 0;
    for(int i=0; i<=m; i++)
        LCSuff[i][0] = 0;

    for(int i=1; i<=m; i++){
        for(int j=1; j<=n; j++){
            if(x[i-1] == y[j-1])
                LCSuff[i][j] = LCSuff[i-1][j-1] + 1;
            else
                LCSuff[i][j] = 0;
        }
    }

    string longest = "";
    for(int i=1; i<=m; i++){
        for(int j=1; j<=n; j++){
            if(LCSuff[i][j] > longest.length())
                longest = x.substr((i-LCSuff[i][j]+1) -1, LCSuff[i][j]);
        }
    }
    return longest;
}

int parseLineForMemory(char* line){
    int i = strlen(line);
    while (*line < '0' || *line > '9') line++;
    line[i-3] = '\0';
    i = atoi(line);
    return i;
}


int getVirtualMemoryUsed(){ //Note: this value is in KB!
    FILE* file = fopen("/proc/self/status", "r");
    int result = -1;
    char line[128];


    while (fgets(line, 128, file) != NULL){
        if (strncmp(line, "VmSize:", 7) == 0){
            result = parseLineForMemory(line);
            break;
        }
    }
    fclose(file);
    return result;
}

int getPeakVirtualMemoryUsed() { //Note: this value is in KB!
    FILE* file = fopen("/proc/self/status", "r");
    int result = -1;
    char line[128];

    while (fgets(line, 128, file) != NULL){
        if (strncmp(line, "VmPeak:", 7) == 0){
            result = parseLineForMemory(line);
            break;
        }
    }
    fclose(file);
    return result;
}

int getPhysicalMemoryUsed(){ //Note: this value is in KB!
    FILE* file = fopen("/proc/self/status", "r");
    int result = -1;
    char line[128];


    while (fgets(line, 128, file) != NULL){
        if (strncmp(line, "VmRSS:", 6) == 0){
            result = parseLineForMemory(line);
            break;
        }
    }
    fclose(file);
    return result;
}

#endif
