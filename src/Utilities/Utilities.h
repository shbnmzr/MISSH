/*---------------------------------------//
Authors: Le Van Vinh
Faculty of Information technology
Hochiminh City University of Technology and Education
and
Faculty of Computer Science and Engineering
Hochiminh City Univeristy of Technology
//---------------------------------------*/

#ifndef	__UTILITIES_H__
#define __UTILITIES_H__ 

#include <sys/stat.h>
#include <fstream>
#include <string>
#include <cstring>
#include <limits>
#include <vector>
#include <regex>

using namespace std;

bool file_exist(string& path);
bool getLines(string& path, vector<string>& line);
void createDirAndSubDir(string path);
void parseLine(string line, vector<string>& lineV, vector<string> delimiter);
string LCSubstr(const string& x, const string& y);

int parseLineForMemory(char* line);
int getVirtualMemoryUsed();
int getPeakVirtualMemoryUsed();
int getPhysicalMemoryUsed();
#endif
