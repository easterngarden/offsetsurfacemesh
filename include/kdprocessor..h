#pragma once
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>

struct Pos {
	double x;
	double y;
	double z;
};

class kdprocessor {
public:
	kdprocessor()
	: pkdtree(0)
	{
	}
	~kdprocessor();

	int buildtree(const std::string & file, double &dOffsetBValue);
	double findOffset(double x, double y, double z);

private:
	std::map< unsigned int, Pos >  nodeMap;	//node id <--> point coords
	void *pkdtree;
	std::map< unsigned int, double > n2tMap;	//node id <--> thickness
};
