/*
This file is part of ``offset_surface'', a library for point based surface mesh offset.
Copyright (C) 2019 Bill He <github.com/easterngarden>

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
3. The name of the author may not be used to endorse or promote products
   derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY
OF SUCH DAMAGE.
*/
// Author : Bill He

#include "kdprocessor..h"
#include <algorithm>
#include <assert.h> 
#include "kdtree.h"
#include <string.h>

#define	MAX_LINE	256

double stringToFloat(std::string str)
{
	std::string::iterator it = str.begin();
	for (; it != str.end(); ) {
		if ((*it) == ' ')
			++it;
		else
			break;
	}
	++it;
	for (; it != str.end(); ++it) {
		if ((*it) == '-')
			break;
	}
	if (it != str.end()) {
		if (isdigit(*(--it))) {
			++it;
			str.insert(it, 1, 'e');
		}
	}
	return atof(str.c_str());
}

kdprocessor::~kdprocessor() {
	kd_free((struct kdtree *)pkdtree);
}

int kdprocessor::buildtree(const std::string & file, double &dOffsetBValue)
{
	/* create a k-d tree for 3-dimensional points */
	pkdtree = kd_create(3);

	std::ifstream ifs(file.c_str());
	assert(ifs.good());
	std::map< unsigned int, std::vector<double> > ntsMap;
	std::string line;
	while (std::getline(ifs, line))
	{
		// skip empty lines:
		if (line.empty())
			continue;

		std::string sstr;
		if (strncmp(line.c_str(), "GRID ", 5) == 0)
		{
			sstr.assign(line);
			unsigned int id = atoi((sstr.substr(8, 8)).c_str());
			Pos point;
			point.x = stringToFloat(sstr.substr(24, 8));
			point.y = stringToFloat(sstr.substr(32, 8));
			point.z = stringToFloat(sstr.substr(40, 8));
			nodeMap.insert(std::make_pair(id, point));
		}
		else if (strncmp(line.c_str(), "CTRIA3 ", 7) == 0)
		{
			std::vector<unsigned int> vec(3);
			sstr.assign(line);
			vec[0] = atoi((sstr.substr(24, 8)).c_str());
			vec[1] = atoi((sstr.substr(32, 8)).c_str());
			vec[2] = atoi((sstr.substr(40, 8)).c_str());
			std::vector<double> thickness(3);
			if (std::getline(ifs, line))
			{
				sstr.assign(line);
				thickness[0] = stringToFloat(sstr.substr(24, 8));
				thickness[1] = stringToFloat(sstr.substr(32, 8));
				thickness[2] = stringToFloat(sstr.substr(40, 8));
				for (int m = 0; m < 3; m++)
				{
					unsigned int key = vec[m];
					if (ntsMap.find(key) != ntsMap.end())
					{
						n2tMap[key] = thickness[m];
						std::vector<double>&  v = ntsMap[key];
						v.push_back(thickness[m]);
					}
					else
					{
						n2tMap[key] = thickness[m];
						std::vector<double> v;
						v.push_back(thickness[m]);
						ntsMap.insert(std::make_pair(key, v));
					}
				}
			}
			else
			{
				assert(false);
				return 1;
			}
		}
		else if (strncmp(line.c_str(), "CTRIA6 ", 7) == 0) 
		{
			std::vector<unsigned int> vec(6);
			sstr.assign(line);
			vec[0] = atoi((sstr.substr(24, 8)).c_str());
			vec[1] = atoi((sstr.substr(32, 8)).c_str());
			vec[2] = atoi((sstr.substr(40, 8)).c_str());
			vec[3] = atoi((sstr.substr(48, 8)).c_str());
			vec[4] = atoi((sstr.substr(56, 8)).c_str());
			vec[5] = atoi((sstr.substr(64, 8)).c_str());
			std::vector<double> thickness(6);
			if (std::getline(ifs, line))
			{
				sstr.assign(line);
				thickness[0] = stringToFloat(sstr.substr(24, 8));
				thickness[1] = stringToFloat(sstr.substr(32, 8));
				thickness[2] = stringToFloat(sstr.substr(40, 8));
				thickness[3] = thickness[0];
				thickness[4] = thickness[1];
				thickness[5] = thickness[2];
				for (int m = 0; m < 6; m++)
				{
					unsigned int key = vec[m];
					if (ntsMap.find(key) != ntsMap.end())
					{
						n2tMap[key] = thickness[m];
						std::vector<double>&  v = ntsMap[key];
						v.push_back(thickness[m]);
					}
					else
					{
						n2tMap[key] = thickness[m];
						std::vector<double> v;
						v.push_back(thickness[m]);
						ntsMap.insert(std::make_pair(key, v));
					}
				}
			}
			else
			{
				assert(false);
				return 1;
			}
		}
	}
	ifs.close();

	double factor = dOffsetBValue > 0 ? 1.0 : -1.0;
	bool undefined = true;
	std::map< unsigned int, double >::const_iterator ntIter;
	for (ntIter = n2tMap.begin(); ntIter != n2tMap.end(); ntIter++)
	{
		unsigned int key = ntIter->first;
		if (ntsMap.find(key) != ntsMap.end())
		{
			std::vector<double>&  v = ntsMap[key];
			int siz = v.size();
			if (siz > 0)
			{
				double dThickness = 0.0;
				for (int m = 0; m < siz; m++)
				{
					dThickness += v[m];
				}
				n2tMap[key] = dThickness / (double)siz;
				if (undefined)
				{
					dOffsetBValue = factor * n2tMap[key];
					undefined = false;
				}
				else
				{
					if (dOffsetBValue < factor * n2tMap[key])
					{
						dOffsetBValue = factor * n2tMap[key];
					}
				}
			}
		}
		std::map< unsigned int, Pos >::const_iterator iter;
		iter = nodeMap.find(key);
		if (iter != nodeMap.end())
		{
			Pos point = iter->second;
			assert(0 == kd_insert3((struct kdtree *)pkdtree, point.x, point.y, point.z, &n2tMap[key]));
		}
		else
		{
			assert(false);
			return 1;
		}
	}
	return 0;
}

double kdprocessor::findOffset(double x, double y, double z)
{
	struct kdres *presult = kd_nearest3((struct kdtree *)pkdtree, x, y, z);
	if (presult != 0)
	{
		// get the data and position of the current result item 
		double pos[3];
		pos[0] = x;
		pos[1] = y;
		pos[2] = z;
		double *ptdata = (double*)kd_res_item(presult, pos);
		if (ptdata)
		{
			return *ptdata;
		}
		else
		{
			assert(false);
		}
	}
	else
	{
		assert(false);
	}
	return 0;
}
