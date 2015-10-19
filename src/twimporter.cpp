/*
 * tickdata09importer.cpp
 *
 *  Created on: Mar 26, 2013
 *      Author: martin
 */

#include "hfd.hpp"
#include <sys/stat.h>
#include <sstream>


void tickdatatwimporter::importstock(string& stock, vector<string>& markets,  vector<char>& marketids)
{
	if(type!=usstock)
	{
		cerr << "Only US stocks importable from tickdata by the current version " << endl;
		throw 1;
	}
	stringstream t;
    stringstream q;

	q << folder << "/QUOTES/" << stock << ".csv";
	t << folder << "/TRADES/" << stock << ".csv";

	importusfile(stock,markets,marketids,q.str(),t.str());
}



