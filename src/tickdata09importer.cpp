/*
 * tickdata09importer.cpp
 *
 *  Created on: Mar 26, 2013
 *      Author: martin
 */

#include "hfd.hpp"
#include <sys/stat.h>
#include <sstream>


void tickdata09importer::importstock(string& stock, vector<string>& markets,  vector<char>& marketids)
{
	if(type!=usstock)
	{
		cerr << "Only US stocks importable from tickdata by the current version " << endl;
		throw 1;
	}
	int year=start.y;
	int month=start.m;
	int day=start.d;

	for(;;)
	{
		stringstream s;
		s << folder << "/US/" << year << "/";
		if(month < 10)
			s << "0";
		s << month;

		stringstream q;
		q << s.str();
		s << "/QUOTES/" << year <<  "_";
		q << "/TRADE/" << year <<  "_";
		if(month < 10)
		{
			q << "0";
			s << "0";
		}
		s << month << "_";
		q << month << "_";

		stringstream unpacked[2];
		bool zipmissing = false;
		for(int k=0; k<2 && !zipmissing ; k++)
		{
			stringstream ps;
			ps << (k==0 ? s.str() : q.str());
			if(day<10)
				ps << "0";
			ps << day;
			if(k==0)
				ps << "_X_Q_" << stock[0];
			else
				ps << "_X_T";
			ps << ".zip";
			struct stat buffer;
			if ( !stat( ps.str().c_str(), &buffer) )
			{
		// now unpack it
				string fz = zippath() + " e " + ps.str() + " -y -o" + hfd::tempdir();
				if(system(fz.c_str()))
				{
					cout << "Error executing " << fz << endl;
					throw 1;
				}
			}
			else
			{
				cerr << ps.str() << " not found." << endl;
				zipmissing = true;
				break;
			}
			unpacked[k] <<  hfd::tempdir() << "/" << stock << "_" << year << "_";
			if(month < 10)
				unpacked[k] << "0";
			unpacked[k] << month << "_";
			if(day < 10)
				unpacked[k] << "0";
			unpacked[k] << day << (k==0 ? "_X_Q" : "_X_T");
            string unpzip = unpacked[k].str() + ".zip" ;
			string zzz = zippath() + " e " + unpzip + " -y -o" + hfd::tempdir();
			if(system(zzz.c_str()))
			{
				cout << "Error executing " << zzz << endl;
				cout << "when importing " << stock << endl;
				throw 1;
			}
			int r = remove(unpzip.c_str());
			if(r != 0)
                cerr << "Cannote remove " << unpzip << endl;
            else
                cout << unpzip << " removed." << endl;

		}
		if(!zipmissing)
		{
            stringstream u1;
            u1 << unpacked[0].str();
            u1 << ".asc";
            stringstream u2;
            u2 << unpacked[1].str();
            u2 << ".asc";

			importusfile(stock,markets,marketids,u1.str(),u2.str());
			int r1 = remove(u1.str().c_str());
			if(r1 != 0)
                cerr << "Cannote remove " << u1.str() << endl;
            else
                cout << u1.str() << " removed." << endl;
			int r2 = remove(u2.str().c_str());
			if(r2 != 0)
                cerr << "Cannote remove " << u2.str() << endl;
            else
                cout << u2.str() << " removed." << endl;
		}

		if(++day > 31)
		{
			day = 1;
			if(++month > 12)
			{
				month=1;
				year++;
			}
		}
		date newdate(year,month,day);
		if(!(newdate <= end))
			break;
	}
}



