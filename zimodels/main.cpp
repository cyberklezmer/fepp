/*
 * main.cpp
 *
 *  Created on: 10 Feb 2010
 *      Author: martin
 */

#include "zimodels.hpp"
#include <iomanip>

#include "epp.hpp"


using namespace epp;


int main(int argc, char ** argv)
{
    const char *ms[] =
    {
        "NASDAQ OMX BX",
        "NSE" ,
        "Chicago",
        "NYSE",
        "ARCA",
        "NASDAQ T",
        "CBOE",
        "BATS" ,
        "ISE",
		0
    };

	const char *ss[] =
	{
 //"JCOM",
         "MKTX",
        "JCOM",
        "PNY",
        "ACU",
        "ARL",
        "FCCY",
 0,
        "XOM",
        "MSFT",
        "GE",
        "MKTX",
        "JCOM",
        "PNY",
        "ACU",
        "ARL",
        "FCCY",
        0
    };


	int todo;

    if(argc <= 2)
        todo = -1;
    else
    {
        switch(argv[1][0])
        {
            case 'I':
                todo = 0;
                break;
            case 'S':
                todo = 1;
                break;
            case 'E':
                todo = 2;
                break;
            case 'z':
                todo = 10;
                break;
            case 'Z':
                todo = 20;
                break;
            case 'g':
                todo = 11;
                break;
            case 'G':
                todo = 21;
                break;
            case 'n':
                todo = 12;
                break;
            case 'N':
                todo = 22;
                break;

            default:
                todo = -1;
        }
    }

    if(todo == -1 || (todo == 0 && argc < 4))
    {
        cout << "zimodels action homefolder [importfolder]" << endl;
        return 1;
    }

    bool onestock = argc >= 4;
    bool onemarket = argc >= 5;

    for(int i=1; i<argc; i++)
    {
        for(int j=0; argv[i][j]; j++)
            if(argv[i][j]=='_')
                argv[i][j] = ' ';
    }

    hfd h(argv[2]);
 	date start(2009,3,1);
	date end(2009,4,30);

    if(todo == 0)
    {
        string folder = argv[3];
 //       importer09::import(ss,ms,start,end,folder);
        importertw::import(ss,ms,start,end,folder);
          return 0;
    }

    vector<smpair> pairs;

	for(int s=0; ss[s]; s++)
		for(int m=0; ms[m]; m++)
		{
            if(onestock && strcmp(ss[s],argv[3]))
                continue;
            if(onemarket && strcmp(ms[m],argv[4]))
                continue;
			smpair p(ss[s],ms[m],9.5,15.5);
			pairs.push_back(p);
		}
	// to do  - add other pairs

    if(pairs.size()==0)
    {
        if(onestock && !onemarket)
            cout << "Cannot find stock " << argv[3] << endl;
        else if(onestock && onemarket)
            cout << "Cannot find pair " << argv[3] << " " << argv[4] << endl;
        else
            cout << "Internal error" << endl;
        return 1;
    }

    ostringstream name;
    bool less = argv[1][1]=='-';
    if(todo >= 10)
    {
        if(less)
            name << "1000";
        else
            name << "5000";
        if(todo < 20)
            name << "kappa";
        else
            name << "phi";
        if(onestock)
            name << argv[3];
        if(onemarket)
            name << argv[4];
    }

	switch(todo)
	{
		case 1:
		{
			zisummary s(2008,2009);
			program pgm(pairs,s,start,end);
			pgm.process();
		}
		break;
        case 2:
        {
            tbtexporter t;
            date end(2009,3,1);
            program pgm(pairs,t,start,end);
            pgm.process();
        }
        break;
        case 10:
		case 20:
        {
            name << "z";
            zianalysis e(name.str());
            e.phipar = todo==10 ? false: true;

            e.modeltype = zianalysis::tail;
            e.maxmletime = 2500;
            e.unitvolume = true;
            e.includenu = false;
            e.includegamma = false;
            e.includezeta = false;
            e.includeeta = false;
            e.setsample(less ? 1000 : 5000,false);
            e.resample = false;
            e.firstn = 1;
            e.maxn = less ? 2 : 3;

            e.twodimestimation = true;
            e.extendedlogging = false;

            program pgm(pairs,e,start,end);
            pgm.process();
            return 0;
        }

		case 11:
		case 21:
		{
            name << "g";
            zianalysis e(name.str());
            e.phipar = todo==11 ? false: true;
            e.modeltype = zianalysis::tail;
            e.maxmletime = 2500;
            e.unitvolume = false;
            e.includenu = false;
            e.includegamma = false;
            e.includezeta = false;
            e.includeeta = true;
            e.setsample(less ? 1000 : 5000,true);
            e.firstn = 1;
            e.maxn = less ? 2 : 3;
            e.twodimestimation = true;
            e.extendedlogging = false;

            program pgm(pairs,e,start,end);
            pgm.process();
            return 0;
		}
		case 12:
		case 22:
        {
            name << "n";
            zianalysis e(name.str());
            e.phipar = todo==12 ? false: true;

            e.modeltype = zianalysis::tail;
            e.maxmletime = 2500;
            e.unitvolume = true;
            e.includenu = false;
            e.includegamma = false;
            e.includezeta = false;
            e.includeeta = false;
            e.setsample(less ? 1000 : 5000,true);
            e.resample = false;
            e.firstn = 1;
            e.maxn = less ? 2 : 3;
            e.twodimestimation = true;
            e.extendedlogging = true;

            program pgm(pairs,e,start,end);
            pgm.process();
            return 0;
        }

		break;
	}
	cout << "Finished!" << endl;
	return 0;
}


