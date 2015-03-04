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
        "XOM",
        "MSFT",
        "GE",
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
            case 'Z':
                todo = 11;
                break;
            case 'G':
                todo = 12;
                break;
            case 'X':
                todo = 20;
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
        importer09::import(ss,ms,start,end,folder);
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
		case 11:
        {
            zianalysis e("zi");
            e.modeltype = zianalysis::tail;
            e.maxmletime = 600;
            e.unitvolume = false;
            e.includenu = false;
            e.includegamma = false;
            e.includezeta = false;
            e.includeeta = false;
            e.setsample(1000,false);
            e.firstn = 1;
            e.maxn = 2;
            e.twodimestimation = true;
            e.extendedlogging = false;

            program pgm(pairs,e,start,end);
            pgm.process();
            return 0;
        }
		case 12:
		{
            zianalysis e("g");
            e.modeltype = zianalysis::tail;
            e.maxmletime = 600;
            e.unitvolume = false;
            e.includenu = false;
            e.includegamma = false;
            e.includezeta = false;
            e.includeeta = false;
            e.setsample(1000,true);
            e.firstn = 1;
            e.maxn = 2;
            e.twodimestimation = true;
            e.extendedlogging = false;

            program pgm(pairs,e,start,end);
            pgm.process();
            return 0;
		}
        case 20:
        {
            zinparmodel m(1,true,false,false,false);
            zianalysis e("zi");

            vector<double> pars(3);
            pars[0]=0.06;
            pars[1]=0.14;
            pars[2]=-1.64;
            e.simulateby = &m;
            m.setinitparams(pars);
            return 0;
        }

		break;
	}
	cout << "Hotovo!" << endl;
	return 0;
}


/*    double d=1;
    for(int i=-1; i<2; i++)
    {
        for(int j=-1; j<2; j++)
        {
            vector<double> grad(zero_vector(2));
            vector<double> g(2);
            vector<double> x(2);
            x[0] = p[0].initial + i*d;
            x[1] = p[1].initial + j*d;
            double s=0;
            for(int n=0; n<9; n++)
            {
                s+= m.evallogdensity(n,x,g);
                grad = grad + g;
            }
            cout << x[0] << " " << x[1] << " = " << s
                 << " (" << grad[0] << "," << grad[1] << ")" << endl;
        }
        cout << endl;
    }*/
