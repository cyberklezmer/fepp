/*
 * analysis.cpp
 *
 *  Created on: 22 Apr 2010
 *      Author: martin
 */


#include "hfd.hpp"

void analysis::start()
{
	if(!latex)
	{
		cout << "cannot open the latex file " << gettexfilename(name) << endl;
		throw 1;
	}
	else
        cout << "Writing latex to " << gettexfilename(name) << endl;

	latex << "\\tiny" << endl;

	latex.setf(ios::fixed,ios::floatfield);

	if(!csv)
	{
		cout << "cannot open the csv file " << getcsvfilename(name) << endl;
		throw 1;
	}
	else
        cout << "Writing csv to " << getcsvfilename(name) << endl;


	onstart();
	column = 0;
	row = 0;
}

void analysis::startpair(smpair& apair)
{
	pair = apair;
	if(column==0)
	{
		latex << "\\begin{tabular}{";
		for(int i=0; i<numcols; i++)
		{
			latex << "c";
			if(i < numcols-1)
				latex << "|";
		}
		latex << "} \\hline " << endl;
	}
	else
		latex << "&";
	latex << "\\begin{minipage}" <<  "{4.5cm}" << endl;
	latex << "\\vspace{3mm}" << endl;

	onstartpair(pair);
}


void analysis::record(hfdrecord& rec, int recn)
{
	onrecord(rec,recn);
	history.push_back(rec);

	if(data::isqouterecord(rec))
		qhistory.push_back(rec);
}


void analysis::finishtable(bool anewpage)
{
	latex << "\\end{tabular}" << endl;
	if(anewpage)
		latex << "\\newpage" << endl;
}

void analysis::endpair(bool aexceptionthrown, const char *aexctext)
{
    bool exc = aexceptionthrown;
    // being static in order not to allocate memory in enpair when it is called inside a chatch block
    static char etxt[1000];

    const char* etext = aexceptionthrown? aexctext : 0;
    if(!exc)
        try
        {
            try
            {
                onendpair(pair);
            }
            catch(std::exception& e)
            {
                strncpy(etxt,e.what(),sizeof(etxt)-1);
                etext = etxt;
                throw;
             }
        }
        catch(...)
        {
            exc = true;
        }
    if(exc)
    {
        latex << "Exception occured " << endl;
        if(etext)
            latex << etext << endl;
    }
	latex << "\\end{minipage}" << endl;
	if(column==numcols - 1)
	{
		bool newpage = (row++==numrows-1);
		finishtable(newpage);
		column=0;
		if(newpage)
			row=0;
	}
	else
		column++;
}

void analysis::end()
{
	if(column != 0)
		finishtable(false);
	onend();
}
