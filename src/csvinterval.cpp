/*
 * csvinterval.cpp
 *
 *  Created on: Jun 19, 2013
 *      Author: martin
 */


#include "hfd.hpp"

void csvinterval::onstartpair(smpair& togo)
{
	ofs = new ofstream(csvfn().c_str());
	if(!(*ofs))
	{
		cout << "Cannot open " << csvfn() << endl;
		throw 1;
	}
	ofs->setf(ios::fixed);
	ofs->precision(3);
	*ofs << "Day,Month,Year,Time,actA,actB,nfSnapA,nfSnapB,avgA,avgB,dNA,dNB,N" << endl;
}

void csvinterval::onsnapshot(interval::rec& r)
{
	*ofs << r.day << "," << r.month << "," << r.year << ","
		<< r.time << ","
	    << r.actA << "," << r.actB << ","
	    << r.nfSnapA << "," << r.nfSnapB  << ","
	    << r.avgA << "," << r.avgB << ","
	    << r.dNA << "," << r.dNB  << "," << r.N
	    << endl;
}

void csvinterval::onendpair(smpair& togo)
{
	delete ofs;
	ofs = 0;
	interval::onendpair(togo);
}

