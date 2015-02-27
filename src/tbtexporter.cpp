/*
 * xiexporter.cpp
 *
 *  Created on: 27 Apr 2010
 *      Author: martin
 */

#include "hfd.hpp"

void tbtexporter::onstartpair(smpair& togo)
{
	ofs = new ofstream(csvfn().c_str());
	if(!(*ofs))
	{
		cout << "Cannot open " << csvfn() << endl;
		throw 1;
	}
	ofs->setf(ios::fixed);
	ofs->precision(3);
	*ofs << "D,M,Y,T,B,A,BN,AN,Q,P,QT" << endl;
}


void tbtexporter::onrecord(hfdrecord& rec, int recn)
{
    *ofs << getday() << ","
         << getmonth() << ","
         << getyear() << ","
         << rec.t << ","
         << rec.b << ","
         << rec.a << ","
         << rec.bn << ","
         << rec.an << ","
         << rec.q << ","
         << rec.p << ","
         << rec.qt << endl;
}

void tbtexporter::onendpair(smpair& togo)
{
	delete ofs;
	ofs = 0;
}
