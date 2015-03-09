/*
 * data.cpp
 *
 *  Created on: 10 Mar 2010
 *      Author: martin
 */

#include "hfd.hpp"


void data::writerecord(ofstream& ofs, hfdrecord& rec)
{
	ofs.write((char*) (&rec),sizeof(rec));
	//		ofs.precision(10);
	//ofs << rec.t << ", " << rec.a << ", "  << rec.b << ", "<< rec.an << ", "  << rec.bn << ", "  << rec.q << ", " << rec.p << ", " << rec.qt << endl;
}

string data::getfilename(string& astock, string& amarket, date d)
{
    #define ID(s,i,y,m,d) s+'_'+i+'_'+(y)+(m)+(d)

	char mm[3];
	if(d.m < 10)
	{
		mm[0] = '0';
		mm[1] = d.m+'0';
		mm[2] = 0;
	}
	else
	{
		mm[0] = '1';
		mm[1] = d.m-10+'0';
		mm[2] = 0;
	}
	char dd[3];
	dd[0] = d.d / 10 +'0';
	dd[1] = d.d % 10 +'0';
	dd[2] = 0;

	stringstream yyyy;
	yyyy << d.y;

	return datadir()+"/"+ID(astock,amarket,yyyy.str(),mm,dd)+".bin";
}


void data::openfile(date d)
{
	if(file)
	{
		delete file;
		file = 0;
	}

	string s=getfilename(pair.stock, pair.market, d);
	file = new ifstream(s.c_str());
	if(!*file)
	{
		cerr << "Error opening " << s << endl;
		throw 1;
	}
	filedate = d;
}

bool data::firstunfilteredrecord(hfdrecord& r, date& d)
{
	if(!findfirstday(d))
		return false;
	openfile(d);
	bool foo, phoo, hoo;
	nextunfilteredrecord(r, d, foo, hoo, phoo);
	return true;
}


bool data::nextunfilteredrecord(hfdrecord& r, date& d, bool& adaybreak,
		bool& amonthbreak, bool& ayearbreak)
{
	adaybreak = false;
	for(;;)
	{
		file->read((char*) (&r),sizeof(r));
		if(file->eof() || r.t >= pair.closinghour * 3600.00)
		{
			if(!findnextday(d))
				return false;
			ayearbreak = filedate.y != d.y;
			amonthbreak = filedate.m != d.m;
			openfile(d);
			adaybreak = true;
		}
		else if(r.t >= pair.openinghour * 3600.00)
			break;
	}

    d = filedate;
	return true;
}

bool data::firstrecord(hfdrecord& r, date& d)
{
   hfdrecord s;
   date e;
   if(!firstunfilteredrecord(s,e))
	   return false;
   while(!satisfiesfilter(s))
   {
	   bool phoo;
	   if(!nextunfilteredrecord(s,e,phoo,phoo,phoo))
		   return false;
   }
   r=s;
   d=e;
   return true;
}

bool data::nextrecord(hfdrecord& r, date& d, bool& adaybreak,
		bool& amonthbreak, bool& ayearbreak)
{
   hfdrecord s;
   date e;
   bool daybreak=false;
   bool monthbreak = false;
   bool yearbreak = false;
   bool db,mb,yb;
   if(!nextunfilteredrecord(s,e,db,mb,yb))
	   return false;
   for(;;)
   {
	   if(db)
		   daybreak=true;
	   if(mb)
		   monthbreak=true;
	   if(yb)
		   yearbreak=true;
	   if(satisfiesfilter(s))
	   {
		   r=s;
		   d=e;
		   adaybreak=daybreak;
		   amonthbreak=monthbreak;
		   ayearbreak = yearbreak;
		   return true;
	   }
	   if(!nextunfilteredrecord(s,e,db,mb,yb))
		   return false;
   }
   return false; // anyway, should never get here
}

bool data::findfirstday(date& d)
{
	daypointer = 0;
	for(;;)
	{
		date rd(0,0,0);

		if(!findnextday(rd))
			return false;

		if(fd <= rd)
		{
			d = rd;
			return true;
		}
	}

	return true; // anyway never gets here
}

bool data::findnextday(date& d)
{
	if(daypointer == notpointing)
		throw;
	for(;;)
	{
		if(daypointer+1 == hfd::size())
			return false;
		hfd::dayrecord dr = hfd::at(daypointer++);

		if(pair.stock==dr.stock && pair.market==dr.market) // here we rely on ordering of query in \p hfd
		{
			if(dr.daydate <= ld)
			{
				d=dr.daydate;
				return true;
			}
			else
				return false;
		}
	}

	return true; // anyway never gets here
}

data::data(smpair& apair, unsigned int arecordfilter, date afd, date ald) :
	pair(apair),
	fd(afd), ld(ald),
	daypointer(notpointing),
	file(0),
	recordfilter(arecordfilter)
{
}


