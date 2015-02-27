/*
 * hfd.cpp
 *
 *  Created on: Feb 18, 2013
 *      Author: martin
 */




#include "hfd.hpp"

hfd* hfd::self;

#define INDEXFN "index.csv"


string hfd::indexfn()
{
	 return hfddir() + "/" + INDEXFN;
}

hfd::hfd(const string& ahome) :home(ahome)
{
    self=this;

	string s;
	s = indexfn();

	ifstream ind(s.c_str());
	if(!ind)
	{
		cerr << "Cannot open " << s << " (on beginning, zero length one should be created)" << endl;
		throw 1;
	}
	for(;;)
	{
		int y,m,d;
		string stock;
		string market;
		int dt;

		char line[200];
		ind.getline(line,sizeof(line)-1);
		if(ind.eof())
			break;

		char* ptr=line;

		char* yc = ptr;
		for(;*ptr && *ptr!=',' ;ptr++)
			;
		if(!*ptr)
		{
			cerr << "Comma expected after year" << endl;
			throw 1;
		}
		*ptr = 0;
		y = atoi(yc);

		char* mc = ++ptr;
		for(;*ptr && *ptr!=',' ;ptr++)
			;
		if(!*ptr)
		{
			cerr << "Comma expected after month" << endl;
			throw 1;
		}
		*ptr = 0;
		m = atoi(mc);

		char* dc = ++ptr;
		for(;*ptr && *ptr!=',' ;ptr++)
			;
		if(!*ptr)
		{
			cerr << "Comma expected after day" << endl;
			throw 1;
		}
		*ptr = 0;
		d = atoi(dc);

		char* sc = ++ptr;
		for(;*ptr && *ptr!=',' ;ptr++)
			;
		if(!*ptr)
		{
			cerr << "Comma expected after stock" << endl;
			throw 1;
		}
		*ptr = 0;
		stock = sc;

		char* mac = ++ptr;
		for(;*ptr && *ptr!=',' ;ptr++)
			;
		if(!*ptr)
		{
			cerr << "Comma expected after market" << endl;
			throw 1;
		}
		*ptr = 0;
		market = mac;

		char* dtc = ++ptr;
		dt = atoi(dtc);

		date dat(y,m,d);
		dayrecord dr;
		dr.daydate =dat;
		dr.stock = stock;
		dr.market = market;
		dr.type = (datatype) dt;
		if(days.size() > 0 && compare(dr,days[size()-1]) <= 0)
		{
			cerr << "Records in " << indexfn() << " unsorted or duplicated!" << endl;
			throw 1;
		}
		days.push_back(dr);
	}
	cout << days.size() << " records read from " << INDEXFN << endl;
}


int hfd::compare(dayrecord& r1, dayrecord& r2)
{
	if(r1.daydate == r2.daydate )
	{
		if(r1.stock == r2.stock )
		{
			if(r1.market == r2.market)
			{
				if(r1.type == r2.type)
					return 0;
				else if(r1.type < r2.type)
					return -1;
				else
					return 1;
			}
			else if(r1.market < r2.market)
				return -1;
			else
				return 1;
		}
		else if(r1.stock < r2.stock)
			return -1;
		else
			return 1;
	}
	else if(r1.daydate <= r2.daydate )
		return -1;
	else
		return 1;
}

bool hfd::dayintovec(dayrecord& dr)
{
	unsigned int i;
	for(i=0; i < days.size() && compare(dr,days[i])>0; i++)
		;
	if(i < days.size() && compare(dr,days[i])==0)
		return false; // we have found match
	days.resize(days.size()+1);
	for(unsigned int j=days.size()-1;j>i;j--)
		days[j]=days[j-1];
	days[i]=dr;
	return true;
}

bool hfd::addday(const string& stock, const string& market,const date& daydate,
		         datatype type, const string& remark )
{
	dayrecord dr;
	dr.daydate = daydate;
	dr.stock = stock;
	dr.market = market;
	dr.type = type;
//	dr.remark = remark;

	return self->dayintovec(dr);
}

void hfd::writeindex()
{
    ofstream ind(indexfn().c_str());
    if(!ind)
    {
    	cerr << "Cannot create " << indexfn() << endl;
    	throw 1;
    }
    for(unsigned int i=0; i<self->days.size(); i++)
    {
    	dayrecord& dr = self->days[i];
		ind << dr.daydate.y << ",";
		ind << dr.daydate.m << ",";
		ind << dr.daydate.d<< ",";
		ind << dr.stock << ",";
		ind << dr.market << ",";
		ind << dr.type << ",";
		ind << dr.remark;
		ind << endl;
    }
}


