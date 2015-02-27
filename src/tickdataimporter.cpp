/*
 * tickdataimporter.cpp
 *
 *  Created on: Mar 13, 2013
 *      Author: martin
 */



#include <sstream>
#include "hfd.hpp"



tickdataimporter::marketinfo tickdataimporter::kmarkets[] =
{
		{ 'A',"AMEX" },
		{ 'B',"NASDAQ OMX BX" }, // (Boston)
		{ 'C',"NSE" },//National_Stock Exchange (Cincinnati)
		{ 'D',"NASD ADF"}, // (FINRA)
		{'E',"MI"},//Market Independent (SIP - Generated)
		{'I',"ISE"}, //  (International Securities Exchange)
		{'J',"DE A"}, //(DirectEdge A)
		{'K',"DE X"},//(DirectEdge X)
		{'M',"Chicago"},
		{'N',"NYSE"},
//		{'O',"Instinet (Valid only during January and February 1993)"},
		{'P',"ARCA"}, // (formerly Pacific)
		{'S',"CTS"}, //Consolidated Tape System
		{'T',"NASDAQ T"},
		{'Q',"NASDAQ Q"},
		{'W',"CBOE"},
		{'X',"NASDAQ OMX"}, // (Philadelphia)
		{'Y',"BATS Y"},
		{'Z',"BATS" },
		{0,""}
};

void tickdataimporter::import()
{
	vector<string> stocks;
	for(unsigned int i=0; i< pairs.size(); i++)
	{
		const string& stock = pairs[i].stock;
		bool found = false;
		for(unsigned int j=0; j<stocks.size(); j++)
			if(stocks[j]==stock)
			{
				found = true;
				break;
			}
		if(!found)
			stocks.push_back(pairs[i].stock);
	}
	for(unsigned int i=0; i<stocks.size(); i++)
	{
		vector<char> marketids;
		vector<string> markets;
		string& stock = stocks[i];
		for(unsigned int j=0; j<pairs.size(); j++)
		{
			if(pairs[j].stock == stock)
			{
				string& market = pairs[j].market;
				bool found = false;
				for(int k=0; kmarkets[k].id; k++)
				{
					if(kmarkets[k].name==market)
					{
						marketids.push_back(kmarkets[k].id);
						markets.push_back(kmarkets[k].name);
						found = true;
						break;
					}
				}
				if(!found)
				{
					cerr << "Cannot find an id of market " << market << endl;
					throw 1;
				}
			}
		}
		importstock(stock,markets,marketids);
	}
}

void tickdataimporter::importusfile(string& astock, vector<string>& markets,
		vector<char>& marketids,
		string qfn, string tfn)
{
 	ifstream istr(qfn.c_str());
	if(!istr)
	{
	    cerr << "Error opening " << qfn << endl;
		throw 1;
	}

 	ifstream tstr(tfn.c_str());
	if(!tstr)
	{
	    cerr << "Error opening " << tfn << endl;
		throw 1;
	}

	char tradeline[200];
	tstr.getline(tradeline,sizeof(tradeline)-1);

	date olddate(0,0,0);


	vector<hfdrecord> recs[maxmarketsperstock];

	double numzeroquotes[maxmarketsperstock];

	int numunassignedtrades[maxmarketsperstock];

	int nummarkets = marketids.size();

	for(;;)
	{
		bool eof = istr.eof();
		bool our = false;
		char datestr[12]; // jako za starejch časů
		char time[16];
		double t;
		char exchange[4];
		int m;
		date d;
		if(!eof)
		{
			istr.getline(datestr,11,',');
			datestr[2]=datestr[5]=0;

			d = date(atoi(datestr+6),atoi(datestr), atoi(datestr+3));
			istr.getline(time,15,',');
			time[2]=time[5]=0;
			t = atof(time+6)+60.0*atoi(time+3)+3600.0*atoi(time);


			istr.getline(exchange,4,',');

			our = false;
			for(m=0; m< nummarkets; m++)
				if(exchange[0]==marketids[m])
				{
					our = true;
					break;
				}
		}
		if(our || eof)
		{
			if(!(olddate == d) || eof)
			{
				date zerodate(0,0,0);
			   if(!(olddate==zerodate)) // tj máme za sebou sbírání kotací
			   {
				   ofstream* oss[maxmarketsperstock];
				   for(int i=0; i<nummarkets;i++)
					   oss[i] = data::createdayfile(astock,markets[i],olddate);

				   // mejprve spárovat s tradama
				   // v tradeline by měl bejt schovanej prvn řádek, kterej se nás týká



				   vector<tinfo> ti[maxmarketsperstock];

				   bool teof = tstr.eof(); // tj jsme u konce soubur s trady
				   bool recordtooearly = true;
				   bool writetrecord = false;

				   for(;!teof;)
				   {
					   double tt = 0; // čas dalšího "tradu" do kterého musíme "dojet" quoty
					   int p,q, tm;

					   // načíst data do stringu a parsovat
					   tradeline[2]=tradeline[5]=tradeline[10]=0;
					   date tdate(atoi(tradeline+6),atoi(tradeline), atoi(tradeline+3));
					   char *timestr= tradeline + 11;
					   timestr[2]=timestr[5]=timestr[12]=0;
					   tt = atof(timestr+6)+60.0*atoi(timestr+3)+3600.0*atoi(timestr);

					   if(recordtooearly)
						   recordtooearly = (tdate<=olddate && !(tdate==olddate));
					   if(!(tdate==olddate))
						   break; // tohle už je z dalšího dne
					   else if(!recordtooearly)
					   {
						   char* pricestr = timestr+13;
						   char* volumestr = timestr+14;
						   for(;*volumestr && *volumestr !=','; volumestr++)
							   ;
						   if(!*volumestr)
						   {
							   cerr << "Wrong format of trade data" << endl;
							   throw 1;
						   }
						   volumestr++;
						   char* exchange = 0;
						   for(exchange = volumestr+1; *exchange && *exchange !=','; exchange++)
							   ;
						   if(!*exchange)
						   {
							   cerr << "Wrong format of trade data - exchange" << endl;
							   throw 1;
						   }
						   exchange++;

						   tm=0;
						   writetrecord = false;
						   for(tm=0; tm<nummarkets;tm++)
							   if(exchange[0]==marketids[tm])
							   {
								   writetrecord = true;
								   break;
							   }
						   if(writetrecord)
						   {
							   p = atof(pricestr)*knumticks + 0.5;
							   q = atoi(volumestr) / klotsize;

							   tinfo f = {tt,p,q,-1,0};
// zakomentovaný kod slučoval trady se stejným časem
//							   int s = ti[tm].size();
//							   if(s > 0 && ti[tm][s-1].p == p && ti[tm][s-1].tt == tt)
//								   ti[tm][s-1].q+=f.q;
//							   else
								   ti[tm].push_back(f);
						   }
					   }
					   tstr.getline(tradeline,sizeof(tradeline)-1);
					   teof = tstr.eof();
				   }
				   for(int i=0;i < nummarkets; i++)
				   {
					   int qvolume = 0;
					   int unmatchedqvol =0;

					   int pq = 0;
						#define ABS(Q) (Q > 0 ? Q: -Q)
					   for(unsigned int k=0; k< ti[i].size(); k++)
					   {
						   int lastma = -1;
						   int lastmb = -1;
						   tinfo& f = ti[i][k];

//cout << recs[i].size() << endl;
//cout << recs[i][0].t << endl;
						   for(int j=pq;j+1 < ((signed int) recs[i].size())-1 && recs[i][j].t < f.tt; j++)
						   {
							   hfdrecord& r = recs[i][j];
							   if(f.p==r.a)
								   lastma = j;
							   if(f.p==r.b)
								   lastmb = j;
						   }
						   qvolume += ABS(f.q);

						   if(lastma == -1 && lastmb == -1)
						   {
							   numunassignedtrades[i]++;
							   unmatchedqvol += ABS(f.q);
						   }
						   else
						   {
							  bool isb = lastma < lastmb;
							  f.type = isb ? 1 : -1;
							  int lastqoute = isb ? lastmb : lastma;

							  bool matchfound = false;
							  int closestpartialmatch = -1;
							  int j=lastqoute;

							  for(;j>=pq && recs[i][j].quote(isb) == f.p;j--)
							  {
								  int nextvolume = recs[i][j+1].quote(isb) == f.p ? recs[i][j+1].volume(isb) : 0;
								  int thisvolume = recs[i][j].volume(isb);
								  int d = (thisvolume-nextvolume) - ABS(recs[i][j+1].q);
								  if(f.q == d)
								  {
									  recs[i][j+1].q += isb ? f.q : -f.q;
									  recs[i][j+1].p = f.p;
									  recs[i][j+1].qt = f.tt;
									  f.pairedwith = j+1;
									  pq = j;
									  matchfound = true;
									  break;
								  }
								  if(f.q < d)
								  {
									  if(closestpartialmatch < 0)
										  closestpartialmatch = j+1;
								  }
							  }
							  if(!matchfound && closestpartialmatch >= 0)
							  {
								  recs[i][closestpartialmatch].q+=isb ? f.q : -f.q;
								  recs[i][closestpartialmatch].qt = f.tt;
								  recs[i][closestpartialmatch].p = f.p;
								  f.pairedwith = closestpartialmatch;
								  pq = closestpartialmatch;
							  }
						   }
					   }

					   unsigned int uq=0;
					   for(unsigned int j=0; j<recs[i].size(); j++)
					   {
						   hfdrecord& r = recs[i][j];

						   while(uq < ti[i].size() && ti[i][uq].tt < r.t)
						   {
							   tinfo& f=ti[i][uq];
							   if(f.pairedwith == -1 && f.type != 0)
							   {
								   hfdrecord r(f.tt,f.type * f.q, f.p,f.tt);
								   data::writerecord(*oss[i],r);
							   }
							   uq++;
						   }
						   data::writerecord(*oss[i],r);
					   }
					   for(;uq < ti[i].size(); uq++)  // vyflushneme
					   {
						   tinfo& f=ti[i][uq];
						   if(f.pairedwith == -1  && f.type != 0)
						   {
							   hfdrecord r(r.t,f.type ? f.q : -f.q,f.p,f.tt);
							   data::writerecord(*oss[i],r);
						   }
					   }

					   delete oss[i];

					   ostringstream remark;
					   remark << unmatchedqvol << "/" << qvolume << " unmatched";
				       if(numzeroquotes[i] > 0)
						  remark << ", "<< numzeroquotes[i] << " zero quotes.";
				       if(recs[i].size() >0 || ti[i].size() > 0)
				       {
				    	   cout << "Adding " << astock << "," << markets[i] << endl;
				    	   addday(astock,markets[i],olddate,hfd::tradequotetick,remark.str());
				    	   recs[i].resize(0);
				    	   ti[i].resize(0);
				       }
				       else
				       {
				    	   cout << astock << " at " << markets[i] << " found empty!" << endl;
				       }
				   }
			   }
			   if(eof)
				   break;

	           olddate = d;
			   for(int i=0; i<nummarkets; i++)
			   {
				   numzeroquotes[i]=0;
				   numunassignedtrades[i] = 0;
			   }
			}

			bool skipthis = false;

			int b;
			int a;
			int an;
			int bn;
			char bid[20];
			istr.getline(bid,20,',');
			char ask[20];
			istr.getline(ask,20,',');
			char bidn[20];
			istr.getline(bidn,20,',');
			char askn[20];
			istr.getline(askn,20,',');
			b = atof(bid)*knumticks + 0.5;
			a = atof(ask)*knumticks + 0.5;
			bn = atoi(bidn);
			an = atoi(askn);

			if(b==0 || a==0)
			{
//				cout << a << " " << b << endl;
				numzeroquotes[m]++;
			}
			else if(/* lastt[m] && */!skipthis && b < a)
			{
				hfdrecord r(t,b,a,bn,an);
				recs[m].push_back(r);
			}
		}

		char temp[100];
		istr.getline(temp,99);
	  }
}

