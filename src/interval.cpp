/*

 *
 *  Created on: May 22, 2013
 *      Author: martin
 */

#include "hfd.hpp"


void interval::onstartpair(smpair& togo)
{
	sumN = 0;
	totalrecs = 0;
	totalvol = 0;
	totaltrades=0;
	totalspread=0;
}

void interval::resetcounters()
{
	sumNa=sumNb=0;
	suma=sumb=0;
}

void interval::step(bool auselast)
{
	interval::rec r;

	r.day = getday();
	r.month = getmonth();
	r.year = getyear();
	r.time = lasttime;
	r.actA = lasta/100.0;
	r.actB = lastb/100.0;
	r.nfSnapA = lastnfa/100.0;
	r.nfSnapB = lastnfb/100.0;
	if(auselast)
	{
		r.avgA = lasta/100.0;
		r.avgB = lastb/100.0;
	}
	else
	{
		r.avgA = (double) suma / inter / 100.0;
		r.avgB = (double) sumb / inter / 100.0;
	}
	r.dNA = sumNa;
	r.dNB = sumNb;
	r.N = sumN;

	if(lasttime >= (getpair().openinghour + warming)*3600)
	{
#define ABS(x) (x>0?x:-x)
		onsnapshot(r);
		sumN += r.dNB-r.dNA;
		totalrecs++;
//cout << ABS(r.dNB)+ABS(r.dNA);
		totalvol += ABS(r.dNB)+ABS(r.dNA);
//cout << " " << totalvol / totalrecs << endl;
		totalspread +=r.actA-r.actB; // tohel dávalo průměr menší než jedna.
		if(r.actA-r.actB < 0.009999)
			cout << "Spread: " << r.actA-r.actB << endl;
	}

	lasttime += inter;

	resetcounters();
}

void interval::onrecord(hfdrecord& rec, int recn)
{
	if(lasttime == notime )
	{
		if(data::isqouterecord(rec))
		{
			if(rec.t > (getpair().openinghour + warming)*3600)
			{
				cerr << "No records during warming period " << getpair().stock << "/" << getpair().market << endl;
				throw 1;
			}

			lastqrectime = lasttime = (int) (rec.t/inter) * inter; // tady zaokrouhleno na interval.
			lastnfa = lasta;
			lastnfb = lastb;

			lasta = rec.a;
			lastb = rec.b;

			resetcounters();
		}
	}
	else
	{
		bool weareover = rec.t > lasttime + inter;
		if(!weareover)
		{
			switch(data::getrecordtype(rec))
			{
				case data::unpairedquoterecord:
				{
					double efftime = rec.t - lastqrectime;

					suma += lasta * efftime;
					sumb += lastb * efftime;

					lasta = rec.a;
					lastb = rec.b;

					lastqrectime = rec.t;

					break;
				}
				case data::traderecord:
					break;
				case data::pairedquoterecord:
				{
					totaltrades++;

					if(rec.q > 0)
						sumNb += rec.q;
					else
						sumNa += -rec.q;
					break;
				}
				default:
                    throw logic_error("unknown rectype");
			}
		}
		else // !weareover
		{
			double efftime = lasttime + inter - lastqrectime;

			suma += lasta * efftime;
			sumb += lasta * efftime;

			for(;;)
			{
				step(false);

				if(lasttime + inter > rec.t)
				{
					double overtime = rec.t - lasttime;
					suma = overtime * lasta;
					sumb = overtime * lastb;
					lastqrectime = rec.t;
					if(data::getrecordtype(rec) ==  data::unpairedquoterecord)
					{
						lasta = rec.a;
						lastb = rec.b;
					}
					break;
				}
				else
				{
					suma = lasta * inter;
					sumb = lastb * inter;
				}
			}
		} // !weareover

		if(data::isqouterecord(rec))
		{
			lastnfa = rec.a;
			lastnfb = rec.b;
		}

//		if(data::getrecordtype(rec) == data::unpairedquoterecord)
		//		{
		//			lastprice = (double) (rec.a + rec.b) / 2.0;
		//			lastsigma = (double) (rec.a - rec.b);
		//}
	}
}

void interval::onendday(int day)
{
	if(lasttime != notime )
	{
		for(;lasttime <= getpair().closinghour * 3600;)
			step(true);
	}

}

