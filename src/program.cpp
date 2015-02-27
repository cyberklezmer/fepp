/*
 * estimator.cpp
 *
 *  Created on: 10 Mar 2010
 *      Author: martin
 */

#include "hfd.hpp"

void program::process()
{
	anal.setprogram(this);

	anal.start();
	for(unsigned int g=0; g<sms.size(); g++)
	{
		cout << sms[g].stock << " "
			<< sms[g].market << endl;
		smpair& p = sms[g];
		try
		{
            anal.startpair(p);
            data dat(p,anal.getrecordfilter(),startdate, enddate);
            hfdrecord r;
            date d;
            bool noeof = true;
            int numrecs = 0;
            int numdays = 0;
            bool daybreak = false;
            bool yearbreak = false;
            bool monthbreak = false;
            bool firsttime = true;
            bool anyrecord = false;
            for(noeof = dat.firstrecord(r,d);
                     noeof;
                     noeof = dat.nextrecord(r,d,daybreak,monthbreak,yearbreak),numrecs++)
            {
                anyrecord=true;
                if(firsttime)
                {
                    cout << d.m << "/" << d.y << endl;

                    anal.startyear(d.y);
                    anal.startmonth(d.m);
                    anal.startday(d.d);
                    firsttime = false;
                }

                if(daybreak)
                {
                    anal.endday();
                    if(yearbreak)
                        anal.endyear();
                    if(monthbreak)
                        anal.endmonth();
                    if(yearbreak)
                        anal.startyear(d.y);
                    if(monthbreak)
                        anal.startmonth(d.m);
                    anal.startday(++numdays);
                }
                anal.record(r,numrecs++);
            }
            if(anyrecord)
            {
                anal.endday();
                anal.endmonth();
                anal.endyear();
                cout << d.m << "/" << d.y << endl;
            }
            anal.endpair();
        }
        catch(std::exception& e)
        {
            anal.endpair(true, e.what());
        }
        catch(...)
        {
            anal.endpair(true);
        }
	}
	anal.end();
	anal.setprogram(0);

}
