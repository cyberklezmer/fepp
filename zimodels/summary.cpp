/*
 * summary.cpp
 *
 *  Created on: 22 Apr 2010
 *      Author: martin
 */

#include "zimodels.hpp"


void zisummary::onstartpair(smpair& togo)
{
	for(int j=0;j<kmaxrecjmp; j++)
        jumphist[j]=jumpdhist[j]=reldjumps[j]=0;

	for(int j=0;j<kmaxlogqs; j++)
        logqs[j] = 0;

	for(int i=0; i<numyears(); i++)
	{

		for(int j=0; j<12; j++)
		{
			int k= i*12 + j;
			mspreads[k] = mnumbers[k] = 0;
		}
	}
}



void zisummary::onstartmonth(int day)
{
	msumspreads = 0;
	mnumtrans = 0;
	numberdaysinmonth = 0;
	numrecsinmonth = 0;
	msumdeltat = 0;
}

void zisummary::onendmonth(int month)
{
	if(numberdaysinmonth > 0)
	{
		int i = (getyear()-starty)*12+month-1;
		mspreads[i] = msumspreads / numrecsinmonth;
		mnumbers[i] = (double) mnumtrans / numberdaysinmonth;
		mtimes[i] = msumdeltat / mnumtrans;
	}
}

void zisummary::onrecord(hfdrecord& rec, int recn)
{
	static double savet;
	if(getqhistorysize() == 0)
		savet = rec.t;
	else
	{
		hfdrecord h = getqhistory(0);
		if(rec.a != h.a || rec.b != h.b)
		{
			msumspreads+= rec.a-rec.b;
			msumdeltat += rec.t - savet;
			numrecsinmonth++;
			mnumtrans++;
			int da = rec.a - h.a;
			if(da < 0)
			{
                int freespace = (h.a - h.b - 1);
                if(freespace > 0 && rec.b == h.b)
                {
                   double rj = (double) da / (double) (h.a - h.b - 1);
                   reldjumps[(int)(20*(1.0 + rj))]++;
                }

			    if(da < -kmaxrecjmp)
			        da = -kmaxrecjmp;
                jumpdhist[kmaxrecjmp+da]++;
			}
			else if(da > 0)
			{
                if(da > kmaxrecjmp)
			        da = kmaxrecjmp;
                jumphist[da-1]++;
            }
            if(data::getrecordtype(rec) == data::pairedquoterecord
               && rec.q < 0)
            {
                unsigned int l = (unsigned int) log(-rec.q);
                if(l>=kmaxlogqs)
                    l = kmaxlogqs-1;
                logqs[l]++;
            }


/*			int db = rec.b - h.b;
			if(db < 0)
                db = -db;
			if(db > kmaxrecjmp)
			   db = kmaxrecjmp;
			if(db)
                jumphist[db-1]++;*/
//cout << da << " " << db << endl;
			savet = rec.t;
		}
	}
}

void zisummary::onendday(int day)
{
	numberdaysinmonth++;
}


void zisummary::onendpair(smpair& togo)
{
	latex << "{\\normalsize " << togo.stock
		<< " on " << togo.market << "}"
		<< endl << endl;

	latex << "\\begin{tabular}{ccccc}" << endl
			<< "$m/y$ & $\\bar s$  & $\\overline {\\Delta t}$ & $\\#$ \\\\" << endl;
	latex.precision(2);
    for(int i = 0; i<numyears(); i++)
    {
        for(int j = 0; j<12; j++)
        {
            int k = i*12+j;

            if(mnumbers[k])
            {
                latex <<
                        j+1 << "/" << starty + i << " & " << endl <<
                        mspreads[k] << " & " <<
                        mtimes[k] << " & " <<
                        mnumbers[k] / 10000 ;
                latex << "\\\\" << endl;
            }
        }
    }
	latex << endl << "\\end{tabular}" << endl;

    int sumjmps=0;
    int sumdjmps=0;
    int sumrjmps=0;
    for(int i=0; i<kmaxrecjmp; i++)
    {
        sumjmps += jumphist[i];
        sumdjmps += jumpdhist[i];
        sumrjmps += reldjumps[i];
    }

    double sumlogq = 0;
    for(int i=0; i<kmaxlogqs; i++)
        sumlogq += logqs[i];

    double histj[kmaxrecjmp];

    for(int S=/*-1*/1; S<=1; S++)
    {
        int s = S==-1 ? sumdjmps : sumjmps;
        switch(S)
        {
            case -1: s=sumdjmps; break;
            case 0: s=sumrjmps; break;
            case 1: s=sumdjmps; break;
        }
        if(s)
        {
            for(int i=0; i<kmaxrecjmp; i++)
            {
                histj[i] =
                  S < 0 ?
                  (double) jumpdhist[i] / (double) s :
                  (S==0 ? (double) reldjumps[i] / (double) s :
                  (double) jumphist[i] / (double) s);
            }

            string pfn = stockfn() + (S==0 ? "r" : (S<0 ? "d" : "u"));
            gnuplot plot(pfn);

            ofstream ofn(plot.datfn().c_str());
            if(!ofn)
            {
                cout << "Error creating " << plot.datfn() << endl;
                throw 1;
            }
            for(int i=0; i<kmaxrecjmp; i++)
            {
                if(i+1 == kmaxrecjmp)
                    ofn << "\">" << i;
                else
                    ofn << '"' << i+1;
                ofn << "\", " << histj[i] << endl;
            }

            ofn.flush();

            plot.script() << "set yrange [0:1]" << endl;
            plot.script() << "set boxwidth 1 relative" << endl;
            plot.script() << "set style data histograms" << endl;
            plot.script() << "set style fill solid 1.0 border -1" << endl;
            plot.script() << "plot '" << plot.datfn() << "' using 2:xticlabels(1) notitle" << endl;

            plot.process();

            latex << "\\vspace{-0.2cm}" << endl;
            latex << "\\begin{center}" << endl;
            latex << endl << "\\includegraphics[width=2cm, angle = 270]{"
                    << pfn
                    << "}" << endl;
            latex << "\\end{center}  " << endl;
            latex << "\\vspace{-0.2cm}" << endl;
        }
    }
    double hist[kmaxlogqs];
    for(int i=0; i<kmaxrecjmp; i++)
        hist[i] = sumlogq ? (double) logqs[i] / sumlogq : 0;

    string pfn = stockfn() + "lq";
    gnuplot plot(pfn);

    ofstream ofn(plot.datfn().c_str());
    if(!ofn)
    {
        cout << "Error creating " << plot.datfn() << endl;
        throw 1;
    }
    for(int i=0; i<kmaxlogqs; i++)
    {
        if(i+1 == kmaxlogqs)
            ofn << "\">" << i-1;
        else
            ofn << '"' << i;
        ofn << "\", " << hist[i] << endl;
    }

    ofn.flush();

    plot.script() << "set yrange [0:1]" << endl;
    plot.script() << "set boxwidth 1 relative" << endl;
    plot.script() << "set style data histograms" << endl;
    plot.script() << "set style fill solid 1.0 border -1" << endl;
    plot.script() << "plot '" << plot.datfn() << "' using 2:xticlabels(1) notitle" << endl;

    plot.process();

    latex << "\\vspace{-0.2cm}" << endl;
    latex << "\\begin{center}" << endl;
    latex << endl << "\\includegraphics[width=2cm, angle = 270]{"
            << pfn
            << "}" << endl;
    latex << "\\end{center}  " << endl;
    latex << "\\vspace{-0.2cm}" << endl;

}

