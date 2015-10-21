/*
 *      Author: martin
 */

#ifndef HFD_HPP_
#define HFD_HPP_

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <stdexcept>


#include <assert.h>
#include <stdlib.h>
#include <vector>

using namespace std;

/// \version 0.5

/** \mainpage

\section Prerequisities

7z installed in path

\p hfd folder present elsewhere

\p data, \p latex, \p inter/gnuplot

\p inter/csv folders present in the hfd folder

file index.csv (of zero length in case of first time import) present in the hfd folder
 7
tbd gnuplot gretl 7z

*/


class date
{
public:
	int y; int m; int d;
	date(int ay, int am, int ad) : y(ay), m(am), d(ad) {}
	date() : y(0), m(0), d(0) {}
	bool operator <=( date& d2)
	{
		if(this->y < d2.y)
			return true;
		else if(this->y > d2.y)
			return false;
		else if(this->m < d2.m)
			return true;
		else if(this->m > d2.m)
			return false;
		else if(this->d > d2.d)
			return false;
		else
			return true;
	}
	bool operator ==( date& d2)
	{
        return d==d2.d && m==d2.m && y ==d2.y;
	}
 };


struct smpair
{
	string stock;
	string market;
	double openinghour;
	double closinghour;
    void* additional;
	smpair(string as, string am, double aopening=0,
	       double aclosing=24, void* aadditional = 0)
	 : stock(as), market(am), openinghour(aopening),
	     closinghour(aclosing), additional(aadditional)
	     {}
	smpair(): openinghour(0), closinghour(24), additional(0)
	     {}
};
// todo brát opening a closing v úvahu v datech


const string khfd = "/hfd";

class hfd
{
	static hfd* self;
	friend class data;
	friend class importer;

	string home;


public:
	enum datatype { tradequotetick=1 };
private:
	struct dayrecord
	{
		date daydate;
		string stock;
		string market;
		datatype type;
		string remark;
	};
	vector<dayrecord> days;
	static dayrecord& at(int i) { return self->days[i]; }
	static int size() { return self->days.size(); }
	static int compare(dayrecord& r1, dayrecord& r2);
	bool dayintovec(dayrecord& dr);
	static void writeindex();
	static bool addday(const string& stock, const string& market,const date& daydate,
			            datatype type, const string& remark);
public:
	static string gethome() { return self->home; }
	static string hfddir() { return gethome()+khfd; }
	static string interdir() { return hfddir()+"/inter"; }
	static string tempdir() { return hfddir()+"/temp"; }
	static string indexfn();
	hfd(const string& ahome);
	~hfd() { self=0; }
};



struct hfdrecord
{
	/// time since start of the day
	double t;
    /// bid (in cents)
	int b;
	/// ask (in cents)
	int a;
	/// volume at bid (in pieces)
	int bn;
	/// volume at ask (in pieces)
	int an;
	/// volume of a trade (if the record corresponds to a trade), negative if it is a purchase (reaction to a sell order)
	int q;
	/// price of a trade
	int p;
	/// time of the trade
	double qt;
	hfdrecord(): t(0), b(0),a(0),bn(0),an(0),q(0), p(0), qt(0) {}
	hfdrecord(double at, int ab, int aa, int abn, int aan):
		t(at), b(ab),a(aa),bn(abn),an(aan),q(0),p(0),qt(0) {}
	hfdrecord(double at, int aq, int ap, double aqt):
		t(at), b(0),a(0),bn(0),an(0),q(aq),p(ap), qt(aqt) {}
	int& quote(bool isb) { return isb ? b : a; }
	int& volume(bool isb) { return isb ? bn : an; }
};


class data
{
private:

	smpair pair;
	static string datadir() { return hfd::hfddir()+ "/data"; }
	static string dayspdir() { return hfd::hfddir()+"/days"; }

	date fd;
	date ld;

	int daypointer;
	ifstream* file;
	date filedate;

	void openfile(date d);

	static const int notpointing = -1;
	bool findfirstday(date& d);
	bool findnextday(date& d);
public:
	static ofstream* createdayfile(string& astock, string& amarket, date d)
	{
	    string ofn = getfilename(astock,amarket,d);
	    ofstream* result = new ofstream(ofn.c_str());
		if(!result || !*result)
		{
		    cerr << "Error creating " << ofn << endl;
			throw 1;
		}
		return result;
	}
	static void writerecord(ofstream& ofs, hfdrecord& rec);
protected:
	static string getfilename(string& astock, string& amarket, date d);

public:
	enum recordtype { traderecord, unpairedquoterecord, pairedquoterecord, numrecortypes };
	enum recordmask { filtertrade = 2<<traderecord, filterunpairedquote = 2<<unpairedquoterecord,
	                  filterpairedquote = 2 << pairedquoterecord,
	                  filterquote = filterpairedquote | filterunpairedquote  };
protected:
	bool satisfiesfilter(hfdrecord& rec)
	{
		recordtype t=getrecordtype(rec);
		return (2<<t) & recordfilter;
	}
private:
	unsigned int recordfilter;
	bool firstunfilteredrecord(hfdrecord& r, date& d);
	bool nextunfilteredrecord(hfdrecord& r, date& d,
			          bool& adaybreak, bool& amonthbreak, bool& ayearbreak);

public:
	static recordtype getrecordtype(hfdrecord& rec)
	{
		if(rec.q != 0)
		{
			return rec.a != 0 ? pairedquoterecord : traderecord;
		}
		else
			return unpairedquoterecord;
	}
	static bool isqouterecord(hfdrecord& rec)
	{
		recordtype rt = getrecordtype(rec);
		return rt == pairedquoterecord || rt == unpairedquoterecord;
	}
	/// \p d is an output parameter
	bool firstrecord(hfdrecord& r, date& d);
	bool nextrecord(hfdrecord& r, date& d,
			          bool& adaybreak, bool& amonthbreak, bool& ayearbreak);
	data(smpair& apair, unsigned int arecordfilter, date afd, date ald);
	~data()
	{
		delete file;
	}
};



typedef vector<smpair> smpairs;

class program;

class analysis
{
private:

	const program* pgm;
	void setprogram(program* p) { pgm=p; }
	friend class program;

	// temporary "state" variables
	int year;
	int month;
	int day;
	smpair pair;

	void start();
protected:
	void startpair(smpair& apair);
	void startyear(int ayear)
	{ onstartyear(year=ayear); };
	void startday(int aday)
	{
		history.clear();
		qhistory.clear();
		onstartday(day=aday);
	};
	void startmonth(int amonth)
	{ onstartmonth(month=amonth); };

	void record(hfdrecord& rec, int recn);

	void endpair(bool aexeptionthown = false, const char *aexctext = 0 );
	void endday() { onendday(day); };
	void endmonth() { onendmonth(month); };
	void endyear() { onendyear(year); };
	void end();
protected:
	string name;
	ofstream latex;
	ofstream csv;
	unsigned int recordfilter;
	int numcols;
	int numrows;
	int getyear() { return year; }
	int getday() { return day; } // poradove cislo v roce ?
	int getmonth() { return month; }
	const smpair& getpair() { return pair; }
	string gettexfilename(const string& aname)
	{
		return latexdir() + "/" + aname + ".tex";
	}
	string getcsvfilename(const string& aname)
	{
		return csvdir() + "/" + aname + ".csv";
	}
	string pairid()
	{
		const string und = "_";
		string res;
		for(unsigned int i=0;i<pair.stock.length(); i++ )
		{
			char c=pair.stock[i];
			res += c==' ' ? '_' : c;
		}
		res+= und;
		for(unsigned int i=0;i<pair.market.length(); i++ )
		{
			char c = pair.market[i];
			res += c==' ' ? '_' : c;
		}
		return res;
	}

private:
	int column;
	int row;

	vector<hfdrecord> history;
	vector<hfdrecord> qhistory;
	void finishtable(bool anewpage);
protected:
	hfdrecord& gethistory(unsigned int offset)
	{
		assert(offset < history.size());
		return history[history.size()-1-offset];
	}
//	int numinhistory;
	int gethistorysize()
	{
		return history.size();
	}
	hfdrecord& getqhistory(unsigned int offset)
	{
		assert(offset < qhistory.size());
		return qhistory[qhistory.size()-1-offset];
	}
	int getqhistorysize()
	{
		return qhistory.size();
	}
	const program* getprogram() { return pgm; }
public:
	unsigned int getrecordfilter() { return recordfilter; }
	analysis(const string& aname, unsigned int arecordfilter, int anumcols, int anumrows=1): //, bool affafterline = false):
			name(aname),
			latex(gettexfilename(aname).c_str()),
			csv(getcsvfilename(aname).c_str()),
			recordfilter(arecordfilter),
			numcols(anumcols),
			numrows(anumrows)
		{
		}
	virtual ~analysis() {}
	virtual void onstart() {}
	virtual void onstartpair(smpair& togo) {}
	virtual void onstartyear(int year) {}
	virtual void onstartmonth(int month) {}
	virtual void onstartday(int day) {}
	virtual void onrecord(hfdrecord& rec, int recn) {}
	virtual void onendday(int day) {};
	virtual void onendmonth(int month) {};
	virtual void onendyear(int year) {};
	virtual void onendpair(smpair& togo) {}
	virtual void onend() {}

    static string latexdir()
    {
        return hfd::hfddir()+"/latex";
    }
    static string csvdir()
    {
        return hfd::interdir()+"/csv";
    }
};



class program
{
private:
	smpairs& sms;
	analysis& anal;
	date startdate;
	date enddate;
public:
	program(smpairs& g, analysis& a, date& astartdate, date& aenddate)
	  :  sms(g), anal(a), startdate(astartdate), enddate(aenddate)
	{
	}
	void process();
};


const string kgnuplotpath = "gnuplot";

class gnuplot
{
	string name;
	string fn;
	ofstream sc;
public:
	ofstream& script() { return sc; }
	const string datfn() { return fn+".dat"; }
	const string scriptfn() { return fn+".sc"; }
	const string getfn() { return fn; }
	gnuplot(const string& aname) :
	  name(aname),
	  fn(hfd::interdir() + "/gnuplot/" + aname),
	  sc(scriptfn().c_str())
      {
		 if(!sc) //fixme - should not be done in constructor
		 {
			 cout << "Error opening script " << scriptfn() << endl;
			 throw 1;
		 }
		 sc << "set terminal postscript landscape solid " << endl
				<< "set output '" << analysis::latexdir() + "/" + aname + ".eps" << "'" << endl;
	  }
	void process()
	{
		sc.flush();

		string cmd = kgnuplotpath + " " + scriptfn();
		if(system(cmd.c_str()))
		{
			cerr << "Error executing " << cmd << endl;
			throw 1;
		}
	}

};


class gretl
{
	static string tmpfn() { return hfd::tempdir()+"/_gretl_tmp"; }
	string scriptresult;
public:
	static string scdir()
	{
		return hfd::hfddir() + "/gretl";
	}
	void runscript(const string& scfn)
	{
		string cmd = "gretlcli -b \"" + scfn + "\" > " + tmpfn();
		cout << cmd << endl;
		if(system(cmd.c_str()))
		{
			cerr << "Failed to run " << cmd << endl;
			throw 1;
		}
		ifstream p(tmpfn().c_str());
		if(!p)
		{
			cerr << "Error opening " << tmpfn() << endl;
			throw 1;
		}

		std::getline(p,scriptresult,'\0');
	}
	string getscriptoutput()
	{
		ifstream s(tmpfn().c_str());
		if(!s)
		{
			cerr << "Error opening " << s << endl;
			throw 1;
		}
		string str;
		std::getline(s,str);
		return str;
	}
	double findpar(const string& parname) // assumes that "? parname" appears in the latex
	{
		string label = "? " + parname;
		unsigned int pos = 0;
		for(;;)
		{
            pos = scriptresult.find(label, pos);
			if(pos==string::npos)
			{
				cerr << "Cannot find parameter " << parname << " in gretl script output " << endl;
				throw 1;
			}
			pos += label.length();
			if(scriptresult.c_str()[pos++] == '\n')
				break;
		}
		return atof(scriptresult.c_str()+pos);
	}
};


class importer
{
protected:
	date start;
	date end;
	smpairs& pairs;
	importer(date astart, date aend, smpairs& apairs) : start(astart), end(aend), pairs(apairs) {}
	virtual ~importer() {};
	virtual void import() = 0;
	static bool addday(const string& stock, const string& market,const date& daydate, hfd::datatype type,
						const string& remark)
	{
		return hfd::addday(stock,market,daydate,type,remark);
	}
public:
	void process()
	{
		import();
		hfd::writeindex();
	}
};

class tickdataimporter : public importer
{
public:
	enum datatype { usstock };
	struct marketinfo
	{
	    char id;
		const char* name;
	};
private:
	static marketinfo kmarkets[];
public:
	static const marketinfo& getmarketinfo(int i) // the last one had id=0
	{
	    return kmarkets[i];
	};
public:
	int numticks;
	static const int klotsize = 100;
protected:
	static string zippath() { return "7z"; }
	datatype type;
    bool tick16;

	void import();
	static const int maxmarketsperstock = 30;
	tickdataimporter(date astart, date aend, smpairs& apairs,
	     datatype atype, bool atick16=true):
		importer(astart, aend, apairs),
		numticks(atick16 ? 16 : 100), type(atype) {}
	virtual void importstock(string& stock, vector<string>& markets, vector<char>& marketids) = 0;

	struct tinfo {double tt; int p; int q; int pairedwith; int type; };
	  // type = 1 if price == bid, -1 if price = ask, 0 otherwise

	void importusfile(string& astock, vector<string>& markets, vector<char>& marketids, string qfn, string tfn);

};

class tickdatatwimporter : public tickdataimporter
{
	string folder;
public:
	tickdatatwimporter(date astart, date aend, smpairs& apairs,
	         datatype atype, string& afolder, bool tick16):
		tickdataimporter(astart, aend, apairs, atype, tick16), folder(afolder) {}
protected:
	void importstock(string& stock, vector<string>& markets, vector<char>& marketids);
};

class tickdata09importer : public tickdataimporter
{
	string folder;
public:
	tickdata09importer(date astart, date aend, smpairs& apairs, datatype atype, string& afolder):
		tickdataimporter(astart, aend, apairs, atype), folder(afolder) {}
protected:
	void importstock(string& stock, vector<string>& markets, vector<char>& marketids);
};

class tbtexporter : public analysis
{
    ofstream *ofs;
protected:
	string csvfn()
	{
		stringstream s;
		s << pairid();
		const string kcsvdir = "/csv/";
		const string kcsvsuff = ".csv";
		return hfd::interdir() + kcsvdir + s.str() + kcsvsuff;
	}

public:
	tbtexporter() :
	   analysis("export", data::filterquote, 3),
	   ofs(0)
	         {}

	void onstartpair(smpair& togo);
	void onendpair(smpair& togo);
	void onrecord(hfdrecord& rec, int recn);
};

class interval : public analysis
{
public:
	struct rec
	{
	    int day;
	    int month;
	    int year;
	    double time;
	    double actA;
	    double actB;
	    double nfSnapA;
	    double nfSnapB;
	    double avgA;
	    double avgB;
	    double dNA;
	    double dNB;
	    double N;
	};
private:
	double inter;

	static constexpr double notime = 0;
	static constexpr double noprice = -1;

	double lasttime; // cas posledniho snimku
	/// snapshot of the price, excluding price movements after trades
	double lastqrectime; // cas posledniho recordu s kotama
	int sumNa;
	int sumNb;
	int lasta;
	int lastb;
	int lastnfa;
	int lastnfb;
	double suma;
	double sumb;
	double sumN;
	void resetcounters();
	int totalrecs;
	int totalvol;
	int totaltrades;
	double totalspread;
	double warming;
	void step(bool auselast);
protected:
	double getintervalsize() { return inter; }
	int gettotalrecs() { return totalrecs; }
	int gettotalvol() { return totalvol; }
	int gettotaltrades() { return totaltrades; }
	double getavgspread() { return (double) totalspread / (double) totalrecs; }
	virtual void onsnapshot( interval::rec& r) {}
public:
	/// \p warming is a time for which initial values are gathered. If no record occurs during warming, exception is thrown
	interval(const string& aname, double aint, int anumrows,double awarming = 0.5) : analysis(aname, data::filterquote, 3,anumrows),
	         inter(aint), warming(awarming) {}

	void onstartpair(smpair& togo);
	void onrecord(hfdrecord& rec, int recn);
	void onstartday(int aday)
	{
		lasttime=notime;
	}
	void onendday(int day);
};


class csvinterval: public interval
{
	static const int numcolumns = 13;
	ofstream* ofs;
	string csvfn()
	{
		stringstream s;
		s << pairid() << "_i" << getintervalsize();
		const string kcsvdir = "/csv/";
		return hfd::interdir() + kcsvdir + s.str();
	}

public:
	csvinterval(const string& aname, double aint, int anumrows) : interval(aname, aint, anumrows), ofs(0) {}
	virtual ~csvinterval()
	{
		delete ofs;
	}

	void onendday(int )
	{
		for(unsigned int i=0; i<numcolumns-1; i++)
			*ofs << ",";
		*ofs << endl;
	}
	void onstartpair(smpair& togo);
	void onendpair(smpair& togo);
	void onsnapshot( interval::rec& r);
};


class importer09
{
public:
    static void  import(const char** stocks, const char** markets,
            date& start, date& end, string& folder);
};


class importertw
{
public:
    static void  import(const char** stocks, const char** markets,
            date& start, date& end, string& folder, bool tick16 = false);
};

#endif /* HFD_HPP_ */
