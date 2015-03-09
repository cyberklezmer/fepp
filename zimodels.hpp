/*
 * zimodels.hpp
 *
 *  Created on: Mar 13, 2013
 *      Author: martin
 */

#ifndef ZIMODELS_HPP_
#define ZIMODELS_HPP_

#include "epp.hpp"
#include "hfd.hpp"
#include <sstream>
#include <boost/random/mersenne_twister.hpp>

using boost::mt19937;
using namespace epp;


class zisummary : public analysis
{
    static const int kmaxrecjmp = 20;
	int starty;
    int endy;
	const int numyears() { return endy-starty+1;}
	static const int kmaxyears = 10;

	double msumspreads;
	double msumdeltat;

	int mnumtrans;
	int numberdaysinmonth;
	int numrecsinmonth;

	double mspreads[kmaxyears*12];
	double mnumbers[kmaxyears*12];
	double mtimes[kmaxyears*12];
    int jumphist[kmaxrecjmp];
    int jumpdhist[kmaxrecjmp];
    int reldjumps[kmaxrecjmp];

    static const int kmaxlogqs = 10;
    int logqs[kmaxlogqs];

	string stockfn()
	{
		const string sum = "sum_";
		string result = sum + getpair().stock + getpair().market;
		for(unsigned int i=0; i<result.length(); i++)
		{
			if(result[i] == ' ')
				result[i] = '_';
		}
		return result;
	}
public:
	zisummary(int astarty, int aendy):
	    analysis("summary", data::filterquote, 3, 5),
		starty(astarty), endy(aendy)  {}

	void onstart() {}
	void onstartpair(smpair& togo);
	void onstartmonth(int month);
	void onrecord(hfdrecord& rec, int recn);
	void onendmonth(int month);
	void onendday(int day);
	void onendpair(smpair& togo);
};

struct zirec
{
    double t; // time of firs apperance
    int a;
    int b;
    double q; // last observed volume with quotes a and b
    int s; // volume traded
    int m; // difference of volume wrt last jump of a
    bool firstinday;
};


class zianalysis;


class zimodel
{
protected:
    vector <paraminfo> params;
    vector<double> phoograd;
public:
    virtual string getname() = 0;

    virtual double getrho(int a, int b, int p,
             const vector<double>& theta, vector<double>& grad) = 0;
     virtual double getphi(int a, int b, int p,
             const vector<double>& theta, vector<double>& grad) = 0;
    virtual double getkappa(int a, int b, int p,
             const vector<double>& theta, vector<double>& grad)
    {
        int q = theta.size();
        vector<double> pg(q);
        vector<double> rg(q);
        double rho = getrho(a,b,p,theta,rg);
        double phi = getphi(a,b,p,theta,pg);
        grad = pg * rho + rg * phi;

        return phi * rho;
    }

    double getkappa(int a,int b, int p, const vector<double>& theta)
    {
        return getkappa(a,b,p > a ? p : a+1,theta, phoograd);
    }
    double getrho(int a, int b, int p, const vector<double>& theta)
    {
        return getrho(a,b,p > a ? p : a+1,theta, phoograd);
    }

    virtual double getphi(int a, int b, int p, const vector<double>& theta)
    {
        return getphi(a,b,p,theta, phoograd);
    }

    virtual double gettheta(int a, int b, int p, const vector<double>& theta)
    {
        return getkappa(a,b,p,theta);
    }

    double getlambda(int a,int b, int p, const vector<double>& theta)
    {
        return getkappa(b,2*a-b,b-(p-a),theta, phoograd);
    }
    double getsigma(int a, int b, int p, const vector<double>& theta)
    {
        return getrho(b,2*a-b,b-(p-a),theta, phoograd);
    }

    virtual double getpsi(int a, int b, int p, const vector<double>& theta)
    {
        return getphi(a,2*a-b,b-(p-a),theta, phoograd);
    }

    virtual double getvartheta(int a, int b, int p, const vector<double>& theta)
    {
        return gettheta(a,2*a-b,b-(p-a),theta);
    }

    virtual double getnu(const vector<double>& theta, vector<double>& grad)
    {
        grad = zero_vector(theta.size());
        return 1;
    }
    double getnu(const vector<double>& theta)
    {
        return getnu(theta,phoograd);;
    }
    virtual double getgamma(const vector<double>& theta, vector<double>& grad)
    {
        grad = zero_vector(theta.size());
        return 1;
    }

    virtual double geteta(const vector<double>& theta, vector<double>& grad)
    {
        grad = zero_vector(theta.size());
        return 1;
    }

    virtual double getzeta(const vector<double>& theta, vector<double>& grad)
    {
        grad = zero_vector(theta.size());
        return 0;
    }


    void setinitparams(const vector<double>& p)
    {
        if(p.size() != params.size())
            throw runtime_error("zimodel - size mismatch");
        for(unsigned int i=0;i<p.size(); i++)
            params[i].initial = p[i];
    }

    vector<double> getinitparams()
    {
        vector<double> p(params.size());
        for(unsigned int i=0; i<p.size(); i++)
            p[i] = params[i].initial;
        return p;
    }

    vector<paraminfo>& getparams() { return params; };
    zimodel(int anumpars)
      : params(anumpars), phoograd(anumpars)
    {}
    virtual ~zimodel() {}
};

class ziextmodel : public zimodel
{
protected:
    bool isnu;
    int nuindex;
    bool isgamma;
    int gammaindex;
    bool iseta;
    int etaindex;
    bool iszeta;
    int zetaindex;

    virtual double getnu(const vector<double>& theta, vector<double>& grad)
    {
        grad = zero_vector(theta.size());
        if(isnu)
        {
            double nuu = exp(theta[nuindex]);
            grad[nuindex] = nuu;
            return nuu;
        }
        else
        {
            return 1;
        }
    }
    virtual double getgamma(const vector<double>& theta, vector<double>& grad)
    {
        grad = zero_vector(theta.size());
        if(isgamma)
        {
            double gamma = theta[gammaindex];
            grad[gammaindex] = 1;
            return gamma;
        }
        else
        {
            return 0;
        }
    }
    virtual double geteta(const vector<double>& theta, vector<double>& grad)
    {
        grad = zero_vector(theta.size());
        if(iseta)
        {
            double eta = theta[etaindex];
            grad[etaindex] = 1;
            return eta;
        }
        else
        {
            return 1;
        }
    }
    virtual double getzeta(const vector<double>& theta, vector<double>& grad)
    {
        grad = zero_vector(theta.size());
        if(iszeta)
        {
            double zeta = theta[zetaindex];
            grad[zetaindex] = 1;
            return zeta;
        }
        else
        {
            return 0;
        }
    }

    ziextmodel(int anumpars, bool aisnu, bool aisgamma,
              bool aiseta, bool aiszeta)
      : zimodel(anumpars+(aisnu ? 1 : 0)+(aisgamma ? 1 : 0)
                      +(aiseta ? 1 : 0) + (aiszeta ? 1 : 0)),
        isnu(aisnu), isgamma(aisgamma), iseta(aiseta), iszeta(aiszeta)
    {
        int i=anumpars;
        if(aisnu)
        {
            nuindex = i++;
            params[nuindex].name = "\\nu";
            params[nuindex].initial = 0;
        }
        if(aisgamma)
        {
            gammaindex = i++;
            params[gammaindex].name = "\\gamma";
            params[gammaindex].initial = 0.001;
            params[gammaindex].lower = 0;
        }
        if(aiseta)
        {
            etaindex = i++;
            params[etaindex].name = "\\eta";
            params[etaindex].initial = 1;
            params[etaindex].lower = 0;
            params[etaindex].upper = 1;
        }
        if(aiszeta)
        {
            zetaindex = i++;
            params[zetaindex].name = "\\zeta";
            params[zetaindex].initial = 0;
            params[zetaindex].lower = 0;
            params[zetaindex].upper = 0.49999999999;
        }

    }
};

class aqestimator
{
public:
    static int debug;
protected:
    zimodel& model;
    zianalysis& analysis;

    class hrec
    {
public:
        double point;
        double csum;
        double A;
        vector<double> gpoint;
        vector<double> gcsum;
        hrec(): point(0), csum(0), A(0) {}
    };
private:
    bool snapshots;
    vector<hrec> snaps;
    vector<int> snaptimes;
    vector<double> params;

    bool issnap(int p, int i, int j);
    const hrec& getsnap(int p, int i, int j);
    void setsnap(int p, int i,int j, const hrec& r);
    void resetsnaps();

    vector<double> parameters;
    /// assumes correct size of grads
    void getomega
       (int i, int p,
        double& point,vector<double>& gpoint,
        double& bin,vector<double>& gbin,
        double& bip, vector<double>& gbip);
public:
    double getomega(int i, int p, vector<double>& g);

    aqestimator(zimodel& amodel, zianalysis& aanalysis);
    void setparameters(const vector<double>& p)
    {
        params = p;
        resetsnaps();
    }
    void addpo(vector<double>& dist,
        vector<vector<double> >& grads,
        double point, const vector<double>& gpoint);
    void addbi(vector<double>& dist,
            vector<vector<double> >& grads,
            int bin,
            double bip, const vector<double>& gbip);


    double getlogdensity
       (int i, int a, int m, int v, vector<double>& grad);
    double getlogdensity
       (int i, int a,  vector<double>& grad);


    double logdensity(int i, vector<double>& grad);

    void resetparameters()
    {
        params = vector<double>();
    }
    virtual double estimate(vector<paramresult>& r ) = 0;
    int getsamplesize();
    void testgrads(int start, int num);
};

class zianalysis : public analysis
{
public:
    enum emodeltype {individual, tail};
    emodeltype modeltype;
    int firstn;
    int maxn;
    double warmuptime;
    bool unitvolume;
    bool resample;
    bool twodimestimation;
    bool extendedlogging;
    bool extendedoutput;
    bool includenu;
    bool includegamma;
    bool includeeta;
    bool includezeta;
	double maxmletime;
protected:
	bool onlytrades;
 	int maxajumps;
	int maxqchanges;

protected:
    empdistn* qdist;
    /// records quote changes
    zirec* qchanges;
    int M;

    int* ajumps;
	int N;
	int Noff;


    // filled by onenepair
    int mina;
    int maxa;

    // "local" in onrecord
    double firstt;
    int sp;
    double sda;
    int nsda;
    double ss;
    int nss;

    static const int offsampleratio = 10;
public:
    int getN() { return N; }
    int getmina() { return mina; }
    int getmaxa() { return maxa; }
    int numprices()
    {
        return mina > maxa ? 0 :  getmaxa()-getmina()+1;
    }
    int getjumps(int i) { return ajumps[i]; }

    /// index of the jth item of history relative to itj (j should be positive)
    /// it is guqranteed that \p X(i, 0).a != X(i, 1).a
    zirec& X(int i, int j)
    {
        return qchanges[ajumps[i]-j];
    }

    bool iskp(int i, int j, int p, bool& kpiszero)
    {
        zirec& r = X(i,j);
        kpiszero = r.firstinday;
        return kpiszero || r.a >= p;
    }

    void setsample(int samplesize, bool aonlytrades)
    {
    	maxajumps = samplesize * (offsampleratio +1 ) / offsampleratio;
	    maxqchanges = maxajumps * 10000;
        onlytrades = aonlytrades;
    }

	zianalysis(const string& aname = "zianalysis") :
		  analysis(aname, data::filterquote,3,3),
		  modeltype(tail),
		  firstn(1),
		  maxn(4),
		  warmuptime(10*60),
		  unitvolume(false),
		  resample(false),
		  twodimestimation(true),
		  extendedlogging(false),
		  extendedoutput(false),
		  includenu(false),
		  includegamma(false),
		  includeeta(false),
		  includezeta(false),
	      maxmletime(10*60),
		  qchanges(0),
		  ajumps(0)
	{
        setsample(500,true);
        if(resample && !unitvolume)
            throw logic_error("Resample implemented only for unit volume");
	}

	~zianalysis ()
	{
        delete qchanges;
        delete qdist;
        delete ajumps;
    }

	void onstartpair(smpair& togo);
	void onrecord(hfdrecord& rec, int recn);
	void onendpair(smpair& togo);

    void estimate(zimodel& model, vector<paramresult>& res,
          bool catchex, bool& failed,
          bool& significant, double& ll,
          string& failmsg);
    double evaluate(zimodel* modelused, const vector<paramresult>& pars,
                                                    int from, int len, double avgdaplus);

};

class zimleestimator : public aqestimator, public mle
{
public:
    zimleestimator(zimodel& amodel, zianalysis& aanalysis) :
       aqestimator(amodel,aanalysis),
        mle(amodel.getparams())
       {}
	virtual int getN() { return analysis.getN(); };
    virtual void beforeloglikeval(vector<double> aparams)
    {
/*        debug = 1;
        vector<double> ph(2);
        ph[0]=2;
        ph[1]=1;
        setparameters(ph);
        vector<double> phoog(2);
        cout << "ld=" << getlogdensity(2,6726,2,phoog) << endl;
        debug = 0;*/
//        throw;
        setparameters(aparams);
    }

	virtual double evallogdensity
	  (int i, const vector<double>&, vector<double>& g)
	  {
//cout << i << endl;
        double r = logdensity(i, g);;
//cout << "r= " << r << "g=";
//for(int i=0; i< g.size(); i++)
    //cout << g[i] << " ";
//cout << endl;
        return r;
	  }
	virtual void afterloglikeval(const vector<double>& aParams,
	                     const vector<double>& grad, double loglik)
            { /*resetparameters(); */}
    virtual double estimate(vector<paramresult>& r )
    {
        return mle::estimate(r,true);
    }
};



class zinparmodel : public ziextmodel
{
    int n;
public:
    int getn() { return n; }
    virtual string getname()
    {
        ostringstream s;
        s <<  "IP(" << n << ")";
        return s.str();
    };

    virtual double getrho(int a,int b, int p,
             const vector<double>& theta, vector<double>& g)
    {
        for(unsigned int j=0; j<g.size(); j++)
            g[j] = 0;

        int i = p-a;
        if(i<=0)
        {
            if(n > 1)
                throw runtime_error("getrho called for insread rho");
            else
                i=1;
        }
        else if(i>n)
            i=n;


        g[n+i-1] = 1;
        return theta[n+i-1];
    };



    virtual double getphi(int a, int b, int p,
             const vector<double>& theta, vector<double>& g)
    {
        for(unsigned int j=0; j<g.size(); j++)
            g[j] = 0;

        int i = p-a;
        if(i<=0)
           i=1;
        else if(i>n)
            i=n;

        g[i-1] = 1;
        return theta[i-1];
    };

    zinparmodel(int an, bool aisnu, bool aisgamma, bool aiseta, bool aiszeta) :
       ziextmodel(2*an, aisnu, aisgamma, aiseta, aiszeta), n(an)
    {

        for(int i=0; i<n; i++)
        {
            ostringstream str;
            str << i;
            params[i].name = "\\phi_" + str.str();
            params[i].initial = 2;
            params[i].lower = 0.00001;
            params[n+i].name = "\\rho_" + str.str();
            params[n+i].initial = 1;
            params[n+i].lower = 0.00001;
        }
    }

};


class zinptmodel : public ziextmodel
{
    int n;
public:
    int getn() { return n; }

    int firstphi() { return 0; }
    int alphaphi() { return n; }
    int firstrho() { return n+1; }
    int alpharho() { return 2*n+1; }
    virtual string getname()
    {
        ostringstream s;
        s <<  "IPT(" << n << ")";
        return s.str();
    };


    virtual double getrho(int a,int b, int p,
             const vector<double>& theta, vector<double>& g)
    {
        for(unsigned int j=0; j<g.size(); j++)
            g[j] = 0;

        int i = p-a;
        if(i<=0)
        {
            throw runtime_error("getrho called for insread rho");
        }
        if(i<=n)
        {
            g[firstrho()+i-1] = 1;
            return theta[firstrho()+i-1];
        }
        else
        {
            double c = theta[firstrho()+n-1];
            double alpha = theta[alpharho()];
            double pw = pow(i-n+1,alpha);
            double rho = c * pw;
            g[firstrho()+n-1] = pw;
            g[alpharho()] = rho * log(i-n+1);
            return rho;
        }
    };
    virtual double getphi(int a,int b, int p,
             const vector<double>& theta, vector<double>& g)
    {
        for(unsigned int j=0; j<g.size(); j++)
            g[j] = 0;

        int i = p-a;
        if(i<=0)
           i=1;
        if(i<=n)
        {
            g[firstphi()+i-1] = 1;
            return theta[firstphi()+i-1];
        }
        else
        {
            double k = theta[firstphi()+n-1];
            double alpha = theta[alphaphi()];
            double pw = pow(i-n+1,alpha);
            double phi = k * pw;
            g[firstphi()+n-1] = pw;
            g[alphaphi()] = phi * log(i-n+1);
            return phi;
        }
    };


    zinptmodel(int an, bool aisnu, bool aisgamma, bool aiseta, bool aiszeta) :
       ziextmodel(2*an+2,  aisnu, aisgamma, aiseta, aiszeta), n(an)
    {

        for(int i=0; i<n; i++)
        {
            ostringstream str;
            str << i;
            params[firstphi()+i].name = "\\phi_" + str.str();
            params[firstphi()+i].initial = 2;
            params[firstphi()+i].lower = 0.00001;
            params[firstrho()+i].name = "\\rho_" + str.str();
            params[firstrho()+i].initial = 1;
            params[firstrho()+i].lower = 0.00001;
        }
        params[alphaphi()].name = "\\alpha_\\phi";
        params[alphaphi()].initial = 0;
        params[alphaphi()].upper = 0;
        params[alpharho()].name = "\\alpha_\\rho";
        params[alpharho()].initial = 0;
        params[alpharho()].upper = 0;
    }

};


class cdasimulator
{
    zimodel& model;
    int N;
    mt19937 gen;
    int warmupn;
public:
    int getpoisson(double lambda);
    int getbi(int n, double p);
    cdasimulator(zimodel& amodel,
          int amaxprice, int awarmupn):
          model(amodel),
          N(amaxprice), warmupn(awarmupn)
    {
    }
    void simulate(int an, zirec* az, int& resnum);
};

#endif /* ZIMODELS_HPP_ */
