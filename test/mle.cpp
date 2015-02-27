#include "epp.hpp"


using namespace epp;


double sample[] = {1,2,5,3,2,6,1,8,1};

class testmle : public mle
{
	virtual int getN()
	{
        return sizeof(sample)/sizeof(sample[0]);
    }
public:
	virtual double evallogdensity(int i, const vector<double>& x,
	                                              vector<double>& g)
    {
        double s = (sample[i]-x[0])*(sample[i]-x[0]);
        g[0] = (sample[i]-x[0]) / x[1] ;
        g[1] = - 1.0 / x[1] / 2 + s / (x[1] * x[1])/2.0;
        return -log(2*3.1415*x[1]) / 2.0 - s / x[1] / 2.0;
    }
public:


	testmle(const vector<paraminfo>& aparams):
		mle(aparams)
      {}
};


int main()
{
    vector<paraminfo> p(2);
    p[0].name = "m";
    p[0].initial = 1;
    p[1].name = "V";
    p[1].lower = 0.00001;
    p[1].initial = 1;

    testmle m(p);

    vector<paramresult> r(2);

    m.estimate(r,true);
    cout.precision(15);
    cout << r[0] << endl << r[1] << endl;
}

