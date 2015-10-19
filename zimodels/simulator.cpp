

#include "zimodels.hpp"

#include <boost/random/mersenne_twister.hpp>
  using boost::mt19937;
#include <boost/random/poisson_distribution.hpp>
  using boost::poisson_distribution;
#include <boost/random/binomial_distribution.hpp>
  using boost::binomial_distribution;
#include <boost/random/exponential_distribution.hpp>
  using boost::exponential_distribution;
#include <boost/random/discrete_distribution.hpp>
#include <boost/random/variate_generator.hpp>
  using boost::variate_generator;
using namespace boost::random;



int cdasimulator::getpoisson(double lambda)
{
    if(lambda == 0)
        return 0;
    boost::poisson_distribution<int> pdist(lambda);
    boost::variate_generator< boost::mt19937&, boost::poisson_distribution<int> > rvt(gen, pdist);
    return rvt();
}

int cdasimulator::getbi(int n, double p)
{
    boost::binomial_distribution<int> pdist(n,p);
    boost::variate_generator< boost::mt19937&, boost::binomial_distribution<int> > rvt(gen, pdist);
    return rvt();
}



//
// We can repeat the above example, with other distributions:
//

void cdasimulator::simulate(int an, zirec* az, int& resnum)
{
    cout << "Simulation: " << an << " required" << endl;
    vector<double> theta = model.getinitparams();

    gen.seed(0);

    int nuinv = (int) (1.0 / model.getnu(theta) + 0.5);
    int z=0;
    int wn=0;
    vector<double> A(N);

    if(N<4)
        throw runtime_error("too litle n");

    int b=N/2;
    int a=b+1;

    for(int j=0; j<b; j++)
        A[j]=-1;

    for(int j=a-1; j<N; j++)
        A[j]=1;

    double t=0;
    double tg=0;

    bool firstinday = true;
    for(int i=0; ; i++)
    {
        double sint = 0;
        int d=a-b;
        int ps = 2*d+4;
        vector<double> p(ps);
        for(int j=0; j<d; j++)
        {
            sint +=
              (p[j]=model.getkappa(a,b,b+j+1,theta))+
              (p[d+j]=model.getlambda(a,b,b+j,theta));
        }
        sint +=
           (p[2*d]=A[a-1]*model.getrho(a,b,a,theta))
           +(p[2*d+1]=-A[b-1]*model.getsigma(a,b,b,theta));
        sint +=
           (p[2*d+2] = model.gettheta(a,b,a,theta))
           +(p[2*d+3] = model.getvartheta(a,b,b,theta));

        for(int j=0; j<ps; j++)
            p[j] /= sint;


        boost::exponential_distribution<double> edist(sint);
        boost::variate_generator< boost::mt19937&, boost::exponential_distribution<double> > egen(gen, edist);

        boost::random::discrete_distribution<int> dist(p);
        boost::variate_generator< boost::mt19937&, boost::random::discrete_distribution<int> > vgen(gen, dist);
        int v = vgen();
        double c = egen();
        int outs = 0; // signalizes a jump out of the spread
        bool ins = false;
        bool bmo = false;
        int newa = a;
        int newb = b;
        if(v==2*d+2 || v==2*d)
        {
            if(a<=N)
            {
                if(A[a-1] > 0)
                    A[a-1]--;
                if(A[a-1]==0)
                    outs = 1;
                if(v==2*d+2)
                    bmo = true;
            }
        }
        else if(v == 2*d+1 || v==2*d+3)
        {
            if(b>0)
            {
                if(A[b-1] < 0)
                    A[b-1]++;
                if(A[b-1]==0)
                    outs = -1;
            }
        }
        else if(v<d)
        {
            int p=b+v+1;
            if(p <= N)
            {
                A[p-1]++;
                if(p!=a)
                {
                    newa=p;
                    ins = true;
                }
            }
        }
        else
        {
            if(v>=2*d)
                throw runtime_error("bad v");
            int p = b+v-d;
            if(p>=1)
            {
                A[p-1]--;
                if(p != b)
                {
                    newb=p;
                    ins = true;
                }
            }
        }
        t += c;
        bool regenerate = ins || outs;
        if(regenerate)
        {
//cout << "generating" << endl;
            double dt=t-tg;

            for(int j=1; j<b; j++)
            {
                double e = exp(-model.getsigma(a,b,j,theta)*dt);
                int bi = getbi(-A[j-1],e);
                double lambda = model.getpsi(a,b,j,theta)
                    * (1.0 - e);
                int po = getpoisson(lambda);
//cout << "e=" << e << " lambda=" << lambda << " po=" << po << " bi=" << bi << endl;
                A[j-1] = -po - bi;
            }
            for(int j=a+1; j<=N; j++)
            {
                double e = exp(-model.getrho(a,b,j,theta)*dt);
                int bi = getbi(A[j-1],e);
                double lambda = model.getphi(a,b,j,theta)
                    * (1.0 - e);
                int po = getpoisson(lambda);
                A[j-1] = po + bi;
            }
            tg = t;
        }
        if(outs < 0)
        {
            int j=b-1;
            for(; j>=1; j--)
            {
                if(A[j-1] < 0)
                {
                    newb = j;
                    break;
                }
            }
            if(j==0)
                break;
        }
        else if(outs > 0)
        {
            int j=a+1;
            for(; j<=N; j++)
            {
                if(A[j-1] > 0)
                {
                    newa = j;
                    break;
                }
            }
            if(j==N+1)
                break;
        }
        a = newa;
        b = newb;
//cout << a << "..." << b << endl;
        bool newzi = ins || outs || bmo;
        if(newzi)
        {
            if(wn++ > warmupn)
            {
                zirec& r=az[z++];
                r.t = t;
                r.a = a;
                r.b = b;
                r.firstinday = firstinday;
                r.q = A[a-1]* nuinv;
                r.s = (bmo ? 1 : 0) * nuinv;
//                if(z%100==0)
//                    cout << z << ": t=" << t << " a=" << a << " b=" << b << " q=" << r.q << " s=" << r.s << endl;

                firstinday = false;
                if(z==an)
                    break;
            }
        }

    }
    cout << z << " obtained." << endl;
    resnum = z;
}

