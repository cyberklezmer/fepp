/*
 * estimator.cpp
 *
 *  Created on: 23 Juand 2015
 *      Author: martin
 */

#include "zimodels.hpp"
#include <boost/math/special_functions/gamma.hpp>

int ziestimator::debug = 0;


void dumpparams(const vector<double>& params)
{
    for(unsigned int i=0; i<params.size(); i++)
    {
        cout << "#" << i+1 << " " << params[i] << ", ";
    }
    cout << endl;
}

struct pairinfo
{
    const char* excludefor;
    double dataextfactor;

    pairinfo(): excludefor(0), dataextfactor(1) {}
};

ziestimator::ziestimator(zimodel& amodel, zianalysis& aanalysis) :
    model(amodel), analysis(aanalysis),
    snapshots(false), parameters(0)
{
    resetsnaps();
}

bool ziestimator::issnap(int p, int i, int j)
{
//cout << "issnap " << p << "," << i << "," << j << endl;
    if(!snapshots)
        return false;
    int k=p-analysis.getmina();
    return snaptimes[k] == analysis.getjumps(i)-j;
}

const ziestimator::hrec& ziestimator::getsnap(int p, int i, int j)
{
//cout << "getsnap " << p << "," << i << "," << j << endl;
    if(!issnap(p,i,j))
        throw runtime_error("snapshot not present here");
    int k=p-analysis.getmina();
    if(k<0 || (unsigned)k>snaps.size())
        throw runtime_error("bounds exceeded in getsnap");
    return snaps[k];
}

void ziestimator::setsnap(int p, int i, int j,
                          const ziestimator::hrec& r )
{
//cout << "setsnap " << p << "," << i << "," << j << endl;
//cout << "shaps.size()=" << snaps.size() << endl;
//cout << "(" << analysis.getmaxa()
//     << "," << analysis.getmina() << ")" << endl;
    if(snapshots)
    {
        int k=p-analysis.getmina();
        if(k<0 || (unsigned)k>snaps.size())
            throw runtime_error("bounds exceeded in getsnap");
        snaptimes[k] = analysis.getjumps(i)-j;
        snaps[k]=r;
    }
//cout << "finished." << endl;
}

void ziestimator::resetsnaps()
{
    unsigned int n = analysis.numprices();
    if(snaptimes.size() != n)
    {
        snaps.resize(n);
        snaptimes.resize(n);
    }
    for(unsigned int i=0; i<snaptimes.size(); i++)
        snaptimes[i] = -1;
}


void ziestimator::getomega
(int i, int p,
 double& point,vector<double>& gpoint,
 double& bin,vector<double>& gbin,
 double& bip, vector<double>& gbip)
{
    int q = params.size();

    point = 0;
    gpoint = zero_vector(q);

    double csum=0;
    vector<double> gcsum(zero_vector(q));

    int A=0;

    int j = 1;

    if(debug)
        cout << "getomega(" << i << "," << p << ")" << endl;
    for(;;)
    {
        double eckminus1 = exp(-csum);

        zirec& kth = analysis.X(i,j); // k=i-j;
        double deltat = analysis.X(i,j-1).t-kth.t;
        if(debug>1)
            cout << "kth: a=" << kth.a << ", t=" << kth.t << " deltat=" << deltat << endl;
        vector<double> gro(q);
        double dc = deltat * model.getrho(kth.a,kth.b,p,params,gro);
        vector<double> gdc = deltat * gro;

        vector<double> gphi(q);
        double phi=model.getphi(kth.a,kth.b,p,params, gphi);

        double eck = exp(-csum-dc);

        point += phi * (eckminus1 - eck);

        gpoint = gpoint +  gphi * (eckminus1 - eck)
                 + phi * ( eck * (gcsum+gdc) - eckminus1 * gcsum);

        csum += dc;
        gcsum = gcsum + gdc;

        j++;

        bool kpiszero;
        if(analysis.iskp(i,j,p,kpiszero))
        {
            if(kpiszero || analysis.X(i,j).a != p)
            {
                if(!kpiszero)
                {
                    vector<double> ggamma(q);
                    double gamma = model.getgamma(params,ggamma);
                    point += eck * gamma;
                    gpoint = gpoint + eck * ggamma - gamma * eck * gcsum;
                };

                bin = 0;
                gbin = zero_vector(q);
                bip = 0;
                gbip = zero_vector(q);
            }
            else
            {
                A = analysis.X(i,j).q;
                vector<double> gnu(q);
                double nu=model.getnu(params,gnu);
                vector<double> geta(q);
                double eta=model.geteta(params,geta);
                bin = A*nu;
                gbin = A*gnu;
                bip = eck * eta;
                gbip = eck * (-1 * gcsum * eta + geta);
            }
            break;
        }
    }

    if(debug)
        cout << "t=" << analysis.X(i,0).t-analysis.X(i,j).t << " point="
             << point << endl;
}


double ziestimator::getomega(int i, int p, vector<double>& g)
{
    int q = params.size();
    double point;
    vector<double> gpoint(q);
    double bin;
    vector<double> gbin(q);
    double bip;
    vector<double> gbip(q);

    getomega(i, p, point,gpoint,bin,gbin,bip,gbip);

    if(bip==1)
        throw runtime_error("minus infinite omega");

    double l=log(1-bip);
    g = gpoint - l * gbin + bin / (1-bip) * gbip;
    return point - bin * l;
}


double ziestimator::getlogdensity
(int i, int a, vector<double>& grad)
{
    if(debug)
        cout << "getlogdensity(" << i << "," << a << ")"
             << " t=" << analysis.X(i,0).t << endl;
    int q = params.size();
    zirec& last=analysis.X(i,1);
    if(a > last.a)
    {
        grad = zero_vector(q);
        double s = 0;

        for(int p=last.a+1; p < a; p++)
        {
            vector<double> g(q);
            s -= getomega(i, p, g);
            if(debug)
                cout << "p=" << p << " s=" << s << endl;
            grad = grad - g;
        }

        vector<double> go(q);
        double o = getomega(i, a, go);

        vector<double> gzeta(q);
        double zeta = model.getzeta(params, gzeta);

        double prm1 = pow(zeta,a-last.a-1);
        double prm0 = prm1 * zeta;
        double prshifta = prm0 ;
        vector<double> gprshifta
         // =((p-a) * prma * (1-zeta)+ pow(zeta,p-a)) / (1-zeta)^2
           = (a-last.a) * prm1 * gzeta;
        double prnotshorga
          = 1.0 - (zeta - pow(zeta,a-last.a+1)) / (1.0-zeta);
        vector<double> gprnotshorga
          // =- (1.0-zeta+zeta) / (1-zeta)^2
           = - ((1.0 - (a-last.a+1) * prm0) * (1.0 - zeta)
                + zeta - prm0 ) / (1.0 - zeta)/ (1.0 - zeta)
                * gzeta;

        double pr = prshifta +
                    prnotshorga * (1-exp(-o));


        vector<double> gpr  //= gprwasshift
                      // - (1-exp(-o)) * gprwasshift
                      //+ (1-prwasshift)* exp(-o) * go
            = gprshifta + (1-exp(-o)) * gprnotshorga
                  + prnotshorga * exp(-o) * go;

//grad = grad + go / (exp(o)-1);
        grad = grad + 1.0 / pr * gpr;

        if(0 && debug)
            cout << "o=" << o << " dens=" << exp(s)*(1.0-exp(-o))
                 << "=" << exp(s + log(1.0-exp(-o)))
                 << endl;

        return  s + log(pr); // log(1.0-exp(-o));
    }
    else if(model.isinspread() && a < last.a && a> last.b)
    {
        double ksum = 0;
        vector<double> gsum(zero_vector(q));

        double ka;
        vector<double> ga(q);

        for(int p=last.b+1; p<last.a; p++)
        {
            vector<double> g(q);

            double k = model.getkappa(last.a, last.b, p, params, g);
            if(p==a)
            {
                ka=k;
                ga=g;
            }

            ksum += k;
            gsum = gsum + g;
        }

        for(int j=0; j<q; j++)
            grad[j] = ga[j] / ka - gsum[j] / ksum;
        return log(ka) - log(ksum);
    }
    else
    {
        grad = zero_vector(q);
        return 0;
    }
}


void ziestimator::addpo(vector<double>& dist,
                        vector<vector<double> >& gdist,
                        double point, const vector<double>& gpoint)
{
    int q = gpoint.size();
    int m=dist.size();

    vector<double> res = zero_vector(m);
    vector<vector<double> > gres(m);
    for(int i=0; i<m; i++)
        gres[i]=zero_vector(q);
    double ppo = exp(-point);
    vector<double> gppo = -ppo * gpoint;

    for(int i=0;; )
    {
        for(int j=i; j<m; j++)
        {
            res[j]+= ppo*dist[j-i];
            gres[j] = gres[j] + dist[j-i]*gppo
                      + ppo*gdist[j-i];
        }
        i++;
        if(i>=m)
            break;
        gppo = 1.0 / (double) i * (ppo * gpoint + point * gppo);
        ppo *= point / (double) i;
    }
    dist = res;
    for(int k=0; k<m; k++)
        gdist[k] = gres[k];


}

void ziestimator::addbi(vector<double>& dist,
                        vector<vector<double> >& gdist,
                        int bin,
                        double bip, const vector<double>& gbip)
{
    int q = gbip.size();
    double pbi = pow(1.0-bip,bin);
    vector<double> gpbi = - bin / (1.0-bip) * pbi * gbip ;
    int m=dist.size();

    vector<double> res = zero_vector(m);
    vector<vector<double> > gres(m);
    for(int i=0; i<m; i++)
        gres[i]=zero_vector(q);

    for(int i=0;; )
    {
        for(int j=i; j<m; j++)
        {
            res[j]+= pbi*dist[j-i];
            gres[j] = gres[j] + dist[j-i]*gpbi
                      + pbi*gdist[j-i];
        }
        i++;
        if(i>=m || i>bin)
            break;
        double f1 = 1.0 / (1.0 / bip - 1.0);
        double f2 = (double) (bin-i+1) / i;
        vector<double> gf1 = 1.0 / (1.0-bip)/(1.0-bip) * gbip;

////cout << "gf1 = (" <<gf1[0]<<","<<gf1[1]<<","<<gf1[2]<<")"<<endl;
//cout << "gf2 = (" <<gf2[0]<<","<<gf2[1]<<","<<gf2[2]<<")"<<endl;
//cout << "gbin = (" <<gbin[0]<<","<<gbin[1]<<","<<gbin[2]<<")"<<endl;
//cout << "gbip = (" <<gbip[0]<<","<<gbip[1]<<","<<gbip[2]<<")"<<endl;
        gpbi = (f2*gf1)*pbi + f1*f2 * gpbi;
        pbi *= f1*f2;
    }
    dist = res;
    for(int k=0; k<m; k++)
        gdist[k] = gres[k];
}

double ziestimator::getlogdensity
(int i, int a, int m, int v, /*bool istrade,*/ vector<double>& grad)
{
    if(debug)
        cout << "getlogdensity(" << i << "," << a
             << "," << m << "," << v << ")"
             << " t=" << analysis.X(i,0).t << endl;
    int q = params.size();
    zirec& last=analysis.X(i,1);
    if(a > last.a)
    {
        if(m < 0)
            throw runtime_error("m<0 in getlogdensity");
        if(v < 1)
            throw runtime_error("v<1 in getlogdensity");

        vector<double> gnu(q);

        double nu = model.getnu(params,gnu);
        for(int k=0; k<q; k++)
            if(gnu[k] != 0)
                throw runtime_error("Parameters have to be constant in nu.");

        int mnu = (int) m*nu;

        int vnu = (int) v*nu;
        if(vnu == 0)
            vnu = 1;

        double apoint = 0;
        vector<double> gapoint = zero_vector(q);

        vector<double> adist = zero_vector(mnu+1);
        adist[0] = 1;
        vector<vector<double> > gadist(mnu+1);
        for(int k=0; k<mnu+1; k++)
            gadist[k] = zero_vector(q);

        int bn = mnu+vnu;
        vector<double> bdist = zero_vector(bn+1);

        bdist[0] = 1;
        vector<vector<double> > gbdist(bn+1);
        for(int k=0; k<=bn; k++)
            gbdist[k] = zero_vector(q);

        for(int p=last.a+1;; p++)
        {
            vector<double> g(q);

            double point;
            vector<double> gpoint(q);
            double bin;
            vector<double> gbin(q);
            double bip;
            vector<double> gbip(q);

            getomega(i, p, point, gpoint, bin, gbin, bip, gbip);

            int n = (int) (bin + 0.5);

            if(p<a)
            {
                apoint += point;
                gapoint = gapoint + gpoint;
                if(n>0)
                {
                    addbi(adist, gadist, n, bip, gbip);
                    if(debug)
                    {
                        cout.precision(2);
                        cout << endl << "addbi a mnu=" << mnu << " p=" << p
                             << " prob=" << bip
                             << " n=" << n << " adist=(";

                        for(unsigned int k=0; k<adist.size(); k++)
                            cout << adist[k] << " ";
                        cout << ")" << endl << "grad[0]=(";
                        for(int k=0; k<q; k++)
                            cout << gadist[0][k] << " ";
                        cout << ")" << endl << "gbip=(";
                        for(int k=0; k<q; k++)
                            cout << gbip[k] << " ";
                        cout << ")" << endl;

                    }
                }
            }
            else
            {
                addpo(adist, gadist, apoint, gapoint);
                if(debug)
                {
                    cout.precision(2);
                    cout << endl << "addpo a apoint=" << apoint
                         << " p=" << p << " adist=(";
                    for(unsigned int k=0; k<adist.size(); k++)
                        cout << adist[k] << " ";
                    cout << ")" << endl << "gadist[0]=(";
                    for(int k=0; k<q; k++)
                        cout << gadist[0][k] << " ";
                    cout << ")" << endl << "gapoint=(";
                    for(int k=0; k<q; k++)
                        cout << gapoint[k] << " ";

                    cout << ")" << endl;
                }

                addpo(bdist, gbdist, point, gpoint);
                if(debug)
                {
                    cout.precision(2);
                    cout << endl << "addpo b bpoint=" << point
                         << " p=" << p << " dist=(";
                    for(unsigned int k=0; k<bdist.size(); k++)
                        cout << bdist[k] << " ";
                    cout << ")" << endl << "gbdist[0]=(";
                    for(int k=0; k<q; k++)
                        cout << gbdist[0][k] << " ";
                    cout << ")" << endl;
                }

                if(n>0)
                {
                    addbi(bdist, gbdist,n,bip,gbip);
                    if(debug)
                    {
                        cout.precision(2);
                        cout << endl << "addbi b mnu=" << mnu << " p=" << p
                             << " pbi=" << bip
                             << " n=" << n << " dist=(";
                        for(unsigned int k=0; k<bdist.size(); k++)
                            cout << bdist[k] << " ";
                        cout << ")" << endl << "gbdist[0]=(";
                        for(int k=0; k<q; k++)
                            cout << gbdist[0][k] << " ";
                        cout << ")" << endl << endl;
                    }
                }

                break;
            }
        }

        double prob = 0;
        vector<double> gprob = zero_vector(q);
        for(int k=0; k<= mnu; k++)
        {
            prob += adist[k] * bdist[bn-k];
            gprob = gprob + adist[k] * gbdist[bn-k]
                    + bdist[bn-k] * gadist[k];
        }

        if(prob<=1e-40)
        {
            grad = zero_vector(q);
            return -0.001;
        }
        else
        {
            grad =1.0 / prob * gprob;
            return log(prob);
        }
    }
    else if(model.isinspread() && a < last.a && a> last.b)
    {
        double ksum = 0;
        vector<double> gsum(zero_vector(q));

        double ka;
        vector<double> ga(q);

        for(int p=last.b+1; p<last.a; p++)
        {
            vector<double> g(q);

            double k = model.getkappa(last.a, last.b, p, params, g);
            if(p==a)
            {
                ka=k;
                ga=g;
            }

            ksum += k;
            gsum = gsum + g;
        }

        for(int j=0; j<q; j++)
            grad[j] = ga[j] / ka - gsum[j] / ksum;
        return log(ka) - log(ksum);
    }
    else
    {
        grad = zero_vector(q);
        return 0;
    }
}

double ziestimator::logdensity(int i, vector<double>& grad)
{
    if(!params.size())
        throw runtime_error("params not set!");
    const zirec& rec = analysis.X(i,0);

    double res = analysis.twodimestimation
                 ? getlogdensity(i, rec.a, rec.m, rec.q, grad)
                 : getlogdensity(i, rec.a, grad);

    return res;
}



int ziestimator::getsamplesize()
{
    return analysis.getN();
}




void ziestimator::testgrads(int start, int num)
{
    try
    {
        cout << "Testing gradients, pars=";
        vector<double> saveparams = params;

        for(unsigned int h=0; h<params.size(); h++)
        {
            cout << params[h] << " ";
        }
        cout << endl;
        for(int i=start; i < start + num && i < analysis.getN(); i++)
        {
            vector<double> g(params.size());
            setparameters(saveparams);

//            debug = 1;
            double d = logdensity(i,g);
//            debug = 0;

            cout << "{f=" << d << " g=(";
            for(unsigned int h=0; h<g.size(); h++)
            {
                cout << g[h] << " ";
            }
            cout << ")}" << endl;
            for(unsigned int k=0; k<params.size(); k++)
            {
                vector<double> dg(params.size());

                cout << "#" << k+1 << ": ";
                for(double c = 0.0001; c >= 0.0000001; c /= 10.0)
                {
                    setparameters(saveparams);
                    params[k] -= c;
                    double x = logdensity(i,dg);
                    params[k] += 2*c;
                    double y = logdensity(i,dg);

                    cout << "(" << (d-x)/c -g[k] << "," << (y-d)/c -g[k] << ") ";
                }
                cout << endl;
            }
        }
    }
    catch(...)
    {
        cout << "error testing grads!" << endl;
    }
}
