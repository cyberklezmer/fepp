/*
 * analysis.cpp
 *
 *  Created on: 14 Jun 2010
 *      Author: martin
 */

#include "zimodels.hpp"
#include <ctime>
#include <iostream>



void zianalysis::onstartpair(smpair& togo)
{
    cout << togo.stock << " at " << togo.market << endl;

    qchanges = new zirec[maxqchanges];
    ajumps = new int[maxajumps];

    qdist = new empdistn;
    M=0;

    mina = 1000000000L;
    maxa = 0;

    firstt = -1;
    sp=0;
    sda = 0;
    nsda=0;
    ss = 0;
    nss = 0;

}

void zianalysis::onrecord(hfdrecord& rec, int recn)
{
    if(M < maxqchanges && sp < maxajumps)
    {
        int hs = getqhistorysize();

        int s = rec.q >= 0 ? 0 : (unitvolume ? 1 : -rec.q);
//        qdist->add(s);

        bool firsttoday = (hs == 0);
        if(M==0 && !firsttoday)
            throw logic_error("M==0 && !firsttoday");
        bool sameaslast;
        if(firsttoday)
            sameaslast = false;
        else
        {
            if(M==0)
                throw logic_error("M==0 (1)");
            sameaslast = (rec.a == qchanges[M-1].a && rec.b == qchanges[M-1].b);
        }
        if(sameaslast) // if trade records were merged then the previous folume would be unknown
        {
            if(M==0)
                throw logic_error("M==0 (1)");
            qchanges[M-1].q = rec.an;

            qchanges[M-1].s = s;
        }
        else
        {
            bool sametime;
            if(firsttoday)
                sametime = false;
            else
            {
                if(M==0)
                    throw logic_error("M==0 (2)");
                sametime = (qchanges[M-1].t == rec.t);
            }
            if(sametime)
            {
               if(M==0)
                    throw logic_error("M==0 (2)");

                zirec& r = qchanges[M-1];
                r.a = rec.a;
                r.b = rec.b;
                r.t = rec.t;
                r.q = rec.an;
                r.s = s;
            }
            else
            {
                zirec& r = qchanges[M];
                r.firstinday = firsttoday;
                r.a = rec.a;
                r.b = rec.b;
                r.t = rec.t;
                r.q = rec.an;
                r.s = s;
                r.m = -1;


                if(r.firstinday)
                    firstt = r.t;
                if(firstt<0)
                    throw runtime_error("first qchange not firstinday");
                if( !r.firstinday && r.a != qchanges[M-1].a  && r.t - firstt > warmuptime )
                {
                    if(!onlytrades || r.s > 0) // tbd check of tjere are recprds vopůagom s-q >=0
                    {
                        r.m = r.s - qchanges[M-1].q;
                        r.m = r.m >= 0 ? r.m : 0;
                        int da = r.a-qchanges[M-1].a;
                        if(da > 0)
                        {
                            sda += da;
                            nsda++;
                            ajumps[sp++] = M;
                        }
                        if(r.s > 0)
                        {
                            ss += r.s;
                            nss++;
                        }

                        if(r.a < mina)
                            mina = r.a;
                        if(r.a > maxa)
                            maxa = r.a;

                    }
                }

                M++;
            }
        }
    }
}

void logrealtime()
{
    time_t t = time(0);   // get time now
    struct tm * now = localtime( & t );
    cout << (now->tm_hour) << ':'
         << (now->tm_min)
         << endl;
}

void zianalysis::estimate(zimodel& model, vector<paramresult>& res,
                          bool catchex, bool& failed,
                          bool& significant, double& ll, string& emsg)
{
    cout << "Estimating model " << model.getname() << endl;
    zimleestimator e(model,*this);
    e.setmaxtime(maxmletime);
    e.setlogging(extendedlogging ? object::extendedlogging :
                 object::basiclogging);

    if(0)
    {
        vector<double> pars(2);
        vector<double> g(2);
        pars[0] = 2.09518;
        pars[1] = 0.0984748;
        e.setparameters(pars);
        e.testgrads(274,5);
        throw;
    }

    failed = false;
    emsg = "";
    bool grad;

    logrealtime();
    try
    {
        try
        {
            e.setgradient(grad = true);
            ll = e.estimate(res);
        }
        catch(const runtime_error& e)
        {
            emsg = "nograd: ";
            emsg += e.what();
            cout << emsg << endl;
            throw;
        }
    }
    catch(...)
    {
        try
        {
            e.setgradient(grad = false);
            e.setxtolrel(10e-5);
            ll = e.estimate(res);
        }
        catch(...)
        {

            cout << "Estimation failed. Gradient = " << grad
                 << endl;
            cout << "Grads:" << endl;
            e.testgrads(1,5);
            if(!catchex)
            {
                latex << "Estmation failed" << endl;
                throw;
            }
            failed = true;
            if(emsg == "")
                emsg = "Failed";
        }
    }
    logrealtime();
    significant = true;
    for(unsigned int i=0; i<res.size(); i++)
    {
        if(res[i].z() != na)
        {
            double q = fabs(res[i].z());
            if(q<1.645)
            {
                significant = false;
                break;
            }
        }
        else
        {
            significant = false;
            break;
        }
    }
}


double zianalysis::evaluate(zimodel* modelused, const vector<paramresult>& res,
                            int from, int len, double avgdaplus)
{
    double result = na;
    try
    {
        try
        {
            latex.precision(3);
            vector<double> pars(res.size());
            for(unsigned int i=0; i<pars.size(); i++)
                pars[i] = res[i].value;

            zimleestimator e2(*modelused,*this);
            e2.setparameters(pars);
            double SSE = 0;
            double MAE = 0;
            double S = 0;
            double S2 = 0;
            vector<double> smpl;
            int nobs = 0;

            csv << "err,dt,da,t,s" << endl;

            clock_t start = clock();

            for(int i=from; i<from+len; i++)
            {
                clock_t now = clock();
                if((now-start) / (double)CLOCKS_PER_SEC > maxmletime)
                {
                    cout << "Evaluation time exceede - only "
                         << i-from << " observation used" << endl;
                    break;
                }
                const zirec& last=X(i,1);
                const zirec& r=X(i,0);
                double da=r.a - last.a;
                if(da > 0)
                {
                    smpl.push_back(da);
                    S += da;
                    S2 += da*da;

                    double E = 0;
                    double sump=0.0;
                    vector<double> phoo(pars.size());
                    int pd=1;
                    double sumomega=0;
                    for(; sump<0.99; pd++)
                    {
                        int p = last.a+pd;
                        double omega =  e2.getomega(i,p,phoo);

                        double pi = exp(-sumomega)*(1.0 - exp(-omega));

                        sumomega += omega;
                        // exp(e2.getlogdensity(i,p,r.m,phoo));
                        E += pi*pd;
                        sump+=pi;
                        if(pd == 1 && pi==0)
                            throw runtime_error("zero first probability");
                        if(p>4*maxa)
                        {
                            cout << "Warning: nobs=" << nobs
                                 << " cum. probability=" << sump << endl;
                            break;
                        }
                    }
                    SSE += (da-E)*(da-E);
                    MAE += fabs(da-E);
                    if(r.t-last.t > 0)
                        csv << fabs(da-E) << "," << r.t - last.t  << ","
                            << da << "," << r.t << "," << r.s << endl;
                    nobs++;
                }
            }
            e2.resetparameters();

            double SST = S2 - S*S / (double) nobs;
            double P2 = 1-SSE/SST;
            double MAT = 0;

            double avg = avgdaplus == 0
              ? (double) S / (double) nobs : avgdaplus;
            for(unsigned int k=0; k<smpl.size(); k++)
                MAT += fabs(smpl[k]-avg);

            double P1 = 1 - MAE / MAT;
            cout << "P1=" << P1 << " P2=" << P2 << endl;
            //latex << "$\\bar a=" << S/(double) nobs << "$," << endl;

            if(extendedoutput)
            {
                latex << "$P_1 = " << P1 << "$";
                latex << ", $P_2=" << P2 << "$";
            }
            result = P1;
        }
        catch(const runtime_error& e)
        {
            cout << e.what() << endl;
            latex << e.what() << endl;
        }
    }
    catch(...)
    {
        cout << "Exception occured when evaluating" << endl;
        latex << "Exception occured when evaluating" << endl;
    }
    return result;
}



void zianalysis::onendpair(smpair& togo)
{

    N = sp * offsampleratio / (offsampleratio + 1);
    Noff = sp-N;

    latex << "{ \\normalsize "
          << togo.stock << " at " << togo.market
          << "}" << endl ;

    latex << "\\\\" << endl << endl << "\\smallskip" << endl;
    if(extendedoutput)
        latex << "$N/M/maxq=" << getN() << "/" << M
          << "/" << maxqchanges << "$";
    else
        latex << "$N=" << getN() << "$, " << endl;

    csv << togo.stock << " at " << togo.market << endl;

    double avgdaplus = nsda ? sda / nsda : 0;
    latex << "$\\overline{a}^+=" << avgdaplus
          << "$, $\\overline q=" << (nss ? ss / nss : 0)
          << "$\\\\" << endl << endl << "\\smallskip" << endl;

    vector<paramresult> result;
    double loglik = -HUGE_VAL;
    zimodel* modelused = 0;
    int nofmodelused = 0;
    bool tailmodelused;
    bool anysignificant = false;

    cout << "Onlytrades = " << onlytrades << ", unitvolume = " << unitvolume
         << endl;

    int n=firstn;
    for(;;)
    {
        cout << "Estimating: n=" << n <<  endl;

        zimodel* model=0;
        try
        {
            switch(modeltype)
            {
            case individual:
                model = new zinparmodel(n,phipar, includenu,includegamma,
                                        includeeta, includezeta);
                break;
            case tail:
                if(n==1)
                    model = new zinparmodel(n,phipar, includenu,includegamma,
                                            includeeta, includezeta);
                else
                {
                    zinptmodel* npmod;
                    model = npmod = new zinptmodel(n-1,phipar,includenu,includegamma,
                                                   includeeta, includezeta);
                    if(modelused)
                    {
                        zinptmodel* mod = (zinptmodel*) modelused;
                        vector<double> p(model->getparams().size());
                        for(unsigned int i=0; i<p.size(); i++)
                            p[i] = model->getparams()[i].initial;
                        if(nofmodelused==1)
                        {
                            for(int k=0; k<npmod->getn(); k++)
                            {
                                p[npmod->firstkappaorphi()+k] = result[0].value;
                                p[npmod->firstrho()+k] = result[1].value;
                            }
                        }
                        else
                        {
                            int m=0;
                            for(; m<nofmodelused-1; m++)
                                p[npmod->firstkappaorphi()+m]
                                    = result[mod->firstkappaorphi()+m].value;
                            p[npmod->firstkappaorphi()+m]
                                = result[mod->firstkappaorphi()+m-1].value;
                            int k=0;
                            for(; k<nofmodelused-1; k++)
                                p[npmod->firstrho()+k]
                                    = result[mod->firstkappaorphi()+k].value;
                            p[npmod->firstkappaorphi()+k]
                                = result[mod->firstkappaorphi()+k-1].value;
                            p[npmod->alphakappaorphi()]
                                = result[mod->alphakappaorphi()].value;
                            p[npmod->alpharho()]
                                = result[mod->alpharho()].value;
                        }
                        model->setinitparams(p);
                    }
                }
                break;
            }

            double ll;
            bool failed;
            bool significant;
            vector<paramresult> r;
            string emsg;
            estimate(*model, r, true, failed, significant, ll, emsg);

            if(extendedoutput)
                latex << model->getname();
            if(failed)
            {
                latex << ": " << emsg << "\\\\" << endl;
                delete model;
                model = 0;
                ll = loglik;
            }
            else
            {
                if(significant)
                    anysignificant = true;
                if(n==1 || extendedoutput)
                {
                    if(significant)
                    {
                        latex << "$P_0="
                            << evaluate(model,r, N, Noff, avgdaplus)
                            << "$ ";
                    }
                }
            }
            double ratio = 2*(ll-loglik);
            cout << "ratio=" << ratio << endl;

            if(ratio < 5.991)// 3.841
            {
                if(modelused == 0)
                    throw runtime_error("internal error: ratio of first model");
                else
                {
                    delete model;
                    break;
                }
            }
            else
            {
                if(significant)
                {
                    result = r;
                    loglik = ll;
                    if(modelused)
                        delete modelused;
                    modelused = model;
                    nofmodelused = n;
                    tailmodelused = n>1 ? true : false;
                }
                if(anysignificant && !significant)
                    break;
                if(n==maxn)
                    break;
                else
                    n++;
            }
            loglik = ll;
        }
        catch(...)
        {
            if(model)
                delete model;
            throw;
        }
    }

//    cout << "Evaluating" << endl;

//    evaluate(modelused, result);

    if(modelused)
    {
        if(!extendedoutput && tailmodelused)
        {
            cout << "Evaluating..." << endl;
            if(tailmodelused)
            {
                if(avgdaplus < 1.02 && !onlytrades)
                    latex << ", $P_=n/a$" << endl;
                else
                    latex << "$P="
                      << evaluate(modelused,result, N, Noff, avgdaplus)
                      << "$" << endl;
            }
        }
        if(resample)
        {
            vector<double> pars(result.size());
            cout << "Simulating pars: ";
            for(unsigned int k=0; k<result.size(); k++)
            {
                pars[k] = result[k].value;
                cout << pars[k] << " ";
            }
            cout << endl;

            modelused->setinitparams(pars);

            int resnum;
            cdasimulator s(*modelused,Noff*10,Noff);
            s.simulate(Noff,qchanges,resnum);
            int m=0;
            for(int i=0; i<resnum; i++)
            {
                zirec& r = qchanges[i];
                if( !r.firstinday && r.a != qchanges[i-1].a )
                {
                   r.m = r.s - qchanges[M-1].q;
                   r.m = r.m >= 0 ? r.m : 0;
                   int da = r.a-qchanges[M-1].a;
                   if(da > 0)
                   {
                      ajumps[m++] = i;
                      if(m == Noff)
                        break;
                   }
                }
            }

            cout << "maxqchanges=" << maxqchanges << " successfully got="
                << resnum << endl;
            latex << ", $P_{\\mathrm{sim}}="
                << evaluate(modelused,result, 0, m, 0)
                << "$";
        }

        cout << "Writing latex..." << endl;
        latex << "\\\\" << endl << endl
              << " \\smallskip" << endl;

        if(extendedoutput)
            latex << modelused->getname() << "\\\\" << endl;

        for(unsigned int i=0; i< result.size(); i++)
        {
            latex.precision(4);
            result[i].output(latex,true);
            //            if((i%2))
            latex << "\\\\" ;
        }
        latex << "\\smallskip" << endl;
        if(tailmodelused && phipar) // i.e. the model is with tails
        {
            zinptmodel* m = (zinptmodel*) modelused;
            double alphak = result[m->alphakappaorphi()].value
                              + result[m->alpharho()].value;
            double s = sqrt(result[m->alphakappaorphi()].std
                              + result[m->alpharho()].std);
            latex << endl << "\\smallskip" << endl << endl
                  << "$\\alpha_{\\kappa}=" << alphak
                  << "(" << s << ")^{"
                  << paramresult::stars(alphak/s)
                  << "}$\\\\"<< endl;
            latex << "\\smallskip" << endl;
        }

        delete modelused;
    }
}




