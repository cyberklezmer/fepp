#include "hfd.hpp"


void importertw::import(const char** stocks, const char** markets,
        date& start, date& end, string& folder)
{
    vector<smpair> pairs;

    for(int s=0; stocks[s]; s++)
        for(int m=0; markets[m]; m++)
        {
            smpair p(stocks[s],markets[m],9.5,15.5);
            pairs.push_back(p);
        }
    tickdatatwimporter i(start,end,pairs,tickdataimporter::usstock,folder);
    i.process();
}

