//
// Created by Chuang on 2020/12/10.
//

#include "testFuncs.h"
#include "random"

int main(){
    try {
        string target = "tdexpand.txt";
        double qts[] = {300,900,1500,2100,2700,3300,3900,4500};
        double seglens[] = {600,900,1500,2100,2700,3300,3900,4500};
        cerr<<"seglen: ";
        for(auto len:seglens){cerr<<len<<" ";}
        cerr<<endl;
        xStore x(target, testFileName(target), true);
        for(auto qt:qts) {
            cerr<<"qt is " << qt<<endl;
            vector<xTrajectory> queries;
            fillQuerySet(queries,x,qt);
            {
                MTQ q;
                q.prepareTrees(&x, [](auto x) { return buildTBTreeWP(x); });
                q.appendQueries(queries);
                std::cerr << q.runQueries().toString();
            }
            {
                MTQ q;
                q.prepareTrees(&x, [](auto x) { return buildSTRTreeWP(x); });
                q.appendQueries(queries);
                std::cerr << q.runQueries().toString();
            }
            for (auto len:seglens) {
                MTQ q;
                q.prepareTrees(&x, [&len](auto x) { return buildMBCRTreeWP(x, xTrajectory::ISS, len); });
                q.appendQueries(queries);
                std::cerr << q.runQueries().toString();
            }
        }
        cerr<<"mission complete.\n";
    }catch (Tools::Exception &e) {
        cerr << "******ERROR******" << endl;
        std::string s = e.what();
        cerr << s << endl;
        return -1;
    }
    return 0;
}