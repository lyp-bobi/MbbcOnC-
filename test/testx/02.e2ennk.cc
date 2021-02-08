//
// Created by Chuang on 2020/12/10.
//
#include "testFuncs.h"
#include "random"

int main(){
    try {
        string target = "tdfilter.txt";
        double qt = 3600;
        double seglens[] = {600,900,1200,1800,2700};
        cerr<<"seglen: ";
        for(auto len:seglens){cerr<<len<<" ";}
        cerr<<endl;
        xStore x(target, testFileName(target), true);
        vector<xTrajectory> queries;
        fillQuerySet(queries,x,qt);
        for(int k = 5;k<200;k+=10) {
            cerr<<"k is " << k<<endl;
            {
                MTQ q;
                q.prepareTrees(&x, [](auto x) { return buildTBTreeWP(x); });
                q.appendQueries(queries,k);
                std::cerr << q.runQueries().toString();
            }
            {
                MTQ q;
                q.prepareTrees(&x, [](auto x) { return buildSTRTreeWP(x); });
                q.appendQueries(queries,k);
                std::cerr << q.runQueries().toString();
            }
            {
                MTQ q;
                q.prepareTrees(&x, [](auto x) { return buildMBCRTreeWP(x, xTrajectory::ISS, 1200); });
                q.appendQueries(queries,k);
                std::cerr << q.runQueries().toString();
            }
            {
                MTQ q;
                SBBFMAP lens;
                lens[make_pair(0,1500)]=600;
                lens[make_pair(1500,3000)]=900;
                lens[make_pair(3000,4500)]=1200;
                lens[make_pair(4500,1e300)]=1800;
                q.prepareForest(&x,lens);
                q.appendQueries(queries,k);
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