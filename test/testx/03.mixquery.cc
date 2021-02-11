//
// Created by Chuang on 2020/12/10.
//
#include "testFuncs.h"
#include "random"

int main(){
    try {
        string target = "tdexpand.txt";
        double qts[] = {300,1800,3600,7200,10800};
        cerr<<"03,mix, with TB,STR,SBB1200,SBBF(600,900,1200,1800)"<<endl;
        xStore x(target, testFileName(target), true);
        vector<xTrajectory> queries;
        fillQuerySetRand(queries,x);
        {
            MTQ q;
            q.prepareTrees(&x, [](auto x) { return buildTBTreeWP(x); });
            q.appendQueries(queries,0);
            std::cerr << q.runQueries().toString();
        }
        {
            MTQ q;
            q.prepareTrees(&x, [](auto x) { return buildSTRTreeWP(x); });
            q.appendQueries(queries,0);
            std::cerr << q.runQueries().toString();
        }
        {
            MTQ q;
            q.prepareTrees(&x, [](auto x) { return buildMBCRTreeWP(x, xTrajectory::OPTS, 1800); });
            q.appendQueries(queries,0);
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
            q.appendQueries(queries,0);
            std::cerr << q.runQueries().toString();
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