//
// Created by Chuang on 2020/12/10.
//
#include "testFuncs.h"
#include "random"

int main(int argc,char *argv[]){
    try {
        string target = "tdfilter.txt";
        double qt = 3600;
        cerr<<"02, e2ennk, TB,STR, and SBB1800";
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
                q.prepareTrees(&x, [](auto x) { return buildMBCRTreeWP(x, xTrajectory::OPTS, 1800); });
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