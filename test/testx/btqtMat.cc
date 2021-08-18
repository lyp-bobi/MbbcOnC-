//
// Created by Chuang on 2020/12/10.
//
#include "testFuncs.h"
#include "random"

int main(int argc,char *argv[]){
    try {
        vector<double> seglens;
        string target;
        if(argc==1) {
            target = "tdexpand.datas";
            seglens = {600,900,1200,1800,2700,3600,5400,7200,9000};
        }else {
            target = "glexpand.datas";
            seglens = {100,200,300,600,900,1200,1800,2700};
        }
//        testtime=400;
        cerr<<"seglen: ";
        for(auto len:seglens){cerr<<len<<" ";}
        cerr<<"TB STR ";
        cerr<<endl;
        xStore x(target, testFileName(target), true);
        for(double qt=300;qt<=10800;qt+=500) {
            cerr<<"qt is " << qt<<endl;
            vector<xTrajectory> queries;
//            xTrajectory tj;
//            tj.loadFromString("116.502520,40.007630,3921202.717629 116.502520,40.007631,3921502.717629");
//            queries.emplace_back(tj);
            fillQuerySet(queries,x,qt);

            for (auto len:seglens) {
                MTQ q;
                q.prepareTrees(&x, [&len](auto x) { return buildMBCRTreeWP(x, xTrajectory::OPTS, len); });
                q.appendQueries(queries);
                std::cerr << q.runQueries().toString();
            }
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