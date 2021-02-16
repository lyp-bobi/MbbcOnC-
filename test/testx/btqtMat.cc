//
// Created by Chuang on 2020/12/10.
//
#include "testFuncs.h"
#include "random"

int main(int argc,char *argv[]){
    try {
        string target;
        if(argc==1) {
            target = "tdexpand.data";
        }else {
            target = "glexpand.data";
        }
        testtime=400;
        double seglens[] = {300,600,900,1200,1800,2700,3600};
        cerr<<"seglen: ";
        for(auto len:seglens){cerr<<len<<" ";}
        cerr<<endl;
        xStore x(target, testFileName(target), true);
        for(double qt=300;qt<=400;qt+=300) {
            cerr<<"qt is " << qt<<endl;
            vector<xTrajectory> queries;
            xTrajectory tj;
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
                q.prepareTrees(&x, [&len](auto x) { return buildMBCRTreeWP(x, xTrajectory::OPTS, len); });
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