//
// Created by Chuang on 2020/12/10.
//
#include "testFuncs.h"
#include "random"

int main(){
    try {
        string target = "tdexpand.data";
//        string target = "glexpand.data";
        cerr<<"09,mbcknn, SBBs"<<endl;
        double seglens[] = {600,900,1200,1800,2700};
        cerr<<"seglen: ";
        for(auto len:seglens){cerr<<len<<" ";}
        cerr<<endl;
        xStore x(target, testFileName(target), true);
        cerr<<"qt is " << 3600<<endl;
        vector<xTrajectory> queries;
        fillQuerySet(queries,x,3600);
        for (auto len:seglens) {
            MTQ q;
            q.prepareTrees(&x, [&len](auto x) {
                xRTree* r =buildMBRRTreeWP(x, xTrajectory::OPTS, len);
//                r->m_bUsingSBBD=false;
                return r;
            });
            q.appendQueries(queries);
            std::cerr << q.runQueries().toString();
        }
        for (auto len:seglens) {
            MTQ q;
            q.prepareTrees(&x, [&len](auto x) {
                xRTree* r =buildMBCRTreeWP(x, xTrajectory::OPTS, len);
//                r->m_bUsingSBBD=false;
                return r;
            });
            q.appendQueries(queries);
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