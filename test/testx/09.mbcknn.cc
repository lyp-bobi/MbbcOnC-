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
            seglens = {600,900,1200,1800,2700,3600};
        }else {
            target = "glexpand.datas";
            seglens = {100,200,300,600,900,1200,1800,2700};
        }
        cerr<<"09,mbcknn, SBBs"<<endl;
        cerr<<"seglen: ";
        for(auto len:seglens){cerr<<len<<" ";}
        cerr<<endl;
        xStore x(target, testFileName(target), true);
        for(auto &qt:{900,3600}) {
            cerr << "qt is " << qt << endl;
            vector<xTrajectory> queries;
            fillQuerySet(queries, x, qt);
            cerr << "MBR" << endl;
//            for (auto len:seglens) {
//                MTQ q;
//                q.prepareTrees(&x, [&len](auto x) {
//                    xRTree *r = buildMBRRTreeWP(x, xTrajectory::OPTS, len);
////                r->m_bUsingSBBD=false;
//                    return r;
//                });
//                q.appendQueries(queries);
//                std::cerr << q.runQueries().toString();
//            }
            cerr << "MBC" << endl;
            for (auto len:seglens) {
                MTQ q;
                q.prepareTrees(&x, [&len](auto x) {
                    xRTree *r = buildMBCRTreeWP(x, xTrajectory::OPTS, len);
//                r->m_bUsingSBBD=false;
                    return r;
                });
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