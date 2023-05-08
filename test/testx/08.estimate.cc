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
            seglens={100,200,400,800,1600,3200,6400,12800};
        }else {
            target = "glexpand.datas";
            seglens={100,200,400,800,1600,3200,6400,12800};
        }
        cerr<<"08,estimate, SBBs"<<endl;
        cerr<<"seglen: ";
        for(auto len:seglens){cerr<<len<<" ";}
        cerr<<endl;
        xStore x(target, testFileName(target), true);
        cerr<<"qt is " << 3600<<endl;
        vector<xTrajectory> queries;
        fillQuerySet(queries,x,3600);
//        for (auto len:seglens) {
//            MTQ q;
//            q.prepareTrees(&x, [&len](auto x) {
//                xRTree* r =buildMBRRTreeWP(x, xTrajectory::GSS, len);
////                r->m_bUsingSBBD=false;
//                return r;
//            });
//            q.appendQueries(queries);
//            std::cerr << q.runQueries().toString();
//        }
        for (auto len:seglens) {
            MTQ q;
            q.prepareTrees(&x, [&len](auto x) {
                xRTree* r =buildMBCRTreeWP(x, xTrajectory::GSS, len);
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