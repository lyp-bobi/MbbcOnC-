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
            seglens = {1800};
        }else {
            target = "glexpand.datas";
            seglens = {1800};
        }
        cerr<<"debug"<<endl;
        cerr<<"seglen: ";
        for(auto len:seglens){cerr<<len<<" ";}
        cerr<<endl;
        xStore x(target, testFileName(target), true);
        cerr<<"qt is " << 3600<<endl;
        for(int p=0;p<1000;p++) {
            vector<xTrajectory> queries;
            testtime = 1;
            NUMTHREAD = 1;
            fillQuerySet(queries, x, 3600);
//        xTrajectory tj;
//        tj.loadFromString("115.928240,39.716420,3456049.000000 115.928310,39.716470,3456349.000000 115.928290,39.716470,3456649.000000 115.928280,39.716450,3456949.000000 115.928290,39.716470,3457549.000000 115.928300,39.716490,3457849.000000 115.928270,39.716490,3458149.000000 115.928270,39.716450,3458449.000000 115.928280,39.716480,3458749.000000 115.928250,39.716450,3459049.000000 115.928270,39.716480,3459349.000000 115.928290,39.716450,3459649.000000");
//                queries.emplace_back(tj);
            cerr << queries[0].toString()<<"\n";
            partstoreskip = false;
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
            partstoreskip = true;
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