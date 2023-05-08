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
            seglens = {600};//600,900,1200,1800,2700,3600,5400,7200,9000
        }else {
            target = "glexpand.datas";
            seglens = {600};
        }
//        testtime=400;
        cerr<<"seglen: ";
        for(auto len:seglens){cerr<<len<<" ";}
        cerr<<"TB STR ";
        cerr<<endl;
        xStore x(target, testFileName(target), true);
        for(current_distance = IED; current_distance <= RMDTW; current_distance = supported_distance(current_distance + 1)) {
            for (double qt = 1200; qt <= 7200; qt += 600) {
                cerr << "qt is " << qt << endl;
                vector<xTrajectory> queries;
                fillQuerySet(queries, x, qt);

                for (auto len:seglens) {
                    MTQ q;
                    q.prepareTrees(&x, [&len](auto x) {
                        xRTree* r = buildMBCRTreeWP(x, xTrajectory::GSS, len);
                        r->m_bUsingLoadleaf = false;
                        return r;
                    });
                    q.appendQueries(queries);
                    std::cerr << q.runQueries().toString();
                }

                for (auto len:seglens) {
                    MTQ q;
                    q.prepareTrees(&x, [&len](auto x) {
                        return buildMBCRTreeWP(x, xTrajectory::GSS, len);
                    });
                    q.appendQueries(queries);
                    std::cerr << q.runQueries().toString();
                }
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