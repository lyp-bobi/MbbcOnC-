//
// Created by Chuang on 2021/2/17.
//

#include "testFuncs.h"
#include "random"

int main(int argc,char *argv[]){
    try {
        string target;
        vector<double> seglens;
        if(argc==1) {
            target = "tdexpand.datas";
            seglens = {1200};
        }else {
            target = "glexpand.datas";
            seglens = {1200};
        }
        cerr<<"seglen: ";
        for(auto len:seglens){cerr<<len<<" ";}
        cerr<<endl;
        cerr<<"12,sbbdknn"<<target<<endl;
        xStore x(target, testFileName(target), true);
        for(current_distance = IED; current_distance <= RMDTW; current_distance = supported_distance(current_distance + 1)) {
            for (double qt = 300; qt <= 5300; qt += 600) {
                cerr << "qt is " << qt << endl;
                vector<xTrajectory> queries;
                fillQuerySet(queries, x, qt);
                for (auto len:seglens) {
                    MTQ q;
                    q.prepareTrees(&x, [&len](auto x) {
                                       xRTree *r = buildMBCRTreeWP(x, xTrajectory::GSS, len);
                                       r->m_bUsingSBBD = false;
                                       return r;
                                   }
                    );
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
                for (auto len:seglens) {
                    MTQ q;
                    q.prepareTrees(&x, [&len](auto x) {
                                       xRTree *r = buildMBCRTreeWP(x, xTrajectory::GSS, len);
                                       r->m_bUsingLoadleaf = false;
                                       return r;
                                   }
                    );
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