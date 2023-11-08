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
            seglens = {600};
        }else {
            target = "glexpand.datas";
            seglens = {600};
        }
        double qt = 3600;
        cerr<<"02, e2ennk, TB,STR, and SBB 600";
        xStore x(target, testFileName(target), true);
        vector<xTrajectory> queries;
        fillQuerySet(queries,x,qt);
        for(current_distance = IED; current_distance <= RMDTW; current_distance = supported_distance(current_distance + 1)) {
            for (int k = 10; k < 100; k += 20) {
                cerr << "k is " << k << endl;
                for (auto len:seglens) {
                    MTQ q;
                    q.prepareTrees(&x, [&len](auto x) {
                        return buildMBCRTreeWP(x, xTrajectory::GSS, len);
                    });
                    q.appendQueries(queries, k);
                    std::cerr << q.runQueries().toString();
                }
                {
                    MTQ q;
                    q.prepareTrees(&x,
                                   [](auto x) { return buildSTRTreeWP(x); });
                    q.appendQueries(queries, k);
                    std::cerr << q.runQueries().toString();
                }
                {
                    MTQ q;
                    q.prepareTrees(&x, [](auto x) { return buildTBTreeWP(x); });
                    q.appendQueries(queries, k);
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