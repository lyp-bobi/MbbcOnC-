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
            seglens = {600};
        }
        double qt =3600;
        std::cerr<<testFileName(target)<<endl;
        xStore x(target, testFileName(target), true);
        for(auto bt:seglens) {
            std::cerr<<"bt is "<<bt <<endl;
            vector<xTrajectory> queries;
            fillQuerySet(queries, x, qt);
            {
                MTQ q;
                q.prepareTrees(&x, [&bt](auto x) { return buildMBCRTreeWP(x, xTrajectory::ISS, bt, "ISS"); });
                q.appendQueries(queries);
                std::cerr << q.runQueries().toString();
            }
            {
                MTQ q;
                q.prepareTrees(&x, [&bt](auto x) { return buildMBCRTreeWP(x, xTrajectory::GSS, bt, "GSS"); });
                q.appendQueries(queries);
                std::cerr << q.runQueries().toString();
            }
            {
                MTQ q;
                q.prepareTrees(&x, [&bt](auto x) { return buildMBCRTreeWP(x, xTrajectory::OPTS, bt, "OPTS"); });
                q.appendQueries(queries);
                std::cerr << q.runQueries().toString();
            }
            {
                MTQ q;
                q.prepareTrees(&x,
                               [&bt](auto x) { return buildMBCRTreeWP(x, xTrajectory::FP, bt / tjstat->tl, "FP"); });
                q.appendQueries(queries);
                std::cerr << q.runQueries().toString();
            }
            {
                MTQ q;
                q.prepareTrees(&x,
                               [&bt](auto x) { return buildMBCRTreeWP(x, xTrajectory::RDP, bt, "RDP"); });
                q.appendQueries(queries);
                std::cerr << q.runQueries().toString();
            }
        }
    }catch (Tools::Exception &e) {
        cerr << "******ERROR******" << endl;
        std::string s = e.what();
        cerr << s << endl;
        return -1;
    }
    return 0;
}