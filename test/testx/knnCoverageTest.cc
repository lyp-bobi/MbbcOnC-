//
// Created by Chuang on 2020/11/24.
//

#include "testFuncs.h"
int main(){
    xStore x("test", "/root/dumpedtraj.txt",true,true);

    double segLens[]={ 800, 1600, 10000};
    double queryLen[]={1000};
    vector<vector<xTrajectory>> querySets;
    for(auto ql:queryLen)
    {
        querySets.emplace_back(vector<xTrajectory>());
        for(int i=0;i<10;i++){
            querySets.back().emplace_back(x.randomSubtraj(ql));
        }
    }
    MyVisitor vis;
    try {
        for(auto sl:segLens) {
            CUTFUNC f=[](auto tj){return xTrajectory::OPTS(tj,tjstat->bt);};
            for (auto qs:querySets) {
                bUsingSBBD = false;
                for (auto querylen:segLens) {
                    xRP r(buildMBCRTreeWP(&x, xTrajectory::ISS,sl));
                    kNNQueryBatch(r, qs, &x, 5);
                    vis.clear();
                    xRP r2(buildMBRRTreeWP(&x, xTrajectory::ISS,sl));

                    kNNQueryBatch(r2, qs, &x, 5);
                    vis.clear();
                }
                {
                    xRP r(buildTBTreeWP(&x));

                    kNNQueryBatch(r, qs, &x, 5);
                    vis.clear();
                }
                {
                    xRP r(buildSTRTreeWP(&x));

                    kNNQueryBatch(r, qs, &x, 5);
                    vis.clear();
                }
                bUsingSBBD = true;
                for (auto querylen:segLens) {
                    tjstat->bt = querylen;
                    xRP r(buildMBCRTreeWP(&x, xTrajectory::ISS,sl));

                    kNNQueryBatch(r, qs, &x, 5);
                    vis.clear();
                    xRP r2(buildMBRRTreeWP(&x, xTrajectory::ISS,sl));

                    kNNQueryBatch(r2, qs, &x, 5);
                    vis.clear();
                }
                {
                    xRP r(buildTBTreeWP(&x));

                    kNNQueryBatch(r, qs, &x, 5);
                    vis.clear();
                }
                {
                    xRP r(buildSTRTreeWP(&x));

                    kNNQueryBatch(r, qs, &x, 5);
                    vis.clear();
                }
            }
        }
    }
    catch (Tools::Exception &e) {
        cerr << "******ERROR******" << endl;
        std::string s = e.what();
        cerr << s << endl;
        return -1;
    }
}