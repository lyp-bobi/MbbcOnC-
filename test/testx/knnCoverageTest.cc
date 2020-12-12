//
// Created by Chuang on 2020/11/24.
//

#include "testFuncs.h"
int main(){
    xStore x("test", "D://TRI-framework/dumpedtraj.txt",true,true);

    double segLens[]={50, 100, 200, 400, 800, 1600, 10000};
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
                    auto r = buildMBCRTreeWP(&x, f);
                    kNNQueryBatch(r, qs, &x, 5);
                    vis.clear();
                    auto r2 = buildMBRRTreeWP(&x, f);

                    kNNQueryBatch(r2, qs, &x, 5);
                    vis.clear();
                }
                {
                    auto r = buildTBTreeWP(&x);

                    kNNQueryBatch(r, qs, &x, 5);
                    vis.clear();
                }
                {
                    auto r = buildSTRTreeWP(&x);

                    kNNQueryBatch(r, qs, &x, 5);
                    vis.clear();
                }
                bUsingSBBD = true;
                for (auto querylen:segLens) {
                    tjstat->bt = querylen;
                    auto r = buildMBCRTreeWP(&x, f);

                    kNNQueryBatch(r, qs, &x, 5);
                    vis.clear();
                    auto r2 = buildMBRRTreeWP(&x, f);

                    kNNQueryBatch(r2, qs, &x, 5);
                    vis.clear();
                }
                {
                    auto r = buildTBTreeWP(&x);

                    kNNQueryBatch(r, qs, &x, 5);
                    vis.clear();
                }
                {
                    auto r = buildSTRTreeWP(&x);

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