//
// Created by Chuang on 2020/11/30.
//

#include "testFuncs.h"
using namespace std;
int main(){
    try {
        double queryLen[] = {1000};
        vector<xTrajectory> queries;
        xStore x("test", "D://TRI-framework/dumpedtraj.txt", true);
        MTQ qmt;
        qmt.prepareTrees(&x,
                [](IStorageManager* r){return buildTBTreeWP(r);});
        for (int i = 0; i < testtime; i++) {
            queries.emplace_back(qmt.m_stores[0]->randomSubtraj(queryLen[0]));
        }
        qmt.appendQueries(queries);
        cerr<<qmt.runQueries().toString();
    }catch (Tools::Exception &e) {
        cerr << "******ERROR******" << endl;
        std::string s = e.what();
        cerr << s << endl;
        return -1;
    }
    return 0;
//    for(int i=0;i<NUMTHREAD;i++){
//        queryInput q;
//        q.tree = trees[i];
//        thread th(kNNQueryBatchThread,q, &res[i]);
//    }
}