//
// Created by Chuang on 2020/11/30.
//

#include "testFuncs.h"
using namespace std;
int main(){
    vector<xStore*> stores;
    vector<xRTree*> trees;
    vector<queryInput> queries;
    vector<queryRet> res;
    vector<thread> ths;
    try {
        double queryLen[] = {1000};
        vector<xTrajectory> queries;
        MTQ qmt;
        qmt.prepareTrees([](){return new xStore("test", "D://TRI-framework/dumpedtraj.txt", true);},
                [](IStorageManager* r){return buildTBTreeWP(r);});
        for (int i = 0; i < 100; i++) {
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