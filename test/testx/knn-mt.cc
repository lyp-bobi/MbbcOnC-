//
// Created by Chuang on 2020/11/30.
//

#include "testFuncsMultithread.h"
using namespace std;
int main(){
    vector<xStore*> stores;
    vector<xRTree*> trees;
    vector<queryInp> queries;
    vector<queryRet> res;
    vector<thread> ths;
    try {
        xStore x("test", "D://TRI-framework/dumpedtraj.txt", true);
        //buildTBTreeWP(&x);
        x.flush();
        double queryLen[] = {1000};
        for(int i=0;i<NUMTHREAD;i++){
            stores.emplace_back(new xStore("test", "D://TRI-framework/dumpedtraj.txt", true));
            trees.emplace_back(buildTBTreeWP(stores.back()));
            res.emplace_back(queryRet());
            queryInp q;
            q.tree = trees.back();
            for (int i = 0; i < 10; i++) {
                q.knn_queries.emplace_back(x.randomSubtraj(queryLen[0]));
            }
            queries.emplace_back(q);
        }
        for(int i=0;i<NUMTHREAD;i++) {
            ths.emplace_back(thread(kNNQueryBatchThread, queries[i], &res[i]));
        }
        for(int i=0;i<NUMTHREAD;i++) {
            ths[i].join();
            cerr<<res[i].toString();
        }
        cerr<<average(res).toString();
    }catch (Tools::Exception &e) {
        cerr << "******ERROR******" << endl;
        std::string s = e.what();
        cerr << s << endl;
        return -1;
    }
    return 0;
//    for(int i=0;i<NUMTHREAD;i++){
//        queryInp q;
//        q.tree = trees[i];
//        thread th(kNNQueryBatchThread,q, &res[i]);
//    }
}