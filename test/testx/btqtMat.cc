//
// Created by Chuang on 2020/12/10.
//

#include "testFuncs.h"
#include "random"

int main(){
    try {
        string target = "tdfilter.txt";
        double qts[] = {300,900,1500,2100,2700,3300,3900,4500};
        for(auto qt:qts) {
            xStore x(target, testFileName(target), true);
            auto queryLen = qt;
            vector<xTrajectory> queries;
            for (int i = 0; i < testtime * 10; i++) {
                queries.emplace_back(x.randomSubtraj(queryLen));
            }
            tjstat->bt = qt;
            {
                MTQ q;
                q.prepareTrees(&x, [](auto x) { return buildTBTreeWP(x); });
                q.appendQueries(queries);
                std::cerr << q.runQueries().toString();
            }
            {
                MTQ q;
                q.prepareTrees(&x, [](auto x) { return buildSTRTreeWP(x); });
                q.appendQueries(queries);
                std::cerr << q.runQueries().toString();
            }
            double seglens[] = {300,900,1500,2100,2700,3300,3900,4500};
            for (auto len:seglens) {
                MTQ q;
                q.prepareTrees(&x, [&len](auto x) { return buildMBCRTreeWP(x, xTrajectory::OPTS, len); });
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