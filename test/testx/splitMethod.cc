//
// Created by Chuang on 2020/12/10.
//

#include "testFuncs.h"
#include "random"

int main(){
    try {
        string target = "tdfilter.txt";
        double avgQL = 1800;
        std::cerr<<testFileName(target)<<endl;
        xStore x(target, testFileName(target), true);
        default_random_engine e;
        auto queryLen =normal_distribution<double>(avgQL,500.0);
        vector<xTrajectory> queries;
        for (int i = 0; i < testtime; i++) {
            queries.emplace_back(x.randomSubtraj(queryLen(e)));
        }
        double len =1800;
        {
            MTQ q;
            q.prepareTrees(&x,[&len](auto x){return buildMBCRTreeWP(x,xTrajectory::ISS, len, "ISS");});
            q.appendQueries(queries);
            std::cerr<<q.runQueries().toString();
        }
        {
            MTQ q;
            q.prepareTrees(&x,[&len](auto x){return buildMBCRTreeWP(x,xTrajectory::GSS, len, "GSS");});
            q.appendQueries(queries);
            std::cerr<<q.runQueries().toString();
        }
        {
            MTQ q;
            q.prepareTrees(&x,[&len](auto x){return buildMBCRTreeWP(x,xTrajectory::OPTS, , "OPTS");});
            q.appendQueries(queries);
            std::cerr<<q.runQueries().toString();
        }
        {
            MTQ q;
            q.prepareTrees(&x,[&len](auto x){return buildMBCRTreeWP(x,xTrajectory::FP, len/tjstat->tl, "FP");});
            q.appendQueries(queries);
            std::cerr<<q.runQueries().toString();
        }
    }catch (Tools::Exception &e) {
        cerr << "******ERROR******" << endl;
        std::string s = e.what();
        cerr << s << endl;
        return -1;
    }
    return 0;
}