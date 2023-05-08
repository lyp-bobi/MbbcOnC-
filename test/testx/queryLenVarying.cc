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
        vector<xTrajectory> queries;
        fillQuerySet(queries,x,avgQL,500.0);
        tjstat->bt = avgQL;
        {
            MTQ q;
            q.prepareTrees(&x,[](auto x){return buildTBTreeWP(x);});
            q.appendQueries(queries);
            std::cerr<<q.runQueries().toString();
        }
        {
            MTQ q;
            q.prepareTrees(&x,[](auto x){return buildSTRTreeWP(x);});
            q.appendQueries(queries);
            std::cerr<<q.runQueries().toString();
        }
        double seglens[] ={1200,1500,1800,2100,2400,3000};
        for(auto len:seglens)
        {
            MTQ q;
            q.prepareTrees(&x,[&len](auto x){return buildMBCRTreeWP(x,xTrajectory::GSS, len);});
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