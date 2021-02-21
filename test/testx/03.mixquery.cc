//
// Created by Chuang on 2020/12/10.
//
#include "testFuncs.h"
#include "random"

int main(int argc,char *argv[]){
    if(argc==1) {
        string target = "tdexpand.data";
        double qts[] = {300,1800,3600,7200,10800};
        cerr<<"03,mix, with TB,STR,SBB1800,SBBF(600,900,1200,1800)"<<endl;
        xStore x(target, testFileName(target), true);
        vector<xTrajectory> queries;
        fillQuerySetRand(queries,x);
        vector<int> nnks;
        for(int i=0;i<testtime;i++){
            nnks.emplace_back(random(6,201));
        }
        {
            MTQ q;
            q.prepareTrees(&x, [](auto x) { return buildTBTreeWP(x); });
            q.appendQueries(queries,nnks);
            std::cerr << q.runQueries().toString();
        }
        {
            MTQ q;
            q.prepareTrees(&x, [](auto x) { return buildSTRTreeWP(x); });
            q.appendQueries(queries,nnks);
            std::cerr << q.runQueries().toString();
        }
        {
            MTQ q;
            q.prepareTrees(&x, [](auto x) { return buildMBCRTreeWP(x, xTrajectory::OPTS, 1800); });
            q.appendQueries(queries,nnks);
            std::cerr << q.runQueries().toString();
        }
        {
            MTQ q;
            SBBFMAP lens;
            lens[make_pair(0,900)]=600;
            lens[make_pair(900,3000)]=1200;
            lens[make_pair(3000,1e300)]=2700;
            q.prepareForest(&x,lens);
            q.appendQueries(queries,nnks);
            std::cerr << q.runQueries().toString();
        }
        cerr<<"mission complete.\n";
    }else{
        string target = "glexpand.data";
        cerr<<"03,mix, with TB,STR,SBB600,SBBF(200,)"<<endl;
        xStore x(target, testFileName(target), true);
        vector<xTrajectory> queries;
        fillQuerySetRand(queries,x);
        vector<int> nnks;
        for(int i=0;i<testtime;i++){
            nnks.emplace_back(random(6,201));
        }
        {
            MTQ q;
            q.prepareTrees(&x, [](auto x) { return buildTBTreeWP(x); });
            q.appendQueries(queries,nnks);
            std::cerr << q.runQueries().toString();
        }
        {
            MTQ q;
            q.prepareTrees(&x, [](auto x) { return buildSTRTreeWP(x); });
            q.appendQueries(queries,nnks);
            std::cerr << q.runQueries().toString();
        }
        {
            MTQ q;
            q.prepareTrees(&x, [](auto x) { return buildMBCRTreeWP(x, xTrajectory::OPTS, 1800); });
            q.appendQueries(queries,nnks);
            std::cerr << q.runQueries().toString();
        }
        {
            MTQ q;
            SBBFMAP lens;
            lens[make_pair(0,500)]=200;
            lens[make_pair(500,1e300)]=600;
            q.prepareForest(&x,lens);
            q.appendQueries(queries,nnks);
            std::cerr << q.runQueries().toString();
        }
        cerr<<"mission complete.\n";
    }
    return 0;
}
// td
//300 - 600
//600-600
//900-900
//1200-900
//1300 - 1200
//1800-1200
//2100-1800
//2300 - 1800
//3300 - 1800
//4300-1800
//5300-1800
//6300-1800
//7300-2700
//8300-2700
//9300-2700
//10300-2700