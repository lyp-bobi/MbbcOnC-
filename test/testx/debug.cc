//
// Created by Chuang on 2020/12/10.
//
#include "testFuncs.h"
#include "random"

int main(int argc,char *argv[]){
    {
        string target = "tdexpand.datas";
        double qts[] = {300,1800,3600,7200,10800};
        cerr<<"03,mix, with TB,STR,SBB1800,SBBF(600,900,1200,1800)"<<endl;
        xStore x(target, testFileName(target), true);
        vector<xTrajectory> queries;
        default_random_engine e;
        auto queryLen1 =uniform_real_distribution<double>(60,300);
        auto queryLen2 =uniform_real_distribution<double>(3600,7200);
        auto queryLen3 = uniform_real_distribution<double>(86400);
        double count1=0,count2=0,count3=0;
        testtime*=2;
        for (int i = 0; i < testtime; i++) {
            int c = uniform_int_distribution<int>(1,3)(e);
            switch (c) {
                case 1:
                    queries.emplace_back(x.randomSubtraj(queryLen1(e)));
                    count1++;
                    break;
                case 2:
                    queries.emplace_back(x.randomSubtraj(queryLen2(e)));
                    count2++;
                    break;
                case 3:
                    queries.emplace_back(x.randomSubtraj(queryLen3(e)));
                    count3++;
                    break;
            }
        }
        cerr<<count1<<"\t"<<count2<<"\t"<<count3<<"\t";
        vector<int> nnks;
        for(int i=0;i<testtime;i++){
            nnks.emplace_back(random(6,201));
        }
//        {
//            MTQ q;
//            q.prepareTrees(&x, [](auto x) { return buildTBTreeWP(x); });
//            q.appendQueries(queries,nnks);
//            std::cerr << q.runQueries().toString();
//        }
        {
            MTQ q;
            q.prepareTrees(&x, [](auto x) { return buildSTRTreeWP(x); });
            q.appendQueries(queries,nnks);
            std::cerr << q.runQueries().toString();
        }
//        {
//            MTQ q;
//            q.prepareTrees(&x, [](auto x) { return buildMBCRTreeWP(x, xTrajectory::OPTS, 1800); });
//            q.appendQueries(queries,nnks);
//            std::cerr << q.runQueries().toString();
//        }
//        {
//            MTQ q;
//            SBBFMAP lens;
//            lens[make_pair(0,500)]=300;
//            lens[make_pair(500,20000)]=1800;
//            lens[make_pair(20000,1e300)]=10800;
//            q.prepareForest(&x,lens,4000);
//            q.appendQueries(queries,nnks);
//            std::cerr << q.runQueries().toString();
//        }
        cerr<<"mission complete.\n";
    }
    cerr<<"mission complete.\n";
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