//
// Created by Chuang on 2021/2/9.

#include "testFuncs.h"

int main(int argc,char *argv[]){
    double seglens[] = {20}; // 5,10,20,40,80,160
    cerr<<"seglen: ";
    for(auto len:seglens){cerr<<len<<" ";}
    cerr<<endl;
    vector<string> files;
    for(int i=0;i<10;i++){
        files.emplace_back(fileFolder+to_string(i)+".out");
    }
    xStore x("od","od");
    for(int i=0;i<6;i++){
        cerr<<"size is "<<i<<endl;
        //don't load it if exist
        if(!x.m_property.contains("TBWP"+to_string(i)))
            x.loadFile(files[i]);
        string idxnum = to_string(i);
        vector<xTrajectory> queries;
        fillQuerySet(queries,x,1000);
        {
            MTQ q;
            q.prepareTrees(&x, [&idxnum](auto x) { return buildTBTreeWP(x,idxnum); });
            q.appendQueries(queries);
            std::cerr << q.runQueries().toString();
        }
//        {
//            MTQ q;
//            q.prepareTrees(&x, [&idxnum](auto x) { return buildSTRTreeWP(x,idxnum); });
//            q.appendQueries(queries);
//            std::cerr << q.runQueries().toString();
//        }
        for (auto len:seglens) {
            {
                MTQ q;
                q.prepareTrees(&x, [&len, &idxnum](auto x) {
                    return buildMBCRTreeWP(x, xTrajectory::OPTS, len, idxnum);
                });
                q.appendQueries(queries);
                std::cerr << q.runQueries().toString();
            }
//            {
//                MTQ q;
//                q.prepareTrees(&x, [&len, &idxnum](auto x) {
//                    return buildMBRRTreeWP(x, xTrajectory::OPTS, len, idxnum);
//                });
//                q.appendQueries(queries);
//                std::cerr << q.runQueries().toString();
//            }
        }
    }
}