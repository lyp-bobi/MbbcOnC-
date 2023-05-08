//
// Created by Chuang on 2021/1/20.
//

#include "testFuncs.h"

int main(int argc,char *argv[]){
    vector<double> seglens;
    string target;
    if(argc==1) {
        target = "tdexpand.datas";
    }else {
        target = "glexpand.datas";
    }
    xStore x(target, testFileName(target), true);
    vector<id_type> ids;
//    for(auto &s:*(x.m_trajIdx)) ids.emplace_back(s.first);
    std::random_shuffle(ids.begin(),ids.end());
    ofstream outFile(testFileName(target)+"s", ios::out);
    xTrajectory tj;
    for(auto id:ids){
        x.loadTraj(tj,xStoreEntry(id,0,1e8));
        outFile<<id<<"\n"<<tj.toString()<<"\n";
    }
    outFile.flush();
    outFile.close();
    cerr<<"mission complete";
    return 0;
}
