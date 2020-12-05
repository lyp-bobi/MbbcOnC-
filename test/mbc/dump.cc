//
// Created by Chuang on 2020/9/24.
//

#include "testFuncs.h"

int main() {

    tjstat->usedata("td");
//    for(int j = 1000;j<5000;j+=500)
//        std::cerr<<j<<"\t"<<biSearchMax(5,j,50,true)<<"\n";
    biSearchMax(5,3600,50,true,-1,300,3000);
//    vector<pair<id_type, Trajectory> > trajs = loadGLToTrajs("D://simp.csv");
//    dumpToFile(trajs,"dumpedtraj.txt",500);
    return 0;
}