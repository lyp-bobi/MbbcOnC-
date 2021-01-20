//
// Created by Chuang on 2020/11/24.
//

#include "testFuncs.h"

int main() {
    auto trajs = loadDumpedFiledToTrajs("D://TRI-framework/dumpedtraj.txt");
    dumpToFile(trajs,"dumpedtraj.txt",50);
//    affine_transform(trajs,xPoint(116.3972282409668,39.90960456049752,0),
//            0.2, xPoint(0.2,0.1,0));
    dumpToFile_append(trajs,"dumpedtraj.txt",50);
//std::cerr<<getLastId("dumpedtraj.txt");
    return 0;
}