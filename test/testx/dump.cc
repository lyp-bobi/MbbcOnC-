//
// Created by Chuang on 2020/11/24.
//

#include "testFuncs.h"

int main() {
    auto trajs = loadGLToTrajs("D://simp.csv");
    dumpToFile(trajs,"dumpedtraj.txt",500);
    affine_transform(trajs,xPoint(116.3972282409668,39.90960456049752,0),
            0.2, xPoint(0.2,0.1,0));
    dumpToFile_append(trajs,"dumpedtraj.txt",500);
    return 0;
}