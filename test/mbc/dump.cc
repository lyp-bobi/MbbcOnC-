//
// Created by Chuang on 2020/9/24.
//

#include "testFuncs.h"

int main() {
    vector<pair<id_type, Trajectory> > trajs = loadGLToTrajs("D://simp.csv");
    dumpToFile(trajs,"dumpedtraj.txt",500);
    return 0;
}