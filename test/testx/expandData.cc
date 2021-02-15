//
// Created by Chuang on 2021/1/20.
//

#include "testFuncs.h"

int main() {
    auto trajs = loadDumpedFiledToTrajs(testFileName("tdfilter.txt"));
    dumpToFile(trajs,testFileName("tdexpand.data"));
    for(int i=0;i<10;i++) {
        auto repli = trajs;
        affine_transform(repli, xPoint(116.3972282409668, 39.90960456049752, 0),
                         0.2*i, xPoint(0, 0, 604800*i));
        dumpToFile_append(repli, testFileName("tdexpand.data"));
    }
    return 0;
}

//int main() {
//    auto trajs = loadDumpedFiledToTrajs(testFileName("glfilter.txt"));
//    dumpToFile(trajs,testFileName("glexpand.data"));
//    for(int i=0;i<10;i++) {
//        auto repli = trajs;
//        affine_transform(repli, xPoint(116.3972282409668, 39.90960456049752, 0),
//                         0.2*i, xPoint(0, 0, 86400*i));
//        dumpToFile_append(repli, testFileName("glexpand.data"));
//    }
//    return 0;
//}