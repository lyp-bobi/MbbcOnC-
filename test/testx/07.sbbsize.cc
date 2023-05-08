//
// Created by Chuang on 2021/2/9.
//

#include "testFuncs.h"

int main(int argc,char *argv[]){
    string target = "tdfilter.txt";
    cerr<<"07,sbbsize";
    xStore x(target, testFileName(target), true);
    xTrajectory tj;
    for(double bt=300;bt<3600;bt+=300) {
        double sum1, sum2;
//        for (auto &s:*(x.m_trajIdx)) {
//            x.loadTraj(tj, xStoreEntry(s.first, 0, 10000));
//            auto bbs = xTrajectory::GSS(tj, bt);
//            while (!bbs.empty()) {
//                auto b = bbs.front();
//                sum1 += b.second.br.getArea();
//                sum2 += b.second.bc.getArea();
//            }
//        }
        cerr << "br" << sum1 / tjstat->M << endl << sum2 / tjstat->M << endl;
    }
}