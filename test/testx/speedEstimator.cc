//
// Created by Chuang on 2021/2/1.
//
#include "testFuncs.h"

int main(int argc,char *argv[]){
    string target;
    if(argc==1) {
        target = "tdfilter.txt";
    }else {
        target = "glfilter.txt";
    }
    std::cerr<<testFileName(target)<<endl;
    xStore x(target, testFileName(target), true);
    for (int stLen = 100; stLen <= 10000; stLen += 200) {
        double sum =0, count =0, v1=0,v2=0;
        for(int i=0;i<testtime*10;i++) {
            auto r = x.randomSubtraj(100000);
            auto sbbs = xTrajectory::OPTS(r,stLen);
            while(!sbbs.empty()){
                xSBB b=sbbs.front().second;
                sbbs.pop();
                double d = sqrt(sq(b.bc.m_pe.m_x - b.bc.m_ps.m_x)
                                + sq(b.bc.m_pe.m_y - b.bc.m_ps.m_y));
                double t = stLen;//b.m_endTime-b.m_startTime;
                double v = d / t;
                v1 += b.br.getArea() / t;
                v2 += b.bc.getArea() / t;
                sum += v;
                count += 1;
            }
        }
        cerr<<"len "<<stLen<<"\t v "<<sum/count<<"\t"<<v1/count<<"\t"<<v2/count<<endl;
    }
    xTrajectory tj;
    double vmax = 0;
//    for(auto &s:*(x.m_trajIdx)){
//        x.loadTraj(tj,xStoreEntry(s.first,0,100000));
//        if (0.1 < tj.maxSpeed()) {
//            cerr<<s.first <<"\t"<<vmax<<endl;
//            vmax = tj.maxSpeed();
//        }
//    }
    return 0;
}