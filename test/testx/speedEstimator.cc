//
// Created by Chuang on 2021/2/1.
//
#include "testFuncs.h"

int main() {
    string target = "tdfilter.txt";
    std::cerr<<testFileName(target)<<endl;
    xStore x(target, testFileName(target), true);
    for (int stLen = 100; stLen <= 10000; stLen += 200) {
        double sum =0, count =0;
        for(int i=0;i<testtime;i++) {
            auto r = x.randomSubtraj(stLen);
            double d=sqrt(sq(r.m_points.back().m_x- r.m_points.front().m_x)
                    + sq(r.m_points.back().m_y- r.m_points.front().m_y));
            double t =r.m_points.back().m_t- r.m_points.front().m_t;
            double v = d/t;
            sum +=v;
            count+=1;
        }
        cerr<<"len "<<stLen<<"\t v "<<sum/count<<endl;
    }
    return 0;
}