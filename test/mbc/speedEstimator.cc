//
// Created by Chuang on 2019/8/19.
//

#include "testFuncs.h"

int main(){
    calcuTime[0] = 0;
    srand((int) time(NULL));
//        vector<pair<id_type, Trajectory> > trajs = loadGTFolder();
    auto stat=trajStat::instance();
    int maxseg = 0;
    double avgSegLen=100;
    vector<pair<id_type, Trajectory> > trajs = loadGLToTrajs("/root/TD.csv");
    for (int queryLen=100;queryLen<=5000;queryLen+=100) {
        double dist=0;
        for (int i = 0; i < 10000; i++) {
            auto ori = &trajs[(int(random(0, trajs.size()))) % trajs.size()].second;
            Trajectory *concate = new Trajectory();
            double ts = std::max(ori->m_startTime(), random(ori->m_startTime(), ori->m_endTime() - queryLen));
            ori->getPartialTrajectory(ts, ts + queryLen, *concate);
            if (!concate->m_points.empty()){
                dist+=concate->m_points.front().getMinimumDistance(concate->m_points.back())/(concate->m_endTime()-concate->m_startTime());
            }
        }
        std::cerr<<queryLen<<"\t"<<dist/10000<<"\n";
    }
    return 0;
}