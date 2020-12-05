//
// Created by Chuang on 2019/8/19.
//

#include "testFuncs.h"

int main() {
    calcuTime[0] = 0;
    srand(0);
//        vector<pair<id_type, Trajectory> > trajs = loadGTFolder();

    int maxseg = 0;
    double avgSegLen = 100;
    vector<pair<id_type, Trajectory> > trajs = loadDumpedFiledToTrajs("/root/tdfilter.txt");
    Region tmpr;
    MBC tmpc;
    for (int stLen = 100; stLen <= 10000; stLen += 200) {
        double dist = 0;
        double vbr=0, vbc=0;
        double t = 0;
        int num = 0;
        for(auto traj:trajs){
            auto subtrajs = traj.second.getSegments(stLen);
            for(auto st:subtrajs){
                st.getMBR(tmpr);
                st.getMBC(tmpc);
                vbr += tmpr.getArea();
                vbc += tmpc.getArea();
                t += st.m_endTime()-st.m_startTime();
            }
        }
        for (int i = 0; i < 20000; i++) {
            auto ori = &trajs[(int(random(0, trajs.size()))) % trajs.size()].second;
            Trajectory *concate = new Trajectory();
            double ts = ori->randomPoint().m_time - stLen/2;
            ori->getPartialTrajectory(ts, ts + stLen, *concate);
            if (!concate->m_points.empty()) {
                dist += concate->m_points.front().getMinimumDistance(concate->m_points.back()) /
                        (concate->m_endTime() - concate->m_startTime());
            }
            delete concate;
        }
        std::cerr << stLen << "\t" << vbr<< "\t" <<vbc<< "\t" <<t <<"\t"<<dist/20000<< "\n";
    }


    trajs = loadDumpedFiledToTrajs("/root/glfilter.txt");
    for (int stLen = 100; stLen <= 10000; stLen += 200) {
        double dist = 0;
        double vbr=0, vbc=0;
        double t = 0;
        int num = 0;
        for(auto traj:trajs){
            auto subtrajs = traj.second.getSegments(stLen);
            for(auto st:subtrajs){
                st.getMBR(tmpr);
                st.getMBC(tmpc);
                vbr += tmpr.getArea();
                vbc += tmpc.getArea();
                t += st.m_endTime()-st.m_startTime();
            }
        }
        for (int i = 0; i < 20000; i++) {
            auto ori = &trajs[(int(random(0, trajs.size()))) % trajs.size()].second;
            Trajectory *concate = new Trajectory();
            double ts = ori->randomPoint().m_time - stLen/2;
            ori->getPartialTrajectory(ts, ts + stLen, *concate);
            if (!concate->m_points.empty()) {
                dist += concate->m_points.front().getMinimumDistance(concate->m_points.back()) /
                        (concate->m_endTime() - concate->m_startTime());
            }
            delete concate;
        }
        std::cerr << stLen << "\t" << vbr<< "\t" <<vbc<< "\t" <<t <<"\t"<<dist/20000<< "\n";
    }

    trajs = loadGTFolder();
    for (int stLen = 10; stLen <= 300; stLen += 10) {
        double dist = 0;
        double vbr=0, vbc=0;
        double t = 0;
        int num = 0;
        for(auto traj:trajs){
            auto subtrajs = traj.second.getSegments(stLen);
            for(auto st:subtrajs){
                st.getMBR(tmpr);
                st.getMBC(tmpc);
                vbr += tmpr.getArea();
                vbc += tmpc.getArea();
                t += st.m_endTime()-st.m_startTime();
            }
        }
        for (int i = 0; i < 20000; i++) {
            auto ori = &trajs[(int(random(0, trajs.size()))) % trajs.size()].second;
            Trajectory *concate = new Trajectory();
            double ts = ori->randomPoint().m_time - stLen/2;
            ori->getPartialTrajectory(ts, ts + stLen, *concate);
            if (!concate->m_points.empty()) {
                dist += concate->m_points.front().getMinimumDistance(concate->m_points.back()) /
                        (concate->m_endTime() - concate->m_startTime());
            }
            delete concate;
        }
        std::cerr << stLen << "\t" << vbr<< "\t" <<vbc<< "\t" <<t <<"\t"<<dist/20000<< "\n";
    }
    return 0;
}