//
// Created by Chuang on 2019/8/19.
//

#include "testFuncs.h"

int main() {
    calcuTime[0] = 0;
    srand((int) time(NULL));
//        vector<pair<id_type, Trajectory> > trajs = loadGTFolder();
    auto stat = trajStat::instance();
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
        std::cerr << stLen << "\t" << vbr<< "\t" <<vbc<< "\t" <<num << "\n";
    }
    return 0;
}