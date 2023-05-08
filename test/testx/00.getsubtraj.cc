//
// Created by Chuang on 2020/12/10.
//
#include "testFuncs.h"
#include "random"

class PrintVisitor : public IVisitor {
public:
    size_t m_indexvisited=0;
    size_t m_leafvisited=0;
    size_t m_resultGet;
    id_type m_lastResult;
    double m_lastDist = 0;
    xTrajectory *m_query;
    xStore *ts = nullptr;

public:
    PrintVisitor() : m_resultGet(0), m_indexvisited(0), m_leafvisited(0) {}

    void visitNode(const INode &n) {
        if (n.isLeaf()) {
            m_leafvisited++;
        }
        else {
            m_indexvisited++;
        }
    }
    void visitData(std::vector<const IData*>& v){m_resultGet++;}
    void visitData(const IData &d) {
        m_resultGet++;
        m_lastResult = d.getIdentifier();
        auto mou=dynamic_cast<const xRTreeNsp::xRTree::simpleData*>(&d);
        if(mou!=nullptr){
            xTrajectory traj, subtraj;
            ts->loadTraj(traj, xStoreEntry(mou->m_id, 0, 1000000));
            m_lastDist=mou->m_dist;
            traj.getPartialxTrajectory(m_query->m_startTime(), m_query->m_endTime(), subtraj);
            cerr <<mou->m_dist << "\t" <<subtraj.toString() << endl;
        }
    }

    void clear(){
        m_indexvisited = m_leafvisited = m_resultGet=0;
    }
};


int main(int argc,char *argv[]){
    try {
        vector<double> seglens;
        string target;
        if(argc==1) {
            target = "tdexpand.datas";
            seglens = {300,600,900,1200,1800,2400,3600};//600,900,1200,1800,2700,3600,5400,7200,9000
        }else {
            target = "glexpand.datas";
            seglens = {300,600,900,1200,1800,2400,3600};
        }
//        testtime=400;
        cerr<<"seglen: ";
        for(auto len:seglens){cerr<<len<<" ";}
        cerr<<"TB STR ";
        cerr<<endl;
        xStore x(target, testFileName(target), true);
        xRTree *r = buildMBCRTreeWP(&x,xTrajectory::OneBox, 10000);
        PrintVisitor vis;
        vis.ts = &x;
        xTrajectory query;
        query.loadFromString("116.282270,39.823860,605368.000000 116.323910,39.829200,605670.000000 116.341250,39.855650,605972.000000 116.374320,39.869290,606274.000000 116.426550,39.869600,606576.000000 116.433110,39.898170,606878.000000 116.428405,39.928466,607168.000000");
        std::cerr<<"";
        vis.m_query = &query;
        cerr << query.toString()<<endl;
        current_distance = SDDTW;
        r->nearestNeighborQuery(6, query, vis);
        current_distance = RMDTW;
        r->nearestNeighborQuery(6, query, vis);
        current_distance = IED;
        r->nearestNeighborQuery(6, query, vis);
        cerr<<"mission complete.\n";
    }catch (Tools::Exception &e) {
        cerr << "******ERROR******" << endl;
        std::string s = e.what();
        cerr << s << endl;
        return -1;
    }
    return 0;
}