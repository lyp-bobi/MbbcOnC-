//
// Created by Chuang on 2020/10/28.
//
#include "storagemanager/xStore.h"
#include "../../src/xrtree/xRTree.h"



class MyVisitor : public IVisitor {
public:
    size_t m_indexvisited;
    size_t m_leafvisited;
    size_t m_resultGet;
    id_type m_lastResult;
    double m_lastDist = 0;
    IShape *m_query;
    xStore *ts = nullptr;

public:
    MyVisitor() : m_resultGet(0), m_indexvisited(0), m_leafvisited(0) {}

    void visitNode(const INode &n) {
        if (n.isLeaf()) {
            m_leafvisited++;
        }
        else {
            m_indexvisited++;
        }
    }
    void visitData(std::vector<const IData*>& v){}
    void visitData(const IData &d) {
        m_resultGet++;
        m_lastResult = d.getIdentifier();
        auto mou=dynamic_cast<const xRTreeNsp::xRTree::simpleData*>(&d);
        if(mou!=nullptr){
            m_lastDist=mou->m_dist;
//            cerr << d.getIdentifier()<<"\t"<<mou->m_dist << endl;
        }
    }

    void clear(){
        m_indexvisited = m_leafvisited = m_resultGet=0;
    }
};



using namespace xRTreeNsp;
int main(){
    xStore x("test", "D://TRI-framework/dumpedtraj.txt",true);
    auto stat = trajStat::instance();
    double querylens[]={50,100,200,400,8000,1600,10000};
    xCylinder query(xPoint(40,116.327,6516),0.0001,6516,9516,2);
    xTrajectory tj;
    x.loadTraj(tj, xStoreEntry(1850,0,1000));
//    x.loadTraj(tj, xStoreEntry(4337,0,1000));
    tj.intersectsxCylinder(query);
    CUTFUNC f=xTrajectory::ISS;
//    cerr<<tj.intersectsxCylinder(query);
    stat->bt=50;
    xSBB s;
    s.loadFromString("0 1 0 39.984438 116.347402 9194.000000 39.998213 116.365543 10194.000000 0.043406 0.000149");
    query.checkRel(s);
    MyVisitor vis;
    try {
        bUsingSBBD = false;
        for (auto querylen:querylens) {
            stat->bt = querylen;
            auto r = buildMBCRTree(&x, f);
            r->intersectsWithQuery(query, vis);
            r->nearestNeighborQuery(5, tj, vis);
            std::cerr << vis.m_resultGet << " " << vis.m_lastResult << endl;
            vis.clear();
            auto r2 = buildMBRRTree(&x, f);
            r2->intersectsWithQuery(query, vis);
            r2->nearestNeighborQuery(5, tj, vis);
            std::cerr << vis.m_resultGet << " " << vis.m_lastResult << endl;
            vis.clear();
        }
        {
            auto r = buildTBTree(&x);
            r->intersectsWithQuery(query, vis);
            r->nearestNeighborQuery(5, tj, vis);
            std::cerr << vis.m_resultGet << " " << vis.m_lastResult << endl;
            vis.clear();
        }
        bUsingSBBD = true;
        for (auto querylen:querylens) {
            stat->bt = querylen;
            auto r = buildMBCRTree(&x, f);
            r->intersectsWithQuery(query, vis);
            r->nearestNeighborQuery(5, tj, vis);
            std::cerr << vis.m_resultGet << " " << vis.m_lastResult << endl;
            vis.clear();
            auto r2 = buildMBRRTree(&x, f);
            r2->intersectsWithQuery(query, vis);
            r2->nearestNeighborQuery(5, tj, vis);
            std::cerr << vis.m_resultGet << " " << vis.m_lastResult << endl;
            vis.clear();
        }
        {
            auto r = buildTBTree(&x);
            r->intersectsWithQuery(query, vis);
            r->nearestNeighborQuery(5, tj, vis);
            std::cerr << vis.m_resultGet << " " << vis.m_lastResult << endl;
            vis.clear();
        }
    }
    catch (Tools::Exception &e) {
        cerr << "******ERROR******" << endl;
        std::string s = e.what();
        cerr << s << endl;
        return -1;
    }
}