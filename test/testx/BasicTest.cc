//
// Created by Chuang on 2020/10/28.
//
#include "storagemanager/xStore.h"
#include "../../src/xrtree/xRTree.h"



class MyVisitor : public IVisitor {
public:
    size_t m_indexIO;
    size_t m_leafIO;
    size_t m_indexvisited;
    size_t m_leafvisited;
    size_t m_resultGet;
    id_type m_lastResult;
    double m_lastDist = 0;
    IShape *m_query;
    xStore *ts = nullptr;

public:
    MyVisitor() : m_indexIO(0), m_leafIO(0), m_resultGet(0), m_indexvisited(0), m_leafvisited(0) {}

    void visitNode(const INode &n) {
//        if (n.isLeaf()) m_leafIO++;
//        else m_indexIO++;
        uint32_t size = n.getIndexByteArraySize();

        if (n.isLeaf()) {
            m_leafvisited++;
            m_leafIO += size;
        }
        else {
            m_indexvisited++;
            m_indexIO += size;
        }
    }

    void visitData(const IData &d) {
        m_resultGet++;
        IShape *pS;
        d.getShape(&pS);
        // do something.
        delete pS;
//        cout<<"data"<<endl;

        // data should be an array of characters representing a Region as a string.
        uint8_t *pData = 0;
        uint32_t cLen = 0;
        d.getData(cLen, &pData);
        m_leafIO += cLen;
        // do something.
//        double *s = reinterpret_cast<double*>(pData);
//        cout << *s << endl;

        m_lastResult = d.getIdentifier();
//        cerr << d.getIdentifier()<<std::endl;
//        auto mou=dynamic_cast<const MBCRTree::MBCRTree::simpleData*>(&d);
//        if(mou!=nullptr){
//            m_lastDist=mou->m_dist;
//            cerr << d.getIdentifier()<<"\t"<<mou->m_dist << endl;
//        }
//        auto mou2=dynamic_cast<const RTree::RTree::simpleData*>(&d);
//        if(mou2!=nullptr){
//            m_lastDist=mou2->m_dist;
//            cerr << d.getIdentifier()<<"\t"<<mou2->m_dist << endl;
//        }
//        //id of the data
//        if(ts== nullptr)
//            cerr << d.getIdentifier() << endl;
//        else
//            cerr << d.getIdentifier() /100<< endl;
//        //the traj
//        Trajectory traj;
//        if(ts== nullptr) {
//            traj.loadFromByteArray(pData);
//            double mindist=m_query->getMinimumDistance(traj);
//            cerr<<"traj dist is"<<mindist<<"\n\n";
//        }
//        else{
//            id_type tid=d.getIdentifier();
//            traj=ts->getTrajByTime(tid,0,1000);
//            auto brs=ts->getMBRsByTime(tid,0,1000);
//            auto bcs=ts->getMBCsByTime(tid,0,1000);
//            cerr<<"traj dist is"<<m_query->getMinimumDistance(traj)<<"\n"
//                    <<m_query->getMinimumDistance(brs)<<"\n"
//                    <<m_query->getMinimumDistance(bcs)<<"\n\n";
//        }

        delete[] pData;
    }

    void visitData(std::vector<const IData *> &v) {
        cout << v[0]->getIdentifier() << " " << v[1]->getIdentifier() << endl;
    }
};



using namespace xRTreeNsp;
int main(){
    xStore x("test", "D://TRI-framework/dumpedtraj.txt");
    auto stat = trajStat::instance();
    double querylens[]={500,1000,1500,100000};
    xCylinder query(xPoint(40,116.327,6516),0.0001,6516,9516,2);
    xTrajectory tj;
    x.loadTraj(tj, xStoreEntry(1850,0,1000));
    CUTFUNC f=xTrajectory::ISS;
    stat->bt=500;
    f(tj);
    for(auto querylen:querylens){
        stat->bt = querylen;
        MyVisitor vis;
        vis.ts = &x;
        auto r = buildMBRRTree(&x,f);
//        r->intersectsWithQuery(query,vis);
        r->nearestNeighborQuery(5,tj,vis);
        std::cerr<<vis.m_resultGet<<vis.m_lastResult<<endl;
    }
}