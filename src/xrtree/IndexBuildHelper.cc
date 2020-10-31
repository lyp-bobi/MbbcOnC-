////
//// Created by Chuang on 2019/6/12.
////
//
#include <spatialindex/SpatialIndex.h>
#include "xRTree.h"
#include "BulkLoader.h"
#include "storagemanager/xStore.h"
using namespace SpatialIndex::xRTreeNsp;

SpatialIndex::xRTreeNsp::xRTree * SpatialIndex::xRTreeNsp::createNewxRTree(IStorageManager *store, long indexfan, long leaffan) {
    Tools::Variant var;
    Tools::PropertySet ps;

    var.m_varType = Tools::VT_DOUBLE;
    var.m_val.dblVal = 0.99;
    ps.setProperty("FillFactor", var);

    var.m_varType = Tools::VT_ULONG;
    var.m_val.ulVal = indexfan;
    ps.setProperty("IndexCapacity", var);

    var.m_varType = Tools::VT_ULONG;
    var.m_val.ulVal = leaffan;
    ps.setProperty("LeafCapacity", var);

    var.m_varType = Tools::VT_ULONG;
    var.m_val.ulVal = 3;
    ps.setProperty("Dimension", var);

    var.m_varType = Tools::VT_LONG;
    var.m_val.lVal = SpatialIndex::xRTreeNsp::RV_RSTAR;
    ps.setProperty("TreeVariant", var);

    SpatialIndex::xRTreeNsp::xRTree *ret = new SpatialIndex::xRTreeNsp::xRTree(*store, ps);

    ret->m_ts= static_cast<xStore *>(store);
//    var.m_varType = Tools::VT_LONGLONG;
//    var = ps.getProperty("IndexIdentifier");
//    indexIdentifier = var.m_val.llVal;

    return ret;
}

inline int idsize(){return 8;}
inline int mbrsize(){return 48;}
inline int mbcsize(){
    if(bCompactMBC) return 56;
    else return 64;
}
inline int linesize(){return 48;}
inline int pointersize(){return 16;}
inline int entrysize(){return 32;}

xRTree * xRTreeNsp::buildMBRRTree(IStorageManager *mng,
                                  const CUTFUNC &f) {
    auto store=static_cast<xStore*>(mng);
    auto stream = new xSBBStream(store,f);
    int bindex= PageSizeDefault/(idsize()+mbrsize()),
            bleaf=PageSizeDefault/(idsize()+mbrsize()+pointersize()+entrysize());
    xRTree* r = createNewxRTree(store,bindex,bleaf);
    r->m_bUsingMBR=true;
    BulkLoader bl;
    bl.bulkLoadUsingSTR(r,*stream,bindex,bleaf,PageSizeDefault*10,100000);
    return r;
}




//ISpatialIndex* SpatialIndex::xRTree::createAndBulkLoadNewRTreeWithxStore(IStorageManager *tsm,
//                                                                                 uint32_t pageSize, uint32_t dimension,
//                                                                                 id_type &indexIdentifier) {
//    std::cerr<<"start bulkloading with xStore\n";
//    auto start = std::chrono::system_clock::now();
//    xStore *ts= static_cast<xStore*>(tsm);
//    int indexCapacity=std::floor((pageSize-56)/56);
//    int leafCapacity=std::floor((pageSize-56)/88);
//    ISpatialIndex *tree;
//    if (ts->m_stream!= nullptr) {
//        tree = createAndBulkLoadNewxRTree(SpatialIndex::xRTree::BLM_STR, *ts->m_stream, *ts, 0.9,
//                                            indexCapacity, leafCapacity, dimension,
//                                            SpatialIndex::xRTree::RV_RSTAR, indexIdentifier,true);
//    } else {
//        auto dataStream = new xRTreeSegmentStream(ts);
//        tree = createAndBulkLoadNewxRTree(SpatialIndex::xRTree::BLM_STR, *dataStream, *ts, 0.9,
//                                            indexCapacity, leafCapacity, dimension,
//                                            SpatialIndex::xRTree::RV_RSTAR, indexIdentifier,true);
//        delete dataStream;
//    }
//    xRTree* r= static_cast<xRTree*>(tree);
//    auto end = std::chrono::system_clock::now();
//    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
//    double time = double(duration.count()) * std::chrono::microseconds::period::num / std::chrono::microseconds::period::den;
//    std::cerr<<"Tree built with time\t"<<time<<"\n";
//    ts->releaseTmp();
////    std::cerr<<"leaf time period"<<calcuTime[0]/calcuTime[1]<<"\n";
////    calcuTime[0]=calcuTime[1]=0;
//    return r;
//}