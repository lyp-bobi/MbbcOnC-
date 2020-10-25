//
// Created by Chuang on 2019/6/12.
//

#include <spatialindex/SpatialIndex.h>
#include "xRTree.h"
#include "storagemanager/xStore.h"

class xRTreeSegmentStream:public baseSegmentStream{
public:
    xRTreeSegmentStream(xStore *ts): baseSegmentStream(ts){}
    IData* constructData(id_type id,xMBR mbr,xMBC xMBC) override {
        xRTree::Data* d=new xRTree::Data(0, nullptr, xMBC,mbr, id);
        return d;
    }
};


ISpatialIndex* SpatialIndex::xRTree::createAndBulkLoadNewxRTreeWithxStore(IStorageManager *tsm,
                                                                           uint32_t pageSize, uint32_t dimension,
                                                                           id_type &indexIdentifier) {
    std::cerr<<"start bulkloading with xStore\n";
    auto start = std::chrono::system_clock::now();
    xStore *ts = static_cast<xStore *>(tsm);
    int indexCapacity = std::floor((pageSize - 56) / 56);
    int leafCapacity = std::floor((pageSize - 56) / 104);
    ISpatialIndex *tree;
    if (ts->m_stream!= nullptr) {
        tree = createAndBulkLoadNewxRTree(SpatialIndex::xRTree::BLM_STR, *ts->m_stream, *ts, 0.9,
                                                           indexCapacity, leafCapacity, dimension,
                                                           SpatialIndex::xRTree::RV_RSTAR, indexIdentifier);
    } else {
        auto dataStream = new xRTreeSegmentStream(ts);
        tree = createAndBulkLoadNewxRTree(SpatialIndex::xRTree::BLM_STR, *dataStream, *ts, 0.9,
                                                           indexCapacity, leafCapacity, dimension,
                                                           SpatialIndex::xRTree::RV_RSTAR, indexIdentifier);
        delete dataStream;
    }
    xRTree* r= static_cast<xRTree*>(tree);
    auto end = std::chrono::system_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    double time = double(duration.count()) * std::chrono::microseconds::period::num / std::chrono::microseconds::period::den;
    std::cerr<<"Tree built with time\t"<<time<<"\n";
    ts->releaseTmp();
//    std::cerr<<"leaf time period"<<calcuTime[0]/calcuTime[1]<<"\n";
    calcuTime[0]=calcuTime[1]=0;

    return r;
}

ISpatialIndex* SpatialIndex::xRTree::createAndBulkLoadNewRTreeWithxStore(IStorageManager *tsm,
                                                                                 uint32_t pageSize, uint32_t dimension,
                                                                                 id_type &indexIdentifier) {
    std::cerr<<"start bulkloading with xStore\n";
    auto start = std::chrono::system_clock::now();
    xStore *ts= static_cast<xStore*>(tsm);
    int indexCapacity=std::floor((pageSize-56)/56);
    int leafCapacity=std::floor((pageSize-56)/88);
    ISpatialIndex *tree;
    if (ts->m_stream!= nullptr) {
        tree = createAndBulkLoadNewxRTree(SpatialIndex::xRTree::BLM_STR, *ts->m_stream, *ts, 0.9,
                                            indexCapacity, leafCapacity, dimension,
                                            SpatialIndex::xRTree::RV_RSTAR, indexIdentifier,true);
    } else {
        auto dataStream = new xRTreeSegmentStream(ts);
        tree = createAndBulkLoadNewxRTree(SpatialIndex::xRTree::BLM_STR, *dataStream, *ts, 0.9,
                                            indexCapacity, leafCapacity, dimension,
                                            SpatialIndex::xRTree::RV_RSTAR, indexIdentifier,true);
        delete dataStream;
    }
    xRTree* r= static_cast<xRTree*>(tree);
    auto end = std::chrono::system_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    double time = double(duration.count()) * std::chrono::microseconds::period::num / std::chrono::microseconds::period::den;
    std::cerr<<"Tree built with time\t"<<time<<"\n";
    ts->releaseTmp();
//    std::cerr<<"leaf time period"<<calcuTime[0]/calcuTime[1]<<"\n";
//    calcuTime[0]=calcuTime[1]=0;
    return r;
}