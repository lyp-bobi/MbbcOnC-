//
// Created by Chuang on 2019/6/12.
//

#include <spatialindex/SpatialIndex.h>
#include "MBCRTree.h"
#include "storagemanager/TrajStore.h"

class mbcrtreeSegmentStream:public baseSegmentStream{
public:
    mbcrtreeSegmentStream(TrajStore *ts): baseSegmentStream(ts){}
    IData* constructData(id_type id,Region mbr,MBC mbc) override {
        MBCRTree::Data* d=new MBCRTree::Data(0, nullptr, mbc,mbr, id);
        return d;
    }
};


ISpatialIndex* SpatialIndex::MBCRTree::createAndBulkLoadNewMBCRTreeWithTrajStore(IStorageManager *tsm,
                                                                           uint32_t pageSize, uint32_t dimension,
                                                                           id_type &indexIdentifier) {
    std::cerr<<"start bulkloading with TrajStore\n";
    TrajStore *ts = static_cast<TrajStore *>(tsm);
    int indexCapacity = std::floor((pageSize - 56) / 56);
    int leafCapacity = std::floor((pageSize - 56) / 104);
    ISpatialIndex *tree;
    if (ts->m_stream!= nullptr) {
        tree = createAndBulkLoadNewMBCRTree(SpatialIndex::MBCRTree::BLM_STR, *ts->m_stream, *ts, 0.9,
                                                           indexCapacity, leafCapacity, dimension,
                                                           SpatialIndex::MBCRTree::RV_RSTAR, indexIdentifier);
    } else {
        auto dataStream = new mbcrtreeSegmentStream(ts);
        tree = createAndBulkLoadNewMBCRTree(SpatialIndex::MBCRTree::BLM_STR, *dataStream, *ts, 0.9,
                                                           indexCapacity, leafCapacity, dimension,
                                                           SpatialIndex::MBCRTree::RV_RSTAR, indexIdentifier);
        delete dataStream;
    }
    MBCRTree* r= static_cast<MBCRTree*>(tree);
    std::cerr<<"Tree built.\n";
    ts->releaseTmp();
//    std::cerr<<"leaf time period"<<calcuTime[0]/calcuTime[1]<<"\n";
    calcuTime[0]=calcuTime[1]=0;

    return r;
}

ISpatialIndex* SpatialIndex::MBCRTree::createAndBulkLoadNewRTreeWithTrajStore(IStorageManager *tsm,
                                                                                 uint32_t pageSize, uint32_t dimension,
                                                                                 id_type &indexIdentifier) {
    std::cerr<<"start bulkloading with TrajStore\n";
    TrajStore *ts= static_cast<TrajStore*>(tsm);
    int indexCapacity=std::floor((pageSize-56)/56);
    int leafCapacity=std::floor((pageSize-56)/88);
    ISpatialIndex *tree;
    if (ts->m_stream!= nullptr) {
        tree = createAndBulkLoadNewMBCRTree(SpatialIndex::MBCRTree::BLM_STR, *ts->m_stream, *ts, 0.9,
                                            indexCapacity, leafCapacity, dimension,
                                            SpatialIndex::MBCRTree::RV_RSTAR, indexIdentifier,true);
    } else {
        auto dataStream = new mbcrtreeSegmentStream(ts);
        tree = createAndBulkLoadNewMBCRTree(SpatialIndex::MBCRTree::BLM_STR, *dataStream, *ts, 0.9,
                                            indexCapacity, leafCapacity, dimension,
                                            SpatialIndex::MBCRTree::RV_RSTAR, indexIdentifier,true);
        delete dataStream;
    }
    MBCRTree* r= static_cast<MBCRTree*>(tree);
    std::cerr<<"Tree built.\n";
    ts->releaseTmp();
//    std::cerr<<"leaf time period"<<calcuTime[0]/calcuTime[1]<<"\n";
//    calcuTime[0]=calcuTime[1]=0;
    return r;
}