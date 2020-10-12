//
// Created by Chuang on 2020/10/9.
//

#include "testFuncs.h"

int main() {
//    vector<pair<id_type, Trajectory> > trajs = loadDumpedFiledToTrajs("D://TRI-framework/dumpedtraj.txt");
    auto stat= trajStat::instance();
    vector<pair<id_type, Trajectory> > trajs = loadDumpedFiledToTrajs("/root/tdfilter.txt");
    double segLenParas[] = {300, 600, 900, 1200, 1500, 1800, 2100, 2400, 2700, 3000, 4000};
    int maxseg;
    for (double segLen:segLenParas) {
        vector<pair<id_type, vector<Trajectory>>> segs;
        vector<pair<id_type, vector<Trajectory>>> emptyseg;
        maxseg = 300;
        segs.clear();
        int totallen = 0, totalseg = 0;
        for (const auto &traj:trajs) {
            totallen += traj.second.m_points.size();
            auto seg = traj.second.getSegments(segLen);
            totalseg += seg.size();
            maxseg = std::max(int(seg.size()), maxseg);
            segs.emplace_back(make_pair(traj.first, seg));
        }
        string name0 = "name0", name1 = "name1", name2 = "name2";
        id_type indexIdentifier0, indexIdentifier1, indexIdentifier2;
        IStorageManager *diskfile1 = StorageManager::createNewDiskStorageManager(name1, 4096);
        // Create a new storage manager with the provided base name and a 4K page size.
        StorageManager::IBuffer
                *file1 = StorageManager::createNewRandomEvictionsBuffer(*diskfile1, 10, false);

        TrajStore *ts1 = new TrajStore(name1, diskfile1, 4096, maxseg + 1);
        ts1->loadSegments(segs, true);
        ISpatialIndex *r = MBCRTree::createAndBulkLoadNewRTreeWithTrajStore(ts1, 4096, 3, indexIdentifier1);
        delete r;
        delete ts1;
        delete file1;
        delete diskfile1;
    }
    trajs.clear();
    stat->init();
    trajs = loadDumpedFiledToTrajs("/root/glfilter.txt");
    for (double segLen:segLenParas) {
        vector<pair<id_type, vector<Trajectory>>> segs;
        vector<pair<id_type, vector<Trajectory>>> emptyseg;
        maxseg = 300;
        segs.clear();
        int totallen = 0, totalseg = 0;
        for (const auto &traj:trajs) {
            totallen += traj.second.m_points.size();
            auto seg = traj.second.getSegments(segLen);
            totalseg += seg.size();
            maxseg = std::max(int(seg.size()), maxseg);
            segs.emplace_back(make_pair(traj.first, seg));
        }
        string name0 = "name0", name1 = "name1", name2 = "name2";
        id_type indexIdentifier0, indexIdentifier1, indexIdentifier2;
        IStorageManager *diskfile1 = StorageManager::createNewDiskStorageManager(name1, 4096);
        // Create a new storage manager with the provided base name and a 4K page size.
        StorageManager::IBuffer
                *file1 = StorageManager::createNewRandomEvictionsBuffer(*diskfile1, 10, false);

        TrajStore *ts1 = new TrajStore(name1, diskfile1, 4096, maxseg + 1);
        ts1->loadSegments(segs, true);
        ISpatialIndex *r = MBCRTree::createAndBulkLoadNewRTreeWithTrajStore(ts1, 4096, 3, indexIdentifier1);

        delete r;
        delete ts1;
        delete file1;
        delete diskfile1;
    }
    trajs.clear();
    stat->init();
    trajs = loadGTFolder();
    double segLenParas2[] = {10,30,50,70,90,110,130,150};
    for (double segLen:segLenParas2) {
        vector<pair<id_type, vector<Trajectory>>> segs;
        vector<pair<id_type, vector<Trajectory>>> emptyseg;
        maxseg = 300;
        segs.clear();
        int totallen = 0, totalseg = 0;
        for (const auto &traj:trajs) {
            totallen += traj.second.m_points.size();
            auto seg = traj.second.getSegments(segLen);
            totalseg += seg.size();
            maxseg = std::max(int(seg.size()), maxseg);
            segs.emplace_back(make_pair(traj.first, seg));
        }
        string name0 = "name0", name1 = "name1", name2 = "name2";
        id_type indexIdentifier0, indexIdentifier1, indexIdentifier2;
        IStorageManager *diskfile1 = StorageManager::createNewDiskStorageManager(name1, 4096);
        // Create a new storage manager with the provided base name and a 4K page size.
        StorageManager::IBuffer
                *file1 = StorageManager::createNewRandomEvictionsBuffer(*diskfile1, 10, false);

        TrajStore *ts1 = new TrajStore(name1, diskfile1, 4096, maxseg + 1);
        ts1->loadSegments(segs, true);
        ISpatialIndex *r = MBCRTree::createAndBulkLoadNewRTreeWithTrajStore(ts1, 4096, 3, indexIdentifier1);

        delete r;
        delete ts1;
        delete file1;
        delete diskfile1;
    }

}