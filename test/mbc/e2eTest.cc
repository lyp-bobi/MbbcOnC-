//
// Created by Chuang on 2019/9/1.
//

#include "testFuncs.h"


int main() {
    try {
        calcuTime[0] = 0;
        srand((int) time(NULL));
        vector<pair<id_type, Trajectory> > trajs = loadDumpedFiledToTrajs("/root/tdfilter.txt");
        vector<pair<id_type, vector<Trajectory>>> segs1, segs2;
        vector<pair<id_type, vector<Trajectory>>> emptyseg;
        int maxseg = 0;
        rsimpli = false;
        bUsingSimp = true;
        for (double queryLen = 100; queryLen <= 3600; queryLen += 200) {
            maxseg = 0;
            segs1.clear();
            segs2.clear();
            emptyseg.clear();
            for (auto &traj:trajs) {
                auto seg = traj.second.getFixedSegments();
                maxseg = std::max(int(seg.size()), maxseg);
                segs1.emplace_back(make_pair(traj.first, seg));
            }
            double segpara1 = biSearchMax(5, queryLen, 40, false, 0.012);
            std::cerr << "query len:" << queryLen << ",partial traj len:" << segpara1 << "\n";
            for (auto &traj:trajs) {
                auto seg = traj.second.getSegments(segpara1);
                maxseg = std::max(int(seg.size()), maxseg);
                segs2.emplace_back(make_pair(traj.first, seg));
            }

            string name0 = "name0", name1 = "name1", name2 = "name2";
            id_type indexIdentifier0, indexIdentifier1, indexIdentifier2;
            IStorageManager *diskfile0 = StorageManager::createNewDiskStorageManager(name0, 4096),
                    *diskfile1 = StorageManager::createNewDiskStorageManager(name1, 4096),
                    *diskfile2 = StorageManager::createNewDiskStorageManager(name2, 4096);
            // Create a new storage manager with the provided base name and a 4K page size.
            StorageManager::IBuffer *file0 = StorageManager::createNewRandomEvictionsBuffer(*diskfile0, 10, false),
                    *file1 = StorageManager::createNewRandomEvictionsBuffer(*diskfile1, 10, false),
                    *file2 = StorageManager::createNewRandomEvictionsBuffer(*diskfile2, 10, false);

            TrajStore *ts1 = new TrajStore(name1, file1, 4096, maxseg + 1);
            ts1->loadSegments(segs1);
            segs1.clear();
            segs1.swap(emptyseg);
            ISpatialIndex *r = MBCRTree::createAndBulkLoadNewRTreeWithTrajStore(ts1, 4096, 3, indexIdentifier1);

            TrajStore *ts2 = new TrajStore(name2, file2, 4096, maxseg + 1);
            ts2->loadSegments(segs2, true);
            segs2.clear();
            segs2.swap(emptyseg);
            ISpatialIndex *rc = MBCRTree::createAndBulkLoadNewMBCRTreeWithTrajStore(ts2, 4096, 3, indexIdentifier2);

            //kNN
            emptyseg.clear();
            vector<IShape *> queries;
            for (int i = 0; i < 500; i++) {
                auto ori = &trajs[(int(random(0, trajs.size()))) % trajs.size()].second;
                Trajectory *concate = new Trajectory();
                double ts = ori->randomPoint().m_time - queryLen/2;
                ori->getPartialTrajectory(ts, ts + queryLen, *concate);
                if (!concate->m_points.empty())
                    queries.emplace_back(concate);
            }
            bUsingSimp = bUsingSBBD = false;
            kNNQueryBatch(r, queries, ts1);
            bUsingSimp = bUsingSBBD = true;
            kNNQueryBatch(rc, queries, ts2);
            cerr << "\n";
            for (auto &tras:queries) {
                delete tras;
            }
            queries.clear();
            delete r;
            delete ts1;
            delete rc;
            delete ts2;
            delete file0;
            delete file1;
            delete file2;
            delete diskfile0;
            delete diskfile1;
            delete diskfile2;
        }
    }
    catch (Tools::Exception &e) {
        cerr << "******ERROR******" << endl;
        std::string s = e.what();
        cerr << s << endl;
        return -1;
    }
    return 0;
}