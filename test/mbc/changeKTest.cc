//
// Created by Chuang on 2019/11/24.
//

#include "testFuncs.h"


int main() {
    try {
        calcuTime[0] = 0;
        srand(0);
//        vector<pair<id_type, Trajectory> > trajs = loadGTFolder(1);
        vector<pair<id_type, Trajectory> > trajs = loadDumpedFiledToTrajs("/root/glfilter.txt");
        vector<pair<id_type, vector<Trajectory>>> segs1, segs2;
        vector<pair<id_type, vector<Trajectory>>> emptyseg;
        int maxseg = 0;
        rsimpli = false;
        bUsingSimp = true;
        vector<IShape *> queries;
        auto stat = trajStat::instance();
        double queryLen = 1800;
        for (int i = 0; i < testtime; i++) {
            auto ori = &trajs[(int(random(0, trajs.size()))) % trajs.size()].second;
            Trajectory *concate = new Trajectory();
            double ts = ori->randomPoint().m_time - queryLen/2;
            ori->getPartialTrajectory(ts, ts + queryLen, *concate);
            if (!concate->m_points.empty())
                queries.emplace_back(concate);
        }
        for (int nnk = 10; nnk <= 100; nnk += 10) {
            std::cerr << "nnk is" << nnk << "\n";
            maxseg = 0;
            segs1.clear();
            segs2.clear();
            emptyseg.clear();
            for (auto &traj:trajs) {
                auto seg = traj.second.getFixedSegments();
                maxseg = std::max(int(seg.size()), maxseg);
                segs1.emplace_back(make_pair(traj.first, seg));
            }
            double segpara1 = biSearchMax(nnk, 1800, 40, false, 0.012);
            std::cerr << "segpara is estimated as " << segpara1;
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

            TrajStore *ts1 = new TrajStore(name1, diskfile1, 4096, maxseg + 1);
            ts1->loadSegments(segs1);
            segs1.clear();
            segs1.swap(emptyseg);
            ISpatialIndex *r = MBCRTree::createAndBulkLoadNewRTreeWithTrajStore(ts1, 4096, 3, indexIdentifier1);

            TrajStore *ts2 = new TrajStore(name2, diskfile2, 4096, maxseg + 1);
            ts2->loadSegments(segs2, true);
            segs2.clear();
            segs2.swap(emptyseg);
            ISpatialIndex *rc = MBCRTree::createAndBulkLoadNewMBCRTreeWithTrajStore(ts2, 4096, 3, indexIdentifier2);

            bUsingSimp = bUsingSBBD = false;
            kNNQueryBatch(r, queries, ts1, nnk);
            bUsingSimp = bUsingSBBD = true;
            kNNQueryBatch(rc, queries, ts2, nnk);
            cerr << "\n";
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




    try {
        calcuTime[0] = 0;
        srand(0);
//        vector<pair<id_type, Trajectory> > trajs = loadGTFolder(1);
        vector<pair<id_type, Trajectory> > trajs = loadDumpedFiledToTrajs("/root/tdfilter.txt");
        vector<pair<id_type, vector<Trajectory>>> segs1, segs2;
        vector<pair<id_type, vector<Trajectory>>> emptyseg;
        int maxseg = 0;
        rsimpli = false;
        bUsingSimp = true;
        vector<IShape *> queries;
        auto stat = trajStat::instance();
        double queryLen = 1800;
        for (int i = 0; i < testtime; i++) {
            auto ori = &trajs[(int(random(0, trajs.size()))) % trajs.size()].second;
            Trajectory *concate = new Trajectory();
            double ts = ori->randomPoint().m_time - queryLen/2;
            ori->getPartialTrajectory(ts, ts + queryLen, *concate);
            if (!concate->m_points.empty())
                queries.emplace_back(concate);
        }
        for (int nnk = 10; nnk <= 100; nnk += 10) {
            std::cerr << "nnk is" << nnk << "\n";
            maxseg = 0;
            segs1.clear();
            segs2.clear();
            emptyseg.clear();
            for (auto &traj:trajs) {
                auto seg = traj.second.getFixedSegments();
                maxseg = std::max(int(seg.size()), maxseg);
                segs1.emplace_back(make_pair(traj.first, seg));
            }
            double segpara1 = biSearchMax(nnk, 1800, 40, false, 0.012);
            std::cerr << "segpara is estimated as " << segpara1;
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

            TrajStore *ts1 = new TrajStore(name1, diskfile1, 4096, maxseg + 1);
            ts1->loadSegments(segs1);
            segs1.clear();
            segs1.swap(emptyseg);
            ISpatialIndex *r = MBCRTree::createAndBulkLoadNewRTreeWithTrajStore(ts1, 4096, 3, indexIdentifier1);

            TrajStore *ts2 = new TrajStore(name2, diskfile2, 4096, maxseg + 1);
            ts2->loadSegments(segs2, true);
            segs2.clear();
            segs2.swap(emptyseg);
            ISpatialIndex *rc = MBCRTree::createAndBulkLoadNewMBCRTreeWithTrajStore(ts2, 4096, 3, indexIdentifier2);

            bUsingSimp = bUsingSBBD = false;
            kNNQueryBatch(r, queries, ts1, nnk);
            bUsingSimp = bUsingSBBD = true;
            kNNQueryBatch(rc, queries, ts2, nnk);
            cerr << "\n";
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