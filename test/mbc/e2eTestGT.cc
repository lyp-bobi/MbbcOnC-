//
// Created by Chuang on 2019/9/2.
//

//
// Created by Chuang on 2019/9/1.
//

#include "testFuncs.h"
#include "dirent.h"


int main() {
    try {
        calcuTime[0] = 0;
        srand(0);
        bUsingSimp = true;
        string name0 = "name0", name1 = "name1", name2 = "name2";
        id_type indexIdentifier0, indexIdentifier1, indexIdentifier2;
        IStorageManager *diskfile0 = StorageManager::createNewDiskStorageManager(name0, 4096),
                *diskfile1 = StorageManager::createNewDiskStorageManager(name1, 4096),
                *diskfile2 = StorageManager::createNewDiskStorageManager(name2, 4096);
        TrajStore *ts1 = new TrajStore(name1, diskfile1, 4096, 500);
        TrajStore *ts2 = new TrajStore(name2, diskfile2, 4096, 500);
        vector<pair<id_type, Trajectory> > trajs;
        vector<IShape *> queries;
        for (int dataSize = 1; dataSize <= 50; dataSize+=5) {
            trajs = loadGTFolder(dataSize);
            vector<pair<id_type, vector<Trajectory>>> segs1;
            vector<pair<id_type, vector<Trajectory>>> segs2;
            for (auto &traj:trajs) {
                auto seg = traj.second.getFixedSegments();
                segs1.emplace_back(make_pair(traj.first, seg));
            }
            double segpara1 = biSearchMax(5, 122, 40, false);
            segpara1 = std::max(segpara1, 20.0);
            std::cerr << "query len:" << 122 << ",partial traj len:" << segpara1 << "\n";
            for (auto &traj:trajs) {
                auto seg = traj.second.getSegments(segpara1);
                segs2.emplace_back(make_pair(traj.first, seg));
            }
            ts1->loadSegments(segs1, true);
            ts2->loadSegments(segs2, true);
            segs1.clear();
            segs2.clear();

            long count = ts1->m_entries.size();
            ISpatialIndex *r = MBCRTree::createAndBulkLoadNewRTreeWithTrajStore(ts1, 4096, 3, indexIdentifier1);
            ISpatialIndex *rc = MBCRTree::createAndBulkLoadNewMBCRTreeWithTrajStore(ts2, 4096, 3, indexIdentifier2);
            ts1->flush();
            ts2->flush();
            drop_cache(1);
            if (dataSize == 1) {
                for (int i = 0; i < testtime; i++) {
                    Trajectory *t = new Trajectory();
                    *t= trajs[random(0,trajs.size())].second;
                    queries.emplace_back(t);
                }
            }
            std::cerr << dataSize << endl;
            bUsingSBBD=bUsingSimp=false;
            kNNQueryBatch(r, queries, ts1);
            bUsingSBBD=bUsingSimp=true;
            kNNQueryBatch(rc, queries, ts2);
            delete r;
            delete rc;
        }
        for (auto &tras:queries) {
            delete tras;
        }
        queries.clear();
    }
    catch (Tools::Exception &e) {
        cerr << "******ERROR******" << endl;
        std::string s = e.what();
        cerr << s << endl;
        return -1;
    }
    return 0;
}