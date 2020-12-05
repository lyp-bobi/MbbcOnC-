//
// Created by Chuang on 2019/8/19.
//

#include "testFuncs.h"

int main() {
    try {
        calcuTime[0] = 0;
        srand(0);
//        vector<pair<id_type, Trajectory> > trajs = loadGLToTrajs("/root/TD.csv");
        vector<pair<id_type, Trajectory> > trajs = loadDumpedFiledToTrajs("/root/tdfilter.txt");
//        vector<pair<id_type, Trajectory> > trajs = loadDumpedFiledToTrajs("D://geolifeCleaner/tdfilter.txt");

        int maxseg = 0;

        maxseg = 3000;
        string name0 = "name0";
        id_type indexIdentifier0, indexIdentifier1, indexIdentifier2;
        IStorageManager *diskfile0 = StorageManager::createNewDiskStorageManager(name0, 4096);
        // Create a new storage manager with the provided base name and a 4K page size.
        StorageManager::IBuffer *file0 = StorageManager::createNewRandomEvictionsBuffer(*diskfile0, 10, false);

        maxseg = std::max(3000, Trajectory::cutTrajsIntoFile(trajs, 2400));
        std::cerr << "maxseg is " << maxseg << "\n";

        TrajStore *ts0 = new TrajStore(name0, diskfile0, 4096, maxseg + 1);

        ts0->loadSegments(subTrajFile, true);
        ISpatialIndex *r = MBCRTree::createAndBulkLoadNewRTreeWithTrajStore(ts0, 4096, 3, indexIdentifier1);


        for(double qr = 1; qr>0.00001;qr=qr/2) {
            vector<IShape *> queries;
            for (int i = 0; i < 10000; i++) {
                auto p = trajs[int(random(0, trajs.size() - 1))].second.randomPoint();
                Cylinder *rg = new Cylinder(p.m_pCoords, qr, p.m_time, p.m_time, 2);
                queries.emplace_back(rg);
            }
            std::cerr<<qr<<"\n";
            rangeQueryBatch(r,queries,ts0);
        }
        delete r;
        delete ts0;
        delete file0;
        delete diskfile0;
    }
    catch (Tools::Exception &e) {
        cerr << "******ERROR******" << endl;
        std::string s = e.what();
        cerr << s << endl;
        return -1;
    }
    catch (...) {

    }
    return 0;
}