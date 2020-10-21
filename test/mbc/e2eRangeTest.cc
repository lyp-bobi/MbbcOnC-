//
// Created by Chuang on 2020/5/27.
//
#include "testFuncs.h"



int main() {
    bool splitfile = true;

    try {
        calcuTime[0] = 0;
        srand((int) time(NULL));
        vector<pair<id_type, Trajectory> > trajs = loadDumpedFiledToTrajs("/root/tdfilter.txt");
//        vector<pair<id_type, Trajectory> > trajs = loadDumpedFiledToTrajs("D://TRI-framework/dumpedtraj.txt");

        auto stat = trajStat::instance();
        int maxseg = 0;
        double queryLen = 3600;
        double queryRadParas[] = {0.02, 0.04, 0.06, 0.08, 0.1};
        std::cerr << "\n";
//        double seglenParas[]={20,50,80,120,160,200,300,500};// for GeoLife
        double seglenParas[]={100,150,200,250,300,400,500,600,700,800,900,1000};//for td
        if(splitfile) {
            maxseg = std::max(maxseg, Trajectory::cutTrajsIntoFile(trajs, 1, 1, "tdline.stj"));
            maxseg = std::max(maxseg, Trajectory::cutTrajsIntoFile(trajs, 170, 1, "tdtb.stj"));
            for (double seglen:seglenParas) {
                maxseg = std::max(maxseg, Trajectory::cutTrajsIntoFile(trajs, seglen, 0,
                                                                       "tdopts" + to_string(seglen) + ".stj"));
            }
        }
        string name0 = "name0", name1 = "name1", name2 = "name2";
        id_type indexIdentifier0, indexIdentifier1, indexIdentifier2;
        IStorageManager
                *diskfile1 = StorageManager::createNewDiskStorageManager(name1, 4096),
                *diskfile2 = StorageManager::createNewDiskStorageManager(name2, 4096);
        std::cerr << "maxseg is " << maxseg << "\n";


        vector<vector<IShape *>> querySet;
        for (double radlen:queryRadParas) {
            vector<IShape *> queries;
            for (int i = 0; i < 2*testtime; i++) {
                auto p =trajs[int(random(0,trajs.size()-1))].second.randomPoint();
                Cylinder *rg = new Cylinder(p.m_pCoords, radlen, p.m_time - queryLen/2, p.m_time + queryLen/2, 2);
                queries.emplace_back(rg);
                queries.emplace_back(rg);
            }
            querySet.emplace_back(queries);
            queries.clear();
        }
        bUsingSBBD = true;
        for(double seglen:seglenParas){
            stat->bt=seglen;
            IStorageManager *diskfile0 = StorageManager::createNewDiskStorageManager(name0, 4096);
            TrajStore *ts0 = new TrajStore(name0, diskfile0, 4096, maxseg);
            ts0->loadSegments("./tdopts"+to_string(seglen)+".stj",true);
            ISpatialIndex *rc = MBCRTree::createAndBulkLoadNewMBCRTreeWithTrajStore(ts0, 4096, 3, indexIdentifier0);
            ts0->flush();
            for(auto qs: querySet){
                rangeQueryBatch(rc, qs, ts0);
                ts0->flush();
            }
            delete rc;
            delete ts0;
            delete diskfile0;
        }
        TrajStore *ts1 = new TrajStore(name1, diskfile1, 4096, maxseg);
        TrajStore *ts2 = new TrajStore(name2, diskfile2, 4096, maxseg);

        stat->bt = stat->tl;
        ts1->loadSegments("./tdline.stj", true);
        ISpatialIndex *r1 = MBCRTree::createAndBulkLoadNewMBCRTreeWithTrajStore(ts1, 4096, 3, indexIdentifier1);
        stat->bt = 170*stat->tl;
        ts2->loadSegments("./tdtb.stj", true);
        ISpatialIndex *r2 = MBCRTree::createAndBulkLoadNewRTreeWithTrajStore(ts2, 4096, 3, indexIdentifier2);
        ts1->flush();
        ts2->flush();
        bUsingSBBD = false;
        for (auto qs: querySet) {
            rangeQueryBatch(r1, qs, ts1);
            ts1->flush();
            rangeQueryBatch(r2, qs, ts2);
            ts2->flush();
        }
        for (auto &qs:querySet) {
            for (auto &shape:qs) {
                delete (shape);
            }
        }
        delete r1;
        delete r2;
        delete ts1;
        delete ts2;
        delete diskfile1;
        delete diskfile2;
    }
    catch (Tools::Exception &e) {
        cerr << "******ERROR******" << endl;
        std::string s = e.what();
        cerr << s << endl;
        return -1;
    }
    catch (...) {

    }

    try {
        calcuTime[0] = 0;
        srand((int) time(NULL));
        vector<pair<id_type, Trajectory> > trajs = loadDumpedFiledToTrajs("/root/glfilter.txt");
        auto stat = trajStat::instance();
        int maxseg = 0;
        double queryLen = 3600;
        double queryRadParas[] = {0.02, 0.04, 0.06, 0.08, 0.1};
        std::cerr << "\n";
        double seglenParas[]={20,50,80,120,160,200,300,500};// for GeoLife
    //    double seglenParas[]={100,150,200,250,300,400,500,600,700,800,900,1000};//for td
        if(splitfile) {
            maxseg = std::max(maxseg, Trajectory::cutTrajsIntoFile(trajs, 1, 1, "glline.stj"));
            maxseg = std::max(maxseg, Trajectory::cutTrajsIntoFile(trajs, 170, 1, "gltb.stj"));
            for (double seglen:seglenParas) {
                maxseg = std::max(maxseg,
                                  Trajectory::cutTrajsIntoFile(trajs, seglen, 0, "glopts" + to_string(seglen) + ".stj"));
            }
        }
        string name0 = "name0", name1 = "name1", name2 = "name2";
        id_type indexIdentifier0, indexIdentifier1, indexIdentifier2;
        IStorageManager
                *diskfile1 = StorageManager::createNewDiskStorageManager(name1, 4096),
                *diskfile2 = StorageManager::createNewDiskStorageManager(name2, 4096);
        std::cerr << "maxseg is " << maxseg << "\n";


        vector<vector<IShape *>> querySet;
        for (double radlen:queryRadParas) {
            vector<IShape *> queries;
            for (int i = 0; i < 2*testtime; i++) {
                auto p =trajs[int(random(0,trajs.size()-1))].second.randomPoint();
                Cylinder *rg = new Cylinder(p.m_pCoords, radlen, p.m_time - queryLen/2, p.m_time + queryLen/2, 2);
                queries.emplace_back(rg);
            }
            querySet.emplace_back(queries);
            queries.clear();
        }

        bUsingSBBD = true;
        for(double seglen:seglenParas){
            stat->bt = seglen;
            IStorageManager *diskfile0 = StorageManager::createNewDiskStorageManager(name0, 4096);
            TrajStore *ts0 = new TrajStore(name0, diskfile0, 4096, maxseg);
            ts0->loadSegments("./glopts"+to_string(seglen)+".stj",true);
            ISpatialIndex *rc = MBCRTree::createAndBulkLoadNewMBCRTreeWithTrajStore(ts0, 4096, 3, indexIdentifier0);
            ts0->flush();
            for(auto qs: querySet){
                rangeQueryBatch(rc, qs, ts0);
                ts0->flush();
            }
            delete rc;
            delete ts0;
            delete diskfile0;
        }

        TrajStore *ts1 = new TrajStore(name1, diskfile1, 4096, maxseg);
        TrajStore *ts2 = new TrajStore(name2, diskfile2, 4096, maxseg);
        stat->bt=stat->tl;
        ts1->loadSegments("./glline.stj", true);
        ISpatialIndex *r1 = MBCRTree::createAndBulkLoadNewMBCRTreeWithTrajStore(ts1, 4096, 3, indexIdentifier1);
        stat->bt = 170*stat->tl;
        ts2->loadSegments("./gltb.stj", true);
        ISpatialIndex *r2 = MBCRTree::createAndBulkLoadNewRTreeWithTrajStore(ts2, 4096, 3, indexIdentifier2);
        ts1->flush();
        ts2->flush();
        bUsingSBBD=false;
        for (auto qs: querySet) {
            rangeQueryBatch(r1, qs, ts1);
            ts1->flush();
            rangeQueryBatch(r2, qs, ts2);
            ts2->flush();
        }
        for (auto &qs:querySet) {
            for (auto &shape:qs) {
                delete (shape);
            }
        }
        delete r1;
        delete r2;
        delete ts1;
        delete ts2;
        delete diskfile1;
        delete diskfile2;
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