//
// Created by Chuang on 2020/5/27.
//
#include "testFuncs.h"

int main(){
    try {
        calcuTime[0] = 0;
        srand((int) time(NULL));
        //vector<pair<id_type, Trajectory> > trajs = loadGLToTrajs("/root/TD.csv");
        vector<pair<id_type, Trajectory> > trajs = loadGLToTrajs("D://simp.csv");
//        vector<pair<id_type, Trajectory> > trajs = loadGTFolder();
        auto stat=trajStat::instance();
        int maxseg = 0;
        double queryLen=1800;
        double queryRadParas[]={0.02,0.04,0.06,0.08,0.1,
                                0.12,0.14,0.16,0.18,0.2,
                                0.22,0.24,0.26,0.28,0.3,
                                0.32,0.34,0.36,0.38,0.4,
                                0.42,0.44,0.46,0.48,0.5};
        std::cerr<<"\n";
//        double seglenParas[]={20,50,80,120,160,200,300,500};// for GeoLife
//        double seglenParas[]={100,150,200,250,300,400,500,600,700,800,900,1000};//for td
//        maxseg=std::max(maxseg,Trajectory::cutTrajsIntoFile(trajs,1,1,"line.stj"));
//        maxseg=std::max(maxseg,Trajectory::cutTrajsIntoFile(trajs,170,1,"tb.stj"));
//        for(double seglen:seglenParas){
//            maxseg=std::max(maxseg,Trajectory::cutTrajsIntoFile(trajs,seglen,0,"opts"+to_string(seglen)+".stj"));
//        }

        string name0 ="name0", name1 ="name1", name2 = "name2";
        id_type indexIdentifier0, indexIdentifier1, indexIdentifier2;
        IStorageManager
                *diskfile1 = StorageManager::createNewDiskStorageManager(name1, 4096),
                *diskfile2 = StorageManager::createNewDiskStorageManager(name2, 4096);
        std::cerr<<"maxseg is "<<maxseg<<"\n";


        vector<vector<IShape *>> querySet;
        for(double radlen:queryRadParas) {
            vector<IShape *> queries;
            for (int i = 0; i < 10000; i++) {
                double t = int(random(stat->mint, stat->maxt - queryLen));
                double pLow[3] = {random(stat->minx, stat->maxx), random(stat->miny, stat->maxy), t};
//                double pHigh[3] = {pLow[0] + random(0.05,0.1),
//                                   pLow[1] + random(0.05,0.1), t + queryLen};
//                Region *rg = new Region(pLow, pHigh, 3);
                Cylinder *rg = new Cylinder(pLow, radlen, t, t + queryLen, 2);
                queries.emplace_back(rg);
            }
            querySet.emplace_back(queries);
            queries.clear();
        }

//        for(double seglen:seglenParas){
//            IStorageManager *diskfile0 = StorageManager::createNewDiskStorageManager(name0, 4096);
//            TrajStore *ts0 = new TrajStore(name0, diskfile0, 4096, maxseg);
//            ts0->loadSegments("./opts"+to_string(seglen)+".stj",true);
//            ISpatialIndex *rc = MBCRTree::createAndBulkLoadNewMBCRTreeWithTrajStore(ts0, 4096, 3, indexIdentifier0);
//            ts0->flush();
//            for(auto qs: querySet){
//                rangeQueryBatch(rc, qs, ts0);
//                ts0->flush();
//            }
//            delete rc;
//            delete ts0;
//            delete diskfile0;
//        }
//        for (double radlen:queryRadParas) {
//            vector<IShape *> queries;
//            for (int i = 0; i < 5000; i++) {
//                double t = int(random(stat->mint, stat->maxt - queryLen));
//                double pLow[3] = {random(stat->minx, stat->maxx), random(stat->miny, stat->maxy), t};
////                double pHigh[3] = {pLow[0] + random(0.05,0.1),
////                                   pLow[1] + random(0.05,0.1), t + queryLen};
////                Region *rg = new Region(pLow, pHigh, 3);
//                Cylinder *rg = new Cylinder(pLow, radlen, t, t + queryLen, 2);
//                queries.emplace_back(rg);
//            }
//            std::cerr << "rad len:" << radlen << "\n";
////            rangeQueryBatch(r1, queries, ts1);
////            rangeQueryBatch(r2, queries, ts2);
//            for(double seglen:seglenParas){
//                IStorageManager *diskfile0 = StorageManager::createNewDiskStorageManager(name0, 4096);
//                TrajStore *ts0 = new TrajStore(name0, diskfile0, 4096, maxseg);
//                ts0->loadSegments("./opts"+to_string(seglen)+".stj",true);
//                ISpatialIndex *rc = MBCRTree::createAndBulkLoadNewMBCRTreeWithTrajStore(ts0, 4096, 3, indexIdentifier0);
//                ts0->flush();
//                rangeQueryBatch(rc, queries, ts0);
//                delete rc;
//                delete ts0;
//                delete diskfile0;
//            }
//            for (auto s:queries) {
//                delete (s);
//            }
//        }
//            std:cerr<<*r<<*rc<<endl;
        TrajStore *ts1 = new TrajStore(name1, diskfile1, 4096, maxseg);
        TrajStore *ts2 = new TrajStore(name2, diskfile2, 4096, maxseg);

        ts1->loadSegments("./line.stj",true);
        ISpatialIndex *r1 = MBCRTree::createAndBulkLoadNewMBCRTreeWithTrajStore(ts1, 4096, 3, indexIdentifier1);

        ts2->loadSegments("./tb.stj", true);
        ISpatialIndex *r2 = MBCRTree::createAndBulkLoadNewRTreeWithTrajStore(ts2, 4096, 3, indexIdentifier2);
        ts1->flush();
        ts2->flush();
        for(auto qs: querySet){
            rangeQueryBatch(r1, qs, ts1);
            ts1->flush();
            rangeQueryBatch(r2, qs, ts2);
            ts2->flush();
        }
        for(auto &qs:querySet){
            for(auto &shape:qs){
                delete(shape);
            }
        }
        delete r1;
        delete r2;
        delete ts1;
        delete ts2;
        delete diskfile1;
        delete diskfile2;
    }
    catch (Tools::Exception& e)
    {
        cerr << "******ERROR******" << endl;
        std::string s = e.what();
        cerr << s << endl;
        return -1;
    }
    catch(...){

    }
    return 0;
}