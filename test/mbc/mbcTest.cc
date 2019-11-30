#include "testFuncs.h"

int main(){
    try {
        calcuTime[0]=0;
        srand((int) time(NULL));
//        srand(21);
//        vector<pair<id_type, Trajectory> > trajs = loadGTToTrajs();
        vector<pair<id_type, Trajectory> > trajs = loadGLToTrajs("D://simp.csv");
        vector<pair<id_type, vector<Trajectory>>> segs;
        vector<pair<id_type, Trajectory> > empty1;
        int totallen=0,totalseg=0;
        int maxseg=1000;
        for(auto &traj:trajs){
            totallen+=traj.second.m_points.size();
            auto seg= traj.second.getSegments(3000);
            totalseg+= seg.size();
            maxseg=std::max(int(seg.size()),maxseg);
            segs.emplace_back(make_pair(traj.first,seg));
        }
        std::cerr<<"segments' average length is "<<totallen*1.0/totalseg<<"\n";
        vector<IShape *> queries;
        cerr<<"generating queries\n";
        int realtesttime=(QueryType==2)?testtime:(100*testtime);
        double queryLen=10800;
        auto stat=trajStat::instance();
//        double mp[2]={38.021442825891902828,114.06330331672108969};
//        Cylinder cy(mp,0.041490676595355088785,55307,66107,2);
//        queries.emplace_back(&cy);
//        for (int i = 0; i < 10000; i++) {
//            double t = int(random(stat->mint, stat->maxt - queryLen));
//            double pLow[3] = {random(stat->minx, stat->maxx), random(stat->miny, stat->maxy), t};
////                double pHigh[3] = {pLow[0] + random(0.05,0.1),
////                                   pLow[1] + random(0.05,0.1), t + queryLen};
////                Region *rg = new Region(pLow, pHigh, 3);
//            Cylinder *rg = new Cylinder(pLow,random(0.025,0.05),t,t+queryLen,2);
//            queries.emplace_back(rg);
//        }
        for (int i = 0; i < realtesttime; i++) {
            if (QueryType == 1) {
                double t = int(random(0, 1000));
                double pLow[3] = {random(0, 25000), random(0, 30000), t};
                double pHigh[3] = {pLow[0] + random(500, 2000), pLow[1] + random(500, 2000), t+20};
                Region *rg = new Region(pLow, pHigh, 3);
                queries.emplace_back(rg);
            }
            else if(QueryType==2){
                auto ori = &trajs[0].second;
//                auto ori = &trajs[(int(random(0, trajs.size()))) % trajs.size()].second;
                Trajectory *concate = new Trajectory();
                double ts=ori->m_startTime();
//                double ts = std::max(ori->m_startTime(),random(ori->m_startTime(), ori->m_endTime() - queryLen));
                ori->getPartialTrajectory(ts, ts + queryLen, *concate);
                if (!concate->m_points.empty())
                    queries.emplace_back(concate);
            }
        }

        cerr<<"queries generated\n";
        trajs.swap(empty1);
        string name0 = "name0", name1 = "name1", name2 = "name2";
        id_type indexIdentifier0, indexIdentifier1, indexIdentifier2;
        IStorageManager *diskfile0 = StorageManager::createNewDiskStorageManager(name0, 4096),
                *diskfile1 = StorageManager::createNewDiskStorageManager(name1, 4096),
                *diskfile2 = StorageManager::createNewDiskStorageManager(name2, 4096);
        // Create a new storage manager with the provided base name and a 4K page size.
        StorageManager::IBuffer *file0 = StorageManager::createNewRandomEvictionsBuffer(*diskfile0, 10, false),
                *file1 = StorageManager::createNewRandomEvictionsBuffer(*diskfile1, 10, false),
                *file2 = StorageManager::createNewRandomEvictionsBuffer(*diskfile2, 10, false);

//        StorageManager::DiskStorageManager *dsm1,*dsm2;
//        dsm1= dynamic_cast<StorageManager::DiskStorageManager*>(diskfile1);
//        dsm2= dynamic_cast<StorageManager::DiskStorageManager*>(diskfile2);

        TrajStore* ts1=new TrajStore(file1,4096,maxseg);
        ts1->loadSegments(segs);

        ISpatialIndex *r=RTree::createAndBulkLoadNewRTreeWithTrajStore(ts1,4096,3,indexIdentifier1);

//        TreeQueryBatch(r, queries,ts1);
//        delete r;
//        delete ts1;


        TrajStore* ts2=new TrajStore(file2,4096,maxseg);
        ts2->loadSegments(segs,true);
        ISpatialIndex *rc=MBCRTree::createAndBulkLoadNewMBCRTreeWithTrajStore(ts2,4096,3,indexIdentifier2);
//        TreeQueryBatch(rc, queries,ts2);
//        delete rc;
//        delete ts2;

        segs.clear();

//        TrajMbrStream ds1;
//        ds1.feedTraj(&trajs);
//        ds1.rewind();
//        ISpatialIndex* real = RTree::createAndBulkLoadNewRTree(
//                RTree::BulkLoadMethod::BLM_STR, ds1, *file0, 0.9, 10000,100000, 3,RTree::RV_RSTAR, indexIdentifier0);
//        real->m_DataType=TrajectoryType;





//        id_type that=242800;
//        cout<<ts1.getTraj(that);
//        cout<<trajs[0].second;
//        for(int i=11;i<12;i++) {
//            id_type that = i*100;
//            auto brs = ts1.getMBRsByTime(that,0,1000);
//            auto bcs = ts1.getMBCsByTime(that,0,1000);
////        cout<<trajs[10].second;
////            cout << trajs[71].second << endl;
////            cout << brs << endl;
////            cout << bcs << endl;
//            cout<< i<<endl;
//            double a,b,c;
//            a=queries[0]->getMinimumDistance(trajs[i].second);
//            b=queries[0]->getMinimumDistance(brs);
//            c=queries[0]->getMinimumDistance(bcs);
//            cout<<a<<" "<<b<<" "<<c<<endl;
//            if((b>a||c>a||std::isnan(b)||std::isnan(c))&&a!=std::numeric_limits<double>::max()) {
//                cout << queries[0]->getMinimumDistance(trajs[i].second) << endl;
//                cout << queries[0]->getMinimumDistance(brs) << endl;
//                cout << queries[0]->getMinimumDistance(bcs) << endl;
//            }
//        }




        cerr << "start query!" << endl << endl << endl;

//        TreeQueryBatch(real,queries);

        cerr.precision(20);
        double aa,bb,oo;
//        for(int j=0;j<queries.size();j++){
        for(int j=0;j<queries.size();j++){
            auto q=queries[j];
//            oo=TreeQuery(real,q);
            Trajectory *qtraj= dynamic_cast<Trajectory*>(q);
//            std::cout<<qtraj->m_startTime()<<" "<<qtraj->m_endTime()<<"\n";
            aa=TreeQuery(r,q,ts1);
            bb=TreeQuery(rc,q,ts2);
            if(aa!=bb){
//                Cylinder *qcy= dynamic_cast<Cylinder*>(q);
                cerr<<"error"<<j<<endl;
                cerr<<aa<<" "<<bb<<" "<<*qtraj<<endl;//<<*qcy<<endl;
                return 1;
//                cerr<<*qtraj<<endl;
            }
        }



//        std::cerr<<"index IO:"<<ts1.m_indexIO<<" "<<ts2.m_indexIO<<endl;
//        std::cerr<<"traj IO:"<<ts1.m_trajIO<<" "<<ts2.m_trajIO<<endl;
//        std::cerr<<"total IO:"<<ts1.m_indexIO+ts1.m_trajIO<<" "<<ts2.m_indexIO+ts2.m_trajIO<<endl;
//        std::cerr<<"bounding IO:"<<ts1.m_boundingVisited<<" "<<ts2.m_boundingVisited<<endl;
//        std::cerr<<"IO time:"<<ts1.m_IOtime<<" "<<ts2.m_IOtime<<"\n";
//        std::cerr<<"calculation time"<<calcuTime[0]<<" "<<calcuTime[1]<<"\n";
        delete file0;delete file1;delete file2;
        delete diskfile0;delete diskfile1;delete diskfile2;
    }
    catch (Tools::Exception& e)
    {
        cerr << "******ERROR******" << endl;
        std::string s = e.what();
        cerr << s << endl;
        return -1;
    }
    return 0;
}