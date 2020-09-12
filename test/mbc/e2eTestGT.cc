//
// Created by Chuang on 2019/9/2.
//

//
// Created by Chuang on 2019/9/1.
//

#include "testFuncs.h"
#include "dirent.h"


int main(){
    try {
        vector<string> files;
#if !WIN32
        DIR *dir;
        string PATH = fileFolder;
        struct dirent **namelist;
        int n;
        n = scandir(PATH.c_str(),&namelist,0,alphasort);
        int index=-1;
        while(index<n&&files.size()<=20)
        {
            index++;
            if(namelist[index]->d_name[0] == '.')
                continue;
            //cout << ptr->d_name << endl;
            files.emplace_back(PATH+namelist[index]->d_name);
        }
#endif
//        files.emplace_back("D://00.txt");
//        files.emplace_back("D://01.txt");
        calcuTime[0] = 0;
        srand((int) time(NULL));
        rsimpli=false;
        bUsingSimp=true;
        string name0 ="name0", name1 ="name1", name2 = "name2";
        id_type indexIdentifier0, indexIdentifier1, indexIdentifier2;
        IStorageManager *diskfile0 = StorageManager::createNewDiskStorageManager(name0, 4096),
                *diskfile1 = StorageManager::createNewDiskStorageManager(name1, 4096),
                *diskfile2 = StorageManager::createNewDiskStorageManager(name2, 4096);
        // Create a new storage manager with the provided base name and a 4K page size.
        StorageManager::IBuffer *file0 = StorageManager::createNewRandomEvictionsBuffer(*diskfile0, 10, false),
                *file1 = StorageManager::createNewRandomEvictionsBuffer(*diskfile1, 10, false),
                *file2 = StorageManager::createNewRandomEvictionsBuffer(*diskfile2, 10, false);
        TrajStore *ts1 = new TrajStore(name1, file1, 4096, 500);
        TrajStore *ts2 = new TrajStore(name2, file2, 4096, 500);
        vector<pair<id_type, Trajectory> > trajs;
        vector<IShape*> queries;
        for(int dataSize=1;dataSize<=20;dataSize++){
            vector<pair<id_type, Trajectory> > tmptrajs = loadGTToTrajs(files[dataSize-1]);
            trajs.insert(trajs.begin(),tmptrajs.begin(),tmptrajs.end());
            vector<pair<id_type, vector<Trajectory>>> segs1;
            vector<pair<id_type, vector<Trajectory>>> segs2;
            for (auto &traj:trajs) {
                auto seg = traj.second.getFixedSegments();
                segs1.emplace_back(make_pair(traj.first, seg));
            }
            double segpara1=biSearchMax(5,177,40,false);
            segpara1=std::max(segpara1,20.0);
            std::cerr<<"query len:"<<177<<",partial traj len:"<<segpara1<<"\n";
            for (auto &traj:trajs) {
                auto seg = traj.second.getSegments(segpara1);
                segs2.emplace_back(make_pair(traj.first, seg));
            }
            ts1->loadSegments(segs1,true);
            ts2->loadSegments(segs2,true);
            segs1.clear();
            segs2.clear();
            auto stat=trajStat::instance();
            long count=ts1->m_entries.size();
            ISpatialIndex *r = RTree::createAndBulkLoadNewRTreeWithTrajStore(ts1, 4096, 3, indexIdentifier1);
            ISpatialIndex *rc = MBCRTree::createAndBulkLoadNewMBCRTreeWithTrajStore(ts2, 4096, 3, indexIdentifier2);
            if(dataSize==1) {
                for (int i = 0; i < 200; i++) {
                    id_type randId = long(random(0, count - 1));
                    auto it = ts1->m_entries.begin();
                    for (long j = 0; j < randId; j++) it++;
                    auto segId = it->first;
                    auto ori = ts1->getTrajByTime(segId, 0, 5000);
                    Trajectory *concate = new Trajectory(ori);
                    queries.emplace_back(concate);
                }
            }
            std::cerr<<dataSize<<endl;
            kNNQueryBatch(r, queries, ts1);
            kNNQueryBatch(rc, queries, ts2);
            delete r;
            delete rc;
        }
        for(auto &tras:queries){
            delete tras;
        }
        queries.clear();
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