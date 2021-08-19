//
// Created by Chuang on 2019/6/10.
//
#include "storagemanager/xStore.h"
#include "DiskStorageManager.h"
#include <cmath>
#include <cstring>
#include <sys/stat.h>
#include <spatialindex/tools/rand48.h>

#define fp (int(m_pageSize/sizeof(prexp)/3))

bool CheckFilesExists(const string &s) {
    struct stat stats;
    int ret = 1;
    std::ostringstream os;
    std::string name;
    do {

        os << s << ".idx"; //diskfileidx
        name = os.str();
        ret = stat(name.c_str(), &stats);
        if (ret != 0) break;

        os.str("");
        os << s << ".dat";//diskfiledat
        name = os.str();
        ret = stat(name.c_str(), &stats);
        if (ret != 0) break;

        os.str("");
        os << s << ".sqlite";//sqlite file
        name = os.str();
        ret = stat(name.c_str(), &stats);
        if (ret != 0) break;

        os.str("");
        os << s << ".property"; //json property
        name = os.str();
        ret = stat(name.c_str(), &stats);
        if(ret!=0) break;
    } while(false);

    if (ret != 0){
//        os.str("");
//        os << s << ".sqlite"; //json property
//        name = os.str();
//        remove(name.c_str());
        return false;
    }

    return true;
}

xTrajEntry::xTrajEntry(id_type page, uint32_t len)
        : m_page(page), m_npoint(len) {}

std::string xTrajEntry::xTrajEntry::toString() {
    string s;
    s = std::to_string(m_page) + std::to_string(m_npoint);
    return s;
}

xTrajEntry::xTrajEntry::xTrajEntry(string &s) {
    auto ws = split(s, ',');
    m_page = std::stoll(ws[0]);
    m_npoint = std::stoi(ws[1]);
}

xStore::~xStore() {
    for(auto &it:m_pagecache) delete it.second;
    delete m_pStorageManager;
}

static int getLastId(string s)
{
    string filename = s;
    ifstream fin;

    fin.open(filename);
    if(fin.is_open()) {
        fin.seekg(0,ios_base::end);                // go to one spot before the EOF
        for(int i=0;i<2;i++) {
            fin.seekg(-2,ios_base::cur);
            bool keepLooping = true;
            while (keepLooping) {
                char ch;
                ch = fin.peek();                            // Get current byte's data
                std::cerr<<ch;
                if ((long long) fin.tellg() <= 1) {             // If the data was at or before the 0th byte
                    fin.seekg(0);                       // The first line is the last line
                    keepLooping = false;                // So stop there
                } else if (ch == '\n') {                   // If the data was a newline
                    keepLooping = false;                // Stop at the current position.
                } else {                                  // If the data was neither a newline nor at the 0 byte
                    fin.seekg(-1,ios_base::cur);
                    // Move to the front of that data, then to the front of the data before it
                }
            }
        }
        string lastLine;
        getline(fin,lastLine);                      // Read the current line
        getline(fin,lastLine);                      // Read the current line
        fin.close();
        return stoll(lastLine);
    }

    return 0;
}

void xStore::loadFile(string file) {
    std::cerr << "start loading " << file << " into " << m_name << "\n";
    string filename = filedirprefix+m_name;
    conn_init(filename+".sqlite");
    if(m_pStorageManager== nullptr){
        if (CheckFilesExists(filename)) {
            m_pStorageManager = loadDiskStorageManager(filename);
        }else{
            m_pStorageManager = createNewDiskStorageManager(filename,m_pageSize);
        }
    }
    ifstream inFile(file, ios::in);
    string lineStr;
    xTrajectory tj;
    xMBR r;
    uint8_t *data;
    uint32_t len;
    int curLine = 0;
    while (getline(inFile, lineStr)) {
        try {
            string str;
            stringstream ss(lineStr);
            getline(ss, str);
            id_type id;
            id = stoll(str);
            if(id>10) {
                id = db_last_trajid()+1;
            }
            getline(inFile, str);
            tj.loadFromString(str);
            //test code
//                std::cerr<<"test id"<<id<<endl;
            if (tj.m_points.size() >= 2) {
                int rem = tj.m_points.size();
                int cur = 0;
                uint8_t *ptr;
                bool isfirst = true;
                id_type firstpage;
                while (rem > 0) {
                    int plen = min(rem, fp);
                    xPoint *p;
                    data = new uint8_t[24 * plen];
                    ptr = data;
                    for (int i = 0; i < plen; i++) {
                        //test code
//                            std::cerr<<" "<<cur<< tj.m_points[cur]<<endl;
                        p = &(tj.m_points[cur++]);
                        p->storeToByteArrayE(&ptr, len);
                        if (i != plen - 1) ptr += len;
                    }
                    ptr = data;
                    id_type newPage = StorageManager::NewPage;
                    m_pStorageManager->storeByteArray(newPage, 24 * plen, data);
                    if (isfirst) {
                        firstpage = newPage;
                        isfirst = false;
                    }
                    rem -= fp;
                    if (m_bSubTraj && rem != 0) {
                        rem++;
                        cur--;
                    }
                    delete[] data;
                }
                db_insert_traj(id,firstpage,tj.m_points.size());
                curLine++;
                tj.getxMBR(r);
                if (r.m_xmax > tjstat->maxx) tjstat->maxx = r.m_xmax;
                if (r.m_xmin < tjstat->minx) tjstat->minx = r.m_xmin;
                if (r.m_ymax > tjstat->maxy) tjstat->maxy = r.m_ymax;
                if (r.m_ymin < tjstat->miny) tjstat->miny = r.m_ymin;
                if (r.m_tmax > tjstat->maxt) tjstat->maxt = r.m_tmax;
                if (r.m_tmin < tjstat->mint) tjstat->mint = r.m_tmin;
                tjstat->dist += tj.m_dist();
                tjstat->lineCount += tj.m_points.size() - 1;
                tjstat->trajCount += 1;
                if (tjstat->vmax < tj.maxSpeed()) {
                    tjstat->vmax = tj.maxSpeed();
                }
                tjstat->M += tj.m_endTime() - tj.m_startTime();
            }
        }
        catch (...) {
            break;
        }
    }
    std::cerr << "load finished\n";
    tjstat->Dx = tjstat->maxx - tjstat->minx;
    tjstat->Dy = tjstat->maxy - tjstat->miny;
    tjstat->Dt = tjstat->maxt - tjstat->mint;
    tjstat->tl = tjstat->M / tjstat->lineCount;
    tjstat->jt = tjstat->M / tjstat->trajCount;
    tjstat->v = tjstat->dist / tjstat->M;
    tjstat->Sr = (tjstat->Dx + tjstat->Dy) / 2;
    tjstat->P = tjstat->Dt;
    std::cerr << file << endl;
    inFile.close();
    if (file.find("tdexpand") != file.npos) tjstat->usedata("tdexpand");
    else if (file.find("glexpand") != file.npos) tjstat->usedata("glexpand");
    else if (file.find("td") != file.npos) tjstat->usedata("td");
    else if (file.find("gl") != file.npos) tjstat->usedata("gl");
    std::cerr << tjstat->toString() << endl;
    m_property["tjstat"] = tjstat->toString();
    flush();
}

xStore::xStore(string myname, string file, bool bsubtraj, bool forceNew)
{
    cerr<<myname<<endl;
    m_name = myname;
    m_pageSize = 4096;
    string filename = filedirprefix+myname;
    conn_init(filename+".sqlite");
    if ((!forceNew) && CheckFilesExists(filename)) {
//        std::cerr << "using existing xStore " << myname << endl;
        m_pStorageManager = loadDiskStorageManager(filename);
        ifstream propFile(filedirprefix+myname + ".property", ios::in);
        propFile >> m_property;
        propFile.close();
        if(tjstat->lineCount==0) { /* only load it when not inited*/
            tjstat->fromString(m_property["tjstat"]);
        }
        m_bSubTraj = m_property["bSubTraj"];
        ifstream trajidxFile(filedirprefix+m_name + ".trajidx", ios::in);
    } else {
        if(m_name=="od") {
            m_bSubTraj = bsubtraj;
            m_property["trajfile"] = file;
            m_property["tjstat"] = tjstat->toString();
            m_property["bSubTraj"] = bsubtraj;
            return;
        }
        m_bSubTraj = bsubtraj;
        std::cerr << "start loading " << file << " into " << myname << "\n";
        m_pStorageManager = createNewDiskStorageManager(filename, m_pageSize);
        m_property["trajfile"] = file;
        tjstat->init();
        ifstream inFile(file, ios::in);
        string lineStr;
        xTrajectory tj;
        xMBR r;
        uint8_t *data;
        uint32_t len;
        int curLine = 0;
        while (getline(inFile, lineStr)) {
            try {
                string str;
                stringstream ss(lineStr);
                getline(ss, str);
                id_type id = stoll(str);
                getline(inFile, str);
                tj.loadFromString(str);
                //test code
//                std::cerr<<"test id"<<id<<endl;
                if (tj.m_points.size() >= 2) {
                    int rem = tj.m_points.size();
                    int cur = 0;
                    uint8_t *ptr;
                    bool isfirst = true;
                    id_type firstpage;
                    while (rem > 0) {
                        int plen = min(rem, fp);
                        xPoint *p;
                        data = new uint8_t[24 * plen];
                        ptr = data;
                        for (int i = 0; i < plen; i++) {
                            //test code
//                            std::cerr<<" "<<cur<< tj.m_points[cur]<<endl;
                            p = &(tj.m_points[cur++]);
                            p->storeToByteArrayE(&ptr, len);
                            if (i != plen - 1) ptr += len;
                        }
                        ptr = data;
                        id_type newPage = StorageManager::NewPage;
                        m_pStorageManager->storeByteArray(newPage, 24 * plen, data);
                        if (isfirst) {
                            firstpage = newPage;
                            isfirst = false;
                        }
                        rem -= fp;
                        if (bsubtraj && rem != 0) {
                            rem++;
                            cur--;
                        }
                        delete[] data;
                    }
                    db_insert_traj(id,firstpage,tj.m_points.size());
                    curLine++;
                    tj.getxMBR(r);
                    if (r.m_xmax > tjstat->maxx) tjstat->maxx = r.m_xmax;
                    if (r.m_xmin < tjstat->minx) tjstat->minx = r.m_xmin;
                    if (r.m_ymax > tjstat->maxy) tjstat->maxy = r.m_ymax;
                    if (r.m_ymin < tjstat->miny) tjstat->miny = r.m_ymin;
                    if (r.m_tmax > tjstat->maxt) tjstat->maxt = r.m_tmax;
                    if (r.m_tmin < tjstat->mint) tjstat->mint = r.m_tmin;
                    tjstat->dist += tj.m_dist();
                    tjstat->lineCount += tj.m_points.size() - 1;
                    tjstat->trajCount += 1;
                    if (tjstat->vmax < tj.maxSpeed()) {
                        tjstat->vmax = tj.maxSpeed();
                    }
                    tjstat->M += tj.m_endTime() - tj.m_startTime();
                }
            }
            catch (...) {
                break;
            }
        }
        std::cerr << "load finished\n";
        tjstat->Dx = tjstat->maxx - tjstat->minx;
        tjstat->Dy = tjstat->maxy - tjstat->miny;
        tjstat->Dt = tjstat->maxt - tjstat->mint;
        tjstat->tl = tjstat->M / tjstat->lineCount;
        tjstat->jt = tjstat->M / tjstat->trajCount;
        tjstat->v = tjstat->dist / tjstat->M;
        tjstat->Sr = (tjstat->Dx + tjstat->Dy) / 2;
        tjstat->P = tjstat->Dt;
        std::cerr << file << endl;
        inFile.close();
        if (file.find("tdexpand") != file.npos) tjstat->usedata("tdexpand");
        else if (file.find("glexpand") != file.npos) tjstat->usedata("glexpand");
        else if (file.find("td") != file.npos) tjstat->usedata("td");
        else if (file.find("gl") != file.npos) tjstat->usedata("gl");
        std::cerr << tjstat->toString() << endl;
        m_property["tjstat"] = tjstat->toString();
        m_property["bSubTraj"] = bsubtraj;
        flush();
    }

}

xStore * xStore::clone()
{
    xStore *res = new xStore();
    assert(!m_isro);
    res->m_name=m_name;
    res->m_pageSize=m_pageSize;
    string filename = filedirprefix+res->m_name;
    conn_init(filename+".sqlite");
    DiskStorageManager* d = dynamic_cast<DiskStorageManager *>(m_pStorageManager);
    res->m_pStorageManager = new DiskStorageManager(*d,filename);
    ifstream propFile(filedirprefix+res->m_name + ".property", ios::in);
    propFile >> res->m_property;
    propFile.close();
    res->m_bSubTraj = res->m_property["bSubTraj"];
    res->m_isro = true;
    return res;
}

// for inner node
void xStore::storeByteArray(id_type &page, const uint32_t len, const uint8_t *const data) {
    m_pStorageManager->storeByteArray(page,len,data);
    auto it = m_pagecache.find(page);
    if(it != m_pagecache.end())
    {
        delete (*it).second;
        m_pagecache.erase(it);
    }
}

void xStore::drop_one_cache() {
    while(true) {
        double random;
        random = drand48();

        uint32_t entry = static_cast<uint32_t>(floor(
                ((double) m_pagecache.size()) * random));

        auto it2 = m_pagecache.begin();
        for (uint32_t cIndex = 0; cIndex < entry; cIndex++) ++it2;
        if (it2->second->m_bUsing == false) {
            delete (*it2).second;
            m_pagecache.erase(it2);
            return;
        }
    }
}

//for inner nodes
void xStore::loadByteArray(const id_type page, uint32_t &len, uint8_t **data) {
    auto start = std::chrono::system_clock::now();
    if(!m_bPageCache)
    {
        m_pStorageManager->loadByteArray(page, len, data);
    }
    else {
        auto it = m_pagecache.find(page);
        if (it != m_pagecache.end() && it->second->m_bDirty == false) {
            it->second->m_bUsing = true;
            len = (*it).second->m_length;
            *data = new uint8_t[len];
            memcpy(*data, (*it).second->m_pData, len);
            it->second->m_bUsing = false;
        } else {
            m_pStorageManager->loadByteArray(page, len, data);
            if (m_pagecache.size() >= m_buffersize) {
                drop_one_cache();
            }
            m_pagecache.insert(std::pair<id_type, PageCahce *>(page,
                                                               new PageCahce(
                                                                       len,
                                                                       *data)));
        }
    }
    auto end = std::chrono::system_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    m_IOtime+=double(duration.count()) * std::chrono::microseconds::period::num
              / std::chrono::microseconds::period::den;
}

void xStore::loadTraj(xTrajectory &out, const xStoreEntry &e) {
    xTrajEntry te = db_load_traj_entry(e.m_id);
    uint32_t ms = min(te.m_npoint - 1, e.m_s), me = min(te.m_npoint - 1, e.m_e);
    //test code
//    std::cerr<<"test id "<<e.m_id<<endl;
    out.m_points.clear();
    out.m_points.reserve(me-ms+1);
    uint8_t *data, *ptr;
    uint32_t len, tmplen;
    if (m_bSubTraj) {
        id_type pages = te.m_page + (ms) / (fp - 1);
        id_type pagee = te.m_page + int(ceil(1.0 * (me) / (fp - 1))) - 1;
        if (pages > pagee) pagee++;
        m_trajIO += pagee - pages + 1;
        id_type cur = ms / (fp - 1) * (fp - 1);
        int ps, pe;
        prexp x, y, t;
        for (auto i = pages; i <= pagee; i++) {
            ps = max(0, int(ms - cur));
            pe = min(fp - 1, int(me - cur));
            if (i != pagee && pe == fp - 1) pe--;
            len = 3 * sizeof(prexp) * (pe + 1);
            m_pStorageManager->loadByteArray(i, tmplen, &data);
            for (ptr = data + 3 * sizeof(prexp) * ps; ptr - data < len;) {
                x = *((double *) ptr);
                ptr += sizeof(prexp);
                y = *((double *) ptr);
                ptr += sizeof(prexp);
                t = *((double *) ptr);
                ptr += sizeof(prexp);
                out.m_points.emplace_back(xPoint(x, y, t));
                //test code
//                std::cerr<< out.m_points.size() << " " << x<<","<<y<<","<<t<<endl;
            }
            delete[] data;
            cur += pe + 1;
        }
        out.m_fakehead = (ms != 0);
        out.m_fakeback = (me == te.m_npoint - 1);
    } else {
        id_type pages = te.m_page + (ms + 1) / fp;
        id_type pagee = te.m_page + int(ceil(1.0 * (me + 1) / fp)) - 1;
        m_trajIO += pagee - pages + 1;
        id_type cur = ms / fp * fp;
        int ps, pe;
        prexp x, y, t;
        for (auto i = pages; i <= pagee; i++) {
            ps = max(0, int(ms - cur));
            pe = min(fp - 1, int(me - cur));
            len = 3 * sizeof(prexp) * (pe + 1);
            m_pStorageManager->loadByteArray(i, tmplen, &data);
            for (ptr = data + 3 * sizeof(prexp) * ps; ptr - data < len;) {
                x = *((double *) ptr);
                ptr += sizeof(prexp);
                y = *((double *) ptr);
                ptr += sizeof(prexp);
                t = *((double *) ptr);
                ptr += sizeof(prexp);
                out.m_points.emplace_back(xPoint(x, y, t));
            }
            cur += fp;
            delete data;
        }
        out.m_fakehead = (ms != 0);
        out.m_fakeback = (me == te.m_npoint - 1);
    }
    //simulate secondary addressing
//    int pagenum = random(0,m_pStorageManager->nextPage() - 1);
//    m_pStorageManager->loadByteArray(pagenum, tmplen, &data);
//    delete[] data;
}

xTrajectory xStore::randomSubtraj(double len) {
    id_type lastid = db_last_trajid();
    int rnd = random(1, lastid);
    int i = 0;
    xTrajectory tj,tj2;
    loadTraj(tj, xStoreEntry(rnd, 0, 1000000));
    double t = tj.m_startTime() + (tj.m_endTime()-tj.m_startTime()) * drand48();
    tj.getPartialxTrajectory(t - len / 2, t + len / 2, tj2);
    return tj2;
}

xPoint xStore::randomPoint() {
    int rnd = random(1, db_last_trajid());
    int rnd2;
    int i = 0;
    xTrajectory tj;
    loadTraj(tj, xStoreEntry(rnd, 0, 1000000));
    rnd2 = random(0, tj.m_points.size() - 1);
    return tj.m_points[rnd2];
}

void xStore::flush() {
    if(m_isro) return;
    m_pStorageManager->flush();
    ofstream propFile(filedirprefix+m_name + ".property", ios::out);
    propFile << m_property;
    propFile.close();
}


ostream & SpatialIndex::StorageManager::operator<<(ostream &os, const xTrajEntry &c) {
    os << c.m_page<<" "<<c.m_npoint;
    return os;
}

istream & SpatialIndex::StorageManager::operator>>(istream &is, xTrajEntry &c) {
    is>>c.m_page>>c.m_npoint;
    return is;
}


xSBBStream::xSBBStream(xStore *p, CUTFUNC f)
        : m_cutFunc(f),m_pstore(p) {
    m_isdb = true;
    m_size = db_last_trajid();
}

bool xSBBStream::hasNext() {
    return !m_buf.empty() || (m_isdb && m_id < m_size);
}

xSBBData *xSBBStream::getNext() {
    bool isFirst = false;
    if (m_buf.empty()) {
        xTrajectory tj;
        if(!m_isdb)
            m_pstore->loadTraj(tj, xStoreEntry(m_it->first, 0, m_it->second.m_npoint));
        else
        {
            m_pstore->loadTraj(tj, xStoreEntry(m_numit, 0, 10000000));
        }
        m_buf = m_cutFunc(tj);
        if(!m_isdb) {
            m_id = m_it->first;
            m_it++;
        }else{
            m_id = m_numit;
            m_numit ++;
        }
        isFirst = true;
    }
    auto b = m_buf.front();
    m_buf.pop();
    return new xSBBData(m_count++,
                        xStoreEntry(m_id, b.first.first, b.first.second), b.second, !isFirst, !m_buf.empty());
}

uint32_t xSBBStream::size() {
    throw Tools::NotSupportedException("xsbbstream has no size");
}

void xSBBStream::rewind() {
    m_id = 0;
    m_count = 0;
}
