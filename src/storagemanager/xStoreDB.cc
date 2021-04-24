//
// Created by Chuang on 2019/6/10.
//
#include "storagemanager/xStore.h"
#include "DiskStorageManager.h"
#include <cmath>
#include <cstring>
#include <sys/stat.h>
#include <storagemanager/tjsql.h>
#define fp (int(m_pageSize/sizeof(prexp)/3)-1)

static bool CheckFilesExists(const string &s) {
    struct stat stats;

    std::ostringstream os;
    std::string name;
    int ret;

    os.str("");
    os << s << ".dbproperty"; //json property
    name = os.str();
    ret = stat(name.c_str(), &stats);
    if (ret != 0) return false;

    return true;
}

static int getLastId(string s)
{
    return 0;
}

void xStoreDB::loadFile(string file) {
    std::cerr << "start loading " << file << " into " << m_name << " DB\n";
    conn_init(m_name);
    ifstream inFile(file, ios::in);
    string lineStr;
    set<id_type> ids;
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
                    uint32_t plen = min(rem, fp);
                    xPoint *p;
                    data = new uint8_t[24 * plen];
                    ptr = data;
                    for (int i = 0; i < plen; i++) {
                        p = &(tj.m_points[cur++]);
                        p->storeToByteArrayE(&ptr, len);
                        if (i != plen - 1) ptr += len;
                    }
                    ptr = data;
                    int64_t newPage = db_last_pageid()+1;
                    db_insert_page(newPage, 24 * plen,data);

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
                ids.insert(id);
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

xStoreDB::xStoreDB(string myname, string file, bool bsubtraj, bool forceNew)
{
    m_needfree = false;
    conn_init(myname);
    cerr<<myname<<endl;
    m_name = myname;
    m_pageSize = 4096;
    string filename = filedirprefix+myname;
    if ((!forceNew) && CheckFilesExists(filename)) {
        ifstream propFile(filedirprefix+myname + ".dbproperty", ios::in);
        propFile >> m_property;
        propFile.close();
        if(tjstat->lineCount==0) { /* only load it when not inited*/
            tjstat->fromString(m_property["tjstat"]);
        }
        m_bSubTraj = m_property["bSubTraj"];
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
        m_property["trajfile"] = file;
        tjstat->init();
        ifstream inFile(file, ios::in);
        string lineStr;
        set<id_type> ids;
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
                        int64_t newPage = db_last_pageid()+1;
                        db_insert_page(newPage, 24 * plen,data);
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
                    ids.insert(id);
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

xStoreDB::xStoreDB() {
    m_needfree = false;
}

xStore * xStoreDB::clone()
{
    xStoreDB* res = new xStoreDB();
    assert(!m_isro);
    conn_init(m_name);
    res->m_name=m_name;
    res->m_pageSize=m_pageSize;
    ifstream propFile(filedirprefix+res->m_name + ".dbproperty", ios::in);
    propFile >> res->m_property;
    propFile.close();
    m_bSubTraj = m_property["bSubTraj"];
    m_isro = true;
    return res;
}

void xStoreDB::loadTraj(xTrajectory &out, const xStoreEntry &e) {
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
            uint32_t tmp;
            db_load_page(i,tmp,&data);
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
            uint32_t tmp;
            db_load_page(i,tmp,&data);
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
        }
        out.m_fakehead = (ms != 0);
        out.m_fakeback = (me == te.m_npoint - 1);
    }
}

xTrajectory xStoreDB::randomSubtraj(double len) {
    id_type lastid = db_last_trajid();
    int rnd = random(1, lastid);
    int i = 0;
    xTrajectory tj,tj2;
    loadTraj(tj, xStoreEntry(rnd, 0, 1000000));
    double t = tj.randomPoint().m_t;
    tj.getPartialxTrajectory(t - len / 2, t + len / 2, tj2);
    return tj2;
}

xPoint xStoreDB::randomPoint() {
    int rnd = random(1, db_last_trajid());
    int rnd2;
    int i = 0;
    xTrajectory tj;
    loadTraj(tj, xStoreEntry(rnd, 0, 1000000));
    rnd2 = random(0, tj.m_points.size() - 1);
    return tj.m_points[rnd2];
    throw Tools::IllegalStateException("...");
}

void xStoreDB::flush() {
    if(m_isro) return;
    ofstream propFile(filedirprefix+m_name + ".dbproperty", ios::out);
    propFile << m_property;
    propFile.close();
}

void xStoreDB::loadByteArray(const id_type page, uint32_t& len, uint8_t** data)
{
    db_load_page(page,len,data);
}
void xStoreDB::storeByteArray(id_type& page, const uint32_t len, const uint8_t* const data){
    db_insert_page(page,len,data);
}
void xStoreDB::deleteByteArray(const id_type page){
    return;
}



xSBBStream::xSBBStream(xStore *p, CUTFUNC f)
        : m_cutFunc(f),m_pstore(p) {
    if(m_pstore->m_trajIdx!=NULL) {
        m_isdb = false;
        m_it = m_pstore->m_trajIdx->begin();
    }else{
        m_isdb = true;
        m_size = db_last_trajid();
    }
}

bool xSBBStream::hasNext() {
    return !m_buf.empty() ||
           (m_isdb && m_id < m_size)
           || (!m_isdb && m_it != m_pstore->m_trajIdx->end());
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
    if(!m_isdb)
        m_it = m_pstore->m_trajIdx->begin();
    else
        m_id = 0;
    m_count = 0;
}
