//
// Created by Chuang on 2019/6/10.
//
#include "storagemanager/xStore.h"
#include <cmath>
#include <cstring>
#include <sys/stat.h>
auto tjstat = trajStat::instance();

bool CheckFilesExists(const string &s)
{
    struct stat stats;

    std::ostringstream os;
    std::string name;
    int ret;
    os << s <<".idx"; //diskfileidx
    name = os.str();
    ret = stat(name.c_str(), &stats);
    if (ret != 0) return false;

    os.str("");
    os << s <<".dat";//diskfiledat
    name = os.str();
    ret = stat(name.c_str(), &stats);
    if (ret != 0) return false;

    os.str("");
    os << s <<".trajidx"; //trajidx(m_trajidx)
    name = os.str();
    ret = stat(name.c_str(), &stats);
    if (ret != 0) return false;

    os.str("");
    os << s <<".property"; //json property
    name = os.str();
    ret = stat(name.c_str(), &stats);
    if (ret != 0) return false;

    return true;
}

xTrajEntry::xTrajEntry(id_type page, uint32_t len)
    :m_page(page),m_npoint(len){}

std::string xTrajEntry::xTrajEntry::toString() {
    string s;
    s=std::to_string(m_page)+std::to_string(m_npoint);
    return s;
}

xTrajEntry::xTrajEntry::xTrajEntry(string &s) {
    auto ws=split(s,',');
    m_page=std::stoll(ws[0]);
    m_npoint=std::stoi(ws[1]);
}

xStore::~xStore() {
    flush();
    for(auto s:m_trajIdx){
        delete s.second;
    }
}

xStore::xStore(string myname,string file,bool forceNew) {
    m_name=myname;
    m_pageSize = 4096;
#define fp (int(m_pageSize/sizeof(prexp)/3))
    if((!forceNew)&&CheckFilesExists(myname)){
        std::cerr<<"using existing xStore "<<myname<<endl;
        m_pStorageManager = loadDiskStorageManager(myname);
        ifstream propFile(myname+".property", ios::in);
        propFile>>m_property;
        propFile.close();
        tjstat->fromString(m_property["stat"]);
        ifstream trajidxFile(m_name+".trajidx", ios::in);
        id_type size,id,page,off;
        trajidxFile>>size;
        for(int i=0;i<size;i++){
            trajidxFile>>id>>page>>off;
            m_trajIdx[id] = new xTrajEntry(page,off);
        }
    }
    else{
        std::cerr<<"start loading "<< file<<" into "<<myname<<"\n";
        m_pStorageManager = createNewDiskStorageManager(myname,m_pageSize);
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
                if (tj.m_points.size() >= 2) {
                    int rem = tj.m_points.size();
                    int cur=0;
                    uint8_t* ptr;
                    bool isfirst = true;
                    id_type firstpage;
                    while (rem>0){
                        int plen = min(rem,fp);
                        xPoint *p;
                        data = new uint8_t[24*plen];
                        ptr = data;
                        for(int i=0;i<plen;i++){
                            p = &(tj.m_points[cur++]);
                            p->storeToByteArrayE(&ptr,len);
                            if(i!= plen-1) ptr += len;
                        }
                        ptr=data;
                        id_type newPage = StorageManager::NewPage;
                        m_pStorageManager->storeByteArray(newPage,24*plen, data);
                        if(isfirst){
                            firstpage = newPage;
                            isfirst=false;
                        }
                        rem-=fp;
                        delete[] data;
                    }
                    m_trajIdx[id] = new xTrajEntry(firstpage, tj.m_points.size());
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
                    if(tjstat->vmax< tj.maxSpeed()){
                        tjstat->vmax=tj.maxSpeed();
                    }
                    tjstat->M += tj.m_endTime() - tj.m_startTime();
                }
            }
            catch (...) {
                break;
            }
        }
        std::cerr<<"load finished\n";
        tjstat->Dx = tjstat->maxx - tjstat->minx;
        tjstat->Dy = tjstat->maxy - tjstat->miny;
        tjstat->Dt = tjstat->maxt - tjstat->mint;
        tjstat->tl = tjstat->M / tjstat->lineCount;
        tjstat->jt = tjstat->M / tjstat->trajCount;
        tjstat->v = tjstat->dist / tjstat->M;
        tjstat->Sr = (tjstat->Dx+tjstat->Dy)/2;
        tjstat->P = tjstat->Dt;
        std::cerr<<file<<endl;
        if (file.find("td")!=file.npos) tjstat->usedata("td");
        if (file.find("gl")!=file.npos) tjstat->usedata("gl");
        std::cerr<<tjstat->toString();
        m_property["stat"]= tjstat->toString();
    }
}


void xStore::loadTraj(xTrajectory &out, const xStoreEntry &e) {
    auto te= m_trajIdx[e.m_id];
    uint32_t ms = min(te->m_npoint-1,e.m_s), me=min(te->m_npoint-1,e.m_e);
    id_type pages = te->m_page + (ms+1)/fp;
    id_type pagee = te->m_page + int(ceil(1.0*(me+1)/fp))-1;
    id_type cur = ms/fp*fp;
    uint8_t *data, *ptr;
    uint32_t len;
    int ps,pe;
    prexp x,y,t;
    out.m_points.clear();
    for(auto i =pages;i<=pagee;i++) {
        ps = max(0, int(ms - cur));
        pe = min(fp-1,  int(me - cur));
        len = 3*sizeof(prexp)*pe;
        m_pStorageManager->loadByteArray(i,len,&data);
        for(ptr = data + 3*sizeof(prexp)*ps;ptr-data<len;){
            x=*((double*)ptr);
            ptr+=sizeof(prexp);
            y=*((double*)ptr);
            ptr+=sizeof(prexp);
            t=*((double*)ptr);
            ptr+=sizeof(prexp);
            out.m_points.emplace_back(xPoint(x, y, t));
        }
        cur+=fp;
        delete data;
    }
    out.m_fakehead = (ms!=0);
    out.m_fakeback = (me==te->m_npoint-1);
}

xPoint xStore::randomPoint() {
    int rnd = random(0,m_trajIdx.size()-1);
    int rnd2;
    int i=0;
    xTrajectory tj;
    for(auto key:m_trajIdx){
        if(i==rnd) {
            rnd2 = random(0,key.second->m_npoint-1);
            loadTraj(tj,xStoreEntry(key.first,0,rnd2));
            return tj.m_points.front();
        }
    }
    throw Tools::IllegalStateException("...");
}

void xStore::flush() {
    m_pStorageManager->flush();
    ofstream propFile(m_name+".property", ios::out);
    propFile<<m_property;
    propFile.close();
    ofstream trajidxFile(m_name+".trajidx", ios::out);
    trajidxFile<<m_trajIdx.size()<<endl;
    for(auto s=m_trajIdx.begin();s!=m_trajIdx.end();s++){
        trajidxFile<<s->first<<" "<<s->second->m_page<<" "<<s->second->m_npoint<<"\n";
    }
    trajidxFile.close();
}

xSBBStream::xSBBStream(xStore *p, CUTFUNC f)
:m_pstore(p), m_cutFunc(f) {
        m_it = m_pstore->m_trajIdx.begin();
}

bool xSBBStream::hasNext() {
    if(!m_buf.empty()||m_it!=m_pstore->m_trajIdx.end()) return true;
    return false;
}
xSBBData * xSBBStream::getNext() {
    bool isFirst=false;
    if(m_buf.empty()){
        xTrajectory tj;
        m_pstore->loadTraj(tj, xStoreEntry(m_it->first,0,m_it->second->m_npoint));
        m_buf = m_cutFunc(tj);
        m_id = m_it->first;
        m_it++;
        isFirst=true;
    }
    auto b = m_buf.front();
    m_buf.pop_front();
    return new xSBBData(m_count++,
                        xStoreEntry(m_id,b.first.first,b.first.second),b.second, !isFirst, !m_buf.empty());
}

uint32_t xSBBStream::size() {
    throw Tools::NotSupportedException("xsbbstream has no size");
}

void xSBBStream::rewind() {
    m_it = m_pstore->m_trajIdx.begin();
    m_count = 0;
}

