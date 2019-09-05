//
// Created by Chuang on 2019/6/10.
//

#include "storagemanager/TrajStore.h"
#include <cmath>
#include <cstring>
//
// tsExternalSorter::Record
//
XZ3Enocder* XZ3Enocder::singleton=nullptr;
XZ3Enocder* XZ3Enocder::instance() {
    if (singleton == nullptr)
    {
        singleton = new XZ3Enocder();
    }
    return singleton;
}
trajStat* trajStat::singleton=nullptr;
trajStat* trajStat::instance() {
    if (singleton == nullptr)
    {
        singleton = new trajStat();
    }
    return singleton;
}
SIDX_DLL std::ostream& SpatialIndex::operator<<(std::ostream& os, const trajStat& r) {
    os<<"Trajectory Statistics:\n";
    os<<"\ttotal time: "<<r.M<<"\n";
    os<<"\tline count: "<<r.lineCount<<"\n";
    os<<"\ttraj count: "<<r.trajCount<<"\n";
    os<<"\tline length: "<<r.tl<<"\n";
    os<<"\ttraj length: "<<r.jt<<"\n";
    os<<"\taverage speed: "<<r.v<<"\n";
    os<<"\texpanding:"<<r.Dx<<"\t"<<r.Dy<<"\t"<<r.Dt<<"\n";
    return os;
}
XZ3Enocder::XZ3Enocder()
    :m_xmin(0),m_xmax(1),m_ymin(0),m_ymax(1),m_zmin(0),m_zmax(1),m_length(8)
    {}

long XZ3Enocder::encode(double x, double y, double z) {
    double xmin=m_xmin,xmax=m_xmax,ymin=m_ymin,ymax=m_ymax,zmin=m_zmin,zmax=m_zmax;
    long cs = 0L;
    int i = 0;
    while (i < m_length) {
        double xCenter = (xmin + xmax) / 2.0;
        double yCenter = (ymin + ymax) / 2.0;
        double zCenter = (zmin + zmax) / 2.0;
        int cmpRes=0;
        if(x<xCenter) cmpRes+=100;
        if(y<yCenter) cmpRes+=10;
        if(z<zCenter) cmpRes+=1;
        int g=10;
        switch (cmpRes){
            case 111:
                cs += 1L;
                xmax = xCenter; ymax = yCenter; zmax = zCenter;
                break;
            case 011:
                cs += 1L + 1L * (long(std::pow(8, g - i)) - 1L) / 7L;
                xmin = xCenter; ymax = yCenter; zmax = zCenter;
                break;
            case 101:
                cs += 1L + 2L * (long(std::pow(8, g - i)) - 1L) / 7L;
                xmax = xCenter; ymin = yCenter; zmax = zCenter;
                break;
            case 001:
                cs += 1L + 3L * (long(std::pow(8, g - i)) - 1L) / 7L;
                xmin = xCenter; ymin = yCenter; zmax = zCenter;
                break;
            case 110:
                cs += 1L + 4L * (long(std::pow(8, g - i)) - 1L) / 7L;
                xmax = xCenter; ymax = yCenter; zmin = zCenter;
                break;
            case 010:
                cs += 1L + 5L * (long(std::pow(8, g - i)) - 1L) / 7L;
                xmin = xCenter; ymax = yCenter; zmin = zCenter;
                break;
            case 100:
                cs += 1L + 6L * (long(std::pow(8, g - i)) - 1L) / 7L;
                xmax = xCenter; ymin = yCenter; zmin = zCenter;
                break;
            case 000:
                cs += 1L + 7L * (long(std::pow(8, g - i)) - 1L) / 7L;
                xmin = xCenter; ymin = yCenter; zmin = zCenter;
                break;
        }
        i += 1;
    }
    return cs;
}
tsExternalSorter::Record::Record()
        : m_pData(0)
{
}

tsExternalSorter::Record::Record(const Region& r, id_type id,id_type pvId,id_type ntId, uint32_t len, uint8_t* pData, uint32_t s,uint32_t level,MBC *mbc)
        :  m_id(id),m_pvId(pvId),m_ntId(ntId), m_len(len), m_pData(pData), m_s(s),m_level(level)
{
    m_r=r;
    if(mbc!= nullptr){
        m_mbc=*mbc;
    }
}


tsExternalSorter::Record::~Record()
{
    delete[] m_pData;
}

bool tsExternalSorter::Record::operator<(const Record& r) const
{
    Point c1,c2;
    m_mbc.getCenter(c1);
    r.m_mbc.getCenter(c2);
    auto encoder=XZ3Enocder::instance();
    return encoder->encode(c1.m_pCoords[0],c1.m_pCoords[1],c1.m_pCoords[2])*4+m_mbc.getOrient()
        <encoder->encode(c2.m_pCoords[0],c2.m_pCoords[1],c2.m_pCoords[2])*4+m_mbc.getOrient();
}

void tsExternalSorter::Record::storeToFile(Tools::TemporaryFile& f)
{
    f.write(static_cast<uint64_t>(m_id));
    f.write(static_cast<uint64_t>(m_pvId));
    f.write(static_cast<uint64_t>(m_ntId));
    f.write(m_r.m_dimension);
    f.write(m_s);
    f.write(m_level);

    for (uint32_t i = 0; i < m_r.m_dimension; ++i)
    {
        f.write(m_r.m_pLow[i]);
        f.write(m_r.m_pHigh[i]);
    }
    if(m_level==0) {
        for (uint32_t i = 0; i < m_mbc.m_dimension; ++i) {
            f.write(m_mbc.m_pLow[i]);
            f.write(m_mbc.m_pHigh[i]);
        }
        f.write(m_mbc.m_startTime);
        f.write(m_mbc.m_endTime);
        f.write(m_mbc.m_rd);
        f.write(m_mbc.m_rv);
    }
    f.write(m_len);
    if (m_len > 0) f.write(m_len, m_pData);
}

void tsExternalSorter::Record::loadFromFile(Tools::TemporaryFile& f)
{
    m_id = static_cast<id_type>(f.readUInt64());
    m_pvId = static_cast<id_type>(f.readUInt64());
    m_ntId = static_cast<id_type>(f.readUInt64());
    uint32_t dim = f.readUInt32();
    m_s = f.readUInt32();
    m_level=f.readUInt32();

    if (dim != m_r.m_dimension)
    {
        delete[] m_r.m_pLow;
        delete[] m_r.m_pHigh;
        m_r.m_dimension = dim;
        m_r.m_pLow = new double[dim];
        m_r.m_pHigh = new double[dim];
    }

    for (uint32_t i = 0; i < m_r.m_dimension; ++i)
    {
        m_r.m_pLow[i] = f.readDouble();
        m_r.m_pHigh[i] = f.readDouble();
    }
    if(m_level==0) {
        m_mbc.makeInfinite(dim - 1);
        for (uint32_t i = 0; i < m_mbc.m_dimension; ++i) {
            m_mbc.m_pLow[i] = f.readDouble();
            m_mbc.m_pHigh[i] = f.readDouble();
        }
        m_mbc.m_startTime = f.readDouble();
        m_mbc.m_endTime = f.readDouble();
        m_mbc.m_rd = f.readDouble();
        m_mbc.m_rv = f.readDouble();
    }
    m_len = f.readUInt32();
    delete[] m_pData; m_pData = 0;
    if (m_len > 0) f.readBytes(m_len, &m_pData);
}

//
// tsExternalSorter
//
tsExternalSorter::tsExternalSorter(uint32_t u32PageSize, uint32_t u32BufferPages)
        : m_bInsertionPhase(true), m_u32PageSize(u32PageSize),
          m_u32BufferPages(u32BufferPages), m_u64TotalEntries(0), m_stI(0)
{
}

tsExternalSorter::~tsExternalSorter()
{
    for (m_stI = 0; m_stI < m_buffer.size(); ++m_stI) delete m_buffer[m_stI];
}

void tsExternalSorter::insert(Record* r)
{
    if (m_bInsertionPhase == false)
        throw Tools::IllegalStateException("tsExternalSorter::insert: Input has already been sorted.");

    m_buffer.emplace_back(r);
    ++m_u64TotalEntries;

    // this will create the initial, sorted buckets before the
    // external merge sort.
    if (m_buffer.size() >= m_u32PageSize * m_u32BufferPages)
    {
        std::sort(m_buffer.begin(), m_buffer.end(), Record::SortAscending());
        Tools::TemporaryFile* tf = new Tools::TemporaryFile();
        for (size_t j = 0; j < m_buffer.size(); ++j)
        {
            m_buffer[j]->storeToFile(*tf);
            delete m_buffer[j];
        }
        m_buffer.clear();
        tf->rewindForReading();
        m_runs.emplace_back(Tools::SmartPointer<Tools::TemporaryFile>(tf));
    }
}

void tsExternalSorter::sort()
{
    if (m_bInsertionPhase == false)
        throw Tools::IllegalStateException("tsExternalSorter::sort: Input has already been sorted.");

    if (m_runs.empty())
    {
        // The data fits in main memory. No need to store to disk.
        std::sort(m_buffer.begin(), m_buffer.end(), Record::SortAscending());
        m_bInsertionPhase = false;
        return;
    }

    if (m_buffer.size() > 0)
    {
        // Whatever remained in the buffer (if not filled) needs to be stored
        // as the final bucket.
        std::sort(m_buffer.begin(), m_buffer.end(), Record::SortAscending());
        Tools::TemporaryFile* tf = new Tools::TemporaryFile();
        for (size_t j = 0; j < m_buffer.size(); ++j)
        {
            m_buffer[j]->storeToFile(*tf);
            delete m_buffer[j];
        }
        m_buffer.clear();
        tf->rewindForReading();
        m_runs.emplace_back(Tools::SmartPointer<Tools::TemporaryFile>(tf));
    }

    if (m_runs.size() == 1)
    {
        m_sortedFile = m_runs.front();
    }
    else
    {
        Record* r = 0;

        while (m_runs.size() > 1)
        {
            Tools::SmartPointer<Tools::TemporaryFile> tf(new Tools::TemporaryFile());
            std::vector<Tools::SmartPointer<Tools::TemporaryFile> > buckets;
            std::vector<std::queue<Record*> > buffers;
            std::priority_queue<PQEntry, std::vector<PQEntry>, PQEntry::SortAscending> pq;

            // initialize buffers and priority queue.
            std::list<Tools::SmartPointer<Tools::TemporaryFile> >::iterator it = m_runs.begin();
            for (uint32_t i = 0; i < (std::min)(static_cast<uint32_t>(m_runs.size()), m_u32BufferPages); ++i)
            {
                buckets.emplace_back(*it);
                buffers.emplace_back(std::queue<Record*>());

                r = new Record();
                r->loadFromFile(**it);
                // a run cannot be empty initially, so this should never fail.
                pq.push(PQEntry(r, i));

                for (uint32_t j = 0; j < m_u32PageSize - 1; ++j)
                {
                    // fill the buffer with the rest of the page of records.
                    try
                    {
                        r = new Record();
                        r->loadFromFile(**it);
                        buffers.back().push(r);
                    }
                    catch (Tools::EndOfStreamException)
                    {
                        delete r;
                        break;
                    }
                }
                ++it;
            }

            // exhaust buckets, buffers, and priority queue.
            while (! pq.empty())
            {
                PQEntry e = pq.top(); pq.pop();
                e.m_r->storeToFile(*tf);
                delete e.m_r;

                if (! buckets[e.m_u32Index]->eof() && buffers[e.m_u32Index].empty())
                {
                    for (uint32_t j = 0; j < m_u32PageSize; ++j)
                    {
                        try
                        {
                            r = new Record();
                            r->loadFromFile(*buckets[e.m_u32Index]);
                            buffers[e.m_u32Index].push(r);
                        }
                        catch (Tools::EndOfStreamException)
                        {
                            delete r;
                            break;
                        }
                    }
                }

                if (! buffers[e.m_u32Index].empty())
                {
                    e.m_r = buffers[e.m_u32Index].front();
                    buffers[e.m_u32Index].pop();
                    pq.push(e);
                }
            }

            tf->rewindForReading();

            // check if another pass is needed.
            uint32_t u32Count = std::min(static_cast<uint32_t>(m_runs.size()), m_u32BufferPages);
            for (uint32_t i = 0; i < u32Count; ++i)
            {
                m_runs.pop_front();
            }

            if (m_runs.size() == 0)
            {
                m_sortedFile = tf;
                break;
            }
            else
            {
                m_runs.emplace_back(tf);
            }
        }
    }

    m_bInsertionPhase = false;
}

tsExternalSorter::Record* tsExternalSorter::getNextRecord()
{
    if (m_bInsertionPhase == true)
        throw Tools::IllegalStateException("tsExternalSorter::getNextRecord: Input has not been sorted yet.");

    Record* ret;

    if (m_sortedFile.get() == 0)
    {
        if (m_stI < m_buffer.size())
        {
            ret = m_buffer[m_stI];
            m_buffer[m_stI] = 0;
            ++m_stI;
        }
        else
            throw Tools::EndOfStreamException("");
    }
    else
    {
        ret = new Record();
        ret->loadFromFile(*m_sortedFile);
    }

    return ret;
}

inline uint64_t tsExternalSorter::getTotalEntries() const
{
    return m_u64TotalEntries;
}


TrajStore::Entry::Entry(id_type page, uint32_t start, uint32_t len, id_type pvId, id_type ntId)
    :m_page(page),m_start(start),m_len(len),m_pvId(pvId),m_ntId(ntId){}


TrajStore::~TrajStore() {
    releaseTmp();
}

TrajStore::TrajStore(IStorageManager *store,uint32_t pageSize,int maxseg)
    :m_pStorageManager(store),m_pageSize(pageSize),m_maxTrajSegs(maxseg){}

struct id3{
    id_type m_id;
    id_type m_pvId;
    id_type m_ntId;
    id3(id_type id,id_type pvId,id_type ntId)
        :m_id(id),m_pvId(pvId),m_ntId(ntId){};
};
void TrajStore::loadSegments(vector<std::pair<id_type, vector<Trajectory>> > &trajs,bool idFisrt){
    if(idFisrt) {
        id_type thisPageId = m_pStorageManager->nextPage();
        int linecount=0;
        double vol1=0,vol2=0;
        for (auto &traj:trajs) {
            auto it=traj.second.begin();
            id_type segid;
            id_type pvId=-1,ntId=-1;
            while(it!=traj.second.end()){
                Trajectory seg=*it;
                vector<id3> ids;
                segid=getSegId(traj.first,it-traj.second.begin());
                pvId=it==traj.second.begin()?-1:segid-1;
                ntId=(it+1)==traj.second.end()?-1:(segid+1);
                ids.emplace_back(id3(segid,pvId,ntId));
                Region thebr;
                MBC thebc;
                it->getMBC(thebc);
                it->getMBRfull(thebr);
                m_entryMbcs[segid] = thebc;
                m_entryMbrs[segid] = thebr;
                vol1 += thebr.getArea();
                vol2 += thebc.getArea();
                double v = (Point(thebc.m_pLow, 2).getMinimumDistance(Point(thebc.m_pHigh, 2))) /
                           (thebc.m_endTime - thebc.m_startTime) + thebc.m_rv;
                if (v > m_maxVelocity) m_maxVelocity = v;
                while((it+1)!=traj.second.end()&&seg.m_points.size()-(seg.m_points.size()/170)*170+(it+1)->m_points.size()<170){
                    it++;
                    seg.linkTrajectory(*it);
                    segid=getSegId(traj.first,it-traj.second.begin());
                    pvId=it==traj.second.begin()?-1:segid-1;
                    ntId=(it+1)==traj.second.end()?-1:(segid+1);
                    ids.emplace_back(id3(segid,pvId,ntId));
                    Region thebr;
                    MBC thebc;
                    it->getMBC(thebc);
                    it->getMBRfull(thebr);
                    m_entryMbcs[segid] = thebc;
                    m_entryMbrs[segid] = thebr;
                    vol1 += thebr.getArea();
                    vol2 += thebc.getArea();
                    double v = (Point(thebc.m_pLow, 2).getMinimumDistance(Point(thebc.m_pHigh, 2))) /
                               (thebc.m_endTime - thebc.m_startTime) + thebc.m_rv;
                    if (v > m_maxVelocity) m_maxVelocity = v;
                }
                m_segCount++;
                m_timeCount += seg.m_endTime() - seg.m_startTime();
                linecount += seg.m_points.size() - 1;
                uint8_t *data;
                uint32_t len;
                seg.storeToByteArray(&data, len);
                for(const auto &s:ids){
                    Entry *sege = new Entry(thisPageId, 0, len, s.m_pvId,s.m_ntId);
                    m_entries[s.m_id] = sege;
                }
                id_type newPage = StorageManager::NewPage;
                m_pStorageManager->storeByteArray(newPage,len, data);
                assert(newPage == thisPageId);
                pvId=segid;
                thisPageId+=std::ceil(double(len)/m_pageSize);
                it++;
            }
        }
        std::cerr << "total segment is" << m_entryMbrs.size() << "=" << m_entryMbcs.size() << "\n";
        std::cerr << "expression effectiveness " << vol1 << " versus " << vol2 << "\n" << vol2 / vol1 << "\n";
    }else{//group by spatial
        Tools::SmartPointer<tsExternalSorter> es = Tools::SmartPointer<tsExternalSorter>(
                new tsExternalSorter(m_pageSize, 200));
        //load segments to sorter's record
        std::cerr << "TrajStore:loading segments\n";
        Region maxbr;
        maxbr.makeInfinite(3);
        double vol1 = 0, vol2 = 0;
        long linecount = 0;
        for (auto &traj:trajs) {
            for (int j = 0; j < traj.second.size(); j++) {//segment
                Trajectory seg = traj.second[j];
                m_segCount++;
                m_timeCount += seg.m_endTime() - seg.m_startTime();
                linecount += seg.m_points.size() - 1;
                uint8_t *data;
                uint32_t len;
                seg.storeToByteArray(&data, len);
                MBC thebc;
                seg.getMBC(thebc);
                double v = (Point(thebc.m_pLow, 2).getMinimumDistance(Point(thebc.m_pHigh, 2))) /
                           (thebc.m_endTime - thebc.m_startTime) + thebc.m_rv;
                if (v > m_maxVelocity) m_maxVelocity = v;
                Region thebr;
                seg.getMBRfull(thebr);
                maxbr.combineRegion(thebr);
                id_type segid = getSegId(traj.first, j),
                        pvId = (j == 0) ? -1 : segid - 1,
                        ntId = (j == traj.second.size() - 1) ? -1 : segid + 1;
                es->insert(new tsExternalSorter::Record(thebr, segid, pvId, ntId, len, data, 0, 0, &thebc));
                m_entryMbcs[segid] = thebc;
                m_entryMbrs[segid] = thebr;
                vol1 += thebr.getArea();
                vol2 += thebc.getArea();
            }
        }
        std::cerr << "total segment is" << m_entryMbrs.size() << "=" << m_entryMbcs.size() << "\n";
        std::cerr << "expression effectiveness " << vol1 << " versus " << vol2 << "\n" << vol2 / vol1 << "\n";
        //adjust the encoder
        auto encoder = XZ3Enocder::instance();
        encoder->m_xmin = maxbr.m_pLow[0];
        encoder->m_xmax = maxbr.m_pHigh[0];
        encoder->m_ymin = maxbr.m_pLow[1];
        encoder->m_ymax = maxbr.m_pHigh[1];
        encoder->m_zmin = maxbr.m_pLow[2];
        encoder->m_zmax = maxbr.m_pHigh[2];

        //sort the data
        std::cerr << "TrajStore:sorting using XZ3 curve\n";
        es->sort();

        //save the data into pages
        id_type thisPageId = m_pStorageManager->nextPage();
#ifndef NDEBUG
        std::cerr << "TrajStore:storing data into pages\n";
#endif
        uint32_t spaceRem = m_pageSize, currentLen = 0;//this two add to m_pagesize
        uint8_t *pageData = new uint8_t[m_pageSize];
        while (true) {
            tsExternalSorter::Record *r;
            try { r = es->getNextRecord(); } catch (Tools::EndOfStreamException) {
                if (currentLen > 0) {
//                m_pStorageManager->storeByteArray(thisPageId, currentLen, pageData);
                    id_type newPage = StorageManager::NewPage;
                    m_pStorageManager->storeByteArray(newPage, currentLen, pageData);
                    assert(newPage == thisPageId);
//                if(newPage!=thisPageId) std::cerr<<newPage<<" "<<thisPageId<<"\n";
                    thisPageId++;
                    spaceRem = m_pageSize;
                    currentLen = 0;
                }
                break;
            }
            //if the traj expand to multiple pages, open a new page and store it
            uint8_t *data = r->m_pData;
            uint32_t len = r->m_len;
            if (len > m_pageSize) {
                //flush the current page
                if (currentLen > 0) {
                    id_type newPage = StorageManager::NewPage;
                    m_pStorageManager->storeByteArray(newPage, currentLen, pageData);
                    assert(newPage == thisPageId);
//                if(newPage!=thisPageId) std::cerr<<newPage<<" "<<thisPageId<<"\n";
                    thisPageId++;
                    spaceRem = m_pageSize;
                    currentLen = 0;
                }
                //record the entry
                Entry *sege = new Entry(thisPageId, 0, len, r->m_pvId, r->m_ntId);
                m_entries[r->m_id] = sege;
                //write the big traj
                id_type newPage = StorageManager::NewPage;
                m_pStorageManager->storeByteArray(newPage, len, data);
                assert(newPage == thisPageId);
//            if(newPage!=thisPageId) std::cerr<<newPage<<" "<<thisPageId<<"\n";
                thisPageId += std::ceil(1.0 * len / m_pageSize);
                spaceRem = m_pageSize;
                currentLen = 0;
            }
                //if this page is full,open a new page and store it
            else if (spaceRem < len) {
                //flush the current page
                if (currentLen > 0) {
                    id_type newPage = StorageManager::NewPage;
                    m_pStorageManager->storeByteArray(newPage, currentLen, pageData);
//                assert(newPage==thisPageId);
                    if (newPage != thisPageId) std::cerr << newPage << " " << thisPageId << "\n";
                    thisPageId++;
                    spaceRem = m_pageSize;
                    currentLen = 0;
                }
                //record the entry
                Entry *sege = new Entry(thisPageId, 0, len, r->m_pvId, r->m_ntId);
                m_entries[r->m_id] = sege;
                //write in the new page
                memcpy(&pageData[currentLen], data, len);
                currentLen += len;
                spaceRem -= len;
            }
                //if the Data could fit in page
            else {
                //record the entry
                Entry *sege = new Entry(thisPageId, currentLen, len, r->m_pvId, r->m_ntId);
                m_entries[r->m_id] = sege;
                //write in the new page
                memcpy(&pageData[currentLen], data, len);
                currentLen += len;
                spaceRem -= len;
            }
            delete r;
        }
        delete[] pageData;
    }
}

const Trajectory TrajStore::getTraj(id_type &id) {
    auto it=m_entries.find(id);
    assert(it!=m_entries.end());
    Entry e=*(it->second);
    uint32_t len=e.m_start+e.m_len;
    uint8_t *load = new uint8_t[len];
    m_pStorageManager->loadByteArray(e.m_page,len,&load);
    m_trajIO+=std::ceil(len/4096.0);
//    std::cerr<<"get Traj "<<id<<" with length"<<len<<"\n";
    uint8_t *data = load+e.m_start;
    Trajectory traj;
    traj.loadFromByteArray(data);
//    std::cerr<<traj;
    delete[](load);
    return traj;
}

const Trajectory TrajStore::getTrajByTime(id_type &id, double tstart, double tend) {
    auto it=m_entries.find(id);
    if(it==m_entries.end()){
        std::cerr<<id<<"\n";
    }
    assert(it!=m_entries.end());
    Entry e=*(it->second);
//    uint32_t len=e.m_start+e.m_len;
//    uint8_t *load = new uint8_t[len];
//    m_pStorageManager->loadByteArray(e.m_page,len,&load);
//    uint8_t *data = load+e.m_start;
    Trajectory traj,tmptraj;
    traj=getTraj(id);
//    traj.loadFromByteArray(data);
    Entry tmpe=e;
    while(tmpe.m_pvId>=0 && traj.m_points.front().m_time>=tstart){
        tmptraj=getTraj(tmpe.m_pvId);
        it=m_entries.find(tmpe.m_pvId);
        tmpe=*(it->second);
        traj.linkTrajectory(tmptraj);
    }
    tmpe=e;
    while(tmpe.m_ntId>=0 && traj.m_points.back().m_time<=tend){
        tmptraj=getTraj(tmpe.m_ntId);
        it=m_entries.find(tmpe.m_ntId);
        tmpe=*(it->second);
        traj.linkTrajectory(tmptraj);
    }
    return traj;
}

const ShapeList TrajStore::getMBCsByTime(id_type &id, double tstart, double tend) {
    auto it=m_entries.find(id);
    auto bc=&((*m_entryMbcs.find(id)).second);
    assert(it!=m_entries.end());
    Entry e=*(it->second);
    ShapeList bcs;
    Entry tmpe=e;
    std::list<MBC*> mbclist;
    mbclist.push_front(bc);
    while(tmpe.m_pvId>=0&&bc->m_startTime>tstart){
        id_type newid=tmpe.m_pvId;
        it=m_entries.find(newid);
        bc=&((*m_entryMbcs.find(newid)).second);
        mbclist.push_front(bc);
        tmpe=*(it->second);
        m_boundingVisited++;
    }
    tmpe=e;
    while(tmpe.m_ntId>=0&&bc->m_endTime<tend){
        id_type newid=tmpe.m_ntId;
        it=m_entries.find(newid);
        bc=&((*m_entryMbcs.find(newid)).second);
        mbclist.emplace_back(bc);
        tmpe=*(it->second);
        m_boundingVisited++;
    }
    for ( const auto &mbc:mbclist) {
        bcs.insertbc(mbc);
    }
    return bcs;
}

const ShapeList TrajStore::getMBRsByTime(id_type &id, double tstart, double tend) {
    auto it=m_entries.find(id);
    auto br=&((*m_entryMbrs.find(id)).second);
    assert(it!=m_entries.end());
    Entry e=*(it->second);
    ShapeList brs;
    Entry tmpe=e;
    std::list<Region*> mbrlist;
    mbrlist.push_front(br);
    while(tmpe.m_pvId>=0&&br->m_pLow[br->m_dimension-1]>tstart){
        id_type newid=tmpe.m_pvId;
        it=m_entries.find(newid);
        br=&((*m_entryMbrs.find(newid)).second);
        mbrlist.push_front(br);
        tmpe=*(it->second);
        m_boundingVisited++;
    }
    tmpe=e;
    while(tmpe.m_ntId>=0&&br->m_pHigh[br->m_dimension-1]<tend){
        id_type newid=tmpe.m_ntId;
        it=m_entries.find(newid);
        br=&((*m_entryMbrs.find(newid)).second);
        mbrlist.emplace_back(br);
        tmpe=*(it->second);
        m_boundingVisited++;
    }
    for ( const auto &mbr:mbrlist) {
        brs.insertbr(mbr);
    }
    return brs;
}

const Region TrajStore::getMBR(id_type &id) {
    return (*m_entryMbrs.find(id)).second;
}

const MBC TrajStore::getMBC(id_type &id) {
    return (*m_entryMbcs.find(id)).second;
}