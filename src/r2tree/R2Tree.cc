//
// Created by chuang on 4/3/19.
//


#include <cstring>
#include <cmath>
#include <limits>

#include <spatialindex/SpatialIndex.h>
#include "Node.h"
#include "Leaf.h"
#include "Index.h"
#include "BulkLoader.h"
#include "R2Tree.h"

using namespace SpatialIndex::R2Tree;
using namespace SpatialIndex;

SpatialIndex::R2Tree::Data::Data(uint32_t len, uint8_t* pData, Mbbc& r, id_type id)
        : m_id(id), m_Mbbc(r), m_pData(0), m_dataLength(len)
{
    if (m_dataLength > 0)
    {
        m_pData = new uint8_t[m_dataLength];
        memcpy(m_pData, pData, m_dataLength);
    }
}


SpatialIndex::R2Tree::Data::~Data()
{
    delete[] m_pData;
}

SpatialIndex::R2Tree::Data* SpatialIndex::R2Tree::Data::clone()
{
    return new Data(m_dataLength, m_pData, m_Mbbc, m_id);
}

id_type SpatialIndex::R2Tree::Data::getIdentifier() const
{
    return m_id;
}

void SpatialIndex::R2Tree::Data::getShape(IShape** out) const
{
    *out = new Mbbc(m_Mbbc);
}

void SpatialIndex::R2Tree::Data::getData(uint32_t& len, uint8_t** data) const
{
    len = m_dataLength;
    *data = 0;

    if (m_dataLength > 0)
    {
        *data = new uint8_t[m_dataLength];
        memcpy(*data, m_pData, m_dataLength);
    }
}

uint32_t SpatialIndex::R2Tree::Data::getByteArraySize()
{
    return
            sizeof(id_type) +
            sizeof(uint32_t) +
            m_dataLength +
            m_Mbbc.getByteArraySize();
}

void SpatialIndex::R2Tree::Data::loadFromByteArray(const uint8_t* ptr)
{
    memcpy(&m_id, ptr, sizeof(id_type));
    ptr += sizeof(id_type);

    delete[] m_pData;
    m_pData = 0;

    memcpy(&m_dataLength, ptr, sizeof(uint32_t));
    ptr += sizeof(uint32_t);

    if (m_dataLength > 0)
    {
        m_pData = new uint8_t[m_dataLength];
        memcpy(m_pData, ptr, m_dataLength);
        ptr += m_dataLength;
    }

    m_Mbbc.loadFromByteArray(ptr);
}

void SpatialIndex::R2Tree::Data::storeToByteArray(uint8_t** data, uint32_t& len)
{
    // it is thread safe this way.
    uint32_t Mbbcsize;
    uint8_t* Mbbcdata = 0;
    m_Mbbc.storeToByteArray(&Mbbcdata, Mbbcsize);

    len = sizeof(id_type) + sizeof(uint32_t) + m_dataLength + Mbbcsize;

    *data = new uint8_t[len];
    uint8_t* ptr = *data;

    memcpy(ptr, &m_id, sizeof(id_type));
    ptr += sizeof(id_type);
    memcpy(ptr, &m_dataLength, sizeof(uint32_t));
    ptr += sizeof(uint32_t);

    if (m_dataLength > 0)
    {
        memcpy(ptr, m_pData, m_dataLength);
        ptr += m_dataLength;
    }

    memcpy(ptr, Mbbcdata, Mbbcsize);
    delete[] Mbbcdata;
    // ptr += Mbbcsize;
}





SpatialIndex::ISpatialIndex* SpatialIndex::R2Tree::returnR2Tree(SpatialIndex::IStorageManager& sm, Tools::PropertySet& ps)
{
    SpatialIndex::ISpatialIndex* si = new SpatialIndex::R2Tree::R2Tree(sm, ps);
    return si;
}

SpatialIndex::ISpatialIndex* SpatialIndex::R2Tree::createNewR2Tree(
        SpatialIndex::IStorageManager& sm,
        double fillFactor,
        uint32_t indexCapacity,
        uint32_t leafCapacity,
        uint32_t dimension,
        id_type& indexIdentifier)
{
    Tools::Variant var;
    Tools::PropertySet ps;

    var.m_varType = Tools::VT_DOUBLE;
    var.m_val.dblVal = fillFactor;
    ps.setProperty("FillFactor", var);

    var.m_varType = Tools::VT_ULONG;
    var.m_val.ulVal = indexCapacity;
    ps.setProperty("IndexCapacity", var);

    var.m_varType = Tools::VT_ULONG;
    var.m_val.ulVal = leafCapacity;
    ps.setProperty("LeafCapacity", var);

    var.m_varType = Tools::VT_ULONG;
    var.m_val.ulVal = dimension;
    ps.setProperty("Dimension", var);


    ISpatialIndex* ret = returnR2Tree(sm, ps);

    var.m_varType = Tools::VT_LONGLONG;
    var = ps.getProperty("IndexIdentifier");
    indexIdentifier = var.m_val.llVal;

    return ret;
}

SpatialIndex::ISpatialIndex* SpatialIndex::R2Tree::createAndBulkLoadNewR2Tree(
        BulkLoadMethod m,
        IDataStream& stream,
        SpatialIndex::IStorageManager& sm,
        double fillFactor,
        uint32_t indexCapacity,
        uint32_t leafCapacity,
        uint32_t dimension,
        id_type& indexIdentifier)
{
    SpatialIndex::ISpatialIndex* tree = createNewR2Tree(sm, fillFactor, indexCapacity, leafCapacity, dimension, indexIdentifier);

    uint32_t bindex = static_cast<uint32_t>(std::floor(static_cast<double>(indexCapacity * fillFactor)));
    uint32_t bleaf = static_cast<uint32_t>(std::floor(static_cast<double>(leafCapacity * fillFactor)));

    SpatialIndex::R2Tree::BulkLoader bl;

    stream.rewind();//rewind for reading
    switch (m)
    {
        case BLM_STR:
            bl.bulkLoadUsingSTR(static_cast<R2Tree*>(tree), stream, bindex, bleaf, 10000, 100);
            break;
        case BLM_STR2:
            bl.bulkLoadUsingSTR2(static_cast<R2Tree*>(tree), stream, bindex, bleaf, 10000, 100);
            break;
        case BLM_STR3:
            bl.bulkLoadUsingSTR3(static_cast<R2Tree*>(tree), stream, bindex, bleaf, 10000, 100);
            break;
        default:
            throw Tools::IllegalArgumentException("createAndBulkLoadNewR2Tree: Unknown bulk load method.");
            break;
    }

    return tree;
}



SpatialIndex::ISpatialIndex* SpatialIndex::R2Tree::loadR2Tree(IStorageManager& sm, id_type indexIdentifier)
{
    Tools::Variant var;
    Tools::PropertySet ps;

    var.m_varType = Tools::VT_LONGLONG;
    var.m_val.llVal = indexIdentifier;
    ps.setProperty("IndexIdentifier", var);

    return returnR2Tree(sm, ps);
}

SpatialIndex::R2Tree::R2Tree::R2Tree(IStorageManager& sm, Tools::PropertySet& ps) :
        m_pStorageManager(&sm),
        m_rootID(StorageManager::NewPage),
        m_headerID(StorageManager::NewPage),
        m_fillFactor(0.7),
        m_indexCapacity(100),
        m_leafCapacity(100),
        m_dimension(2),
        m_bTightMBRs(true),
        m_pointPool(500),
        m_MbbcPool(1000),
        m_indexPool(100),
        m_leafPool(100)
{
#ifdef HAVE_PTHREAD_H
    pthread_mutex_init(&m_lock, NULL);
#endif

    Tools::Variant var = ps.getProperty("IndexIdentifier");
    if (var.m_varType != Tools::VT_EMPTY)
    {
        if (var.m_varType == Tools::VT_LONGLONG) m_headerID = var.m_val.llVal;
        else if (var.m_varType == Tools::VT_LONG) m_headerID = var.m_val.lVal;
            // for backward compatibility only.
        else throw Tools::IllegalArgumentException("R2Tree: Property IndexIdentifier must be Tools::VT_LONGLONG");

        initOld(ps);
    }
    else
    {
        initNew(ps);
        var.m_varType = Tools::VT_LONGLONG;
        var.m_val.llVal = m_headerID;
        ps.setProperty("IndexIdentifier", var);
    }
}

SpatialIndex::R2Tree::R2Tree::~R2Tree()
{
#ifdef HAVE_PTHREAD_H
    pthread_mutex_destroy(&m_lock);
#endif

    storeHeader();
}



//
// ISpatialIndex interface
//

void SpatialIndex::R2Tree::R2Tree::insertData(uint32_t len, const uint8_t* pData, const IShape& shape, id_type id)
{
    throw Tools::NotSupportedException("insertion on R2Tree is not supported now");
}

bool SpatialIndex::R2Tree::R2Tree::deleteData(const SpatialIndex::IShape &shape, SpatialIndex::id_type id) {
    throw Tools::NotSupportedException("deletion on R2Tree is not supported now");
}
void SpatialIndex::R2Tree::R2Tree::containsWhatQuery(const SpatialIndex::IShape &query, SpatialIndex::IVisitor &v) {
    throw Tools::NotSupportedException("not supported now");
}

void SpatialIndex::R2Tree::R2Tree::intersectsWithQuery(const SpatialIndex::IShape &query, SpatialIndex::IVisitor &v) {
    rangeQuery(IntersectionQuery, query, v);
}

void SpatialIndex::R2Tree::R2Tree::pointLocationQuery(const SpatialIndex::Point &query, SpatialIndex::IVisitor &v) {
    throw Tools::NotSupportedException("not supported now");
}

void SpatialIndex::R2Tree::R2Tree::nearestNeighborQuery(uint32_t k, const SpatialIndex::IShape &query,
                                                        SpatialIndex::IVisitor &v) {
    if (query.getDimension() != m_dimension) throw Tools::IllegalArgumentException("nearestNeighborQuery: Shape has the wrong number of dimensions.");
    NNComparator nnc;
    nearestNeighborQuery(k, query, v, nnc);
}
void SpatialIndex::R2Tree::R2Tree::nearestNeighborQuery(uint32_t k, const SpatialIndex::IShape &query,
                                                        SpatialIndex::IVisitor &v,
                                                        SpatialIndex::INearestNeighborComparator &nnc) {
    if (query.getDimension() != m_dimension) throw Tools::IllegalArgumentException("nearestNeighborQuery: Shape has the wrong number of dimensions.");

#ifdef HAVE_PTHREAD_H
    Tools::LockGuard lock(&m_lock);
#endif

    std::priority_queue<NNEntry*, std::vector<NNEntry*>, NNEntry::ascending> queue;

    queue.push(new NNEntry(m_rootID, 0, 0.0));

    uint32_t count = 0;
    double knearest = 0.0;

    while (! queue.empty())
    {
        NNEntry* pFirst = queue.top();

        // report all nearest neighbors with equal greatest distances.
        // (neighbors can be more than k, if many happen to have the same greatest distance).
        if (count >= k && pFirst->m_minDist > knearest)	break;

        queue.pop();

        if (pFirst->m_pEntry == nullptr)
        {
            // n is a leaf or an index.
            NodePtr n = readNode(pFirst->m_id);
            v.visitNode(*n);

            for (uint32_t cChild = 0; cChild < n->m_children; ++cChild)
            {
                if (n->m_level == 0)
                {
                    Data* e = new Data(n->m_pDataLength[cChild], n->m_pData[cChild], *(n->m_ptrMbbc[cChild]), n->m_pIdentifier[cChild]);
                    // we need to compare the query with the actual data entry here, so we call the
                    // appropriate getMinimumDistance method of NearestNeighborComparator.
                    if(m_DataType==TrajectoryType){
                        Trajectory traj;
                        traj.loadFromByteArray(e->m_pData);
                        queue.push(new NNEntry(n->m_pIdentifier[cChild], e, nnc.getMinimumDistance(query,traj)));
                    }else{
                        queue.push(new NNEntry(n->m_pIdentifier[cChild], e, nnc.getMinimumDistance(query, *e)));
                    }
                }
                else
                {
                    queue.push(new NNEntry(n->m_pIdentifier[cChild], 0, nnc.getMinimumDistance(query, *(n->m_ptrMbbc[cChild]))));
                }
            }
        }
        else
        {
            v.visitData(*(static_cast<IData*>(pFirst->m_pEntry)));
            ++(m_stats.m_u64QueryResults);
            ++count;
            knearest = pFirst->m_minDist;
//            std::cout<<"knearest is"<<knearest<<std::endl;
            delete pFirst->m_pEntry;
        }
        delete pFirst;
    }

    while (! queue.empty())
    {
        NNEntry* e = queue.top(); queue.pop();
        if (e->m_pEntry != 0) delete e->m_pEntry;
        delete e;
    }
//    std::cout<<"knearest is"<<knearest<<std::endl;
    m_stats.m_doubleExactQueryResults+=knearest;
}

void SpatialIndex::R2Tree::R2Tree::selfJoinQuery(const SpatialIndex::IShape &s, SpatialIndex::IVisitor &v) {
    throw Tools::NotSupportedException("not supported now");
}

void SpatialIndex::R2Tree::R2Tree::queryStrategy(IQueryStrategy& qs)
{
#ifdef HAVE_PTHREAD_H
    Tools::LockGuard lock(&m_lock);
#endif

    id_type next = m_rootID;
    bool hasNext = true;

    while (hasNext)
    {
        NodePtr n = readNode(next);
        qs.getNextEntry(*n, next, hasNext);
    }
}


void SpatialIndex::R2Tree::R2Tree::getIndexProperties(Tools::PropertySet& out) const
{
    Tools::Variant var;

    // dimension
    var.m_varType = Tools::VT_ULONG;
    var.m_val.ulVal = m_dimension;
    out.setProperty("Dimension", var);

    // index capacity
    var.m_varType = Tools::VT_ULONG;
    var.m_val.ulVal = m_indexCapacity;
    out.setProperty("IndexCapacity", var);

    // leaf capacity
    var.m_varType = Tools::VT_ULONG;
    var.m_val.ulVal = m_leafCapacity;
    out.setProperty("LeafCapacity", var);


    // fill factor
    var.m_varType = Tools::VT_DOUBLE;
    var.m_val.dblVal = m_fillFactor;
    out.setProperty("FillFactor", var);


    // tight MBRs
    var.m_varType = Tools::VT_BOOL;
    var.m_val.blVal = m_bTightMBRs;
    out.setProperty("EnsureTightMBRs", var);

    // index pool capacity
    var.m_varType = Tools::VT_ULONG;
    var.m_val.ulVal = m_indexPool.getCapacity();
    out.setProperty("IndexPoolCapacity", var);

    // leaf pool capacity
    var.m_varType = Tools::VT_ULONG;
    var.m_val.ulVal = m_leafPool.getCapacity();
    out.setProperty("LeafPoolCapacity", var);

    // region pool capacity
    var.m_varType = Tools::VT_ULONG;
    var.m_val.ulVal = m_MbbcPool.getCapacity();
    out.setProperty("RegionPoolCapacity", var);

    // point pool capacity
    var.m_varType = Tools::VT_ULONG;
    var.m_val.ulVal = m_pointPool.getCapacity();
    out.setProperty("PointPoolCapacity", var);
}


void SpatialIndex::R2Tree::R2Tree::addCommand(ICommand* pCommand, CommandType ct)
{
    throw Tools::NotSupportedException("not supported now");
}

bool SpatialIndex::R2Tree::R2Tree::isIndexValid() {
    return true;
}

void SpatialIndex::R2Tree::R2Tree::getStatistics(IStatistics** out) const
{
    *out = new Statistics(m_stats);
}


void SpatialIndex::R2Tree::R2Tree::initNew(Tools::PropertySet& ps)
{
    Tools::Variant var;


    // fill factor
    // it cannot be larger than 50%, since linear and quadratic split algorithms
    // require assigning to both nodes the same number of entries.
    var = ps.getProperty("FillFactor");
    if (var.m_varType != Tools::VT_EMPTY)
    {
        if (var.m_varType != Tools::VT_DOUBLE)
            throw Tools::IllegalArgumentException("initNew: Property FillFactor was not of type Tools::VT_DOUBLE");

        if (var.m_val.dblVal <= 0.0)
            throw Tools::IllegalArgumentException("initNew: Property FillFactor was less than 0.0");

        if ( var.m_val.dblVal >= 1.0)
            throw Tools::IllegalArgumentException(  "initNew: Property FillFactor must be in range "
                                                    "(0.0, 1.0) for RSTAR index type");
        m_fillFactor = var.m_val.dblVal;
    }

    // index capacity
    var = ps.getProperty("IndexCapacity");
    if (var.m_varType != Tools::VT_EMPTY)
    {
        if (var.m_varType != Tools::VT_ULONG || var.m_val.ulVal < 4)
            throw Tools::IllegalArgumentException("initNew: Property IndexCapacity must be Tools::VT_ULONG and >= 4");

        m_indexCapacity = var.m_val.ulVal;
    }

    // leaf capacity
    var = ps.getProperty("LeafCapacity");
    if (var.m_varType != Tools::VT_EMPTY)
    {
        if (var.m_varType != Tools::VT_ULONG )//|| var.m_val.ulVal < 4)
            throw Tools::IllegalArgumentException("initNew: Property LeafCapacity must be Tools::VT_ULONG and >= 4");

        m_leafCapacity = var.m_val.ulVal;
    }


    // dimension
    var = ps.getProperty("Dimension");
    if (var.m_varType != Tools::VT_EMPTY)
    {
        if (var.m_varType != Tools::VT_ULONG)
            throw Tools::IllegalArgumentException("initNew: Property Dimension must be Tools::VT_ULONG");
        if (var.m_val.ulVal <= 1)
            throw Tools::IllegalArgumentException("initNew: Property Dimension must be greater than 1");

        m_dimension = var.m_val.ulVal;
    }

    // tight MBRs
    var = ps.getProperty("EnsureTightMBRs");
    if (var.m_varType != Tools::VT_EMPTY)
    {
        if (var.m_varType != Tools::VT_BOOL)
            throw Tools::IllegalArgumentException("initNew: Property EnsureTightMBRs must be Tools::VT_BOOL");

        m_bTightMBRs = var.m_val.blVal;
    }

    // index pool capacity
    var = ps.getProperty("IndexPoolCapacity");
    if (var.m_varType != Tools::VT_EMPTY)
    {
        if (var.m_varType != Tools::VT_ULONG)
            throw Tools::IllegalArgumentException("initNew: Property IndexPoolCapacity must be Tools::VT_ULONG");

        m_indexPool.setCapacity(var.m_val.ulVal);
    }

    // leaf pool capacity
    var = ps.getProperty("LeafPoolCapacity");
    if (var.m_varType != Tools::VT_EMPTY)
    {
        if (var.m_varType != Tools::VT_ULONG)
            throw Tools::IllegalArgumentException("initNew: Property LeafPoolCapacity must be Tools::VT_ULONG");

        m_leafPool.setCapacity(var.m_val.ulVal);
    }

    // region pool capacity
    var = ps.getProperty("RegionPoolCapacity");
    if (var.m_varType != Tools::VT_EMPTY)
    {
        if (var.m_varType != Tools::VT_ULONG)
            throw Tools::IllegalArgumentException("initNew: Property RegionPoolCapacity must be Tools::VT_ULONG");

        m_MbbcPool.setCapacity(var.m_val.ulVal);
    }

    // point pool capacity
    var = ps.getProperty("PointPoolCapacity");
    if (var.m_varType != Tools::VT_EMPTY)
    {
        if (var.m_varType != Tools::VT_ULONG)
            throw Tools::IllegalArgumentException("initNew: Property PointPoolCapacity must be Tools::VT_ULONG");

        m_pointPool.setCapacity(var.m_val.ulVal);
    }

    m_infiniteMbbc.makeInfinite(m_dimension);

    m_stats.m_u32TreeHeight = 1;
    m_stats.m_nodesInLevel.emplace_back(0);

    Leaf root(this, -1);
    m_rootID = writeNode(&root);

    storeHeader();
}

void SpatialIndex::R2Tree::R2Tree::initOld(Tools::PropertySet& ps)
{
    loadHeader();

    // only some of the properties may be changed.
    // the rest are just ignored.

    Tools::Variant var;

    // tight MBRs
    var = ps.getProperty("EnsureTightMBRs");
    if (var.m_varType != Tools::VT_EMPTY)
    {
        if (var.m_varType != Tools::VT_BOOL) throw Tools::IllegalArgumentException("initOld: Property EnsureTightMBRs must be Tools::VT_BOOL");

        m_bTightMBRs = var.m_val.blVal;
    }

    // index pool capacity
    var = ps.getProperty("IndexPoolCapacity");
    if (var.m_varType != Tools::VT_EMPTY)
    {
        if (var.m_varType != Tools::VT_ULONG) throw Tools::IllegalArgumentException("initOld: Property IndexPoolCapacity must be Tools::VT_ULONG");

        m_indexPool.setCapacity(var.m_val.ulVal);
    }

    // leaf pool capacity
    var = ps.getProperty("LeafPoolCapacity");
    if (var.m_varType != Tools::VT_EMPTY)
    {
        if (var.m_varType != Tools::VT_ULONG) throw Tools::IllegalArgumentException("initOld: Property LeafPoolCapacity must be Tools::VT_ULONG");

        m_leafPool.setCapacity(var.m_val.ulVal);
    }

    // region pool capacity
    var = ps.getProperty("RegionPoolCapacity");
    if (var.m_varType != Tools::VT_EMPTY)
    {
        if (var.m_varType != Tools::VT_ULONG) throw Tools::IllegalArgumentException("initOld: Property RegionPoolCapacity must be Tools::VT_ULONG");

        m_MbbcPool.setCapacity(var.m_val.ulVal);
    }

    // point pool capacity
    var = ps.getProperty("PointPoolCapacity");
    if (var.m_varType != Tools::VT_EMPTY)
    {
        if (var.m_varType != Tools::VT_ULONG) throw Tools::IllegalArgumentException("initOld: Property PointPoolCapacity must be Tools::VT_ULONG");

        m_pointPool.setCapacity(var.m_val.ulVal);
    }

    m_infiniteMbbc.makeInfinite(m_dimension);
}

void SpatialIndex::R2Tree::R2Tree::storeHeader()
{
    const uint32_t headerSize =
            sizeof(id_type) +						// m_rootID
            sizeof(double) +						// m_fillFactor
            sizeof(uint32_t) +						// m_indexCapacity
            sizeof(uint32_t) +						// m_leafCapacity
            sizeof(uint32_t) +						// m_dimension
            sizeof(char) +							// m_bTightMBRs
            sizeof(uint32_t) +						// m_stats.m_nodes
            sizeof(uint64_t) +						// m_stats.m_data
            sizeof(uint32_t) +						// m_stats.m_treeHeight
            m_stats.m_u32TreeHeight * sizeof(uint32_t);	// m_stats.m_nodesInLevel

    uint8_t* header = new uint8_t[headerSize];
    uint8_t* ptr = header;

    memcpy(ptr, &m_rootID, sizeof(id_type));
    ptr += sizeof(id_type);
    memcpy(ptr, &m_fillFactor, sizeof(double));
    ptr += sizeof(double);
    memcpy(ptr, &m_indexCapacity, sizeof(uint32_t));
    ptr += sizeof(uint32_t);
    memcpy(ptr, &m_leafCapacity, sizeof(uint32_t));
    ptr += sizeof(uint32_t);
    memcpy(ptr, &m_dimension, sizeof(uint32_t));
    ptr += sizeof(uint32_t);
    char c = (char) m_bTightMBRs;
    memcpy(ptr, &c, sizeof(char));
    ptr += sizeof(char);
    memcpy(ptr, &(m_stats.m_u32Nodes), sizeof(uint32_t));
    ptr += sizeof(uint32_t);
    memcpy(ptr, &(m_stats.m_u64Data), sizeof(uint64_t));
    ptr += sizeof(uint64_t);
    memcpy(ptr, &(m_stats.m_u32TreeHeight), sizeof(uint32_t));
    ptr += sizeof(uint32_t);

    for (uint32_t cLevel = 0; cLevel < m_stats.m_u32TreeHeight; ++cLevel)
    {
        memcpy(ptr, &(m_stats.m_nodesInLevel[cLevel]), sizeof(uint32_t));
        ptr += sizeof(uint32_t);
    }

    m_pStorageManager->storeByteArray(m_headerID, headerSize, header);

    delete[] header;
}



void SpatialIndex::R2Tree::R2Tree::loadHeader()
{
    uint32_t headerSize;
    uint8_t* header = 0;
    m_pStorageManager->loadByteArray(m_headerID, headerSize, &header);

    uint8_t* ptr = header;

    memcpy(&m_rootID, ptr, sizeof(id_type));
    ptr += sizeof(id_type);
    memcpy(&m_fillFactor, ptr, sizeof(double));
    ptr += sizeof(double);
    memcpy(&m_indexCapacity, ptr, sizeof(uint32_t));
    ptr += sizeof(uint32_t);
    memcpy(&m_leafCapacity, ptr, sizeof(uint32_t));
    ptr += sizeof(uint32_t);
    memcpy(&m_dimension, ptr, sizeof(uint32_t));
    ptr += sizeof(uint32_t);
    char c;
    memcpy(&c, ptr, sizeof(char));
    m_bTightMBRs = (c != 0);
    ptr += sizeof(char);
    memcpy(&(m_stats.m_u32Nodes), ptr, sizeof(uint32_t));
    ptr += sizeof(uint32_t);
    memcpy(&(m_stats.m_u64Data), ptr, sizeof(uint64_t));
    ptr += sizeof(uint64_t);
    memcpy(&(m_stats.m_u32TreeHeight), ptr, sizeof(uint32_t));
    ptr += sizeof(uint32_t);

    for (uint32_t cLevel = 0; cLevel < m_stats.m_u32TreeHeight; ++cLevel)
    {
        uint32_t cNodes;
        memcpy(&cNodes, ptr, sizeof(uint32_t));
        ptr += sizeof(uint32_t);
        m_stats.m_nodesInLevel.emplace_back(cNodes);
    }

    delete[] header;
}


SpatialIndex::R2Tree::NodePtr SpatialIndex::R2Tree::R2Tree::readNode(id_type page)
{
    uint32_t dataLength;
    uint8_t* buffer;

    try
    {
        m_pStorageManager->loadByteArray(page, dataLength, &buffer);
    }
    catch (InvalidPageException& e)
    {
        std::cerr << e.what() << std::endl;
        throw;
    }

    try
    {
        uint32_t nodeType;
        memcpy(&nodeType, buffer, sizeof(uint32_t));

        NodePtr n;

        if (nodeType == PersistentIndex) n = m_indexPool.acquire();
        else if (nodeType == PersistentLeaf) n = m_leafPool.acquire();
        else throw Tools::IllegalStateException("readNode: failed reading the correct node type information");

        if (n.get() == 0)
        {
            if (nodeType == PersistentIndex) n = NodePtr(new Index(this, -1, 0), &m_indexPool);
            else if (nodeType == PersistentLeaf) n = NodePtr(new Leaf(this, -1), &m_leafPool);
        }

        //n->m_pTree = this;
        n->m_identifier = page;
        n->loadFromByteArray(buffer);

        ++(m_stats.m_u64Reads);
        /*
        for (size_t cIndex = 0; cIndex < m_readNodeCommands.size(); ++cIndex)
        {
            m_readNodeCommands[cIndex]->execute(*n);
        }
        */
        delete[] buffer;
        return n;
    }
    catch (...)
    {
        delete[] buffer;
        throw;
    }
}
SpatialIndex::id_type SpatialIndex::R2Tree::R2Tree::writeNode(Node* n)
{
    uint8_t* buffer;
    uint32_t dataLength;
    n->storeToByteArray(&buffer, dataLength);

    id_type page;
    if (n->m_identifier < 0) page = StorageManager::NewPage;
    else page = n->m_identifier;

    try
    {
        m_pStorageManager->storeByteArray(page, dataLength, buffer);
        delete[] buffer;
    }
    catch (InvalidPageException& e)
    {
        delete[] buffer;
        std::cerr << e.what() << std::endl;
        throw;
    }

    if (n->m_identifier < 0)
    {
        n->m_identifier = page;
        ++(m_stats.m_u32Nodes);

#ifndef NDEBUG
        try
        {
            m_stats.m_nodesInLevel[n->m_level] = m_stats.m_nodesInLevel.at(n->m_level) + 1;
        }
        catch(...)
        {
            throw Tools::IllegalStateException("writeNode: writing past the end of m_nodesInLevel.");
        }
#else
        m_stats.m_nodesInLevel[n->m_level] = m_stats.m_nodesInLevel[n->m_level] + 1;
#endif
    }

    ++(m_stats.m_u64Writes);
    /*
    for (size_t cIndex = 0; cIndex < m_writeNodeCommands.size(); ++cIndex)
    {
        m_writeNodeCommands[cIndex]->execute(*n);
    }
    */
    return page;
}
void SpatialIndex::R2Tree::R2Tree::deleteNode(Node* n)
{
    try
    {
        m_pStorageManager->deleteByteArray(n->m_identifier);
    }
    catch (InvalidPageException& e)
    {
        std::cerr << e.what() << std::endl;
        throw;
    }

    --(m_stats.m_u32Nodes);
    m_stats.m_nodesInLevel[n->m_level] = m_stats.m_nodesInLevel[n->m_level] - 1;
    /*
    for (size_t cIndex = 0; cIndex < m_deleteNodeCommands.size(); ++cIndex)
    {
        m_deleteNodeCommands[cIndex]->execute(*n);
    }
     */
}


/*
void SpatialIndex::R2Tree::R2Tree::insertData_impl(uint32_t dataLength, uint8_t* pData, Mbbc& mbbc, id_type id)
{
    assert(mbbc.getDimension() == m_dimension);

    std::stack<id_type> pathBuffer;
    uint8_t* overflowTable = 0;

    try
    {
        NodePtr root = readNode(m_rootID);

        overflowTable = new uint8_t[root->m_level];
        memset(overflowTable, 0, root->m_level);

        NodePtr l = root->chooseSubtree(mbbc, 0, pathBuffer);
        if (l.get() == root.get())
        {
            assert(root.unique());
            root.relinquish();
        }
        l->insertData(dataLength, pData, mbbc, id, pathBuffer, overflowTable);

        delete[] overflowTable;
        ++(m_stats.m_u64Data);
    }
    catch (...)
    {
        delete[] overflowTable;
        throw;
    }
}
 */


std::ostream& SpatialIndex::R2Tree::operator<<(std::ostream& os, const R2Tree& t)
{
    os	<< "Dimension: " << t.m_dimension << std::endl
          << "Fill factor: " << t.m_fillFactor << std::endl
          << "Index capacity: " << t.m_indexCapacity << std::endl
          << "Leaf capacity: " << t.m_leafCapacity << std::endl
          << "Tight MBRs: " << ((t.m_bTightMBRs) ? "enabled" : "disabled") << std::endl;


    if (t.m_stats.getNumberOfNodesInLevel(0) > 0)
        os	<< "Utilization: " << 100 * t.m_stats.getNumberOfData() / (t.m_stats.getNumberOfNodesInLevel(0) * t.m_leafCapacity) << "%" << std::endl
              << t.m_stats;

#ifndef NDEBUG
    os	<< "Leaf pool hits: " << t.m_leafPool.m_hits << std::endl
          << "Leaf pool misses: " << t.m_leafPool.m_misses << std::endl
          << "Index pool hits: " << t.m_indexPool.m_hits << std::endl
          << "Index pool misses: " << t.m_indexPool.m_misses << std::endl
          << "Mbbc pool hits: " << t.m_MbbcPool.m_hits << std::endl
          << "Mbbc pool misses: " << t.m_MbbcPool.m_misses << std::endl
          << "Point pool hits: " << t.m_pointPool.m_hits << std::endl
          << "Point pool misses: " << t.m_pointPool.m_misses << std::endl;
#endif
    return os;
}




void SpatialIndex::R2Tree::R2Tree::rangeQuery(RangeQueryType type, const IShape& query, IVisitor& v)
{
#ifdef HAVE_PTHREAD_H
    Tools::LockGuard lock(&m_lock);
#endif

    std::stack<NodePtr> st;
    NodePtr root = readNode(m_rootID);

    if (root->m_children > 0 && root->m_nodeMbbc.intersectsShape(query)) st.push(root);

    while (! st.empty())
    {
        NodePtr n = st.top(); st.pop();
//        std::cout<<"\n level"<<n->m_level<<"\n node\n"<<n->m_nodeMbbc.toString();
        if (n->m_level == 0)
        {
            v.visitNode(*n);
            for (uint32_t cChild = 0; cChild < n->m_children; ++cChild)
            {
                bool b;
                if (type == ContainmentQuery) b = n->m_ptrMbbc[cChild]->containsShape(query);
                else b = n->m_ptrMbbc[cChild]->intersectsShape(query);

                if (b)
                {
                    Data data = Data(n->m_pDataLength[cChild], n->m_pData[cChild], *(n->m_ptrMbbc[cChild]), n->m_pIdentifier[cChild]);
                    ++(m_stats.m_u64QueryResults);
                    if(m_DataType==TrajectoryType){
                        Trajectory traj;
                        traj.loadFromByteArray(data.m_pData);
                        if(traj.intersectsShape(query)){
                            m_stats.m_doubleExactQueryResults+=1;
                            v.visitData(data);
                        }
                    }else{
                        v.visitData(data);
                    }
                }
//                else{
//                    std::cout<<"ack failed\n"<<query.toString()<<"\n"<<n->m_ptrMbbc[cChild]->toString()<<"\n";
//                }
            }
        }
        else
        {
            v.visitNode(*n);
//            if(n->m_level<3) {
                for (uint32_t cChild = 0; cChild < n->m_children; ++cChild) {
                    if (n->m_ptrMbbc[cChild]->intersectsShape(query)) st.push(readNode(n->m_pIdentifier[cChild]));
                }
//            }else{
//                for (uint32_t cChild = 0; cChild < n->m_children; ++cChild) {
//                    if (n->m_ptrMbbc[cChild]->m_wmbr.intersectsShape(query)) st.push(readNode(n->m_pIdentifier[cChild]));
//                }
//            }
        }
    }
}