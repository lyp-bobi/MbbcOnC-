//
// Created by chuang on 4/3/19.
//


#include <cstring>
#include <cmath>
#include <limits>

#include <spatialindex/SpatialIndex.h>
#include "Node.h"
//#include "Leaf.h"
//#include "Index.h"
#include "BulkLoader.h"
#include "R2Tree.h"

using namespace SpatialIndex::R2Tree;
using namespace SpatialIndex;

SpatialIndex::R2Tree::Data::Data(uint32_t len, byte* pData, Mbbc& r, id_type id)
        : m_id(id), m_Mbbc(r), m_pData(0), m_dataLength(len)
{
    if (m_dataLength > 0)
    {
        m_pData = new byte[m_dataLength];
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

void SpatialIndex::R2Tree::Data::getData(uint32_t& len, byte** data) const
{
    len = m_dataLength;
    *data = 0;

    if (m_dataLength > 0)
    {
        *data = new byte[m_dataLength];
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

void SpatialIndex::R2Tree::Data::loadFromByteArray(const byte* ptr)
{
    memcpy(&m_id, ptr, sizeof(id_type));
    ptr += sizeof(id_type);

    delete[] m_pData;
    m_pData = 0;

    memcpy(&m_dataLength, ptr, sizeof(uint32_t));
    ptr += sizeof(uint32_t);

    if (m_dataLength > 0)
    {
        m_pData = new byte[m_dataLength];
        memcpy(m_pData, ptr, m_dataLength);
        ptr += m_dataLength;
    }

    m_Mbbc.loadFromByteArray(ptr);
}

void SpatialIndex::R2Tree::Data::storeToByteArray(byte** data, uint32_t& len)
{
    // it is thread safe this way.
    uint32_t Mbbcsize;
    byte* Mbbcdata = 0;
    m_Mbbc.storeToByteArray(&Mbbcdata, Mbbcsize);

    len = sizeof(id_type) + sizeof(uint32_t) + m_dataLength + Mbbcsize;

    *data = new byte[len];
    byte* ptr = *data;

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
    SpatialIndex::ISpatialIndex* tree = createNewR2Tree(sm, fillFactor, indexCapacity, leafCapacity, dimension, rv, indexIdentifier);

    uint32_t bindex = static_cast<uint32_t>(std::floor(static_cast<double>(indexCapacity * fillFactor)));
    uint32_t bleaf = static_cast<uint32_t>(std::floor(static_cast<double>(leafCapacity * fillFactor)));

    SpatialIndex::R2Tree::BulkLoader bl;

    switch (m)
    {
        case BLM_STR:
            bl.bulkLoadUsingSTR(static_cast<R2Tree*>(tree), stream, bindex, bleaf, 10000, 100);
            break;
        case BLM_KDT:
            bl.bulkLoadUsingKDT(static_cast<R2Tree*>(tree), stream, bindex, bleaf, 10000, 100);
            break;
        default:
            throw Tools::IllegalArgumentException("createAndBulkLoadNewR2Tree: Unknown bulk load method.");
            break;
    }

    return tree;
}