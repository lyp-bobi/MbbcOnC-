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
//#include "BulkLoader.h"
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
