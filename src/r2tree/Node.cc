//
// Created by chuang on 4/3/19.
//

#include <cstring>
#include <cmath>
#include <limits>

#include <spatialindex/SpatialIndex.h>

#include "Node.h"
#include "R2Tree.h"


using namespace SpatialIndex;
using namespace SpatialIndex::R2Tree;


//
// Tools::IObject interface
//
Tools::IObject* Node::clone()
{
    throw Tools::NotSupportedException("IObject::clone should never be called.");
}


//
// Tools::ISerializable interface
//

uint32_t Node::getByteArraySize(){
    return
            (sizeof(uint32_t) +
             sizeof(uint32_t) +
             sizeof(uint32_t) +
             (m_children * (m_nodeMbbc.getByteArraySize() + sizeof(id_type) + sizeof(uint32_t))) +
             m_totalDataLength +
                    m_nodeMbbc.getByteArraySize());
}

void Node::loadFromByteArray(const uint8_t* ptr){

    m_nodeMbbc=m_pTree->m_infiniteMbbc;

    // skip the node type information, it is not needed.
    ptr += sizeof(uint32_t);

    memcpy(&m_level, ptr, sizeof(uint32_t));
    ptr += sizeof(uint32_t);

    memcpy(&m_children, ptr, sizeof(uint32_t));
    ptr += sizeof(uint32_t);

    for (uint32_t u32Child = 0; u32Child < m_children; ++u32Child)
    {
        m_ptrMbbc[u32Child] = m_pTree->m_MbbcPool.acquire();
        *(m_ptrMbbc[u32Child]) = m_pTree->m_infiniteMbbc;

        m_ptrMbbc[u32Child]->loadFromByteArray(ptr);
        ptr += m_ptrMbbc[u32Child]->getByteArraySize();
        memcpy(&(m_pIdentifier[u32Child]), ptr, sizeof(id_type));
        ptr += sizeof(id_type);

        memcpy(&(m_pDataLength[u32Child]), ptr, sizeof(uint32_t));
        ptr += sizeof(uint32_t);

        if (m_pDataLength[u32Child] > 0)
        {
            m_totalDataLength += m_pDataLength[u32Child];
            m_pData[u32Child] = new uint8_t[m_pDataLength[u32Child]];
            memcpy(m_pData[u32Child], ptr, m_pDataLength[u32Child]);
            ptr += m_pDataLength[u32Child];
        }
        else
        {
            m_pData[u32Child] = 0;
        }

        //m_nodeMBR.combineRegion(*(m_ptrMBR[u32Child]));
    }

    m_nodeMbbc.loadFromByteArray(ptr);

}

void Node::storeToByteArray(uint8_t** data, uint32_t& len){
    len = getByteArraySize();

    *data = new uint8_t[len];
    uint8_t* ptr = *data;

    uint32_t nodeType;

    if (m_level == 0) nodeType = PersistentLeaf;
    else nodeType = PersistentIndex;

    memcpy(ptr, &nodeType, sizeof(uint32_t));
    ptr += sizeof(uint32_t);

    memcpy(ptr, &m_level, sizeof(uint32_t));
    ptr += sizeof(uint32_t);

    memcpy(ptr, &m_children, sizeof(uint32_t));
    ptr += sizeof(uint32_t);
    uint8_t* tmpb;
    uint32_t tmplen;
    for (uint32_t u32Child = 0; u32Child < m_children; ++u32Child)
    {

        m_ptrMbbc[u32Child]->storeToByteArray(&tmpb,tmplen);
        memcpy(ptr, tmpb, tmplen);
        ptr += tmplen;

        memcpy(ptr, &(m_pIdentifier[u32Child]), sizeof(id_type));
        ptr += sizeof(id_type);

        memcpy(ptr, &(m_pDataLength[u32Child]), sizeof(uint32_t));
        ptr += sizeof(uint32_t);

        if (m_pDataLength[u32Child] > 0)
        {
            memcpy(ptr, m_pData[u32Child], m_pDataLength[u32Child]);
            ptr += m_pDataLength[u32Child];
        }
    }

    m_nodeMbbc.storeToByteArray(&tmpb,tmplen);
    memcpy(ptr, tmpb, tmplen);
    //ptr += tmplen;
    assert(len == (ptr - *data)+tmplen);

}


//
// SpatialIndex::IEntry interface
//
SpatialIndex::id_type Node::getIdentifier() const
{
    return m_identifier;
}

void Node::getShape(IShape** out) const
{
    *out = new Mbbc(m_nodeMbbc);
}

//
// SpatialIndex::INode interface
//
uint32_t Node::getChildrenCount() const
{
    return m_children;
}

SpatialIndex::id_type Node::getChildIdentifier(uint32_t index) const
{
    if (index >= m_children) throw Tools::IndexOutOfBoundsException(index);

    return m_pIdentifier[index];
}

void Node::getChildShape(uint32_t index, IShape** out) const
{
    if (index >= m_children) throw Tools::IndexOutOfBoundsException(index);

    *out = new Mbbc(*(m_ptrMbbc[index]));
}

void Node::getChildData(uint32_t index, uint32_t& length, uint8_t** data) const
{
    if (index >= m_children) throw Tools::IndexOutOfBoundsException(index);
    if (m_pData[index] == NULL)
    {
        length = 0;
        data = NULL;
    }
    else
    {
        length = m_pDataLength[index];
        *data = m_pData[index];
    }
}

uint32_t Node::getLevel() const
{
    return m_level;
}

bool Node::isLeaf() const
{
    return (m_level == 0);
}

bool Node::isIndex() const
{
    return (m_level != 0);
}


//
// Internal
//

Node::Node() :
        m_pTree(0),
        m_level(0),
        m_identifier(-1),
        m_children(0),
        m_capacity(0),
        m_pData(0),
        m_ptrMbbc(0),
        m_pIdentifier(0),
        m_pDataLength(0),
        m_totalDataLength(0)
{
}

Node::Node(SpatialIndex::R2Tree::R2Tree* pTree, id_type id, uint32_t level, uint32_t capacity) :
        m_pTree(pTree),
        m_level(level),
        m_identifier(id),
        m_children(0),
        m_capacity(capacity),
        m_pData(0),
        m_ptrMbbc(0),
        m_pIdentifier(0),
        m_pDataLength(0),
        m_totalDataLength(0)
{
    m_nodeMbbc.makeInfinite();

    try
    {
        m_pDataLength = new uint32_t[m_capacity + 1];
        m_pData = new uint8_t*[m_capacity + 1];
        m_ptrMbbc = new MbbcPtr[m_capacity + 1];
        m_pIdentifier = new id_type[m_capacity + 1];
    }
    catch (...)
    {
        delete[] m_pDataLength;
        delete[] m_pData;
        delete[] m_ptrMbbc;
        delete[] m_pIdentifier;
        throw;
    }
}


Node::~Node(){
    if (m_pData != 0)
    {
        for (uint32_t u32Child = 0; u32Child < m_children; ++u32Child)
        {
            if (m_pData[u32Child] != 0) delete[] m_pData[u32Child];
        }

        delete[] m_pData;
    }

    delete[] m_pDataLength;
    delete[] m_ptrMbbc;
    delete[] m_pIdentifier;
}

Node& Node::operator=(const Node&){
    throw Tools::IllegalStateException("operator =: This should never be called.");
}


void Node::insertEntry(uint32_t dataLength, uint8_t* pData, Mbbc& mbbc, id_type id)
{
    assert(m_children < m_capacity);

    m_pDataLength[m_children] = dataLength;
    m_pData[m_children] = pData;
    m_ptrMbbc[m_children] = m_pTree->m_MbbcPool.acquire();
    *(m_ptrMbbc[m_children]) = mbbc;
    m_pIdentifier[m_children] = id;

    m_totalDataLength += dataLength;
    ++m_children;

    m_nodeMbbc.combineMbbc(mbbc);
}


/*


bool Node::insertData(uint32_t dataLength, uint8_t* pData, Mbbc& mbbc, id_type id, std::stack<id_type>& pathBuffer, uint8_t* overflowTable)
{
    if (m_children < m_capacity)
    {
        bool adjusted = false;

        // this has to happen before insertEntry modifies m_nodeMBR.
        bool b = m_nodeMbbc.containsMbbc(mbbc);

        insertEntry(dataLength, pData, mbbc, id);
        m_pTree->writeNode(this);

        if ((! b) && (! pathBuffer.empty()))
        {
            id_type cParent = pathBuffer.top(); pathBuffer.pop();
            NodePtr ptrN = m_pTree->readNode(cParent);
            Index* p = static_cast<Index*>(ptrN.get());
            p->adjustTree(this, pathBuffer);
            adjusted = true;
        }

        return adjusted;
    }
    else
    {
        NodePtr n;
        NodePtr nn;
        split(dataLength, pData, mbbc, id, n, nn);

        if (pathBuffer.empty())
        {
            n->m_level = m_level;
            nn->m_level = m_level;
            n->m_identifier = -1;
            nn->m_identifier = -1;
            m_pTree->writeNode(n.get());
            m_pTree->writeNode(nn.get());

            NodePtr ptrR = m_pTree->m_indexPool.acquire();
            if (ptrR.get() == 0)
            {
                ptrR = NodePtr(new Index(m_pTree, m_pTree->m_rootID, m_level + 1), &(m_pTree->m_indexPool));
            }
            else
            {
                //ptrR->m_pTree = m_pTree;
                ptrR->m_identifier = m_pTree->m_rootID;
                ptrR->m_level = m_level + 1;
                ptrR->m_nodeMBR = m_pTree->m_infiniteRegion;
            }

            ptrR->insertEntry(0, 0, n->m_nodeMBR, n->m_identifier);
            ptrR->insertEntry(0, 0, nn->m_nodeMBR, nn->m_identifier);

            m_pTree->writeNode(ptrR.get());

            m_pTree->m_stats.m_nodesInLevel[m_level] = 2;
            m_pTree->m_stats.m_nodesInLevel.emplace_back(1);
            m_pTree->m_stats.m_u32TreeHeight = m_level + 2;
        }
        else
        {
            n->m_level = m_level;
            nn->m_level = m_level;
            n->m_identifier = m_identifier;
            nn->m_identifier = -1;

            m_pTree->writeNode(n.get());
            m_pTree->writeNode(nn.get());

            id_type cParent = pathBuffer.top(); pathBuffer.pop();
            NodePtr ptrN = m_pTree->readNode(cParent);
            Index* p = static_cast<Index*>(ptrN.get());
            p->adjustTree(n.get(), nn.get(), pathBuffer, overflowTable);
        }

        return true;
    }
}
 */