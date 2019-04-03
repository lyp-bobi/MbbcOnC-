//
// Created by chuang on 4/3/19.
//

#include <cstring>
#include <cmath>
#include <limits>

#include <spatialindex/SpatialIndex.h>

#include "Node.h"


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

//todo: implement these
uint32_t Node::getByteArraySize(){
    return 0;
}

void Node::loadFromByteArray(const byte* ptr){}

void Node::storeToByteArray(byte** data, uint32_t& len){}


//
// SpatialIndex::IEntry interface
//
SpatialIndex::id_type Node::getIdentifier() const
{
    return m_identifier;
}

void Node::getShape(IShape** out) const
{
    *out = new Region(m_nodeMbbc);
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

    *out = new Region(*(m_ptrMbbc[index]));
}

void Node::getChildData(uint32_t index, uint32_t& length, byte** data) const
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

