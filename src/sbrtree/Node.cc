/******************************************************************************
 * Project:  libspatialindex - A C++ library for spatial indexing
 * Author:   Marios Hadjieleftheriou, mhadji@gmail.com
 ******************************************************************************
 * Copyright (c) 2002, Marios Hadjieleftheriou
 *
 * All rights reserved.
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
******************************************************************************/

#include <cstring>
#include <cmath>
#include <limits>

#include <spatialindex/SpatialIndex.h>

#include "SBRTree.h"
#include "Node.h"
#include "Index.h"

using namespace SpatialIndex;
using namespace SpatialIndex::SBRTree;

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
uint32_t Node::getByteArraySize() const
{
    //todo: don't store the dimensions of each sbr/mbc
    if(m_level>0)
        return
                (sizeof(uint32_t) +
                 sizeof(uint32_t) +
                 sizeof(uint32_t) +
                 (m_children * (sizeof(uint32_t)+(m_pTree->m_dimension+1) * sizeof(double) * 2 +sizeof(double) + sizeof(id_type) + sizeof(uint32_t))) +
                 m_totalDataLength +
                 m_nodeSBR.getByteArraySize());
    else
        return
                (sizeof(uint32_t) +
                 sizeof(uint32_t) +
                 sizeof(uint32_t) +
                 (m_children * (sizeof(uint32_t)+(m_pTree->m_dimension+1) * sizeof(double) * 2+ sizeof(double)*2 + sizeof(id_type) + sizeof(uint32_t))) +
                 m_totalDataLength +
                 m_nodeSBR.getByteArraySize());
}

uint32_t Node::getIndexByteArraySize() const
{
    //todo: don't store the dimensions of each sbr/mbc
    if(m_level>0)
        return
                (sizeof(uint32_t) +
                 sizeof(uint32_t) +
                 sizeof(uint32_t) +
                 (m_children * (sizeof(uint32_t)+(m_pTree->m_dimension+1) * sizeof(double) * 2+sizeof(double) + sizeof(id_type) + sizeof(uint32_t))) +
                 m_nodeSBR.getByteArraySize());
    else
        return
                (sizeof(uint32_t) +
                 sizeof(uint32_t) +
                 sizeof(uint32_t) +
                 (m_children * (sizeof(uint32_t)+(m_pTree->m_dimension+1) * sizeof(double) * 2+ sizeof(double)*2 + sizeof(id_type) + sizeof(uint32_t))) +
                 m_nodeSBR.getByteArraySize());
}

void Node::loadFromByteArray(const uint8_t* ptr)
{
    m_nodeSBR=m_pTree->m_infiniteSBR;

    // skip the node type information, it is not needed.
    ptr += sizeof(uint32_t);

    memcpy(&m_level, ptr, sizeof(uint32_t));
    ptr += sizeof(uint32_t);

    memcpy(&m_children, ptr, sizeof(uint32_t));
    ptr += sizeof(uint32_t);
    if(m_level>0)
        m_ptrSBR = new SBRPtr[m_children + 1];
    else
        m_ptrMBC = new MBCPtr[m_children + 1];
    for (uint32_t u32Child = 0; u32Child < m_children; ++u32Child) {
        if(m_level>0) {
            m_ptrSBR[u32Child] = m_pTree->m_SBRPool.acquire();
            *(m_ptrSBR[u32Child]) = m_pTree->m_infiniteSBR;
            m_ptrSBR[u32Child]->loadFromByteArray(ptr);
            ptr += m_ptrSBR[u32Child]->getByteArraySize();
        }
        else{
            m_ptrMBC[u32Child] = m_pTree->m_mbcPool.acquire();
            MBC bc=MBC();bc.makeInfinite(m_pTree->m_dimension);
            *(m_ptrMBC[u32Child]) = bc;
            m_ptrMBC[u32Child]->loadFromByteArray(ptr);
            ptr += m_ptrMBC[u32Child]->getByteArraySize();
        }
        memcpy(&(m_pIdentifier[u32Child]), ptr, sizeof(id_type));
        ptr += sizeof(id_type);

        memcpy(&(m_pDataLength[u32Child]), ptr, sizeof(uint32_t));
        ptr += sizeof(uint32_t);

        if (m_pDataLength[u32Child] > 0) {
            m_totalDataLength += m_pDataLength[u32Child];
            m_pData[u32Child] = new uint8_t[m_pDataLength[u32Child]];
            memcpy(m_pData[u32Child], ptr, m_pDataLength[u32Child]);
            ptr += m_pDataLength[u32Child];
        } else {
            m_pData[u32Child] = 0;
        }

        //m_nodeSBR.combineSBR(*(m_ptrSBR[u32Child]));
    }
    m_nodeSBR.loadFromByteArray(ptr);
}

void Node::storeToByteArray(uint8_t** data, uint32_t& len)
{
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

    for (uint32_t u32Child = 0; u32Child < m_children; ++u32Child) {
        if(m_level>0)
            m_ptrSBR[u32Child]->storeToByteArray(&tmpb, tmplen);
        else
            m_ptrMBC[u32Child]->storeToByteArray(&tmpb, tmplen);
        memcpy(ptr, tmpb, tmplen);
        ptr += tmplen;

        memcpy(ptr, &(m_pIdentifier[u32Child]), sizeof(id_type));
        ptr += sizeof(id_type);

        memcpy(ptr, &(m_pDataLength[u32Child]), sizeof(uint32_t));
        ptr += sizeof(uint32_t);

        if (m_pDataLength[u32Child] > 0) {
            memcpy(ptr, m_pData[u32Child], m_pDataLength[u32Child]);
            ptr += m_pDataLength[u32Child];
        }
    }
    m_nodeSBR.storeToByteArray(&tmpb,tmplen);
    memcpy(ptr, tmpb, tmplen);
    //ptr += tmplen;
//    std::cerr<<len<<" "<<(ptr - *data)+tmplen<<" "<<m_children<<"\n";
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
	*out = new SBR(m_nodeSBR);
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

	*out = new SBR(*(m_ptrSBR[index]));
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
	m_ptrSBR(0),
	m_pIdentifier(0),
	m_pDataLength(0),
	m_totalDataLength(0)
{
}

Node::Node(SpatialIndex::SBRTree::SBRTree* pTree, id_type id, uint32_t level, uint32_t capacity) :
	m_pTree(pTree),
	m_level(level),
	m_identifier(id),
	m_children(0),
	m_capacity(capacity),
	m_pData(0),
	m_ptrSBR(0),
    m_ptrMBC(0),
	m_pIdentifier(0),
	m_pDataLength(0),
	m_totalDataLength(0)
{
	m_nodeSBR.makeInfinite(m_pTree->m_dimension);

	try
	{
		m_pDataLength = new uint32_t[m_capacity + 1];
		m_pData = new uint8_t*[m_capacity + 1];
		if(m_level==0)
            m_ptrMBC = new MBCPtr[m_capacity + 1];
		else
            m_ptrSBR = new SBRPtr[m_capacity + 1];
		m_pIdentifier = new id_type[m_capacity + 1];
	}
	catch (...)
	{
		delete[] m_pDataLength;
		delete[] m_pData;
		delete[] m_ptrSBR;
		delete[] m_pIdentifier;
		throw;
	}
}

Node::~Node()
{
	if (m_pData != 0)
	{
		for (uint32_t u32Child = 0; u32Child < m_children; ++u32Child)
		{
			if (m_pData[u32Child] != 0) delete[] m_pData[u32Child];
		}

		delete[] m_pData;
	}

	delete[] m_pDataLength;
	delete[] m_ptrSBR;
	delete[] m_pIdentifier;
}

Node& Node::operator=(const Node&)
{
	throw Tools::IllegalStateException("operator =: This should never be called.");
}

//void Node::insertEntry(uint32_t dataLength, uint8_t *pData, SpatialIndex::IShape &shape, SpatialIndex::id_type id) {
//    const SBR *pr= dynamic_cast<SBR*>(&shape);
//    if(pr!= nullptr) {
//        SBR newr(*pr);
//        insertEntry(dataLength,pData,newr,id);
//    }
//    const MBC *pbc= dynamic_cast<MBC*>(&shape);
//    if(pbc!= nullptr) {
//        MBC newbc(*pbc);
//        insertEntry(dataLength,pData,newbc,id);
//    }
//}

void Node::insertEntry(uint32_t dataLength, uint8_t* pData, SBR& sbr, id_type id)
{
	assert(m_children < m_capacity);

	m_pDataLength[m_children] = dataLength;
	m_pData[m_children] = pData;
	m_ptrSBR[m_children] = m_pTree->m_SBRPool.acquire();
	*(m_ptrSBR[m_children]) = sbr;
	m_pIdentifier[m_children] = id;

	m_totalDataLength += dataLength;
	++m_children;

	m_nodeSBR.combineSBR(sbr);
}

void Node::insertEntry(uint32_t dataLength, uint8_t *pData,SBR& sbr, MBC &mbc, SpatialIndex::id_type id) {
    assert(m_children < m_capacity);
    m_pDataLength[m_children] = dataLength;
    m_pData[m_children] = pData;
    m_ptrMBC[m_children] = m_pTree->m_mbcPool.acquire();
    *(m_ptrMBC[m_children]) = mbc;

    m_pIdentifier[m_children] = id;

    m_totalDataLength += dataLength;
    ++m_children;

    m_nodeSBR.combineSBR(sbr);
}

void Node::deleteEntry(uint32_t index)
{
	assert(index >= 0 && index < m_children);

	// cache it, since I might need it for "touches" later.
	SBRPtr ptrR = m_ptrSBR[index];

	m_totalDataLength -= m_pDataLength[index];
	if (m_pData[index] != 0) delete[] m_pData[index];

	if (m_children > 1 && index != m_children - 1)
	{
		m_pDataLength[index] = m_pDataLength[m_children - 1];
		m_pData[index] = m_pData[m_children - 1];
		m_ptrSBR[index] = m_ptrSBR[m_children - 1];
		m_pIdentifier[index] = m_pIdentifier[m_children - 1];
	}

	--m_children;

	// WARNING: index has now changed. Do not use it below here.

	if (m_children == 0)
	{
		m_nodeSBR = m_pTree->m_infiniteSBR;
	}
	else if (m_pTree->m_bTightSBRs && m_nodeSBR.touchesSBR(*ptrR))
	{
		for (uint32_t cDim = 0; cDim < m_nodeSBR.m_dimension; ++cDim)
		{
			m_nodeSBR.m_pLow[cDim] = std::numeric_limits<double>::max();
			m_nodeSBR.m_pHigh[cDim] = -std::numeric_limits<double>::max();

			for (uint32_t u32Child = 0; u32Child < m_children; ++u32Child)
			{
				m_nodeSBR.m_pLow[cDim] = std::min(m_nodeSBR.m_pLow[cDim], m_ptrSBR[u32Child]->m_pLow[cDim]);
				m_nodeSBR.m_pHigh[cDim] = std::max(m_nodeSBR.m_pHigh[cDim], m_ptrSBR[u32Child]->m_pHigh[cDim]);
			}
		}
	}
}

bool Node::insertData(uint32_t dataLength, uint8_t* pData, SBR& sbr, id_type id, std::stack<id_type>& pathBuffer, uint8_t* overflowTable)
{
	if (m_children < m_capacity)
	{
		bool adjusted = false;

		// this has to happen before insertEntry modifies m_nodeSBR.
		bool b = m_nodeSBR.containsSBR(sbr);

		insertEntry(dataLength, pData, sbr, id);
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
	else if (m_pTree->m_treeVariant == RV_RSTAR && (! pathBuffer.empty()) && overflowTable[m_level] == 0)
	{
		overflowTable[m_level] = 1;

		std::vector<uint32_t> vReinsert, vKeep;
		reinsertData(dataLength, pData, sbr, id, vReinsert, vKeep);

		uint32_t lReinsert = static_cast<uint32_t>(vReinsert.size());
		uint32_t lKeep = static_cast<uint32_t>(vKeep.size());

		uint8_t** reinsertdata = 0;
		SBRPtr* reinsertsbr = 0;
		id_type* reinsertid = 0;
		uint32_t* reinsertlen = 0;
		uint8_t** keepdata = 0;
		SBRPtr* keepsbr = 0;
		id_type* keepid = 0;
		uint32_t* keeplen = 0;

		try
		{
			reinsertdata = new uint8_t*[lReinsert];
			reinsertsbr = new SBRPtr[lReinsert];
			reinsertid = new id_type[lReinsert];
			reinsertlen = new uint32_t[lReinsert];

			keepdata = new uint8_t*[m_capacity + 1];
			keepsbr = new SBRPtr[m_capacity + 1];
			keepid = new id_type[m_capacity + 1];
			keeplen = new uint32_t[m_capacity + 1];
		}
		catch (...)
		{
			delete[] reinsertdata;
			delete[] reinsertsbr;
			delete[] reinsertid;
			delete[] reinsertlen;
			delete[] keepdata;
			delete[] keepsbr;
			delete[] keepid;
			delete[] keeplen;
			throw;
		}

		uint32_t cIndex;

		for (cIndex = 0; cIndex < lReinsert; ++cIndex)
		{
			reinsertlen[cIndex] = m_pDataLength[vReinsert[cIndex]];
			reinsertdata[cIndex] = m_pData[vReinsert[cIndex]];
			reinsertsbr[cIndex] = m_ptrSBR[vReinsert[cIndex]];
			reinsertid[cIndex] = m_pIdentifier[vReinsert[cIndex]];
		}

		for (cIndex = 0; cIndex < lKeep; ++cIndex)
		{
			keeplen[cIndex] = m_pDataLength[vKeep[cIndex]];
			keepdata[cIndex] = m_pData[vKeep[cIndex]];
			keepsbr[cIndex] = m_ptrSBR[vKeep[cIndex]];
			keepid[cIndex] = m_pIdentifier[vKeep[cIndex]];
		}

		delete[] m_pDataLength;
		delete[] m_pData;
		delete[] m_ptrSBR;
		delete[] m_pIdentifier;

		m_pDataLength = keeplen;
		m_pData = keepdata;
		m_ptrSBR = keepsbr;
		m_pIdentifier = keepid;
		m_children = lKeep;
		m_totalDataLength = 0;

		for (uint32_t u32Child = 0; u32Child < m_children; ++u32Child) m_totalDataLength += m_pDataLength[u32Child];

		for (uint32_t cDim = 0; cDim < m_nodeSBR.m_dimension; ++cDim)
		{
			m_nodeSBR.m_pLow[cDim] = std::numeric_limits<double>::max();
			m_nodeSBR.m_pHigh[cDim] = -std::numeric_limits<double>::max();

			for (uint32_t u32Child = 0; u32Child < m_children; ++u32Child)
			{
				m_nodeSBR.m_pLow[cDim] = std::min(m_nodeSBR.m_pLow[cDim], m_ptrSBR[u32Child]->m_pLow[cDim]);
				m_nodeSBR.m_pHigh[cDim] = std::max(m_nodeSBR.m_pHigh[cDim], m_ptrSBR[u32Child]->m_pHigh[cDim]);
			}
		}

		m_pTree->writeNode(this);

		// Divertion from R*-Tree algorithm here. First adjust
		// the path to the root, then start reinserts, to avoid complicated handling
		// of changes to the same node from multiple insertions.
		id_type cParent = pathBuffer.top(); pathBuffer.pop();
		NodePtr ptrN = m_pTree->readNode(cParent);
		Index* p = static_cast<Index*>(ptrN.get());
		p->adjustTree(this, pathBuffer);

		for (cIndex = 0; cIndex < lReinsert; ++cIndex)
		{
			m_pTree->insertData_impl(
				reinsertlen[cIndex], reinsertdata[cIndex],
				*(reinsertsbr[cIndex]), reinsertid[cIndex],
				m_level, overflowTable);
		}

		delete[] reinsertdata;
		delete[] reinsertsbr;
		delete[] reinsertid;
		delete[] reinsertlen;

		return true;
	}
	else
	{
		NodePtr n;
		NodePtr nn;
		split(dataLength, pData, sbr, id, n, nn);

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
				ptrR->m_nodeSBR = m_pTree->m_infiniteSBR;
			}

			ptrR->insertEntry(0, 0, n->m_nodeSBR, n->m_identifier);
			ptrR->insertEntry(0, 0, nn->m_nodeSBR, nn->m_identifier);

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

void Node::reinsertData(uint32_t dataLength, uint8_t* pData, SBR& sbr, id_type id, std::vector<uint32_t>& reinsert, std::vector<uint32_t>& keep)
{
	ReinsertEntry** v = new ReinsertEntry*[m_capacity + 1];

	m_pDataLength[m_children] = dataLength;
	m_pData[m_children] = pData;
	m_ptrSBR[m_children] = m_pTree->m_SBRPool.acquire();
	*(m_ptrSBR[m_children]) = sbr;
	m_pIdentifier[m_children] = id;

	PointPtr nc = m_pTree->m_pointPool.acquire();
	m_nodeSBR.getCenter(*nc);
	PointPtr c = m_pTree->m_pointPool.acquire();

	for (uint32_t u32Child = 0; u32Child < m_capacity + 1; ++u32Child)
	{
		try
		{
			v[u32Child] = new ReinsertEntry(u32Child, 0.0);
		}
		catch (...)
		{
			for (uint32_t i = 0; i < u32Child; ++i) delete v[i];
			delete[] v;
			throw;
		}

		m_ptrSBR[u32Child]->getCenter(*c);

		// calculate relative distance of every entry from the node SBR (ignore square root.)
		for (uint32_t cDim = 0; cDim < m_nodeSBR.m_dimension; ++cDim)
		{
			double d = nc->m_pCoords[cDim] - c->m_pCoords[cDim];
			v[u32Child]->m_dist += d * d;
		}
	}

	// sort by increasing order of distances.
	::qsort(v, m_capacity + 1, sizeof(ReinsertEntry*), ReinsertEntry::compareReinsertEntry);

	uint32_t cReinsert = static_cast<uint32_t>(std::floor((m_capacity + 1) * m_pTree->m_sbreinsertFactor));

	uint32_t cCount;

	for (cCount = 0; cCount < cReinsert; ++cCount)
	{
		reinsert.emplace_back(v[cCount]->m_index);
		delete v[cCount];
	}

	for (cCount = cReinsert; cCount < m_capacity + 1; ++cCount)
	{
		keep.emplace_back(v[cCount]->m_index);
		delete v[cCount];
	}

	delete[] v;
}

void Node::SBRTreeSplit(uint32_t dataLength, uint8_t* pData, SBR& sbr, id_type id, std::vector<uint32_t>& group1, std::vector<uint32_t>& group2)
{
	uint32_t u32Child;
	uint32_t minimumLoad = static_cast<uint32_t>(std::floor(m_capacity * m_pTree->m_fillFactor));

	// use this mask array for marking visited entries.
	uint8_t* mask = new uint8_t[m_capacity + 1];
	memset(mask, 0, m_capacity + 1);

	// insert new data in the node for easier manipulation. Data arrays are always
	// by one larger than node capacity.
	m_pDataLength[m_capacity] = dataLength;
	m_pData[m_capacity] = pData;
	m_ptrSBR[m_capacity] = m_pTree->m_SBRPool.acquire();
	*(m_ptrSBR[m_capacity]) = sbr;
	m_pIdentifier[m_capacity] = id;
	// m_totalDataLength does not need to be increased here.

	// initialize each group with the seed entries.
	uint32_t seed1, seed2;
	pickSeeds(seed1, seed2);

	group1.emplace_back(seed1);
	group2.emplace_back(seed2);

	mask[seed1] = 1;
	mask[seed2] = 1;

	// find SBR of each group.
	SBRPtr sbr1 = m_pTree->m_SBRPool.acquire();
	*sbr1 = *(m_ptrSBR[seed1]);
	SBRPtr sbr2 = m_pTree->m_SBRPool.acquire();
	*sbr2 = *(m_ptrSBR[seed2]);

	// count how many entries are left unchecked (exclude the seeds here.)
	uint32_t cRemaining = m_capacity + 1 - 2;

	while (cRemaining > 0)
	{
		if (minimumLoad - group1.size() == cRemaining)
		{
			// all remaining entries must be assigned to group1 to comply with minimun load requirement.
			for (u32Child = 0; u32Child < m_capacity + 1; ++u32Child)
			{
				if (mask[u32Child] == 0)
				{
					group1.emplace_back(u32Child);
					mask[u32Child] = 1;
					--cRemaining;
				}
			}
		}
		else if (minimumLoad - group2.size() == cRemaining)
		{
			// all remaining entries must be assigned to group2 to comply with minimun load requirement.
			for (u32Child = 0; u32Child < m_capacity + 1; ++u32Child)
			{
				if (mask[u32Child] == 0)
				{
					group2.emplace_back(u32Child);
					mask[u32Child] = 1;
					--cRemaining;
				}
			}
		}
		else
		{
			// For all remaining entries compute the difference of the cost of grouping an
			// entry in either group. When done, choose the entry that yielded the maximum
			// difference. In case of linear split, select any entry (e.g. the first one.)
			uint32_t sel;
			double md1 = 0.0, md2 = 0.0;
			double m = -std::numeric_limits<double>::max();
			double d1, d2, d;
			double a1 = sbr1->getArea();
			double a2 = sbr2->getArea();

			SBRPtr a = m_pTree->m_SBRPool.acquire();
			SBRPtr b = m_pTree->m_SBRPool.acquire();

			for (u32Child = 0; u32Child < m_capacity + 1; ++u32Child)
			{
				if (mask[u32Child] == 0)
				{
					sbr1->getCombinedSBR(*a, *(m_ptrSBR[u32Child]));
					d1 = a->getArea() - a1;
					sbr2->getCombinedSBR(*b, *(m_ptrSBR[u32Child]));
					d2 = b->getArea() - a2;
					d = std::abs(d1 - d2);

					if (d > m)
					{
						m = d;
						md1 = d1; md2 = d2;
						sel = u32Child;
						if (m_pTree->m_treeVariant== RV_LINEAR || m_pTree->m_treeVariant == RV_RSTAR) break;
					}
				}
			}

			// determine the group where we should add the new entry.
			int32_t group = -1;

			if (md1 < md2)
			{
				group1.emplace_back(sel);
				group = 1;
			}
			else if (md2 < md1)
			{
				group2.emplace_back(sel);
				group = 2;
			}
			else if (a1 < a2)
			{
				group1.emplace_back(sel);
				group = 1;
			}
			else if (a2 < a1)
			{
				group2.emplace_back(sel);
				group = 2;
			}
			else if (group1.size() < group2.size())
			{
				group1.emplace_back(sel);
				group = 1;
			}
			else if (group2.size() < group1.size())
			{
				group2.emplace_back(sel);
				group = 2;
			}
			else
			{
				group1.emplace_back(sel);
				group = 1;
			}
			mask[sel] = 1;
			--cRemaining;
			if (group == 1)
			{
				sbr1->combineSBR(*(m_ptrSBR[sel]));
			}
			else
			{
				sbr2->combineSBR(*(m_ptrSBR[sel]));
			}
		}
	}

	delete[] mask;
}

void Node::rstarSplit(uint32_t dataLength, uint8_t* pData, SBR& sbr, id_type id, std::vector<uint32_t>& group1, std::vector<uint32_t>& group2)
{
	RstarSplitEntry** dataLow = 0;
	RstarSplitEntry** dataHigh = 0;

	try
	{
		dataLow = new RstarSplitEntry*[m_capacity + 1];
		dataHigh = new RstarSplitEntry*[m_capacity + 1];
	}
	catch (...)
	{
		delete[] dataLow;
		throw;
	}

	m_pDataLength[m_capacity] = dataLength;
	m_pData[m_capacity] = pData;
	m_ptrSBR[m_capacity] = m_pTree->m_SBRPool.acquire();
	*(m_ptrSBR[m_capacity]) = sbr;
	m_pIdentifier[m_capacity] = id;
	// m_totalDataLength does not need to be increased here.

	uint32_t nodeSPF = static_cast<uint32_t>(
		std::floor((m_capacity + 1) * m_pTree->m_splitDistributionFactor));
	uint32_t splitDistribution = (m_capacity + 1) - (2 * nodeSPF) + 2;

	uint32_t u32Child = 0, cDim, cIndex;

	for (u32Child = 0; u32Child <= m_capacity; ++u32Child)
	{
		try
		{
			dataLow[u32Child] = new RstarSplitEntry(m_ptrSBR[u32Child].get(), u32Child, 0);
		}
		catch (...)
		{
			for (uint32_t i = 0; i < u32Child; ++i) delete dataLow[i];
			delete[] dataLow;
			delete[] dataHigh;
			throw;
		}

		dataHigh[u32Child] = dataLow[u32Child];
	}

	double minimumMargin = std::numeric_limits<double>::max();
	uint32_t splitAxis = std::numeric_limits<uint32_t>::max();
	uint32_t sortOrder = std::numeric_limits<uint32_t>::max();

	// chooseSplitAxis.
	for (cDim = 0; cDim < m_pTree->m_dimension; ++cDim)
	{
		::qsort(dataLow, m_capacity + 1, sizeof(RstarSplitEntry*), RstarSplitEntry::compareLow);
		::qsort(dataHigh, m_capacity + 1, sizeof(RstarSplitEntry*), RstarSplitEntry::compareHigh);

		// calculate sum of margins and overlap for all distributions.
		double marginl = 0.0;
		double marginh = 0.0;

		SBR bbl1, bbl2, bbh1, bbh2;

		for (u32Child = 1; u32Child <= splitDistribution; ++u32Child)
		{
			uint32_t l = nodeSPF - 1 + u32Child;

			bbl1 = *(dataLow[0]->m_pSBR);
			bbh1 = *(dataHigh[0]->m_pSBR);

			for (cIndex = 1; cIndex < l; ++cIndex)
			{
				bbl1.combineSBR(*(dataLow[cIndex]->m_pSBR));
				bbh1.combineSBR(*(dataHigh[cIndex]->m_pSBR));
			}

			bbl2 = *(dataLow[l]->m_pSBR);
			bbh2 = *(dataHigh[l]->m_pSBR);

			for (cIndex = l + 1; cIndex <= m_capacity; ++cIndex)
			{
				bbl2.combineSBR(*(dataLow[cIndex]->m_pSBR));
				bbh2.combineSBR(*(dataHigh[cIndex]->m_pSBR));
			}

			marginl += bbl1.getMargin() + bbl2.getMargin();
			marginh += bbh1.getMargin() + bbh2.getMargin();
		} // for (u32Child)

		double margin = std::min(marginl, marginh);

		// keep minimum margin as split axis.
		if (margin < minimumMargin)
		{
			minimumMargin = margin;
			splitAxis = cDim;
			sortOrder = (marginl < marginh) ? 0 : 1;
		}

		// increase the dimension according to which the data entries should be sorted.
		for (u32Child = 0; u32Child <= m_capacity; ++u32Child)
		{
			dataLow[u32Child]->m_sortDim = cDim + 1;
		}
	} // for (cDim)

	for (u32Child = 0; u32Child <= m_capacity; ++u32Child)
	{
		dataLow[u32Child]->m_sortDim = splitAxis;
	}

	::qsort(dataLow, m_capacity + 1, sizeof(RstarSplitEntry*), (sortOrder == 0) ? RstarSplitEntry::compareLow : RstarSplitEntry::compareHigh);

	double ma = std::numeric_limits<double>::max();
	double mo = std::numeric_limits<double>::max();
	uint32_t splitPoint = std::numeric_limits<uint32_t>::max();

	SBR bb1, bb2;

	for (u32Child = 1; u32Child <= splitDistribution; ++u32Child)
	{
		uint32_t l = nodeSPF - 1 + u32Child;

		bb1 = *(dataLow[0]->m_pSBR);

		for (cIndex = 1; cIndex < l; ++cIndex)
		{
			bb1.combineSBR(*(dataLow[cIndex]->m_pSBR));
		}

		bb2 = *(dataLow[l]->m_pSBR);

		for (cIndex = l + 1; cIndex <= m_capacity; ++cIndex)
		{
			bb2.combineSBR(*(dataLow[cIndex]->m_pSBR));
		}

		double o = bb1.getIntersectingArea(bb2);

		if (o < mo)
		{
			splitPoint = u32Child;
			mo = o;
			ma = bb1.getArea() + bb2.getArea();
		}
		else if (o == mo)
		{
			double a = bb1.getArea() + bb2.getArea();

			if (a < ma)
			{
				splitPoint = u32Child;
				ma = a;
			}
		}
	} // for (u32Child)

	uint32_t l1 = nodeSPF - 1 + splitPoint;

	for (cIndex = 0; cIndex < l1; ++cIndex)
	{
		group1.emplace_back(dataLow[cIndex]->m_index);
		delete dataLow[cIndex];
	}

	for (cIndex = l1; cIndex <= m_capacity; ++cIndex)
	{
		group2.emplace_back(dataLow[cIndex]->m_index);
		delete dataLow[cIndex];
	}

	delete[] dataLow;
	delete[] dataHigh;
}

void Node::pickSeeds(uint32_t& index1, uint32_t& index2)
{
	double separation = -std::numeric_limits<double>::max();
	double inefficiency = -std::numeric_limits<double>::max();
	uint32_t cDim, u32Child, cIndex;

	switch (m_pTree->m_treeVariant)
	{
		case RV_LINEAR:
		case RV_RSTAR:
			for (cDim = 0; cDim < m_pTree->m_dimension; ++cDim)
			{
				double leastLower = m_ptrSBR[0]->m_pLow[cDim];
				double greatestUpper = m_ptrSBR[0]->m_pHigh[cDim];
				uint32_t greatestLower = 0;
				uint32_t leastUpper = 0;
				double width;

				for (u32Child = 1; u32Child <= m_capacity; ++u32Child)
				{
					if (m_ptrSBR[u32Child]->m_pLow[cDim] > m_ptrSBR[greatestLower]->m_pLow[cDim]) greatestLower = u32Child;
					if (m_ptrSBR[u32Child]->m_pHigh[cDim] < m_ptrSBR[leastUpper]->m_pHigh[cDim]) leastUpper = u32Child;

					leastLower = std::min(m_ptrSBR[u32Child]->m_pLow[cDim], leastLower);
					greatestUpper = std::max(m_ptrSBR[u32Child]->m_pHigh[cDim], greatestUpper);
				}

				width = greatestUpper - leastLower;
				if (width <= 0) width = 1;

				double f = (m_ptrSBR[greatestLower]->m_pLow[cDim] - m_ptrSBR[leastUpper]->m_pHigh[cDim]) / width;

				if (f > separation)
				{
					index1 = leastUpper;
					index2 = greatestLower;
					separation = f;
				}
			}  // for (cDim)

			if (index1 == index2)
			{
				if (index2 == 0) ++index2;
				else --index2;
			}

			break;
		case RV_QUADRATIC:
			// for each pair of SBRs (account for overflow SBR too!)
			for (u32Child = 0; u32Child < m_capacity; ++u32Child)
			{
				double a = m_ptrSBR[u32Child]->getArea();

				for (cIndex = u32Child + 1; cIndex <= m_capacity; ++cIndex)
				{
					// get the combined SBR of those two entries.
					SBR r;
					m_ptrSBR[u32Child]->getCombinedSBR(r, *(m_ptrSBR[cIndex]));

					// find the inefficiency of grouping these entries together.
					double d = r.getArea() - a - m_ptrSBR[cIndex]->getArea();

					if (d > inefficiency)
					{
						inefficiency = d;
						index1 = u32Child;
						index2 = cIndex;
					}
				}  // for (cIndex)
			} // for (u32Child)

			break;
		default:
			throw Tools::NotSupportedException("Node::pickSeeds: Tree variant not supported.");
	}
}

void Node::condenseTree(std::stack<NodePtr>& toReinsert, std::stack<id_type>& pathBuffer, NodePtr& ptrThis)
{
	uint32_t minimumLoad = static_cast<uint32_t>(std::floor(m_capacity * m_pTree->m_fillFactor));

	if (pathBuffer.empty())
	{
		// eliminate root if it has only one child.
		if (m_level != 0 && m_children == 1)
		{
			NodePtr ptrN = m_pTree->readNode(m_pIdentifier[0]);
			m_pTree->deleteNode(ptrN.get());
			ptrN->m_identifier = m_pTree->m_rootID;
			m_pTree->writeNode(ptrN.get());

			m_pTree->m_stats.m_nodesInLevel.pop_back();
			m_pTree->m_stats.m_u32TreeHeight -= 1;
			// HACK: pending deleteNode for deleted child will decrease nodesInLevel, later on.
			m_pTree->m_stats.m_nodesInLevel[m_pTree->m_stats.m_u32TreeHeight - 1] = 2;
		}
	}
	else
	{
		id_type cParent = pathBuffer.top(); pathBuffer.pop();
		NodePtr ptrParent = m_pTree->readNode(cParent);
		Index* p = static_cast<Index*>(ptrParent.get());

		// find the entry in the parent, that points to this node.
		uint32_t child;

		for (child = 0; child != p->m_children; ++child)
		{
			if (p->m_pIdentifier[child] == m_identifier) break;
		}

		if (m_children < minimumLoad)
		{
			// used space less than the minimum
			// 1. eliminate node entry from the parent. deleteEntry will fix the parent's SBR.
			p->deleteEntry(child);
			// 2. add this node to the stack in order to reinsert its entries.
			toReinsert.push(ptrThis);
		}
		else
		{
			// adjust the entry in 'p' to contain the new bounding SBR of this node.
			*(p->m_ptrSBR[child]) = m_nodeSBR;

			// global recalculation necessary since the SBR can only shrink in size,
			// due to data removal.
			if (m_pTree->m_bTightSBRs)
			{
				for (uint32_t cDim = 0; cDim < p->m_nodeSBR.m_dimension; ++cDim)
				{
					p->m_nodeSBR.m_pLow[cDim] = std::numeric_limits<double>::max();
					p->m_nodeSBR.m_pHigh[cDim] = -std::numeric_limits<double>::max();

					for (uint32_t u32Child = 0; u32Child < p->m_children; ++u32Child)
					{
						p->m_nodeSBR.m_pLow[cDim] = std::min(p->m_nodeSBR.m_pLow[cDim], p->m_ptrSBR[u32Child]->m_pLow[cDim]);
						p->m_nodeSBR.m_pHigh[cDim] = std::max(p->m_nodeSBR.m_pHigh[cDim], p->m_ptrSBR[u32Child]->m_pHigh[cDim]);
					}
				}
			}
		}

		// write parent node back to storage.
		m_pTree->writeNode(p);

		p->condenseTree(toReinsert, pathBuffer, ptrParent);
	}
}
