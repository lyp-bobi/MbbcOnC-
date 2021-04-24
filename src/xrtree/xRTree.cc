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
#include "Node.h"
#include "Leaf.h"
#include "Index.h"
#include "BulkLoader.h"
#include "xRTree.h"


using namespace SpatialIndex::xRTreeNsp;
using namespace SpatialIndex;

SpatialIndex::xRTreeNsp::xRTree::xRTree(IStorageManager& sm, Tools::PropertySet& ps) :
	m_pStorageManager(&sm),
	m_rootID(StorageManager::NewPage),
	m_headerID(StorageManager::NewPage),
	m_treeVariant(RV_RSTAR),
	m_fillFactor(0.7),
	m_indexCapacity(100),
	m_leafCapacity(100),
	m_nearMinimumOverlapFactor(32),
	m_splitDistributionFactor(0.4),
	m_reinsertFactor(0.3),
	m_dimension(2),
	m_bTightMBRs(true),
	m_xPointPool(500),
	m_xMBRPool(1000),
    m_xSBBPool(1000),
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
		else throw Tools::IllegalArgumentException("xRTree: Property IndexIdentifier must be Tools::VT_LONGLONG");

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

SpatialIndex::xRTreeNsp::xRTree::~xRTree()
{
#ifdef HAVE_PTHREAD_H
	pthread_mutex_destroy(&m_lock);
#endif

	storeHeader();
}


void SpatialIndex::xRTreeNsp::xRTree::getIndexProperties(Tools::PropertySet& out) const
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

	// R-tree variant
	var.m_varType = Tools::VT_LONG;
	var.m_val.lVal = m_treeVariant;
	out.setProperty("TreeVariant", var);

	// fill factor
	var.m_varType = Tools::VT_DOUBLE;
	var.m_val.dblVal = m_fillFactor;
	out.setProperty("FillFactor", var);

	// near minimum overlap factor
	var.m_varType = Tools::VT_ULONG;
	var.m_val.ulVal = m_nearMinimumOverlapFactor;
	out.setProperty("NearMinimumOverlapFactor", var);

	// split distribution factor
	var.m_varType = Tools::VT_DOUBLE;
	var.m_val.dblVal = m_splitDistributionFactor;
	out.setProperty("SplitDistributionFactor", var);

	// reinsert factor
	var.m_varType = Tools::VT_DOUBLE;
	var.m_val.dblVal = m_reinsertFactor;
	out.setProperty("ReinsertFactor", var);

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

	// xMBR pool capacity
	var.m_varType = Tools::VT_ULONG;
	var.m_val.ulVal = m_xMBRPool.getCapacity();
	out.setProperty("xMBRPoolCapacity", var);

	// xPoint pool capacity
	var.m_varType = Tools::VT_ULONG;
	var.m_val.ulVal = m_xPointPool.getCapacity();
	out.setProperty("xPointPoolCapacity", var);
}

void SpatialIndex::xRTreeNsp::xRTree::getStatistics(IStatistics** out) const
{
	*out = new Statistics(m_stats);
}

void SpatialIndex::xRTreeNsp::xRTree::initNew(Tools::PropertySet& ps)
{
	Tools::Variant var;

	// tree variant
	var = ps.getProperty("TreeVariant");
	if (var.m_varType != Tools::VT_EMPTY)
	{
		if (
			var.m_varType != Tools::VT_LONG ||
			(var.m_val.lVal != RV_LINEAR &&
			var.m_val.lVal != RV_QUADRATIC &&
			var.m_val.lVal != RV_RSTAR))
			throw Tools::IllegalArgumentException("initNew: Property TreeVariant must be Tools::VT_LONG and of xRTreeVariant type");

		m_treeVariant = static_cast<xRTreeVariant>(var.m_val.lVal);
	}

    // bt
    var = ps.getProperty("Bt");
    if (var.m_varType != Tools::VT_EMPTY)
    {
        if (
                var.m_varType != Tools::VT_DOUBLE ||
                var.m_val.dblVal <= 0.0)
            throw Tools::IllegalArgumentException("initNew: Property Bt must be Tools::VT_DOUBLE and in (0.0, inf)");

        m_bt = var.m_val.dblVal;
    }

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

        if (((m_treeVariant == RV_LINEAR || m_treeVariant == RV_QUADRATIC) && var.m_val.dblVal > 0.5))
            throw Tools::IllegalArgumentException(  "initNew: Property FillFactor must be in range "
                                                    "(0.0, 0.5) for LINEAR or QUADRATIC index types");
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

	// near minimum overlap factor
	var = ps.getProperty("NearMinimumOverlapFactor");
	if (var.m_varType != Tools::VT_EMPTY)
	{
		if (
			var.m_varType != Tools::VT_ULONG ||
			var.m_val.ulVal < 1 ||
			var.m_val.ulVal > m_indexCapacity ||
			var.m_val.ulVal > m_leafCapacity)
			throw Tools::IllegalArgumentException("initNew: Property NearMinimumOverlapFactor must be Tools::VT_ULONG and less than both index and leaf capacities");

		m_nearMinimumOverlapFactor = var.m_val.ulVal;
	}

	// split distribution factor
	var = ps.getProperty("SplitDistributionFactor");
	if (var.m_varType != Tools::VT_EMPTY)
	{
		if (
			var.m_varType != Tools::VT_DOUBLE ||
			var.m_val.dblVal <= 0.0 ||
			var.m_val.dblVal >= 1.0)
			throw Tools::IllegalArgumentException("initNew: Property SplitDistributionFactor must be Tools::VT_DOUBLE and in (0.0, 1.0)");

		m_splitDistributionFactor = var.m_val.dblVal;
	}

	// reinsert factor
	var = ps.getProperty("ReinsertFactor");
	if (var.m_varType != Tools::VT_EMPTY)
	{
		if (
			var.m_varType != Tools::VT_DOUBLE ||
			var.m_val.dblVal <= 0.0 ||
			var.m_val.dblVal >= 1.0)
			throw Tools::IllegalArgumentException("initNew: Property ReinsertFactor must be Tools::VT_DOUBLE and in (0.0, 1.0)");

		m_reinsertFactor = var.m_val.dblVal;
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

	// xMBR pool capacity
	var = ps.getProperty("xMBRPoolCapacity");
	if (var.m_varType != Tools::VT_EMPTY)
	{
		if (var.m_varType != Tools::VT_ULONG)
			throw Tools::IllegalArgumentException("initNew: Property xMBRPoolCapacity must be Tools::VT_ULONG");

		m_xMBRPool.setCapacity(var.m_val.ulVal);
	}

	// xPoint pool capacity
	var = ps.getProperty("xPointPoolCapacity");
	if (var.m_varType != Tools::VT_EMPTY)
	{
		if (var.m_varType != Tools::VT_ULONG)
			throw Tools::IllegalArgumentException("initNew: Property xPointPoolCapacity must be Tools::VT_ULONG");

		m_xPointPool.setCapacity(var.m_val.ulVal);
	}

	m_stats.m_u32TreeHeight = 1;
	m_stats.m_nodesInLevel.emplace_back(0);

	Leaf root(this, -1);
	m_rootID = writeNode(static_cast<Node *>(&root));

	storeHeader();
}

void SpatialIndex::xRTreeNsp::xRTree::initOld(Tools::PropertySet& ps)
{
	loadHeader();

	// only some of the properties may be changed.
	// the rest are just ignored.

	Tools::Variant var;

	// tree variant
	var = ps.getProperty("TreeVariant");
	if (var.m_varType != Tools::VT_EMPTY)
	{
		if (
			var.m_varType != Tools::VT_LONG ||
			(var.m_val.lVal != RV_LINEAR &&
			 var.m_val.lVal != RV_QUADRATIC &&
			 var.m_val.lVal != RV_RSTAR))
			throw Tools::IllegalArgumentException("initOld: Property TreeVariant must be Tools::VT_LONG and of xRTreeVariant type");

		m_treeVariant = static_cast<xRTreeVariant>(var.m_val.lVal);
	}

    // bt
    var = ps.getProperty("Bt");
    if (var.m_varType != Tools::VT_EMPTY)
    {
        if (
                var.m_varType != Tools::VT_DOUBLE ||
                var.m_val.dblVal <= 0.0)
            throw Tools::IllegalArgumentException("initNew: Property Bt must be Tools::VT_DOUBLE and in (0.0, inf)");

        m_bt = var.m_val.dblVal;
    }

	// near minimum overlap factor
	var = ps.getProperty("NearMinimumOverlapFactor");
	if (var.m_varType != Tools::VT_EMPTY)
	{
		if (
			var.m_varType != Tools::VT_ULONG ||
			var.m_val.ulVal < 1 ||
			var.m_val.ulVal > m_indexCapacity ||
			var.m_val.ulVal > m_leafCapacity)
			throw Tools::IllegalArgumentException("initOld: Property NearMinimumOverlapFactor must be Tools::VT_ULONG and less than both index and leaf capacities");

		m_nearMinimumOverlapFactor = var.m_val.ulVal;
	}

	// split distribution factor
	var = ps.getProperty("SplitDistributionFactor");
	if (var.m_varType != Tools::VT_EMPTY)
	{
		if (var.m_varType != Tools::VT_DOUBLE || var.m_val.dblVal <= 0.0 || var.m_val.dblVal >= 1.0)
			throw Tools::IllegalArgumentException("initOld: Property SplitDistributionFactor must be Tools::VT_DOUBLE and in (0.0, 1.0)");

		m_splitDistributionFactor = var.m_val.dblVal;
	}

	// reinsert factor
	var = ps.getProperty("ReinsertFactor");
	if (var.m_varType != Tools::VT_EMPTY)
	{
		if (var.m_varType != Tools::VT_DOUBLE || var.m_val.dblVal <= 0.0 || var.m_val.dblVal >= 1.0)
			throw Tools::IllegalArgumentException("initOld: Property ReinsertFactor must be Tools::VT_DOUBLE and in (0.0, 1.0)");

		m_reinsertFactor = var.m_val.dblVal;
	}

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

	// xMBR pool capacity
	var = ps.getProperty("xMBRPoolCapacity");
	if (var.m_varType != Tools::VT_EMPTY)
	{
		if (var.m_varType != Tools::VT_ULONG) throw Tools::IllegalArgumentException("initOld: Property xMBRPoolCapacity must be Tools::VT_ULONG");

		m_xMBRPool.setCapacity(var.m_val.ulVal);
	}

	// xPoint pool capacity
	var = ps.getProperty("xPointPoolCapacity");
	if (var.m_varType != Tools::VT_EMPTY)
	{
		if (var.m_varType != Tools::VT_ULONG) throw Tools::IllegalArgumentException("initOld: Property xPointPoolCapacity must be Tools::VT_ULONG");

		m_xPointPool.setCapacity(var.m_val.ulVal);
	}
}

void SpatialIndex::xRTreeNsp::xRTree::storeHeader()
{
    if(m_pStorageManager->m_isro) return;
	const uint32_t headerSize =
		sizeof(id_type) +						// m_rootID
		sizeof(double)+                         // m_bt
		sizeof(xRTreeVariant) +					// m_treeVariant
		sizeof(double) +						// m_fillFactor
		sizeof(uint32_t) +						// m_indexCapacity
		sizeof(uint32_t) +						// m_leafCapacity
		sizeof(uint32_t) +						// m_nearMinimumOverlapFactor
		sizeof(double) +						// m_splitDistributionFactor
		sizeof(double) +						// m_reinsertFactor
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
    memcpy(ptr, &m_bt, sizeof(double));
    ptr += sizeof(double);
	memcpy(ptr, &m_treeVariant, sizeof(xRTreeVariant));
	ptr += sizeof(xRTreeVariant);
	memcpy(ptr, &m_fillFactor, sizeof(double));
	ptr += sizeof(double);
	memcpy(ptr, &m_indexCapacity, sizeof(uint32_t));
	ptr += sizeof(uint32_t);
	memcpy(ptr, &m_leafCapacity, sizeof(uint32_t));
	ptr += sizeof(uint32_t);
	memcpy(ptr, &m_nearMinimumOverlapFactor, sizeof(uint32_t));
	ptr += sizeof(uint32_t);
	memcpy(ptr, &m_splitDistributionFactor, sizeof(double));
	ptr += sizeof(double);
	memcpy(ptr, &m_reinsertFactor, sizeof(double));
	ptr += sizeof(double);
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

void SpatialIndex::xRTreeNsp::xRTree::loadHeader()
{
	uint32_t headerSize;
	uint8_t* header = 0;
	m_pStorageManager->loadByteArray(m_headerID, headerSize, &header);

	uint8_t* ptr = header;

	memcpy(&m_rootID, ptr, sizeof(id_type));
	ptr += sizeof(id_type);
    memcpy(&m_bt, ptr, sizeof(double));
    ptr += sizeof(double);
	memcpy(&m_treeVariant, ptr, sizeof(xRTreeVariant));
	ptr += sizeof(xRTreeVariant);
	memcpy(&m_fillFactor, ptr, sizeof(double));
	ptr += sizeof(double);
	memcpy(&m_indexCapacity, ptr, sizeof(uint32_t));
	ptr += sizeof(uint32_t);
	memcpy(&m_leafCapacity, ptr, sizeof(uint32_t));
	ptr += sizeof(uint32_t);
	memcpy(&m_nearMinimumOverlapFactor, ptr, sizeof(uint32_t));
	ptr += sizeof(uint32_t);
	memcpy(&m_splitDistributionFactor, ptr, sizeof(double));
	ptr += sizeof(double);
	memcpy(&m_reinsertFactor, ptr, sizeof(double));
	ptr += sizeof(double);
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
    if(m_pStorageManager->m_needfree)
	    delete[] header;
}

SpatialIndex::id_type SpatialIndex::xRTreeNsp::xRTree::writeNode(Node* n)
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
    //test cost
//    std:cerr<<n->toString()<<endl;
	return page;
}

SpatialIndex::xRTreeNsp::NodePtr SpatialIndex::xRTreeNsp::xRTree::readNode(id_type page)
{
    m_ts->m_indexIO++;
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
		uint32_t level;
		memcpy(&level, buffer, sizeof(uint32_t));

		NodePtr n;

		if (level>0) n = m_indexPool.acquire();
		else if (level==0) n = m_leafPool.acquire();
		else
		    throw Tools::IllegalStateException("readNode: failed reading the correct node type information");

		if (n.get() == nullptr)
		{
			if (level>0) n = NodePtr(new Index(this, -1, 0), &m_indexPool);
			else if (level == 0) n = NodePtr(new Leaf(this, -1), &m_leafPool);
		}

//		n->m_pTree = this;
		n->m_identifier = page;

		n->loadFromByteArray(buffer);

		++(m_stats.m_u64Reads);
        if(m_pStorageManager->m_needfree)
		    delete[] buffer;
        //test code
//        std:cerr<<n->toString()<<endl;
		return n;
	}
	catch (...)
	{
        if(m_pStorageManager->m_needfree)
		    delete[] buffer;
		throw;
	}
}

void SpatialIndex::xRTreeNsp::xRTree::deleteNode(Node* n)
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

}



std::ostream& SpatialIndex::xRTreeNsp::operator<<(std::ostream& os, const xRTree& t)
{
	os	<< "Dimension: " << t.m_dimension << std::endl
		<< "Fill factor: " << t.m_fillFactor << std::endl
		<< "Index capacity: " << t.m_indexCapacity << std::endl
		<< "Leaf capacity: " << t.m_leafCapacity << std::endl
		<< "Tight MBRs: " << ((t.m_bTightMBRs) ? "enabled" : "disabled") << std::endl;

	if (t.m_treeVariant == RV_RSTAR)
	{
		os	<< "Near minimum overlap factor: " << t.m_nearMinimumOverlapFactor << std::endl
			<< "Reinsert factor: " << t.m_reinsertFactor << std::endl
			<< "Split distribution factor: " << t.m_splitDistributionFactor << std::endl;
	}

	if (t.m_stats.getNumberOfNodesInLevel(0) > 0)
		os	<< "Utilization: " << 100 * t.m_stats.getNumberOfData() / (t.m_stats.getNumberOfNodesInLevel(0) * t.m_leafCapacity) << "%" << std::endl
			<< t.m_stats;

	#ifndef NDEBUG
	os	<< "Leaf pool hits: " << t.m_leafPool.m_hits << std::endl
		<< "Leaf pool misses: " << t.m_leafPool.m_misses << std::endl
		<< "Index pool hits: " << t.m_indexPool.m_hits << std::endl
		<< "Index pool misses: " << t.m_indexPool.m_misses << std::endl
		<< "xMBR pool hits: " << t.m_xMBRPool.m_hits << std::endl
		<< "xMBR pool misses: " << t.m_xMBRPool.m_misses << std::endl
        << "xPoint pool hits: " << t.m_xPointPool.m_hits << std::endl
        << "xPoint pool misses: " << t.m_xPointPool.m_misses << std::endl;
#endif
    return os;
}



void SpatialIndex::xRTreeNsp::xRTree::insertData_impl(xMBR& mbr, id_type id)
{
    assert(mbr.getDimension() == m_dimension);

    std::stack<id_type> pathBuffer;
    uint8_t* overflowTable = 0;

    try
    {
        NodePtr root = readNode(m_rootID);

        overflowTable = new uint8_t[root->m_level];
        memset(overflowTable, 0, root->m_level);

        NodePtr l = root->chooseSubtree(mbr, 0, pathBuffer);
        if (l.get() == root.get())
        {
            assert(root.unique());
            root.relinquish();
        }
        l->insertData(mbr, id, pathBuffer, overflowTable);

        delete[] overflowTable;
        ++(m_stats.m_u64Data);
    }
    catch (...)
    {
        delete[] overflowTable;
        throw;
    }
}

void SpatialIndex::xRTreeNsp::xRTree::insertData_impl(xMBR& mbr, id_type id, uint32_t level, uint8_t* overflowTable)
{
    assert(mbr.getDimension() == m_dimension);

    std::stack<id_type> pathBuffer;
    NodePtr root = readNode(m_rootID);
    NodePtr n = root->chooseSubtree(mbr, level, pathBuffer);

    assert(n->m_level == level);

    if (n.get() == root.get())
    {
        assert(root.unique());
        root.relinquish();
    }
    n->insertData(mbr, id, pathBuffer, overflowTable);
}

bool SpatialIndex::xRTreeNsp::xRTree::deleteData_impl(const xMBR& mbr, id_type id)
{
    std::stack<id_type> pathBuffer;
    NodePtr root = readNode(m_rootID);
    NodePtr l = root->findLeaf(mbr, id, pathBuffer);
    if (l.get() == root.get())
    {
        assert(root.unique());
        root.relinquish();
    }

    if (l.get() != 0)
    {
        Leaf* pL = static_cast<Leaf*>(l.get());
        pL->deleteData(id, pathBuffer);
        --(m_stats.m_u64Data);
        return true;
    }

    return false;
}


void xRTree::intersectsWithQuery(const xCylinder &query, IVisitor &v) {
    const xCylinder *querycy= &query;
    std::set<id_type > results;
    std::multimap<id_type, xStoreEntry> pending;
    bool isSlice=(querycy->m_startTime==querycy->m_endTime);
    std::stack<NodePtr> st;
    NodePtr root = readNode(m_rootID);
    if (root->m_children > 0 && root->m_nodeMBR.intersectsShape(query)) st.push(root);
    while (!st.empty()) {
        NodePtr n = st.top();
        st.pop();
        if (n->m_level == 0) {
            v.visitNode(*n);
            for (uint32_t cChild = 0; cChild < n->m_children; ++cChild) {
                id_type id = n->m_se[cChild].m_id;
                if(results.count(id)>0) continue;
                int b;
                b = querycy->checkRel(*(n->m_ptrxSBB[cChild]));
                if (b>0) {
                    simpleData data = simpleData(n->m_se[cChild].m_id, 0);
                    if (b==2) {
                        m_stats.m_doubleExactQueryResults += 1;
                        results.insert(n->m_se[cChild].m_id);
                        pending.erase(n->m_se[cChild].m_id);
                        v.visitData(data);
//                        std::cerr<<*querycy << endl<< (n->m_ptrxSBB[cChild])->toString();
                        ++(m_stats.m_u64QueryResults);
                    } else {
                        if(m_bUsingSBBD) {
                            pending.insert(make_pair(id, n->m_se[cChild]));
                        }else{
                            xTrajectory partTraj;
                            m_ts->loadTraj(partTraj,n->m_se[cChild]);
                            if (partTraj.intersectsxCylinder(*querycy)) {
                                m_stats.m_doubleExactQueryResults += 1;
                                results.insert(id);
                                simpleData data = simpleData(id, 0);
                                v.visitData(data);
                                ++(m_stats.m_u64QueryResults);
                            }
                        }
                    }
                }
            }
        } else {
            v.visitNode(*n);
            for (uint32_t cChild = 0; cChild < n->m_children; ++cChild) {
                if (query.intersectsShape(*(n->m_ptrMBR[cChild]))) {
                    st.push(readNode(n->m_pIdentifier[cChild]));
                }
            }
        }
    }
    if(m_bUsingSBBD && !pending.empty()) {
        id_type previd=pending.begin()->first;
        uint32_t lower = pending.begin()->second.m_s, higher=pending.begin()->second.m_e;
        while (!pending.empty()) {
            auto iter = pending.begin();
            id_type id = iter->first;
            xStoreEntry storee = iter->second;
            pending.erase(iter);
            if(id!=previd||pending.empty()){ //check last
                xTrajectory partTraj;
                if(id ==previd){
                    lower = min(lower, storee.m_s);
                    higher = max(higher, storee.m_e);
                }
                m_ts->loadTraj(partTraj,xStoreEntry(previd,lower,higher));
                if (partTraj.intersectsxCylinder(*querycy)) {
                    m_stats.m_doubleExactQueryResults += 1;
                    results.insert(previd);
                    simpleData data(previd, 0);
                    v.visitData(data);
                    ++(m_stats.m_u64QueryResults);
                }
                if(id!=previd && pending.empty()){
                    m_ts->loadTraj(partTraj,storee);
                    if (partTraj.intersectsxCylinder(*querycy)) {
                        m_stats.m_doubleExactQueryResults += 1;
                        results.insert(id);
                        simpleData data(id, 0);
                        v.visitData(data);
                        ++(m_stats.m_u64QueryResults);
                    }
                }
                previd = id;
                lower = storee.m_s, higher=storee.m_e;
            }
            else{//merge
                lower = min(lower, storee.m_s);
                higher = max(higher, storee.m_e);
            }
        }
    }
}

void xRTree::findid(id_type qid) {
    std::set<id_type > results;
    std::multimap<id_type, xStoreEntry> pending;
    std::stack<NodePtr> st;
    NodePtr root = readNode(m_rootID);
    if (root->m_children > 0) st.push(root);
    while (!st.empty()) {
        NodePtr n = st.top();
        st.pop();
        if (n->m_level == 0) {
            for (uint32_t cChild = 0; cChild < n->m_children; ++cChild) {
                id_type id = n->m_se[cChild].m_id;
                if(id==qid){
                    cerr<<n->toString();
                }
            }
        } else {
            for (uint32_t cChild = 0; cChild < n->m_children; ++cChild) {
                st.push(readNode(n->m_pIdentifier[cChild]));
            }
        }
    }
}

void xRTree::nearestNeighborQuery(uint32_t k, const xTrajectory &query, IVisitor &v) {
    if(query.m_points.size()<=2) return;
    const xTrajectory *queryTraj= &query;
    xTrajectory simpleTraj;
    xTrajectory ssTraj;
    double delta=0, ssdelta= 0;
    if(m_bUsingSBBD&& m_bUsingSimp && m_bStoringLinks) {
        vector<xPoint> simpp;
        int segnum = std::floor((queryTraj->m_endTime() - queryTraj->m_startTime()) / (tjstat->bt));
        if (segnum * 2 < queryTraj->m_points.size()){
            simpleTraj = *queryTraj;
        }else {
            vector<vector<xPoint>> simpseg;
            simpseg = xTrajectory::simplifyWithRDPN(queryTraj->m_points,
                                                    std::min(segnum, int(std::sqrt(queryTraj->m_points.size()))));
            for (const auto &s:simpseg) {
                simpp.emplace_back(s.front());
            }
            simpp.emplace_back(simpseg.back().back());
            simpleTraj = xTrajectory(simpp);
        }
        delta = queryTraj->getMinimumDistance(simpleTraj);
        simpp.clear();
        simpp.emplace_back(queryTraj->m_points[0]);
        simpp.emplace_back(queryTraj->m_points[queryTraj->m_points.size()-1]);
        ssTraj=xTrajectory(simpp);
        ssdelta = queryTraj->getMinimumDistance(ssTraj);
    }else{
        simpleTraj=*queryTraj;
    }
#ifdef TJDEBUG
    cerr<<"query is "<< *queryTraj<<endl;
    cerr<<"simped query is"<<simpleTraj<<"\n with delta "<<delta<<endl;
#endif
    double knearest = 0.0;
    int iternum = 0;
    bool btopnode =false;
    /*SBB-Driven*/
    if(m_bUsingSBBD == true && m_bStoringLinks) {
        PartsStore ps(simpleTraj, delta, this,k);
        ps.push(new NNEntry(m_rootID, DISTE(0), 0));

        uint32_t count = 0;

        std::map<id_type, int> insertedTrajId;
        while (!ps.empty()) {
            iternum++;
            NNEntry *pFirst;
            if(btopnode) {
                pFirst = ps.nodetop();
                btopnode = false;
            }
            else pFirst=ps.top();

            // report all nearest neighbors with equal greatest distances.
            // (neighbors can be more than k, if many happen to have the same greatest distance).
            if (count >= k && pFirst->m_dist.opt > knearest) {
//            std::cerr<<"find minDist"<<knearest<<"\n";
                break;
            }
            switch (pFirst->m_type) {
                case 0: {//inner node
                    ps.pop(0);
                    NodePtr n = readNode(pFirst->m_id);
                    v.visitNode(*n);
                    for (uint32_t cChild = 0; cChild < n->m_children; ++cChild) {
                        double pd;
                        if(n->m_level>=2 && m_bUsingSimp){
                            pd = std::max(0.0, ssTraj.nodeDist(*(n->m_ptrMBR[cChild])) - ssdelta);
                        }else{
                            pd = std::max(0.0, simpleTraj.nodeDist(*(n->m_ptrMBR[cChild])) - delta);
                        }
                        if (pd < 1e300) {
                            if (n->m_level == 1)
                                ps.push(new NNEntry(n->m_pIdentifier[cChild], DISTE(pd), 1));
                            else
                                ps.push(new NNEntry(n->m_pIdentifier[cChild], DISTE(pd), 0));
                        }
                    }
                    delete pFirst;
//                n.relinquish();
                    break;
                }
                case 1: {//leaf node
                    ps.pop(1);
                    if (!ps.isLoaded(pFirst->m_id)) {
                        NodePtr n = readNode(pFirst->m_id);
                        m_ts->m_leaf1 += 1;
                        v.visitNode(*n);
                        ps.loadLeaf(*n);
#ifdef TJDEBUG
                        cerr<<iternum<<"\tleaf with "<<pFirst->m_id<<"\t"<<pFirst->m_dist.opt<<"\t"<<n->m_nodeMBR<<endl;
#endif
//                    n.relinquish();
                    }
                    delete pFirst;
                    break;
                }
                case 2: {//incomplete bounding
                    if(m_bStoringLinks) {
                        id_type missing;
                        missing = ps.getOneMissingPart(pFirst->m_id);
#ifdef TJDEBUG
                        cerr<<iternum<<"\tIB with "<<pFirst->m_dist.opt<<"\t"<<ps.explain(pFirst->m_id)<<endl;
#endif
                        NodePtr n = readNode(missing);
                        v.visitNode(*n);
                        ps.loadLeaf(*n, pFirst->m_dist.opt);
                        m_ts->m_leaf2 += 1;
                        break;
                    }else{
                        btopnode = true;
                        break;
                    }
                }
                case 3: {//complete bounding
                    ps.pop(3);
#ifdef TJDEBUG
                    cerr<<iternum<<"\tCB with "<<pFirst->m_dist.opt<<"\t"<<ps.explain(pFirst->m_id)<<endl;
#endif
                    if (pFirst->m_dist.pes < ps.top()->m_dist.opt){
                        // we judge by sbbs instead of subtraj, 1 for error caused by float point numbers
                        ++(m_stats.m_u64QueryResults);
                        ++count;
                        knearest = pFirst->m_dist.opt;
                        simpleData d(pFirst->m_id, pFirst->m_dist.opt);
                        v.visitData(d);
                        delete pFirst;
                        break;
                    }else{
                        xTrajectory traj = ps.getTraj(pFirst->m_id);
//                std::cerr<<"trajIO"<<m_ts->m_trajIO<<"\n";
//                std::cerr<<"getTraj"<<traj<<"\n";
//                xTrajectory traj = m_ts->getTrajByTime(pFirst->m_id, queryTraj->m_startTime(), queryTraj->m_endTime());
                        ps.push(new NNEntry(pFirst->m_id, DISTE(queryTraj->getMinimumDistance(traj)), 4));
                        delete pFirst;
                        break;
                    }

                }
                case 4: {//exact traj
                    ps.pop(4);
                    ++(m_stats.m_u64QueryResults);
                    ++count;
                    knearest = pFirst->m_dist.opt;
                    simpleData d(pFirst->m_id, pFirst->m_dist.opt);
                    v.visitData(d);
                    delete pFirst;
                    break;
                }
                default:
                    throw Tools::IllegalStateException("illegal NNEntry state");
            }
        }
        ps.clean();
    }
    else{/*BFMST*/
        PartsStoreBFMST ps(simpleTraj,0, this,k);
        string str = queryTraj->toString();
        ps.push(new NNEntry(m_rootID, DISTE(0), 0));

        uint32_t count = 0;
        double knearest = 0.0;
        int iternum=0;
        std::map<id_type ,int> insertedTrajId;

        while (! ps.empty()) {
            iternum++;
            NNEntry *pFirst = ps.top();
            //test code
//            std::cerr<<"pfirst\t"<<pFirst->m_dist.opt<<endl;
            // report all nearest neighbors with equal greatest distances.
            // (neighbors can be more than k, if many happen to have the same greatest distance).
            if (count >= k && pFirst->m_dist.opt > knearest) {
//            std::cerr<<"find minDist"<<knearest<<"\n";
                break;
            }
            switch (pFirst->m_type) {
                case 0: {//inner and leaf node
                    ps.pop();
                    NodePtr n = readNode(pFirst->m_id);
                    v.visitNode(*n);
#ifdef TJDEBUG
                    cerr<<iternum<<"\tnode with "<<pFirst->m_id<<"\t"<<pFirst->m_dist.opt<<"\t"<<n->m_nodeMBR<<endl;
#endif
                    for (uint32_t cChild = 0; cChild < n->m_children; ++cChild) {
                        if (n->m_level == 0) {
                            double pd;
                            double ts, te;
                            //test code
//                            if(n->m_pIdentifier[cChild]==788){
//                                std::cerr<<"";
//                            }
                            pd = std::max(0.0, simpleTraj.sbbDistInfer(*n->m_ptrxSBB[cChild], tjstat->vmax).opt);
                            ts = n->m_ptrxSBB[cChild]->m_startTime;
                            te = n->m_ptrxSBB[cChild]->m_endTime;
                            leafInfo *e = new leafInfo();
                            e->m_se = n->m_se[cChild];
                            e->m_hasPrev = (n->m_prevNode[cChild] != -1);
                            e->m_hasNext = (n->m_nextNode[cChild] != -1);
                            e->m_ts = ts;
                            e->m_te = te;
                            ps.push(new NNEntry(n->m_pIdentifier[cChild], e, pd, 1));
                        } else {
                            double pd = std::max(0.0, simpleTraj.nodeDist(*(n->m_ptrMBR[cChild])) - delta);
                            ps.push(new NNEntry(n->m_pIdentifier[cChild], nullptr, pd, 0));
                        }
                    }
                    delete pFirst;
                    break;
                }
                case 1: {//leaf node
                    ps.pop();
#ifdef TJDEBUG
                    cerr<<iternum<<"\tCB with "<<pFirst->m_dist.opt<<"\t"<<ps.explain(pFirst->m_pEntry->m_se.m_id)<<endl;
#endif
                    ps.loadPartTraj(pFirst->m_id, pFirst->m_pEntry,pFirst->m_dist.opt);
                    delete pFirst;
                    break;
                }
                case 3: {
                    ps.pop();
                    ++(m_stats.m_u64QueryResults);
                    ++count;
                    knearest = pFirst->m_dist.opt;
                    simpleData d(pFirst->m_id, pFirst->m_dist.opt);
                    v.visitData(d);
                    delete pFirst;
                    break;
                }

                default:
                    throw Tools::IllegalStateException("illegal NNEntry state" + std::to_string(pFirst->m_type));
            }
        }
        ps.clean();
    }
//    std::cout<<"knearest is"<<knearest<<std::endl;
//    std::cerr<<"iternum is "<<iternum<<"\n";
    m_stats.m_doubleExactQueryResults+=knearest;
}

