//
// Created by chuang on 4/8/19.
//


#include <limits>

#include <spatialindex/SpatialIndex.h>
#include "R2Tree.h"
#include "Node.h"
#include "Leaf.h"
#include "Index.h"

using namespace SpatialIndex;
using namespace SpatialIndex::R2Tree;

Index::~Index()
{
}

Index::Index(SpatialIndex::R2Tree::R2Tree* pTree, id_type id, uint32_t level) : Node(pTree, id, level, pTree->m_indexCapacity)
{
}
/*
NodePtr Index::chooseSubtree(const Mbbc& mbbc, uint32_t insertionLevel, std::stack<id_type>& pathBuffer)
{
    if (m_level == insertionLevel) return NodePtr(this, &(m_pTree->m_indexPool));

    pathBuffer.push(m_identifier);

    uint32_t child = 0;

    child = findLeastEnlargement(mbbc);

    assert(child != std::numeric_limits<uint32_t>::max());

    NodePtr n = m_pTree->readNode(m_pIdentifier[child]);
    NodePtr ret = n->chooseSubtree(mbbc, insertionLevel, pathBuffer);
    assert(n.unique());
    if (ret.get() == n.get()) n.relinquish();

    return ret;
}
*/
NodePtr Index::findLeaf(const Mbbc& mbbc, id_type id, std::stack<id_type>& pathBuffer)
{
    pathBuffer.push(m_identifier);

    for (uint32_t cChild = 0; cChild < m_children; ++cChild)
    {
        if (m_ptrMbbc[cChild]->containsMbbc(mbbc))
        {
            NodePtr n = m_pTree->readNode(m_pIdentifier[cChild]);
            NodePtr l = n->findLeaf(mbbc, id, pathBuffer);
            if (n.get() == l.get()) n.relinquish();
            if (l.get() != 0) return l;
        }
    }

    pathBuffer.pop();

    return NodePtr();
}
/*
void Index::split(uint32_t dataLength, uint8_t* pData, Mbbc& mbbc, id_type id, NodePtr& ptrLeft, NodePtr& ptrRight)
{
    ++(m_pTree->m_stats.m_u64Splits);

    std::vector<uint32_t> g1, g2;

    switch (m_pTree->m_treeVariant)
    {
        case RV_LINEAR:
        case RV_QUADRATIC:
            rtreeSplit(dataLength, pData, mbbc, id, g1, g2);
            break;
        case RV_RSTAR:
            rstarSplit(dataLength, pData, mbbc, id, g1, g2);
            break;
        default:
            throw Tools::NotSupportedException("Index::split: Tree variant not supported.");
    }

    ptrLeft = m_pTree->m_indexPool.acquire();
    ptrRight = m_pTree->m_indexPool.acquire();

    if (ptrLeft.get() == 0) ptrLeft = NodePtr(new Index(m_pTree, m_identifier, m_level), &(m_pTree->m_indexPool));
    if (ptrRight.get() == 0) ptrRight = NodePtr(new Index(m_pTree, -1, m_level), &(m_pTree->m_indexPool));

    ptrLeft->m_nodeMbbc = m_pTree->m_infiniteMbbc;
    ptrRight->m_nodeMbbc = m_pTree->m_infiniteMbbc;

    uint32_t cIndex;

    for (cIndex = 0; cIndex < g1.size(); ++cIndex)
    {
        ptrLeft->insertEntry(0, 0, *(m_ptrMbbc[g1[cIndex]]), m_pIdentifier[g1[cIndex]]);
    }

    for (cIndex = 0; cIndex < g2.size(); ++cIndex)
    {
        ptrRight->insertEntry(0, 0, *(m_ptrMbbc[g2[cIndex]]), m_pIdentifier[g2[cIndex]]);
    }
}
*/
/*
uint32_t Index::findLeastEnlargement(const Mbbc& r) const
{
    double area = std::numeric_limits<double>::max();
    uint32_t best = std::numeric_limits<uint32_t>::max();

    MbbcPtr t = m_pTree->m_MbbcPool.acquire();

    for (uint32_t cChild = 0; cChild < m_children; ++cChild)
    {
        m_ptrMbbc[cChild]->getCombinedMbbc(*t, r);

        double a = m_ptrMbbc[cChild]->getArea();
        double enl = t->getArea() - a;

        if (enl < area)
        {
            area = enl;
            best = cChild;
        }
        else if (enl == area)
        {
            // this will rarely happen, so compute best area on the fly only
            // when necessary.
            if (a < m_ptrMbbc[best]->getArea()) best = cChild;
        }
    }

    return best;
}

uint32_t Index::findLeastOverlap(const Mbbc& r) const
{
    OverlapEntry** entries = new OverlapEntry*[m_children];

    double leastOverlap = std::numeric_limits<double>::max();
    double me = std::numeric_limits<double>::max();
    OverlapEntry* best = 0;

    // find combined region and enlargement of every entry and store it.
    for (uint32_t cChild = 0; cChild < m_children; ++cChild)
    {
        try
        {
            entries[cChild] = new OverlapEntry();
        }
        catch (...)
        {
            for (uint32_t i = 0; i < cChild; ++i) delete entries[i];
            delete[] entries;
            throw;
        }

        entries[cChild]->m_index = cChild;
        entries[cChild]->m_original = m_ptrMbbc[cChild];
        entries[cChild]->m_combined = m_pTree->m_MbbcPool.acquire();
        m_ptrMbbc[cChild]->getCombinedMbbc(*(entries[cChild]->m_combined), r);
        entries[cChild]->m_oa = entries[cChild]->m_original->getArea();
        entries[cChild]->m_ca = entries[cChild]->m_combined->getArea();
        entries[cChild]->m_enlargement = entries[cChild]->m_ca - entries[cChild]->m_oa;

        if (entries[cChild]->m_enlargement < me)
        {
            me = entries[cChild]->m_enlargement;
            best = entries[cChild];
        }
        else if (entries[cChild]->m_enlargement == me && entries[cChild]->m_oa < best->m_oa)
        {
            best = entries[cChild];
        }
    }

    if (me < -std::numeric_limits<double>::epsilon() || me > std::numeric_limits<double>::epsilon())
    {
        uint32_t cIterations;

        if (m_children > m_pTree->m_nearMinimumOverlapFactor)
        {
            // sort entries in increasing order of enlargement.
            ::qsort(entries, m_children,
                    sizeof(OverlapEntry*),
                    OverlapEntry::compareEntries);
            assert(entries[0]->m_enlargement <= entries[m_children - 1]->m_enlargement);

            cIterations = m_pTree->m_nearMinimumOverlapFactor;
        }
        else
        {
            cIterations = m_children;
        }

        // calculate overlap of most important original entries (near minimum overlap cost).
        for (uint32_t cIndex = 0; cIndex < cIterations; ++cIndex)
        {
            double dif = 0.0;
            OverlapEntry* e = entries[cIndex];

            for (uint32_t cChild = 0; cChild < m_children; ++cChild)
            {
                if (e->m_index != cChild)
                {
                    double f = e->m_combined->getIntersectingArea(*(m_ptrMbbc[cChild]));
                    if (f != 0.0) dif += f - e->m_original->getIntersectingArea(*(m_ptrMbbc[cChild]));
                }
            } // for (cChild)

            if (dif < leastOverlap)
            {
                leastOverlap = dif;
                best = entries[cIndex];
            }
            else if (dif == leastOverlap)
            {
                if (e->m_enlargement == best->m_enlargement)
                {
                    // keep the one with least area.
                    if (e->m_original->getArea() < best->m_original->getArea()) best = entries[cIndex];
                }
                else
                {
                    // keep the one with least enlargement.
                    if (e->m_enlargement < best->m_enlargement) best = entries[cIndex];
                }
            }
        } // for (cIndex)
    }

    uint32_t ret = best->m_index;

    for (uint32_t cChild = 0; cChild < m_children; ++cChild)
    {
        delete entries[cChild];
    }
    delete[] entries;

    return ret;
}

void Index::adjustTree(Node* n, std::stack<id_type>& pathBuffer)
{
    ++(m_pTree->m_stats.m_u64Adjustments);

    // find entry pointing to old node;
    uint32_t child;
    for (child = 0; child < m_children; ++child)
    {
        if (m_pIdentifier[child] == n->m_identifier) break;
    }

    // Mbbc needs recalculation if either:
    //   1. the NEW child Mbbc is not contained.
    //   2. the OLD child Mbbc is touching.
    bool bContained = m_nodeMbbc.containsMbbc(n->m_nodeMbbc);
    bool bTouches = m_nodeMbbc.touchesMbbc(*(m_ptrMbbc[child]));
    bool bRecompute = (! bContained || (bTouches && m_pTree->m_bTightMbbcs));

    *(m_ptrMbbc[child]) = n->m_nodeMbbc;

    if (bRecompute)
    {
        for (uint32_t cDim = 0; cDim < m_nodeMbbc.m_dimension; ++cDim)
        {
            m_nodeMbbc.m_pLow[cDim] = std::numeric_limits<double>::max();
            m_nodeMbbc.m_pHigh[cDim] = -std::numeric_limits<double>::max();

            for (uint32_t cChild = 0; cChild < m_children; ++cChild)
            {
                m_nodeMbbc.m_pLow[cDim] = std::min(m_nodeMbbc.m_pLow[cDim], m_ptrMbbc[cChild]->m_pLow[cDim]);
                m_nodeMbbc.m_pHigh[cDim] = std::max(m_nodeMbbc.m_pHigh[cDim], m_ptrMbbc[cChild]->m_pHigh[cDim]);
            }
        }
    }

    m_pTree->writeNode(this);

    if (bRecompute && (! pathBuffer.empty()))
    {
        id_type cParent = pathBuffer.top(); pathBuffer.pop();
        NodePtr ptrN = m_pTree->readNode(cParent);
        Index* p = static_cast<Index*>(ptrN.get());
        p->adjustTree(this, pathBuffer);
    }
}

void Index::adjustTree(Node* n1, Node* n2, std::stack<id_type>& pathBuffer, uint8_t* overflowTable)
{
    ++(m_pTree->m_stats.m_u64Adjustments);

    // find entry pointing to old node;
    uint32_t child;
    for (child = 0; child < m_children; ++child)
    {
        if (m_pIdentifier[child] == n1->m_identifier) break;
    }

    // Mbbc needs recalculation if either:
    //   1. the NEW child Mbbc is not contained.
    //   2. the OLD child Mbbc is touching.
    bool bContained = m_nodeMbbc.containsMbbc(n1->m_nodeMbbc);
    bool bTouches = m_nodeMbbc.touchesMbbc(*(m_ptrMbbc[child]));
    bool bRecompute = (! bContained || (bTouches && m_pTree->m_bTightMbbcs));

    *(m_ptrMbbc[child]) = n1->m_nodeMbbc;

    if (bRecompute)
    {
        for (uint32_t cDim = 0; cDim < m_nodeMbbc.m_dimension; ++cDim)
        {
            m_nodeMbbc.m_pLow[cDim] = std::numeric_limits<double>::max();
            m_nodeMbbc.m_pHigh[cDim] = -std::numeric_limits<double>::max();

            for (uint32_t cChild = 0; cChild < m_children; ++cChild)
            {
                m_nodeMbbc.m_pLow[cDim] = std::min(m_nodeMbbc.m_pLow[cDim], m_ptrMbbc[cChild]->m_pLow[cDim]);
                m_nodeMbbc.m_pHigh[cDim] = std::max(m_nodeMbbc.m_pHigh[cDim], m_ptrMbbc[cChild]->m_pHigh[cDim]);
            }
        }
    }

    // No write necessary here. insertData will write the node if needed.
    //m_pTree->writeNode(this);

    bool bAdjusted = insertData(0, 0, n2->m_nodeMbbc, n2->m_identifier, pathBuffer, overflowTable);

    // if n2 is contained in the node and there was no split or reinsert,
    // we need to adjust only if recalculation took place.
    // In all other cases insertData above took care of adjustment.
    if ((! bAdjusted) && bRecompute && (! pathBuffer.empty()))
    {
        id_type cParent = pathBuffer.top(); pathBuffer.pop();
        NodePtr ptrN = m_pTree->readNode(cParent);
        Index* p = static_cast<Index*>(ptrN.get());
        p->adjustTree(this, pathBuffer);
    }
}
*/