//
// Created by chuang on 4/8/19.
//


#include <cstring>

#include <spatialindex/SpatialIndex.h>

#include "LERTree.h"
#include "Node.h"
#include "Index.h"
#include "L1Index.h"
#include "L1Index.h"

using namespace SpatialIndex;
using namespace SpatialIndex::LERTree;

L1Index::~L1Index()
{
}

L1Index::L1Index(SpatialIndex::LERTree::LERTree* pTree, id_type id): Node(pTree, id, 0, pTree->m_L1indexCapacity)
{
}
/*
NodePtr L1Index::chooseSubtree(const Region&, uint32_t, std::stack<id_type>&)
{
    // should make sure to relinquish other PoolPointer lists that might be pointing to the
    // same leaf.
    return NodePtr(this, &(m_pTree->m_leafPool));
}
*/
//NodePtr L1Index::findL1Index(const Region& mbc, id_type id, std::stack<id_type>&)
//{
//    for (uint32_t cChild = 0; cChild < m_children; ++cChild)
//    {
//        // should make sure to relinquish other PoolPointer lists that might be pointing to the
//        // same leaf.
//        if (m_pIdentifier[cChild] == id && mbc == *(m_ptrMBR[cChild])) return NodePtr(this, &(m_pTree->m_leafPool));
//    }
//
//    return NodePtr();
//}
