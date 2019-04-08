//
// Created by chuang on 4/8/19.
//


#include <cstring>

#include <spatialindex/SpatialIndex.h>

#include "R2Tree.h"
#include "Node.h"
#include "Index.h"
#include "Leaf.h"

using namespace SpatialIndex;
using namespace SpatialIndex::R2Tree;

Leaf::~Leaf()
{
}

Leaf::Leaf(SpatialIndex::R2Tree::R2Tree* pTree, id_type id): Node(pTree, id, 0, pTree->m_leafCapacity)
{
}
/*
NodePtr Leaf::chooseSubtree(const Mbbc&, uint32_t, std::stack<id_type>&)
{
    // should make sure to relinquish other PoolPointer lists that might be pointing to the
    // same leaf.
    return NodePtr(this, &(m_pTree->m_leafPool));
}
*/
NodePtr Leaf::findLeaf(const Mbbc& mbbc, id_type id, std::stack<id_type>&)
{
    for (uint32_t cChild = 0; cChild < m_children; ++cChild)
    {
        // should make sure to relinquish other PoolPointer lists that might be pointing to the
        // same leaf.
        if (m_pIdentifier[cChild] == id && mbbc == *(m_ptrMbbc[cChild])) return NodePtr(this, &(m_pTree->m_leafPool));
    }

    return NodePtr();
}
