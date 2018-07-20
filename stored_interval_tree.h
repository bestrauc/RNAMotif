//
// Created by benjamin on 20.07.18.
//

#ifndef SEQAN_APPS_RNAMOTIF_STORED_INTERVAL_TREE_H
#define SEQAN_APPS_RNAMOTIF_STORED_INTERVAL_TREE_H

#include <seqan/misc/interval_tree.h>

namespace seqan {

template <typename TValue = int, typename TCargo = unsigned int, typename TNodeSpec = StoreIntervals>
class StoredIntervalTree
{
public:
    typedef Graph<Directed<void, WithoutEdgeId> > TGraph;
    typedef IntervalAndCargo<TValue, TCargo> TInterval;
    typedef IntervalTreeNode<TInterval, TNodeSpec> TNode;
    typedef String<TNode> TPropertyMap;

    TGraph g;
    TPropertyMap pm;
    size_t interval_counter;

    StoredIntervalTree()
    {
        interval_counter = 0;
    }

    template <typename TIterator, typename TCargoIterator>
    StoredIntervalTree(TIterator interval_begins,
                       TIterator interval_ends,
                       TCargoIterator interval_cargos,
                       size_t len)
    {
        String<TInterval> intervals;
        resize(intervals, len);
        size_t i = 0;
        while (i < len)
        {
            intervals[i].i1 = value(interval_begins);
            ++interval_begins;
            intervals[i].i2 = value(interval_ends);
            ++interval_ends;
            intervals[i].cargo = value(interval_cargos);
            ++interval_cargos;
            ++i;
        }
        createIntervalTree(*this, intervals);
    }

    template <typename TIterator>
    StoredIntervalTree(TIterator interval_begins,
                       TIterator interval_ends,
                       size_t len)
    {
        String<TInterval> intervals;
        resize(intervals, len);
        size_t i = 0;
        while (i < len)
        {
            intervals[i].i1 = value(interval_begins);
            ++interval_begins;
            intervals[i].i2 = value(interval_ends);
            ++interval_ends;
            intervals[i].cargo = i;
            ++i;
        }
        createIntervalTree(*this, intervals);
    }

    StoredIntervalTree(String<TInterval> intervals)
    {
        createIntervalTree(*this, intervals);
    }

    template <typename TTagSpec>
    StoredIntervalTree(String<TInterval> intervals, Tag<TTagSpec> const tag)
    {
        interval_counter = length(intervals);
        createIntervalTree(g, pm, intervals, tag);
    }

    StoredIntervalTree(String<TInterval> intervals, TValue center)
    {
        interval_counter = length(intervals);
        createIntervalTree(g, pm, intervals, center);
    }

};

// specialize seqan's Interval tree 'findIntervals' function to enable it to
// return IntervalAndCargo to include boundaries of the overlapping intervals
template <typename TValue, typename TCargo, typename TValue2>
inline void
findIntervals(
        String<IntervalAndCargo<TValue, TCargo> > & result,
        StoredIntervalTree<TValue, TCargo> const & it,
        TValue2 query)
{
    findIntervals(result, it.g, it.pm, query);
}

template <typename TSpec, typename TPropertyMap, typename TValue, typename TCargo, typename TValue2>
inline void
findIntervals(
        String< IntervalAndCargo<TValue, TCargo> > & result,
        Graph<TSpec> const & g,
        TPropertyMap const & pm,
		TValue2 query)
{
    typedef Graph<TSpec> const TGraph;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef typename Value<TPropertyMap>::Type TProperty;
    typedef typename Value<TProperty>::Type TPropertyValue;
    typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;

    resize(result, 0);
    if (empty(g))
        return;

    // start at root
    TVertexDescriptor act_knot = 0;
    TProperty act_prop = property(pm, act_knot);
    TProperty next_prop;

    while (true)
    {
        //typename Iterator<Graph<TSpec>, OutEdgeIterator>::Type it7;
        //Iter<Graph<TSpec>, GraphIterator<InternalOutEdgeIterator<OutEdgeIterator> > > it5(g, act_knot);
        //TOutEdgeIterator it4;
        TOutEdgeIterator it(g, act_knot);
        act_prop = property(pm, act_knot);
        if (act_prop.center < (TPropertyValue)query) // look in current node and right subtree
        {
            unsigned int i = 0;
            while (i<length(act_prop.list2) && rightBoundary(value(act_prop.list2, i))>(TPropertyValue) query)
            {
            	auto &val = value(act_prop.list2, i);
            	IntervalAndCargo<TValue, TCargo> tmp(leftBoundary(val), rightBoundary(val), cargo(val));
                appendValue(result, tmp, Generous());
                ++i;
            }
            if (atEnd(it))
                break;

            next_prop = property(pm, targetVertex(it));
            if (next_prop.center <= act_prop.center)
            {
                goNext(it);
                if (atEnd(it))
                    break;
            }
            act_knot = targetVertex(it);
        }
        else
        {
            if ((TPropertyValue)query < act_prop.center) // look in current node and left subtree
            {
                unsigned int i = 0;
                while (i < length(act_prop.list1) && leftBoundary(value(act_prop.list1, i)) <= (TPropertyValue)query)
                {
                	auto &val = value(act_prop.list1, i);
					IntervalAndCargo<TValue, TCargo> tmp(leftBoundary(val), rightBoundary(val), cargo(val));
					appendValue(result, tmp, Generous());
                    ++i;
                }
                if (atEnd(it))
                    break;

                next_prop = property(pm, targetVertex(it));
                if (next_prop.center >= act_prop.center)
                {
                    goNext(it);
                    if (atEnd(it))
                        break;
                }
                act_knot = targetVertex(it);
            }
            else  // look in current node only, as query is center
            {
                for (unsigned int i = 0; i < length(act_prop.list1); ++i){
                	auto &val = value(act_prop.list1, i);
					IntervalAndCargo<TValue, TCargo> tmp(leftBoundary(val), rightBoundary(val), cargo(val));
					appendValue(result, tmp, Generous());
                }
                break;
            }
        }
    }
}

// do an inorder traversal of the tree and report all intervals
template <typename TValue, typename TCargo>
inline void
getAllIntervals(
        String<IntervalAndCargo<TValue, TCargo> > & result,
		StoredIntervalTree<TValue, TCargo, StoreIntervals> const & it)
{
	getAllIntervals(result, it.g, it.pm);
}

template <typename TSpec, typename TPropertyMap, typename TValue, typename TCargo>
inline void
getAllIntervals(
        String<IntervalAndCargo<TValue, TCargo> > & result,
        Graph<TSpec> const & g,
        TPropertyMap const & pm)
{
    typedef Graph<TSpec> const TGraph;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;

	resize(result, 0);
	// start at root
	TVertexDescriptor act_knot = 0;
	getAllIntervals(result, g, pm, act_knot);
}

template <typename TSpec, typename TPropertyMap, typename TValue, typename TCargo, typename TVertexDescriptor>
inline void
getAllIntervals(
        String< IntervalAndCargo<TValue, TCargo> > & result,
        Graph<TSpec> const & g,
        TPropertyMap const & pm,
		TVertexDescriptor & act_knot)
{
	typedef Graph<TSpec> const TGraph;
	typedef typename Value<TPropertyMap>::Type TProperty;
	//typedef typename Value<TProperty>::Type TPropertyValue;
	typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;

	if (empty(g))
		return;

	TProperty act_prop = property(pm, act_knot);
	TProperty next_prop;

	TOutEdgeIterator it(g, act_knot);

	// go left
	if (!atEnd(it)){
		TVertexDescriptor next_knot = targetVertex(it);
		getAllIntervals(result, g, pm, next_knot);
		goNext(it);
	}

	// append center list
	for (unsigned int i = 0; i < length(act_prop.list1); ++i){
		auto &val = value(act_prop.list1, i);
		IntervalAndCargo<TValue, TCargo> tmp(leftBoundary(val), rightBoundary(val), cargo(val));
		appendValue(result, tmp, Generous());
	}

	// go right
	if (!atEnd(it))
	{
		TVertexDescriptor next_knot = targetVertex(it);
		getAllIntervals(result, g, pm, next_knot);
	}
}

}

#endif //SEQAN_APPS_RNAMOTIF_STORED_INTERVAL_TREE_H
