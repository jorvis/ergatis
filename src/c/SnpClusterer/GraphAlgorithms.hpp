#ifndef GRAPH_ALGORITHMS_HPP
#define GRAPH_ALGORITHMS_HPP

#include "Graph.hpp"

#include <set>

template<typename T>
class CGraphAlgorithms
{

public:
    typedef CGraph<T> TGraph;
    typedef std::vector<typename TGraph::TNodePtrs> TNodePtrsContainer;
    typedef std::set<const typename TGraph::CNode *> TVisitedNodes;
    typedef std::set<const typename TGraph::CEdge *> TVisitedEdges;

public:
    static void GetConnectedSubGraphs(const TGraph &graph,
                                      TNodePtrsContainer &subGraphs);
    static void DFS(const TGraph &graph, const typename TGraph::CNode *node,
                    TVisitedNodes &visitedNodes, TVisitedEdges &visitedEdges);
};

template<typename T>
inline
void
CGraphAlgorithms<T>::GetConnectedSubGraphs(const TGraph &graph,
                                           TNodePtrsContainer &subGraphs)
{
    TVisitedNodes allVisitedNodes;
    TVisitedNodes visitedNodes;
    TVisitedEdges visitedEdges;
    typename TGraph::TNodePtrs nodePtrs;
    graph.GetNodePtrs(nodePtrs);
    for (typename TGraph::TNodePtrs::const_iterator i = nodePtrs.begin();
         i != nodePtrs.end(); ++i) {
        if (!allVisitedNodes.count(*i)) {
            visitedNodes.clear();
            visitedEdges.clear();
            DFS(graph, *i, visitedNodes, visitedEdges);
            subGraphs.push_back(typename TNodePtrsContainer::value_type());
            for (typename TVisitedNodes::const_iterator j =
                     visitedNodes.begin();
                 j != visitedNodes.end(); ++j) {
                allVisitedNodes.insert(*j);
                subGraphs.back().push_back(*j);
            }
        }
    }
}

template<typename T>
inline
void
CGraphAlgorithms<T>::DFS(const TGraph &graph,
                         const typename TGraph::CNode *node,
                         TVisitedNodes &visitedNodes,
                         TVisitedEdges &visitedEdges)
{
    visitedNodes.insert(node);
    const typename TGraph::CNode::TEdges &edges = node->GetEdges();
    for (typename TGraph::CNode::TEdges::const_iterator i = edges.begin();
         i != edges.end(); ++i) {
        const typename TGraph::CEdge *edge = *i;
        if (!visitedEdges.count(edge)) {
            visitedEdges.insert(edge);
            const typename TGraph::CNode *otherNode =
                edge->GetFirstNode() != node ?
                edge->GetFirstNode () : edge->GetSecondNode();
            if (!visitedNodes.count(otherNode)) {
                visitedNodes.insert(otherNode);
                DFS(graph, otherNode, visitedNodes, visitedEdges);
            }
        }
    }
}

#endif
