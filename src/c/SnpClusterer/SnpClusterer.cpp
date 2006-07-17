#include "SnpClusterer.hpp"

#include "GraphAlgorithms.hpp"

#include <set>

#include <iostream>

using namespace std;

void
CSnpClusterer::GetClusters(TClusters &clusters)
{
    typedef CGraphAlgorithms<TSnp> TGraphAlgorithms;
    typedef TGraphAlgorithms::TNodePtrsContainer TNodePtrsContainer;
    typedef set<TSnp> TVisitedNodes;
    typedef set<TSnp> TAddedIndels;

    TAddedIndels addedIndels;
    clusters.clear();
    TGraphAlgorithms::TNodePtrsContainer subGraphs;
    TGraphAlgorithms::GetConnectedSubGraphs(m_Graph, subGraphs);
    for (TNodePtrsContainer::const_iterator i = subGraphs.begin();
         i != subGraphs.end(); ++i) {
        clusters.push_back(TSnps());
        const TGraph::TNodePtrs &nodePtrs = *i;
        addedIndels.clear();
        for (TGraph::TNodePtrs::const_iterator j = nodePtrs.begin();
             j != nodePtrs.end(); ++j) {
            const TSnp &snp = (*j)->GetData();
            clusters.back().push_back(snp);
            TIndels::iterator indel = m_Indels.find(snp);
            if (indel != m_Indels.end()) {
                if (!addedIndels.count(indel->second)) {
                    clusters.back().push_back(indel->second);
                    addedIndels.insert(indel->second);
                }
                m_Indels.erase(indel);
            }
        }
    }

    TGraph indelClusters;
    TClusters indelDecomposedClusters;
    TSnps indels;
    for (TIndels::iterator i = m_Indels.begin(); i != m_Indels.end(); ++i) {
        indelClusters.AddData(i->first, i->second);
    }
    subGraphs.clear();
    TGraphAlgorithms::GetConnectedSubGraphs(indelClusters, subGraphs);
    for (TNodePtrsContainer::iterator i = subGraphs.begin();
         i != subGraphs.end(); ++i) {
        indelDecomposedClusters.clear();
        indels.clear();
        indelDecomposedClusters.push_back(TSnps());
        TGraph::TNodePtrs &nodePtrs = *i;
        sort(nodePtrs.begin(), nodePtrs.end(), CSnpComparator());
        for (TGraph::TNodePtrs::iterator j = nodePtrs.begin();
             j != nodePtrs.end(); ++j) {
            TSnp snp = const_cast<TSnp>((*j)->GetData());
            if (snp->GetLength()) {
                x_AddToDecomposedCluster(indelDecomposedClusters, snp);
            }
            else {
                indels.push_back(snp);
            }
        }
        for (TClusters::iterator j = indelDecomposedClusters.begin();
             j != indelDecomposedClusters.end(); ++j) {
            copy(indels.begin(), indels.end(), back_inserter(*j));
            clusters.push_back(TSnps());
            copy(j->begin(), j->end(), back_inserter(clusters.back()));
        }
    }
}

void
CSnpClusterer::x_AddToDecomposedCluster(TClusters &clusters,
                                        CSnpDataType *snp)
{
    TSnps *snps = 0;
    for (TClusters::iterator i = clusters.begin(); i != clusters.end(); ++i) {
        bool found = false;
        for (TSnps::iterator j = i->begin(); j != i->end(); ++j) {
            if ((*j)->GetId() == snp->GetId()) {
                found = true;
                break;
            }
        }
        if (!found) {
            snps = &*i;
            break;
        }
    }
    if (!snps) {
        clusters.push_back(TSnps());
        snps = &clusters.back();
    }
    snps->push_back(snp);
}

void
CSnpClusterer::AddSnps(const TSnpData &snpData)
{
    bool complement1, complement2;
    CSnpDataType *snp1 = x_GetSnp(snpData[0], complement1);
    CSnpDataType *snp2 = x_GetSnp(snpData[1], complement2);
    if (complement1) {
        x_ComplementSnp(snp2);
    }
    if (complement2) {
        x_ComplementSnp(snp1);
    }
    if (snp1->GetLength() && snp2->GetLength()) {
        AddSnp(snp1, snp2);
    }
    else {
        if (snp1->GetLength()) {
            AddIndel(snp1, snp2);
        }
        else {
            AddIndel(snp2, snp1);
        }
    }
}

CSnpDataType*
CSnpClusterer::x_GetSnp(const string &data, bool &complement)
{
    CSnpDataType *snp = new CSnpDataType(data);
    CSnpDataType *&loaded =
        m_Id2Snp[TId2Snp::key_type(snp->GetId(),
                                   TInterval(snp->GetFrom(), snp->GetTo()))];
    if (!loaded) {
        loaded = snp;
        complement = false;
    }
    else {
        complement = snp->IsPlus() ^ loaded->IsPlus();
        delete snp;
    }
    return loaded;
}
