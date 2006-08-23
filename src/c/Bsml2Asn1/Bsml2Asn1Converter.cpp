#include "Bsml2Asn1Converter.hpp"
#include "FastaIndexer.hpp"
#include "Registry.hpp"

#include <stdexcept>
#include <iostream>

#include <xercesc/util/PlatformUtils.hpp>

#include <xercesc/framework/LocalFileInputSource.hpp>

#include <xalanc/XPath/XPathEvaluator.hpp>
#include <xalanc/XPath/XObject.hpp>
#include <xalanc/XPath/NodeRefList.hpp>

#include <xalanc/XalanSourceTree/XalanSourceTreeDOMSupport.hpp>
#include <xalanc/XalanSourceTree/XalanSourceTreeInit.hpp>
#include <xalanc/XalanSourceTree/XalanSourceTreeParserLiaison.hpp>

#include <xalanc/DOMSupport/XalanDocumentPrefixResolver.hpp>

#include <xalanc/XalanDOM/XalanNamedNodeMap.hpp>
#include <xalanc/XalanDOM/XalanNodeList.hpp>

#include <corelib/ncbiobj.hpp>

#include <objects/seqfeat/Seq_feat.hpp>
#include <objects/seqfeat/SeqFeatData.hpp>
#include <objects/seqfeat/Gene_ref.hpp>
#include <objects/seqfeat/RNA_ref.hpp>
#include <objects/seqfeat/Cdregion.hpp>
#include <objects/seqfeat/BioSource.hpp>

#include <objects/seqloc/Seq_loc.hpp>
#include <objects/seqloc/Seq_interval.hpp>

#include <objects/seq/Bioseq.hpp>
#include <objects/seq/Seq_inst.hpp>
#include <objects/seq/Seq_data.hpp>
#include <objects/seq/Seq_descr.hpp>
#include <objects/seq/Seqdesc.hpp>
#include <objects/seq/Pubdesc.hpp>

#include <objects/seqset/Seq_entry.hpp>
#include <objects/seqset/Bioseq_set.hpp>

#include <objects/submit/Seq_submit.hpp>
#include <objects/submit/Submit_block.hpp>
#include <objects/submit/Contact_info.hpp>

#include <objects/biblio/Author.hpp>
#include <objects/biblio/Affil.hpp>
#include <objects/biblio/Cit_sub.hpp>
#include <objects/biblio/Auth_list.hpp>
#include <objects/biblio/Cit_gen.hpp>

#include <objects/pub/Pub_equiv.hpp>
#include <objects/pub/Pub.hpp>

#include <objects/general/Person_id.hpp>
#include <objects/general/Name_std.hpp>
#include <objects/general/Int_fuzz.hpp>

#include <objects/taxon1/taxon1.hpp>

#include <objmgr/scope.hpp>
#include <objmgr/bioseq_ci.hpp>
#include <objmgr/seq_vector.hpp>

#include <objmgr/util/sequence.hpp>

#include <serial/serialbase.hpp>
#include <serial/iterator.hpp>

using namespace std;
XERCES_CPP_NAMESPACE_USE
XALAN_CPP_NAMESPACE_USE
USING_NCBI_SCOPE;
using namespace ncbi::objects;

CBsml2Asn1Converter::CBsml2Asn1Converter()
{
    x_Init();
}

CBsml2Asn1Converter::~CBsml2Asn1Converter()
{
    delete m_XPathEval;
    delete m_Liaison;
    delete m_DomSupport;
    delete m_SrcTreeInit;
    XPathEvaluator::terminate();
    XMLPlatformUtils::Terminate();
    if (m_Taxon) {
        m_Taxon->Fini();
        delete m_Taxon;
    }
}

void
CBsml2Asn1Converter::Convert(const string &bsmlFile, const string &asn1File,
                             istream *contact, bool update, bool genProdSet)
{
    XalanDOMString fn(bsmlFile.c_str());
    LocalFileInputSource input(fn.c_str());
    XalanDocument *doc = m_Liaison->parseXMLStream(input);
    XalanDocumentPrefixResolver resolver(doc);
    XalanNode *contextNode =
        m_XPathEval->selectSingleNode(*m_DomSupport, doc,
                                      XalanDOMString("/Bsml").c_str(),
                                      resolver);
    if (!contextNode) {
        throw runtime_error(bsmlFile + " does not contain a <Bsml> node");
    }
    CRef<CSeq_submit> submit(new CSeq_submit);
    CRef<CSeq_entry> submitEntry(new CSeq_entry);
    submit->SetData().SetEntrys().push_back(submitEntry);
    TIndexes indexes;
    TTitles titles;
    TIds sidGbIds;
    TCompleteness completeness;
    TFeatures features;
    int gcode = 0;
    TWrittenIds writtenGenes;
    TWrittenIds pseudoProts;
    TNodeList results;
    CSubmit_block &submitBlock = submit->SetSub();
    CRef<TSeqDesc> genomeRef;
    submitBlock.SetSubtype(update ? CSubmit_block::eSubtype_update :
                           CSubmit_block::eSubtype_new);

    x_GetNodeList(results, contextNode, "//Genome", resolver);
    for (TNodeList::size_type i = 0; i < results.getLength(); ++i) {
        genomeRef.Reset(x_ProcessGenome(results.item(i), resolver));
        if (!genomeRef.IsNull()) {
            const COrgName &orgName =
                genomeRef->GetSource().GetOrg().GetOrgname();
            if (orgName.IsSetGcode()) {
                gcode = orgName.GetGcode();
            }
        }
    }

    TElement *chromo = x_GetElement(contextNode,
                                    "//Sequence[@class='assembly']/"
                                    "Attribute[@name='chromosome']",
                                    resolver);
    if (chromo) {
        CRef<CSubSource> chromoSrc
            (new CSubSource
             (CSubSource::eSubtype_chromosome,
              x_ConvertToStdStr
              (chromo->getAttribute(XalanDOMString("content")))));
        genomeRef->SetSource().SetSubtype().push_back(chromoSrc);
    }

    x_GetNodeList(results, contextNode, "//Feature", resolver);
    for (TNodeList::size_type i = 0; i < results.getLength(); ++i) {
        TNode *feat = results.item(i);
        string idStr = x_ConvertToStdStr
            (feat->getAttributes()->getNamedItem(XalanDOMString("id"))->
             getNodeValue());
        features.insert(TFeatures::value_type(idStr, feat));
    }

    if (update) {
        x_GetNodeList(results, contextNode, "//Sequence", resolver);
        for (TNodeList::size_type i = 0; i < results.getLength(); ++i) {
            TNode *seq = results.item(i);
            TNode *xref = x_GetNode(seq,
                                    "Cross-reference[@database='Genbank']",
                                    resolver);
            if (!xref) {
                continue;
            }
            TNode *sidNode =
                seq->getAttributes()->getNamedItem(XalanDOMString("id"));
            string sid(x_ConvertToStdStr(sidNode->getNodeValue()));
            TNode *gbNode =
                xref->getAttributes()->getNamedItem(XalanDOMString
                                                    ("identifier"));
            string gbId(x_ConvertToStdStr(gbNode->getNodeValue()));
            sidGbIds.insert(TIds::value_type(sid, gbId));
        }
    }

    TFeatGroups featGroups;
    x_GetNodeList(results, contextNode, "//Sequence[@class='assembly']",
                  resolver);
    for (TNodeList::size_type i = 0; i < results.getLength(); ++i) {
        CRef<CSeq_entry> entry(new CSeq_entry);
        TBioseq *seq = x_ProcessSequence(results.item(i), indexes,
                                         titles, sidGbIds,
                                         completeness, gcode,
                                         features, update, genProdSet,
                                         writtenGenes,
                                         pseudoProts, featGroups, resolver);
        if (seq) {
            entry->SetSeq(*seq);
            if (!featGroups.empty()) {
                submitEntry->SetSet().SetSeq_set().push_back(entry);
            }
            else {
                submitEntry->SetSeq(entry->SetSeq());
            }
        }
    }
    if (!featGroups.empty()) {
        submitEntry->SetSet().SetClass(genProdSet ?
                                       CBioseq_set::eClass_gen_prod_set :
                                       CBioseq_set::eClass_nuc_prot);
        CRef<CScope> scope(new CScope(*m_ObjMgr));
        scope->AddTopLevelSeqEntry(*submitEntry);

        for (TFeatGroups::iterator i = featGroups.begin();
             i != featGroups.end(); ++i) {
            SFeatGroup &featGroup = *i;
            TEntry *entrySeq = genProdSet ?
                entrySeq = x_CreateRnaProtSet(featGroup, *scope) :
                entrySeq = x_CreateProteinSequence(*featGroup.m_Cds,
                                                   *featGroup.m_Rna,
                                                   *scope);
            if (entrySeq) {
                submitEntry->SetSet().SetSeq_set().push_back
                    (CRef<TEntry>(entrySeq));
            }
        }
        if (!genomeRef.IsNull()) {
            submitEntry->SetSet().SetDescr().Set().push_back(genomeRef);
        }
    }
    else {
        if (!genomeRef.IsNull()) {
            submitEntry->SetSeq().SetDescr().Set().push_back(genomeRef);
        }
    }

    if (contact) {
        TRegistry reg(*contact);
        submitBlock.SetContact(*x_ProcessContact(reg));
        submitBlock.SetCit(*x_ProcessCitation(reg));
    }

    else {
        x_GetNodeList(results, contextNode, "//Resource", resolver);
        for (TNodeList::size_type i = 0; i < results.getLength(); ++i) {
            TNode *node = results.item(i);
            TNode *type = x_GetNode(node, "descendant::Type", resolver);
            if (!type) {
                continue;
            }
            TNode *typeValue = type->getFirstChild();
            if (!typeValue) {
                continue;
            }
            string typeValueStr = x_ConvertToStdStr(typeValue->getNodeValue());
            if (typeValueStr == "contact") {
                submitBlock.SetContact(*x_ProcessContact(node, resolver));
            }
            else if (typeValueStr == "citation") {
                submitBlock.SetCit(*x_ProcessCitation(node, resolver));
            }
            else if (typeValueStr == "publication") {
                if (!featGroups.empty()) {
                    submitEntry->SetSet().SetDescr().Set().push_back
                        (CRef<TSeqDesc>(x_ProcessPublication(node, resolver)));
                }
                else {
                    submitEntry->SetSeq().SetDescr().Set().push_back
                        (CRef<TSeqDesc>(x_ProcessPublication(node, resolver)));
                }
            }
        }
    }

    x_SetTitles(*submitEntry);

    x_CleanupIndexes(indexes);
    ofstream ofs(asn1File.c_str());
    ofs << MSerial_AsnText << *submit;
}

CBsml2Asn1Converter::TEntry*
CBsml2Asn1Converter::x_CreateProteinSequence(TSeqFeat &cds, TSeqFeat &rna,
                                             TScope &scope)
{
    string seqStr;
    CCdregion_translate::TranslateCdregion(seqStr, cds, scope, false);
    CSeq_entry *entry = x_CreateSequence(seqStr, cds, false);
    CBioseq &seq = entry->SetSeq();
    if (rna.GetData().GetRna().IsSetExt()) {
        const string &title = rna.GetData().GetRna().GetExt().GetName();
        CRef<CSeq_annot> annot(new CSeq_annot);
        seq.SetAnnot().push_back(annot);
        CRef<CSeq_feat> feat(new CSeq_feat);
        if (cds.IsSetPartial() && cds.GetPartial()) {
            feat->SetPartial(true);
        }
        annot->SetData().SetFtable().push_back(feat);
        feat->SetData().SetProt().SetName().push_back(title);
        CRef<CSeq_loc> loc(new CSeq_loc(cds.SetProduct().SetWhole(), 0,
                                        seqStr.length() - 1));
        if (cds.GetLocation().IsPartialStart(eExtreme_Biological)) {
            loc->SetPartialStart(true, eExtreme_Biological);
        }
        if (cds.GetLocation().IsPartialStop(eExtreme_Biological)) {
            loc->SetPartialStop(true, eExtreme_Biological);
        }
        feat->SetLocation(*loc);
    }
    return entry;
}

CBsml2Asn1Converter::TEntry*
CBsml2Asn1Converter::x_CreateRnaSequence(TSeqFeat &rna, TSeqFeat &gene,
                                         TScope &scope)
{
    string seqStr;
    CSeqVector seqVec(rna.GetLocation(), scope,
                      CBioseq_Handle::eCoding_Iupac);
    seqVec.GetSeqData(seqVec.begin(), seqVec.end(), seqStr);
    CSeq_entry *entry = x_CreateSequence(seqStr, rna, true);
    CBioseq &seq = entry->SetSeq();
    CRef<CSeq_annot> annot(new CSeq_annot);
    seq.SetAnnot().push_back(annot);
    CRef<CSeq_feat> feat(new CSeq_feat);
    annot->SetData().SetFtable().push_back(feat);
    feat->Assign(gene);
    CRef<CSeq_loc> loc(new CSeq_loc(rna.SetProduct().SetWhole(), 0,
                                    seqStr.length() - 1));
    feat->SetLocation(*loc);
    if (gene.GetLocation().IsPartialStart(eExtreme_Biological)) {
        loc->SetPartialStart(true, eExtreme_Biological);
    }
    if (gene.GetLocation().IsPartialStop(eExtreme_Biological)) {
        loc->SetPartialStop(true, eExtreme_Biological);
    }

    CRef<CSeqdesc> titleDesc(new CSeqdesc);
    seq.SetDescr().Set().push_back(titleDesc);
    if (rna.GetData().GetRna().IsSetExt()) {
        titleDesc->SetTitle(rna.GetData().GetRna().GetExt().GetName());
    }
    return entry;
}

CBsml2Asn1Converter::TEntry*
CBsml2Asn1Converter::x_CreateSequence(const std::string &seqStr,
                                       TSeqFeat &feat, bool rna)
{
    CRef<CSeq_entry> entry(new CSeq_entry);
    CBioseq &seq = entry->SetSeq();
    seq.SetId().push_back(CRef<CSeq_id>(&feat.SetProduct().SetWhole()));
    CRef<CSeq_data> data(new CSeq_data(seqStr,
                                       rna ? CSeq_data::e_Iupacna :
                                       CSeq_data::e_Iupacaa));
    seq.SetInst().SetSeq_data(*data);
    seq.SetInst().SetRepr(CSeq_inst::eRepr_raw);
    seq.SetInst().SetMol(rna ? CSeq_inst::eMol_rna : CSeq_inst::eMol_aa);
    seq.SetInst().SetLength(seqStr.length());

    CSeq_descr::Tdata &descs = seq.SetDescr().Set();
    CRef<CSeqdesc> molInfo(new CSeqdesc);
    descs.push_back(molInfo);
    molInfo->SetMolinfo().SetBiomol(rna ? CMolInfo::eBiomol_mRNA :
                                    CMolInfo::eBiomol_peptide);
    if (!rna) {
        molInfo->SetMolinfo().SetTech(CMolInfo::eTech_concept_trans);
    }
    bool partialLeft = feat.GetLocation().IsPartialStart(eExtreme_Biological);
    bool partialRight = feat.GetLocation().IsPartialStop(eExtreme_Biological);
    if (partialLeft && partialRight) {
        molInfo->SetMolinfo().SetCompleteness
            (CMolInfo::eCompleteness_no_ends);
    }
    else if (partialLeft) {
        molInfo->SetMolinfo().SetCompleteness
            (CMolInfo::eCompleteness_no_left);
    }
    else if (partialRight) {
        molInfo->SetMolinfo().SetCompleteness
            (CMolInfo::eCompleteness_no_right);
    }
    return entry.Release();
}

void
CBsml2Asn1Converter::x_Init()
{
    XMLPlatformUtils::Initialize();
    XPathEvaluator::initialize();
    m_SrcTreeInit = new XalanSourceTreeInit;
    m_DomSupport = new XalanSourceTreeDOMSupport;
    m_Liaison = new XalanSourceTreeParserLiaison;
    m_DomSupport->setParserLiaison(m_Liaison);
    m_XPathEval = new XPathEvaluator;
    m_Taxon = 0;
    m_ObjMgr.Reset(CObjectManager::GetInstance());

    m_ValidClasses.insert("GATC_rich_region");
    m_ValidClasses.insert("repeat_region");
    m_ValidClasses.insert("microsatellite");
    m_ValidClasses.insert("transposable_element");
}

void
CBsml2Asn1Converter::x_InitTaxon()
{
    if (!m_Taxon) {
        m_Taxon = new TTaxon;
        STimeout timeout = {5, 0};
        m_Taxon->Init(&timeout);
    }
}

CBsml2Asn1Converter::TNodeList&
CBsml2Asn1Converter::x_GetNodeList(TNodeList &results, TNode *contextNode,
                                   const string &searchPath,
                                   TResolver &resolver)
{
    m_XPathEval->selectNodeList(results, *m_DomSupport, contextNode,
                                XalanDOMString(searchPath.c_str()).c_str(),
                                resolver);
    return results;
}

CBsml2Asn1Converter::TNode*
CBsml2Asn1Converter::x_GetNode(TNode *contextNode, const string &searchPath,
                               TResolver &resolver)
{
    return m_XPathEval->selectSingleNode
        (*m_DomSupport, contextNode,
         XalanDOMString(searchPath.c_str()).c_str(), resolver);
}

CBsml2Asn1Converter::TElement*
CBsml2Asn1Converter::x_GetElement(TNode *contextNode, const string &searchPath,
                                  TResolver &resolver)
{
    return dynamic_cast<TElement *>(x_GetNode(contextNode, searchPath,
                                              resolver));
}

CBsml2Asn1Converter::TBioseq*
CBsml2Asn1Converter::x_ProcessSequence(TNode *seq,
                                       TIndexes &indexes,
                                       TTitles &titles,
                                       TIds &sidGbIds,
                                       TCompleteness &completeness,
                                       int gcode,
                                       TFeatures &features,
                                       bool update,
                                       bool genProdSet,
                                       TWrittenIds &writtenGenes,
                                       TWrittenIds &pseudoProts,
                                       TFeatGroups &featGroups,
                                       TResolver &resolver)
{
    if (pseudoProts.count(x_ConvertToStdStr
                          (seq->getAttributes()->getNamedItem
                           (XalanDOMString("id"))->getNodeValue()))) {

        cout << "Skipping pseudogene protein " << seq->getAttributes()->getNamedItem(XalanDOMString("id"))->getNodeValue() << endl;

        return 0;
    }
    TNode *seqData = x_GetNode(seq, "descendant::Seq-data-import", resolver);
    if (!seqData) {

        cout << "No Seq-data-import found for "
             << seq->getAttributes()->getNamedItem(XalanDOMString("id"))->
            getNodeValue() << endl;

        return 0;
    }

    string seqClass
        (x_ConvertToStdStr(seq->getAttributes()->getNamedItem
                           (XalanDOMString("class"))->getNodeValue()));
    if (seqClass != "assembly" && seqClass != "polypeptide") {
        return 0;
    }
    bool prot = seqClass == "polypeptide";

    TNode *id = 0;
    if (update) {
        if (id = x_GetNode(seq,
                           "descendant::Cross-reference[@database='Genbank']",
                           resolver)) {
            id = id->getAttributes()->
                getNamedItem(XalanDOMString("identifier"));
        }
    }
    if (!id) {
        id = seq->getAttributes()->getNamedItem(XalanDOMString("id"));
    }
    string idStr(x_ConvertToStdStr(id->getNodeValue()));

    CRef<CBioseq> bioseq(new CBioseq);
    bioseq->SetId().push_back(CRef<TSeqId>(x_CreateSeqId(idStr)));
    CSeq_inst &inst = bioseq->SetInst();
    inst.SetRepr(CSeq_inst::eRepr_raw);
    inst.SetMol(prot ? CSeq_inst::eMol_aa : CSeq_inst::eMol_dna);
    string fastaFile
        (x_ConvertToStdStr(seqData->getAttributes()->getNamedItem
                           (XalanDOMString("source"))->getNodeValue()));
    string fastaId
        (x_ConvertToStdStr(seqData->getAttributes()->getNamedItem
                           (XalanDOMString("identifier"))->getNodeValue()));
    TIndexes::iterator idxIter = indexes.find(fastaFile);
    if (idxIter == indexes.end()) {
        idxIter =
            indexes.insert
            (TIndexes::value_type(fastaFile,
                                  new CFastaIndexer(fastaFile))).first;
    }
    string strSeqData;
    idxIter->second->GetSequence(fastaId, strSeqData);
    if (!isalpha(strSeqData[strSeqData.length() - 1])) {
        strSeqData.erase(strSeqData.length() - 1);
    }
    inst.SetLength(strSeqData.size());
    inst.SetSeq_data(*new CSeq_data(strSeqData, prot ? CSeq_data::e_Iupacaa :
                                    CSeq_data::e_Iupacna));

    TNodeList results;
    x_GetNodeList(results, seq, "descendant::Feature-group", resolver);
    CRef<CSeq_annot> annot(new CSeq_annot);
    TFeats &ftable = annot->SetData().SetFtable();
    for (TNodeList::size_type i = 0; i < results.getLength(); ++i) {
        x_ProcessFeatureGroup(results.item(i), id, ftable, titles, 
                              sidGbIds, 
                              completeness, gcode, features, writtenGenes,
                              pseudoProts, featGroups, genProdSet, resolver);
    }

    x_GetNodeList(results, seq, "descendant::Feature", resolver);
    for (TNodeList::size_type i = 0; i < results.getLength(); ++i) {
        x_ProcessMiscFeature(results.item(i), id, ftable, resolver);
    }

    if (!ftable.empty()) {
        bioseq->SetAnnot().push_back(annot);
    }
    CSeq_descr::Tdata &descs = bioseq->SetDescr().Set();
    CRef<CSeqdesc> molInfo(new CSeqdesc);
    descs.push_back(molInfo);
    molInfo->SetMolinfo().SetBiomol(CMolInfo::eBiomol_genomic);
    const CBioseq::TId &seqIds = bioseq->GetId();
    for (CBioseq::TId::const_iterator i = seqIds.begin();
         i != seqIds.end(); ++i) {
        if ((*i)->IdentifyAccession() & CSeq_id::eAcc_wgs) {
            molInfo->SetMolinfo().SetTech(CMolInfo::eTech_wgs);
        }
    }
    return bioseq.Release();
}

void
CBsml2Asn1Converter::x_ProcessFeatureGroup(TNode *featGroup,
                                           TNode *id, TFeats &featRefs,
                                           TTitles &titles,
                                           TIds &sidGbIds,
                                           TCompleteness &completeness,
                                           int gcode,
                                           TFeatures &features,
                                           TWrittenIds &writtenGenes,
                                           TWrittenIds &pseudoProts,
                                           TFeatGroups &featGroups,
                                           bool genProdSet,
                                           TResolver &resolver)
{
    TNode *gene = 0, *rna = 0, *prot = 0, *cds = 0;
    TNodes exons;
    TNodeList feats;
    x_GetNodeList(feats, featGroup,
                  "descendant::Feature-group-member[@feature-type]",
                  resolver);
    for (TNodeList::size_type i = 0; i < feats.getLength(); ++i) {
        TNode *feat = feats.item(i);
        string featType = x_ConvertToStdStr
            (feat->getAttributes()->
             getNamedItem(XalanDOMString("feature-type"))->getNodeValue());
        string featRef =
            x_ConvertToStdStr
            (feat->getAttributes()->
             getNamedItem(XalanDOMString("featref"))->getNodeValue());
        TFeatures::iterator featIter = features.find(featRef);
        if (featIter == features.end()) {
            throw runtime_error("Cannot find feature element for id " +
                                featRef);
        }
        TNode *mappedFeat = featIter->second;
        if (featType == "gene") {
            gene = mappedFeat;
        }
        else if (featType == "transcript" ||
                 featType == "tRNA") {
            rna = mappedFeat;
        }
        else if (featType == "polypeptide") {
            if (!x_IsCanonical(mappedFeat, resolver)) {
                return;
            }
            prot = mappedFeat;
        }
        else if (featType == "exon") {
            exons.push_back(mappedFeat);
        }
        else if (featType == "CDS") {
            cds = mappedFeat;
        }
    }
    CRef<TSeqFeat> geneFeat(x_CreateGeneFeat(id, gene, resolver));
    CRef<TSeqFeat> rnaFeat;
    CRef<TSeqFeat> cdsFeat;

    if (!geneFeat.IsNull()) {
        rnaFeat.Reset(x_CreateRnaFeat(id, rna, exons, sidGbIds, genProdSet,
                                      resolver));
        if (!rnaFeat->GetData().GetRna().IsSetPseudo() ||
            !rnaFeat->GetData().GetRna().GetPseudo()) {
            cdsFeat.Reset(x_CreateCdsFeat(id, prot, cds, rna, exons, gcode,
                                          sidGbIds, resolver));
        }
        else {
            geneFeat->SetData().SetGene().SetPseudo(true);
        }

        string geneLabel;
        geneFeat->GetData().GetGene().GetLabel(&geneLabel);
        if (!writtenGenes.count(geneLabel)) {
            featRefs.push_back(geneFeat);
            writtenGenes.insert(geneLabel);
        }
        if (!rnaFeat.IsNull()) {
            featRefs.push_back(rnaFeat);
        }
        if (!cdsFeat.IsNull()) {
            featRefs.push_back(cdsFeat);
            bool isPartial5p = 
                cdsFeat->GetLocation().IsPartialStart(eExtreme_Positional);
            bool isPartial3p =
                cdsFeat->GetLocation().IsPartialStop(eExtreme_Positional);
            if (isPartial5p || isPartial3p) {
                string protId =
                    cdsFeat->GetProduct().GetWhole().GetSeqIdString(true);
                if (isPartial5p && isPartial3p) {
                    completeness.insert(TCompleteness::value_type
                                        (protId,
                                         CMolInfo::eCompleteness_no_ends));
                }
                else if (isPartial5p) {
                    completeness.insert(TCompleteness::value_type
                                        (protId,
                                         CMolInfo::eCompleteness_no_left));
                }
                else if (isPartial3p) {
                    completeness.insert(TCompleteness::value_type
                                        (protId,
                                         CMolInfo::eCompleteness_no_right));
                }
            }
        }

        if (!rnaFeat.IsNull() && !cdsFeat.IsNull()) {
            featGroups.push_back(SFeatGroup());
            SFeatGroup &featGroup = featGroups.back();
            featGroup.m_Gene = geneFeat.GetPointer();
            featGroup.m_Rna = rnaFeat.GetPointer();
            featGroup.m_Cds = cdsFeat.GetPointer();
        }

    }
    if (rna && prot) {
        TNode *protIdNode = x_GetNode
            (prot, "descendant::Link[@rel='sequence']", resolver);

        if (!protIdNode) {
            string error = "No sequence link for polypeptide " +
                x_ConvertToStdStr(prot->getAttributes()->
                                  getNamedItem(XalanDOMString("id"))->
                                  getNodeValue());
            throw runtime_error(error);
        }

        TNode *titleNode = x_GetNode
            (rna, "descendant::Attribute[@name='gene_product_name']",
             resolver);
        string protId =
            x_ConvertToStdStr(protIdNode->getAttributes()->
                              getNamedItem(XalanDOMString("href"))->
                              getNodeValue()).
            substr(1);
        if (protIdNode && titleNode) {
            string title =
                x_ConvertToStdStr(titleNode->getAttributes()->
                                  getNamedItem(XalanDOMString("content"))->
                                  getNodeValue());
            titles.insert(TTitles::value_type(protId, title));
            if (rnaFeat->GetData().GetRna().IsSetPseudo() &&
                rnaFeat->GetData().GetRna().GetPseudo()) {
                pseudoProts.insert(protId);
            }
        }
    }
}

void
CBsml2Asn1Converter::x_ProcessMiscFeature(TNode *feat, TNode *id,
                                          TFeats &featRefs,
                                          TResolver &resolver)
{
    string classStr = x_ConvertToStdStr
        (feat->getAttributes()->getNamedItem(XalanDOMString("class"))->
         getNodeValue());
    if (!m_ValidClasses.count(classStr)) {
        return;
    }
    CRef<TSeqFeat> miscFeat(new TSeqFeat);
    TNode *interval = x_GetNode(feat, "descendant::Interval-loc", resolver);
    miscFeat->SetLocation(*x_CreateSeqLoc(id, interval));
    if (classStr == "GATC_rich_region") {
        miscFeat->SetData().SetImp().SetKey("misc_feature");
    }
    else if (classStr == "repeat_region" ||
             classStr == "microsatellite" ||
             classStr == "transposable_element") {
        miscFeat->SetData().SetImp().SetKey("repeat_region");
        if (classStr == "microsatellite") {
            miscFeat->AddQualifier("rpt_type", "tandem");
            miscFeat->SetComment("microsatellite");
        }
    }
    TNodeList results;
    x_GetNodeList(results, feat, "descendant::Attribute", resolver);
    for (TNodeList::size_type i = 0; i < results.getLength(); ++i) {
        TNode *attr = results.item(i);
        string name = x_ConvertToStdStr
            (attr->getAttributes()->getNamedItem(XalanDOMString("name"))->
             getNodeValue());
        string content = x_ConvertToStdStr
            (attr->getAttributes()->getNamedItem(XalanDOMString("content"))->
             getNodeValue());
        if (classStr == "transposable_element" &&
            name == "comment") {
            name = "transposon";
        }
        else if (name == "repeat_unit") {
            name = "rpt_unit";
        }
        else if (name == "repeat_family") {
            name = "rpt_family";
        }
        if (name == "comment") {
            miscFeat->SetComment(content);
        }
        else {
            miscFeat->AddQualifier(name, content);
        }
    }
    featRefs.push_back(miscFeat);
}

CBsml2Asn1Converter::TSeqFeat*
CBsml2Asn1Converter::x_CreateGeneFeat(TNode *id, TNode *gene,
                                      TResolver &resolver)
{
    CRef<TSeqFeat> geneFeat(new TSeqFeat);
    TNode *locus = x_GetNode
        (gene, "descendant::Cross-reference[@identifier-type='pub_locus']",
         resolver);

    if (!locus) {
        cout << "No pub_locus found for gene " << gene->getAttributes()->getNamedItem(XalanDOMString("id"))->getNodeValue() << endl;
        return 0;
    }

    TNode *interval = x_GetNode(gene, "descendant::Interval-loc", resolver);
    geneFeat->SetLocation(*x_CreateSeqLoc(id, interval));
    geneFeat->SetData().SetGene().SetLocus_tag
        (x_ConvertToStdStr(locus->getAttributes()->
                           getNamedItem(XalanDOMString("identifier"))->
                           getNodeValue()));
    return geneFeat.Release();
}

CBsml2Asn1Converter::TSeqFeat*
CBsml2Asn1Converter::x_CreateRnaFeat(TNode *id, TNode *rna,
                                     TNodes &exons, TIds &sidGbIds,
                                     bool genProdSet, TResolver &resolver)
{
    CRef<TSeqFeat> rnaFeat(new TSeqFeat);
    TSeqLoc::TPacked_int &ints = rnaFeat->SetLocation().SetPacked_int();
    for (TNodes::iterator i = exons.begin(); i != exons.end(); ++i) {
        CRef<TSeqLoc> exonLoc
            (x_CreateSeqLoc(id,
                            x_GetNode(*i, "descendant::Interval-loc",
                                      resolver)));
        ints.AddInterval(exonLoc->SetInt());
    }
    x_SortPackedInt(ints);
    CRNA_ref::TType rnaType = CRNA_ref::eType_unknown;
    string featClass = x_ConvertToStdStr
        (rna->getAttributes()->getNamedItem(XalanDOMString("class"))->
         getNodeValue());
    if (featClass == "transcript") {
        rnaType = CRNA_ref::eType_mRNA;
    }
    else if (featClass == "tRNA") {
        rnaType = CRNA_ref::eType_tRNA;
    }
    rnaFeat->SetData().SetRna().SetType(rnaType);
    TNode *pseudo =
        x_GetNode(rna,
                  "descendant::Attribute[@name='SO'][@content='pseudogene']",
                  resolver);
    if (pseudo) {
        rnaFeat->SetData().SetRna().SetPseudo(true);
    }
    TNode *title =
        x_GetNode(rna, "descendant::Attribute[@name='gene_product_name']",
                  resolver);
    if (title) {
        rnaFeat->SetData().SetRna().SetExt().SetName
            (x_ConvertToStdStr
             (title->getAttributes()->getNamedItem
              (XalanDOMString("content"))->getNodeValue()));
    }
    if (genProdSet) {
        TElement *pubLocus = dynamic_cast<TElement *>
            (x_GetNode(rna,
                       "descendant::Cross-reference"
                       "[starts-with(@database, 'TIGR')]"
                       "[@identifier-type='pub_locus']", resolver));
        if (!pubLocus) {
            XalanDOMString rnaId = dynamic_cast<TElement *>
                (rna)->getAttribute(XalanDOMString("id"));
            cout << "No pub_locus found for " << rnaId << endl;
            return 0;
        }
        string prodId = "gnl|tigr|mrna." +
            x_ConvertToStdStr
            (pubLocus->getAttribute(XalanDOMString("identifier")));
        rnaFeat->SetProduct().SetWhole(*x_CreateSeqId(prodId));
    }
    return rnaFeat.Release();
}

CBsml2Asn1Converter::TSeqFeat*
CBsml2Asn1Converter::x_CreateCdsFeat(TNode *id, TNode *prot, TNode *cds,
                                     TNode *rna, TNodes &exons, int gcode,
                                     TIds &sidGbIds, 
                                     TResolver &resolver)
{
    CRef<TSeqFeat> cdsFeat(new TSeqFeat);
    TSeqLoc::TPacked_int &ints = cdsFeat->SetLocation().SetPacked_int();
    CRef<TSeqLoc> cdsInterval
        (x_CreateSeqLoc(id,
                        x_GetNode(prot,
                                  "descendant::Interval-loc", resolver)));
    TSeqPos from = cdsInterval->GetStart(eExtreme_Positional);
    TSeqPos to = cdsInterval->GetStop(eExtreme_Positional);
    TNode *partial3p = x_GetNode(cds, "Attribute[@name='three_prime_partial']",
                                 resolver);
    TNode *partial5p = x_GetNode(cds, "Attribute[@name='five_prime_partial']",
                                 resolver);
    bool isPartial3p = partial3p &&
        partial3p->getAttributes()->getNamedItem(XalanDOMString("content"))->
        getNodeValue() == XalanDOMString("1");
    bool isPartial5p = partial5p &&
        partial5p->getAttributes()->getNamedItem(XalanDOMString("content"))->
        getNodeValue() == XalanDOMString("1");
    for (TNodes::iterator i = exons.begin(); i != exons.end(); ++i) {
        CRef<TSeqLoc> exonLoc
            (x_CreateSeqLoc(id,
                            x_GetNode(*i, "descendant::Interval-loc",
                                      resolver)));
        TSeqPos exonFrom = exonLoc->GetStart(eExtreme_Positional);
        TSeqPos exonTo = exonLoc->GetStop(eExtreme_Positional);
        if (exonFrom < from) {
            exonLoc->SetInt().SetFrom(from);
            exonFrom = from;
        }
        if (exonTo > to) {
            exonLoc->SetInt().SetTo(to);
            exonTo = to;
        }

        if (exonFrom == from) {
            if (!exonLoc->IsReverseStrand()) {
                if (isPartial5p) {
                    exonLoc->SetInt().SetFuzz_from().SetLim
                        (CInt_fuzz::eLim_lt);
                }
            }
            else {
                if (isPartial3p) {
                    exonLoc->SetInt().SetFuzz_to().SetLim(CInt_fuzz::eLim_gt);
                }
            }
        }
        if (exonTo == to) {
            if (!exonLoc->IsReverseStrand()) {
                if (isPartial3p) {
                    exonLoc->SetInt().SetFuzz_to().SetLim(CInt_fuzz::eLim_gt);
                }
            }
            else {
                if (isPartial5p) {
                    exonLoc->SetInt().SetFuzz_from().SetLim
                        (CInt_fuzz::eLim_lt);
                }
            }
        }

        if (exonFrom < to && exonTo > from) {
            ints.AddInterval(exonLoc->SetInt());
        }
    }
    x_SortPackedInt(ints);
    CCdregion &cdsRef = cdsFeat->SetData().SetCdregion();
    TNode *frame = x_GetNode(cds, "descendant::Attribute[@name='frame']",
                             resolver);
    if (frame) {
        string frameStr = x_ConvertToStdStr
            (frame->getAttributes()->getNamedItem(XalanDOMString("content"))->
             getNodeValue());
        CCdregion::EFrame frameEnum;
        switch (atoi(frameStr.c_str())) {
        case 1:
            frameEnum = CCdregion::eFrame_one;
            break;
        case 2:
            frameEnum = CCdregion::eFrame_two;
            break;
        case 3:
            frameEnum = CCdregion::eFrame_three;
            break;
        default:
            frameEnum = CCdregion::eFrame_not_set;
            break;
        }
        cdsRef.SetFrame(frameEnum);
    }
    if (gcode) {
        CRef<CGenetic_code::C_E> gcodeRef(new CGenetic_code::C_E);
        gcodeRef->SetId(gcode);
        cdsRef.SetCode().Set().push_back(gcodeRef);
    }
    TNode *seqLink = x_GetNode(prot, "Link[@rel='sequence']", resolver);

    if (!seqLink) {
        cout << "No Link[@rel='sequence'] found for " << prot->getAttributes()->getNamedItem(XalanDOMString("id"))->getNodeValue() << endl;
        return 0;
    }
    TElement *pubLocus = dynamic_cast<TElement *>
        (x_GetNode(rna,
                   "descendant::Cross-reference"
                   "[starts-with(@database, 'TIGR')]"
                   "[@identifier-type='pub_locus']", resolver));
    if (!pubLocus) {
        XalanDOMString rnaId =
            dynamic_cast<TElement *>(rna)->getAttribute(XalanDOMString("id"));
        cout << "No pub_locus found for " << rnaId << endl;
        return 0;
    }
    string prodId = "gnl|tigr|cds." +
        x_ConvertToStdStr
        (pubLocus->getAttribute(XalanDOMString("identifier")));
    cdsFeat->SetProduct().SetWhole(*x_CreateSeqId(prodId));

    return cdsFeat.Release();
}

CBsml2Asn1Converter::TSeqLoc*
CBsml2Asn1Converter::x_CreateSeqLoc(TNode *id, TNode *interval)
{
    if (interval->getNodeName() != XalanDOMString("Interval-loc")) {
        return 0;
    }
    const XalanNamedNodeMap *attrs = interval->getAttributes();
    string idStr = 
        x_ConvertToStdStr(XalanDOMString(id->getNodeValue()));
    CRef<CSeq_id> seqId(x_CreateSeqId(idStr));
    TSeqPos from = atoll
        (x_ConvertToStdStr
         (attrs->getNamedItem(XalanDOMString("startpos"))->getNodeValue()).
         c_str());
    TSeqPos to = atoll
        (x_ConvertToStdStr
         (attrs->getNamedItem(XalanDOMString("endpos"))->getNodeValue()).
         c_str()) - 1;
    ENa_strand strand = attrs->getNamedItem(XalanDOMString("complement"))->
        getNodeValue() == XalanDOMString("0") ?
        eNa_strand_plus : eNa_strand_minus;
    CRef<TSeqLoc> seqLoc(new TSeqLoc(*seqId, from, to, strand));
    return seqLoc.Release();
}

CBsml2Asn1Converter::TSeqId*
CBsml2Asn1Converter::x_CreateSeqId(const string &id)
{
    CRef<TSeqId> seqId(new TSeqId);
    try {
        seqId->Set(id);
    }
    catch (const exception &e) {
        seqId->Set("lcl|" + id);
    }
    return seqId.Release();
}

CBsml2Asn1Converter::TSeqDesc*
CBsml2Asn1Converter::x_ProcessGenome(TNode *genome, TResolver &resolver)
{
    TNode *taxon = x_GetNode(genome,
                             "descendant::Cross-reference[@database='taxon']",
                             resolver);
    if (!taxon) {
        throw runtime_error("No taxid found");
    }
    string taxIdStr = x_ConvertToStdStr
        (taxon->getAttributes()->getNamedItem(XalanDOMString("identifier"))->
         getNodeValue());
    int taxId = atoi(taxIdStr.c_str());
    if (taxId <= 0) {
        return 0;
    }
    x_InitTaxon();
    bool isSpecies, isUncultured;
    string blastName;
    CConstRef<COrg_ref> orgRef(m_Taxon->GetOrgRef(taxId, isSpecies,
                                                  isUncultured, blastName));
    if (orgRef.IsNull()) {
        return 0;
    }

    CRef<CSeqdesc> desc(new CSeqdesc);
    CBioSource &bioSource = desc->SetSource();
    bioSource.SetGenome(CBioSource::eGenome_genomic);
    CRef<COrg_ref> nonConstOrgRef(new COrg_ref);
    nonConstOrgRef->Assign(*orgRef);
    bioSource.SetOrg(*nonConstOrgRef);
    TElement *org = dynamic_cast<TElement *>
        (x_GetNode(genome, "descendant::Organism", resolver));
    if (!org) {
        throw runtime_error("No organism element found");
    }
    string genus =
        x_ConvertToStdStr(org->getAttribute(XalanDOMString("genus")));
    string species =
        x_ConvertToStdStr(org->getAttribute(XalanDOMString("species")));
    if ((genus + " " + species) != orgRef->GetTaxname()) {
        cout << "Organism name in BSML (" << genus << " " << species
             << ") does not match data in taxonomy database ("
             << orgRef->GetTaxname() << ")" << endl;
    }
    return desc.Release();
}

CBsml2Asn1Converter::TContactInfo*
CBsml2Asn1Converter::x_ProcessContact(TNode *contactNode, TResolver &resolver)
{
    CRef<CContact_info> contactInfoRef(new CContact_info);
    TNode *organization = x_GetNode
        (contactNode, "Creator/Organization", resolver);
    TNode *person = x_GetNode(organization, "Person", resolver);
    contactInfoRef->SetContact(*x_CreateAuthor(person, resolver,
                                               organization));
    return contactInfoRef.Release();
}

CBsml2Asn1Converter::TContactInfo*
CBsml2Asn1Converter::x_ProcessContact(TRegistry &info)
{
    TRegistry::TEntries entries;
    info.GetEntriesByName("contact", entries);
    if (entries.empty()) {
        return 0;
    }
    CRef<CContact_info> contact(new CContact_info);
    const TRegistry::TRegistryData *reg = entries.begin()->second;
    contact->SetContact().SetAffil(*x_CreateAffil(*reg));
    contact->SetContact().SetName(*x_CreateName(*reg));
    return contact.Release();
}

CBsml2Asn1Converter::TAffil*
CBsml2Asn1Converter::x_CreateAffil(const TRegistryData &reg)
{
    CRef<TAffil> affilRef(new TAffil);
    CAffil::C_Std &affil = affilRef->SetStd();
    for (TRegistry::TRegistryData::TConstIterator i = reg.begin();
         i != reg.end(); ++i) {
        const std::string &key = i->first;
        const std::string &val = i->second;
        if (key == "affil") {
            affil.SetAffil(val);
        }
        else if (key == "div") {
            affil.SetDiv(val);
        }
        else if (key == "city") {
            affil.SetCity(val);
        }
        else if (key == "sub") {
            affil.SetSub(val);
        }
        else if (key == "country") {
            affil.SetCountry(val);
        }
        else if (key == "street") {
            affil.SetStreet(val);
        }
        else if (key == "email") {
            affil.SetEmail(val);
        }
        else if (key == "fax") {
            affil.SetFax(val);
        }
        else if (key == "phone") {
            affil.SetPhone(val);
        }
        else if (key == "postal-code") {
            affil.SetPostal_code(val);
        }
    }
    return affilRef.Release();
}

CBsml2Asn1Converter::TName*
CBsml2Asn1Converter::x_CreateName(const TRegistryData &reg)
{
    CRef<TName> nameRef(new TName);
    CName_std &name = nameRef->SetName();
    for (TRegistry::TRegistryData::TConstIterator i = reg.begin();
         i != reg.end(); ++i) {
        const std::string &key = i->first;
        const std::string &val = i->second;
        if (key == "last") {
            name.SetLast(val);
        }
        else if (key == "first") {
            name.SetFirst(val);
        }
    }
    return nameRef.Release();
}

CBsml2Asn1Converter::TCitation*
CBsml2Asn1Converter::x_ProcessCitation(TNode *citationNode,
                                       TResolver &resolver)
{
    CRef<CCit_sub> citationRef(new CCit_sub);
    CAuth_list &authList = citationRef->SetAuthors();
    TNode *organization = x_GetNode
        (citationNode, "Creator/Organization", resolver);
    TNodeList results;
    x_GetNodeList(results, organization, "Person", resolver);
    for (TNodeList::size_type i = 0; i < results.getLength(); ++i) {
        authList.SetNames().SetStd().push_back
            (CRef<CAuthor>(x_CreateAuthor(results.item(i), resolver)));
    }
    authList.SetAffil(*x_CreateAffil(organization, resolver));
    return citationRef.Release();
}

CBsml2Asn1Converter::TCitation*
CBsml2Asn1Converter::x_ProcessCitation(const TRegistry &info)
{
    TRegistry::TEntries entries;
    info.GetEntriesByName("contact", entries);
    if (entries.empty()) {
        return 0;
    }
    CRef<TCitation> citation(new TCitation);
    CAuth_list::C_Names::TStd &names =
        citation->SetAuthors().SetNames().SetStd();
    CRef<CAuthor> nameRef(new CAuthor);
    names.push_back(nameRef);
    const TRegistry::TRegistryData *reg = entries.begin()->second;
    citation->SetAuthors().SetAffil(*x_CreateAffil(*reg));
    nameRef->SetName(*x_CreateName(*reg));
    return citation.Release();
}

CBsml2Asn1Converter::TSeqDesc*
CBsml2Asn1Converter::x_ProcessPublication(TNode *publicationNode,
                                          TResolver &resolver)
{
    CRef<TSeqDesc> seqDescRef(new TSeqDesc);
    CPub_equiv &pubs = seqDescRef->SetPub().SetPub();
    CRef<CPub> pub(new CPub);
    pubs.Set().push_back(pub);
    TNode *cit = x_GetNode(publicationNode, "Content[@name='cit']", resolver);
    TNode *title = x_GetNode(publicationNode, "Title", resolver);
    TNodeList results;
    x_GetNodeList(results, publicationNode, "Creator/Person", resolver);
    CAuth_list &authList = pub->SetGen().SetAuthors();
    for (TNodeList::size_type i = 0; i < results.getLength(); ++i) {
        authList.SetNames().SetStd().push_back
            (CRef<CAuthor>(x_CreateAuthor(results.item(i), resolver)));
    }
    if (title) {
        pub->SetGen().SetTitle
            (x_ConvertToStdStr(title->getFirstChild()->getNodeValue()));
    }
    if (cit) {
        pub->SetGen().SetCit
            (x_ConvertToStdStr(cit->getFirstChild()->getNodeValue()));
    }
    return seqDescRef.Release();
}

CBsml2Asn1Converter::TAuthor*
CBsml2Asn1Converter::x_CreateAuthor(TNode *person, TResolver &resolver,
                                    TNode *organization)
{
    CRef<CAuthor> author(new CAuthor);
    TNode *lastName =
        person->getAttributes()->getNamedItem(XalanDOMString("lastname"));
    TNode *firstName =
        person->getAttributes()->getNamedItem(XalanDOMString("firstname"));
    if (lastName) {
        author->SetName().SetName().SetLast(x_ConvertToStdStr
                                            (lastName->getNodeValue()));
    }
    if (firstName) {
        author->SetName().SetName().SetFirst(x_ConvertToStdStr
                                             (firstName->getNodeValue()));
    }
    if (organization) {
        author->SetAffil(*x_CreateAffil(organization, resolver));
    }
    return author.Release();
}

CBsml2Asn1Converter::TAffil*
CBsml2Asn1Converter::x_CreateAffil(TNode *organization, TResolver &resolver)
{
    CRef<CAffil> affilRef(new CAffil);
    TNode *name = organization->getAttributes()->getNamedItem
        (XalanDOMString("name"));
    TNode *street = 0;
    TNode *city = 0;
    TNode *state = 0;
    TNode *country = 0;
    TNode *postalCode = 0;
    TNode *phone = 0;
    TNode *fax = 0;
    TNodeList results;
    x_GetNodeList(results, organization, "descendant::Contact-info", resolver);
    for (TNodeList::size_type i = 0; i < results.getLength(); ++i) {
        TNode *info = results.item(i);
        const XalanNamedNodeMap *infoAttrs = info->getAttributes();
        if (infoAttrs) {
            phone =
                infoAttrs->getNamedItem(XalanDOMString("telephone-number"));
            fax =
                infoAttrs->getNamedItem(XalanDOMString("fax-number"));
        }
        TNode *address = x_GetNode(info, "Postal-address", resolver);
        if (address) {
            const XalanNamedNodeMap *addressAttrs = address->getAttributes();
            if (addressAttrs) {
                street = addressAttrs->getNamedItem
                    (XalanDOMString("street-address1"));
                city = addressAttrs->getNamedItem
                    (XalanDOMString("city"));
                state = addressAttrs->getNamedItem
                    (XalanDOMString("state-province"));
                country = addressAttrs->getNamedItem
                    (XalanDOMString("country"));
                postalCode = addressAttrs->getNamedItem
                    (XalanDOMString("postal-code"));
            }
        }
    }
    if (phone) {
        affilRef->SetStd().SetPhone(x_ConvertToStdStr(phone->getNodeValue()));
    }
    if (fax) {
        affilRef->SetStd().SetFax(x_ConvertToStdStr(fax->getNodeValue()));
    }
    if (street) {
        affilRef->SetStd().SetStreet
            (x_ConvertToStdStr(street->getNodeValue()));
    }
    if (city) {
        affilRef->SetStd().SetCity(x_ConvertToStdStr(city->getNodeValue()));
    }
    if (state) {
        affilRef->SetStd().SetSub(x_ConvertToStdStr(state->getNodeValue()));
    }
    if (country) {
        affilRef->SetStd().SetCountry
            (x_ConvertToStdStr(country->getNodeValue()));
    }
    if (postalCode) {
        affilRef->SetStd().SetPostal_code
            (x_ConvertToStdStr(postalCode->getNodeValue()));
    }
    if (name) {
        affilRef->SetStd().SetAffil(x_ConvertToStdStr(name->getNodeValue()));
    }
    return affilRef.Release();
}

string
CBsml2Asn1Converter::x_ConvertToStdStr(const TXalanString &xalanStr)
{
    XalanDOMString::CharVectorType vec;
    TranscodeToLocalCodePage(xalanStr, vec);
    return string(vec.begin(), vec.end());
}

void
CBsml2Asn1Converter::x_CleanupIndexes(TIndexes &indexes)
{
    for (TIndexes::iterator i = indexes.begin(); i != indexes.end(); ++i) {
        delete i->second;
    }
}

bool
CBsml2Asn1Converter::x_IsCanonical(TNode *feat, TResolver &resolver)
{
    TNode *analysis = x_GetNode
        (feat, "descendant::Link[@rel='analysis'][@role='computed_by']",
         resolver);
    TNode *seqLink = x_GetNode
        (feat, "descendant::Link[@rel='sequence']", resolver);
    return !analysis && seqLink;
}

void
CBsml2Asn1Converter::x_SortPackedInt(TPackedInt &ints)
{
    CSeqIntComparator comp(!ints.IsReverseStrand());
    ints.Set().sort(comp);
}

CBsml2Asn1Converter::TEntry*
CBsml2Asn1Converter::x_CreateRnaProtSet(SFeatGroup &featGroup, TScope &scope)
{
    CRef<TEntry> entry(new TEntry);

    CBioseq_set &seqSet = entry->SetSet();
    seqSet.SetClass(CBioseq_set::eClass_nuc_prot);

    CRef<TEntry> rna(x_CreateRnaSequence(*featGroup.m_Rna,
                                         *featGroup.m_Gene,
                                         scope));
    CRef<TEntry> prot(x_CreateProteinSequence(*featGroup.m_Cds,
                                              *featGroup.m_Rna,
                                              scope));
    seqSet.SetSeq_set().push_back(rna);
    seqSet.SetSeq_set().push_back(prot);
    CRef<CSeq_feat> newCds(new CSeq_feat);
    newCds->Assign(*featGroup.m_Cds);
    CRef<CSeq_loc> newLoc
        (sequence::SourceToProduct(*featGroup.m_Rna,
                                   featGroup.m_Cds->GetLocation()));

    for (CTypeIterator<CSeq_interval> i(Begin(*newLoc)); i; ++i) {
        i->ResetFuzz_from();
        i->ResetFuzz_to();
    }
    if (featGroup.m_Cds->GetLocation().IsPartialStart(eExtreme_Biological)) {
        newLoc->SetPartialStart(true, eExtreme_Biological);
    }
    if (featGroup.m_Cds->GetLocation().IsPartialStop(eExtreme_Biological)) {
        newLoc->SetPartialStop(true, eExtreme_Biological);
    }
    newCds->SetLocation(*newLoc);
    CRef<CSeq_annot> annot(new CSeq_annot);
    annot->SetData().SetFtable().push_back(newCds);
    seqSet.SetAnnot().push_back(annot);

    return entry.Release();
}

void
CBsml2Asn1Converter::x_SetTitles(TEntry &entry)
{
    CRef<CScope> scope(new CScope(*m_ObjMgr));
    scope->AddTopLevelSeqEntry(entry);
    for (CTypeIterator<CBioseq> bioseqIter(Begin(entry));
         bioseqIter; ++bioseqIter) {
        if (bioseqIter->GetInst().GetMol() == CSeq_inst::eMol_dna) {
            continue;
        }
        CSeqdesc *oldTitle = 0;
        for (CTypeIterator<CSeqdesc> descIter(Begin(*bioseqIter));
             descIter; ++descIter) {
            if (descIter->IsTitle()) {
                oldTitle = &*descIter;
                break;
            }
        }
        CBioseq_Handle handle(scope->GetBioseqHandle
                              (*bioseqIter->GetFirstId(),
                               CScope::eGetBioseq_Loaded));
        string newTitle
            (sequence::GetTitle(handle,
                                sequence::fGetTitle_Reconstruct |
                                sequence::fGetTitle_Organism |
                                sequence::fGetTitle_AllProteins));
        if (oldTitle) {
            oldTitle->SetTitle() += " [" + newTitle + "]";
        }
        else {
            CRef<CSeqdesc> titleDesc(new CSeqdesc);
            bioseqIter->SetDescr().Set().push_back(titleDesc);
            titleDesc->SetTitle(newTitle);
        }
    }
}

bool
CBsml2Asn1Converter::CSeqIntComparator::operator()(const TData &int1,
                                                   const TData &int2)
{
    if (m_Forward) {
        return int1->GetFrom() < int2->GetFrom();
    }
    return int1->GetFrom() > int2->GetFrom();
}
