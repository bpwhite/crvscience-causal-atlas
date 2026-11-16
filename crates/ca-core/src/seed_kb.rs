// Copyright (C) 2026 CRV Science Network (https://crvscience.com)
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://gnu.org>.

#![forbid(unsafe_code)]

use crate::graph::{
    CausalEdge, EdgeClass, EdgeSign, EntityKind, EntityNode, KnowledgeGraph, PathwayDef,
    PathwayReactionEdge, PerturbationKind, PerturbationSignature, SubcellularLocation,
};
use std::collections::BTreeMap;

/// Assembles a self-contained open-source starter `KnowledgeGraph` synthesized from curated
/// Reactome, WikiPathways, SIGNOR 3.0, OmniPath / CollecTRI, Gene Ontology, Open Targets,
/// and SigCom LINCS records so `causal-atlas` operates offline out-of-the-box and passes
/// benchmark validation suites immediately prior to remote stream synchronization.
pub fn build_seed_knowledge_graph() -> KnowledgeGraph {
    let mut kg = KnowledgeGraph::new();

    let reg_nodes: &[(&str, &str, EntityKind, SubcellularLocation, &str, &str)] = &[
        // Extracellular ligands, cytokines, growth factors & secreted antagonists
        ("TNF", "Tumor Necrosis Factor", EntityKind::Cytokine, SubcellularLocation::Extracellular, "Cytokine Signaling", "SIGNOR / OmniPath"),
        ("IL1B", "Interleukin 1 Beta", EntityKind::Cytokine, SubcellularLocation::Extracellular, "Cytokine Signaling", "SIGNOR / OmniPath"),
        ("IL6", "Interleukin 6", EntityKind::Cytokine, SubcellularLocation::Extracellular, "Cytokine Signaling", "SIGNOR / OmniPath"),
        ("IFNG", "Interferon Gamma", EntityKind::Cytokine, SubcellularLocation::Extracellular, "Interferon Signaling", "SIGNOR / OmniPath"),
        ("CXCL8", "C-X-C Motif Chemokine Ligand 8", EntityKind::Cytokine, SubcellularLocation::Extracellular, "Chemokine Signaling", "Reactome"),
        ("CCL2", "C-C Motif Chemokine Ligand 2", EntityKind::Cytokine, SubcellularLocation::Extracellular, "Chemokine Signaling", "Reactome"),
        ("CXCL10", "C-X-C Motif Chemokine Ligand 10", EntityKind::Cytokine, SubcellularLocation::Extracellular, "Chemokine Signaling", "Reactome"),
        ("CXCL12", "Stromal Cell-Derived Factor 1", EntityKind::Cytokine, SubcellularLocation::Extracellular, "Chemokine Signaling", "Reactome"),
        ("CSF2", "Granulocyte-Macrophage CSF", EntityKind::Cytokine, SubcellularLocation::Extracellular, "Hematopoietic Cytokines", "Reactome"),
        ("AREG", "Amphiregulin", EntityKind::GrowthFactor, SubcellularLocation::Extracellular, "EGF Family", "SIGNOR"),
        ("BMP2", "Bone Morphogenetic Protein 2", EntityKind::GrowthFactor, SubcellularLocation::Extracellular, "TGF-beta / BMP", "SIGNOR"),
        ("BMP6", "Bone Morphogenetic Protein 6", EntityKind::GrowthFactor, SubcellularLocation::Extracellular, "TGF-beta / BMP", "SIGNOR"),
        ("WNT3A", "Wnt Family Member 3A", EntityKind::GrowthFactor, SubcellularLocation::Extracellular, "Wnt Signaling", "SIGNOR"),
        ("SOST", "Sclerostin", EntityKind::GeneOther, SubcellularLocation::Extracellular, "Wnt Antagonist", "SIGNOR"),
        ("DKK1", "Dickkopf WNT Signaling Inhibitor 1", EntityKind::GeneOther, SubcellularLocation::Extracellular, "Wnt Antagonist", "SIGNOR"),
        ("TFF1", "Trefoil Factor 1", EntityKind::GeneOther, SubcellularLocation::Extracellular, "Estrogen Target", "CollecTRI"),
        ("BGLAP", "Osteocalcin", EntityKind::GeneOther, SubcellularLocation::Extracellular, "Extracellular Matrix / Bone", "WikiPathways"),
        ("COL1A1", "Collagen Type I Alpha 1 Chain", EntityKind::GeneOther, SubcellularLocation::Extracellular, "Extracellular Matrix / Bone", "Reactome"),
        ("MMP9", "Matrix Metallopeptidase 9", EntityKind::Enzyme, SubcellularLocation::Extracellular, "Matrix Remodeling", "SIGNOR"),
        ("MMP13", "Matrix Metallopeptidase 13", EntityKind::Enzyme, SubcellularLocation::Extracellular, "Matrix Remodeling", "SIGNOR"),
        ("LCN2", "Lipocalin 2", EntityKind::GeneOther, SubcellularLocation::Extracellular, "Innate Immunity", "CollecTRI"),

        // Small-molecule ligands, metabolites & open-source ChEMBL/LINCS compounds
        ("LPS", "Lipopolysaccharide (TLR4 Agonist)", EntityKind::ChemicalDrug, SubcellularLocation::Extracellular, "Microbial PAMP", "SIGNOR"),
        ("ESTRADIOL", "17beta-Estradiol", EntityKind::ChemicalDrug, SubcellularLocation::Extracellular, "Steroid Hormone", "ChEMBL / SIGNOR"),
        ("FULVESTRANT", "Fulvestrant (SERD)", EntityKind::ChemicalDrug, SubcellularLocation::Extracellular, "Selective ER Downregulator", "ChEMBL / Open Targets"),
        ("DEXAMETHASONE", "Dexamethasone", EntityKind::ChemicalDrug, SubcellularLocation::Extracellular, "Glucocorticoid Agonist", "ChEMBL / Open Targets"),

        // Plasma membrane receptors & adhesion molecules
        ("TNFRSF1A", "TNF Receptor Superfamily Member 1A", EntityKind::TransmembraneReceptor, SubcellularLocation::PlasmaMembrane, "Death Domain Receptors", "SIGNOR"),
        ("TLR4", "Toll-Like Receptor 4", EntityKind::TransmembraneReceptor, SubcellularLocation::PlasmaMembrane, "Pattern Recognition Receptors", "SIGNOR"),
        ("IFNGR1", "Interferon Gamma Receptor 1", EntityKind::TransmembraneReceptor, SubcellularLocation::PlasmaMembrane, "Cytokine Receptors", "SIGNOR"),
        ("EGFR", "Epidermal Growth Factor Receptor", EntityKind::TransmembraneReceptor, SubcellularLocation::PlasmaMembrane, "Receptor Tyrosine Kinases", "SIGNOR"),
        ("ERBB2", "Erb-B2 Receptor Tyrosine Kinase 2", EntityKind::TransmembraneReceptor, SubcellularLocation::PlasmaMembrane, "Receptor Tyrosine Kinases", "SIGNOR"),
        ("BMPR2", "Bone Morphogenetic Protein Receptor Type 2", EntityKind::TransmembraneReceptor, SubcellularLocation::PlasmaMembrane, "TGF-beta / BMP", "SIGNOR"),
        ("LRP5", "LDL Receptor Related Protein 5", EntityKind::TransmembraneReceptor, SubcellularLocation::PlasmaMembrane, "Wnt Co-Receptors", "SIGNOR"),
        ("ICAM1", "Intercellular Adhesion Molecule 1", EntityKind::TransmembraneReceptor, SubcellularLocation::PlasmaMembrane, "Cell Adhesion", "Reactome"),
        ("VCAM1", "Vascular Cell Adhesion Molecule 1", EntityKind::TransmembraneReceptor, SubcellularLocation::PlasmaMembrane, "Cell Adhesion", "Reactome"),
        ("SELE", "Selectin E", EntityKind::TransmembraneReceptor, SubcellularLocation::PlasmaMembrane, "Cell Adhesion", "Reactome"),
        ("THBD", "Thrombomodulin", EntityKind::TransmembraneReceptor, SubcellularLocation::PlasmaMembrane, "Coagulation Homeostasis", "Reactome"),

        // Cytoplasmic kinases, enzymes & signaling complexes
        ("TRAF2", "TNF Receptor Associated Factor 2", EntityKind::GeneOther, SubcellularLocation::Cytoplasm, "Ubiquitin Signaling Adaptors", "SIGNOR"),
        ("MAP3K7", "TAK1 Mitogen-Activated Kinase", EntityKind::Kinase, SubcellularLocation::Cytoplasm, "MAPK Cascade", "SIGNOR"),
        ("IKBKB", "Inhibitor of Nuclear Factor Kappa B Kinase Beta", EntityKind::Kinase, SubcellularLocation::Cytoplasm, "IKK Complex", "SIGNOR"),
        ("NFKBIA", "NFKB Inhibitor Alpha (IkBa)", EntityKind::GeneOther, SubcellularLocation::Cytoplasm, "NF-kB Feedback", "SIGNOR"),
        ("TNFAIP3", "A20 Deubiquitinase (TNFAIP3)", EntityKind::Enzyme, SubcellularLocation::Cytoplasm, "NF-kB Feedback", "SIGNOR"),
        ("JAK1", "Janus Kinase 1", EntityKind::Kinase, SubcellularLocation::Cytoplasm, "JAK-STAT Signaling", "SIGNOR"),
        ("AKT1", "AKT Serine/Threonine Kinase 1", EntityKind::Kinase, SubcellularLocation::Cytoplasm, "PI3K-AKT Signaling", "SIGNOR"),
        ("GSK3B", "Glycogen Synthase Kinase 3 Beta", EntityKind::Kinase, SubcellularLocation::Cytoplasm, "Wnt Destruction Complex", "SIGNOR"),
        ("AXIN2", "Axin 2", EntityKind::GeneOther, SubcellularLocation::Cytoplasm, "Wnt Feedback", "WikiPathways"),
        ("PTGS2", "Prostaglandin-Endoperoxide Synthase 2 (COX-2)", EntityKind::Enzyme, SubcellularLocation::Cytoplasm, "Eicosanoid Synthesis", "Reactome"),
        ("NOS3", "Endothelial Nitric Oxide Synthase", EntityKind::Enzyme, SubcellularLocation::Cytoplasm, "Vascular Tone", "Reactome"),
        ("ALPL", "Alkaline Phosphatase, Biomineralization Associated", EntityKind::Enzyme, SubcellularLocation::Cytoplasm, "Osteoblast Maturation", "WikiPathways"),
        ("BIRC3", "Baculoviral IAP Repeat Containing 3", EntityKind::Enzyme, SubcellularLocation::Cytoplasm, "Apoptosis Regulation", "Reactome"),
        ("BCL2", "BCL2 Apoptosis Regulator", EntityKind::GeneOther, SubcellularLocation::Organelle, "Mitochondrial Survival", "SIGNOR"),
        ("BAX", "BCL2 Associated X, Apoptosis Regulator", EntityKind::GeneOther, SubcellularLocation::Organelle, "Mitochondrial Apoptosis", "SIGNOR"),
        ("SOD2", "Superoxide Dismutase 2", EntityKind::Enzyme, SubcellularLocation::Organelle, "Redox Homeostasis", "Reactome"),
        ("CTSD", "Cathepsin D", EntityKind::Enzyme, SubcellularLocation::Organelle, "Lysosomal Protease", "CollecTRI"),
        ("GREB1", "Growth Regulating Estrogen Receptor Binding 1", EntityKind::GeneOther, SubcellularLocation::Cytoplasm, "Estrogen Target", "CollecTRI"),

        // Nuclear transcription factors, co-regulators & cell-cycle machinery
        ("NFKB1", "Nuclear Factor Kappa B Subunit 1 (p50/RelA)", EntityKind::TranscriptionRegulator, SubcellularLocation::Nucleus, "NF-kB Family", "CollecTRI / SIGNOR"),
        ("RELA", "RELA Proto-Oncogene, NF-kB Subunit", EntityKind::TranscriptionRegulator, SubcellularLocation::Nucleus, "NF-kB Family", "CollecTRI / SIGNOR"),
        ("STAT1", "Signal Transducer and Activator of Transcription 1", EntityKind::TranscriptionRegulator, SubcellularLocation::Nucleus, "STAT Family", "CollecTRI / SIGNOR"),
        ("STAT3", "Signal Transducer and Activator of Transcription 3", EntityKind::TranscriptionRegulator, SubcellularLocation::Nucleus, "STAT Family", "CollecTRI / SIGNOR"),
        ("IRF1", "Interferon Regulatory Factor 1", EntityKind::TranscriptionRegulator, SubcellularLocation::Nucleus, "IRF Family", "CollecTRI"),
        ("ESR1", "Estrogen Receptor 1 (ER-alpha)", EntityKind::TranscriptionRegulator, SubcellularLocation::Nucleus, "Nuclear Hormone Receptors", "CollecTRI / SIGNOR"),
        ("NCOA3", "Nuclear Receptor Coactivator 3 (SRC-3)", EntityKind::TranscriptionRegulator, SubcellularLocation::Nucleus, "Nuclear Coactivators", "SIGNOR"),
        ("PGR", "Progesterone Receptor", EntityKind::TranscriptionRegulator, SubcellularLocation::Nucleus, "Nuclear Hormone Receptors", "CollecTRI"),
        ("NR3C1", "Glucocorticoid Receptor", EntityKind::TranscriptionRegulator, SubcellularLocation::Nucleus, "Nuclear Hormone Receptors", "CollecTRI / SIGNOR"),
        ("MYC", "MYC Proto-Oncogene, bHLH Transcription Factor", EntityKind::TranscriptionRegulator, SubcellularLocation::Nucleus, "Cell Cycle & Growth", "CollecTRI"),
        ("E2F1", "E2F Transcription Factor 1", EntityKind::TranscriptionRegulator, SubcellularLocation::Nucleus, "Cell Cycle & Growth", "CollecTRI"),
        ("FOXO3", "Forkhead Box O3", EntityKind::TranscriptionRegulator, SubcellularLocation::Nucleus, "Forkhead Transcription Factors", "CollecTRI"),
        ("XBP1", "X-Box Binding Protein 1", EntityKind::TranscriptionRegulator, SubcellularLocation::Nucleus, "Unfolded Protein Response", "CollecTRI"),
        ("CTNNB1", "Catenin Beta 1 (beta-Catenin)", EntityKind::TranscriptionRegulator, SubcellularLocation::Nucleus, "Wnt Transcriptional Effector", "CollecTRI / SIGNOR"),
        ("LEF1", "Lymphoid Enhancer Binding Factor 1", EntityKind::TranscriptionRegulator, SubcellularLocation::Nucleus, "TCF/LEF Family", "CollecTRI"),
        ("SMAD1", "SMAD Family Member 1", EntityKind::TranscriptionRegulator, SubcellularLocation::Nucleus, "BMP/SMAD Effectors", "CollecTRI / SIGNOR"),
        ("RUNX2", "RUNX Family Transcription Factor 2", EntityKind::TranscriptionRegulator, SubcellularLocation::Nucleus, "Master Osteogenic Regulators", "CollecTRI / WikiPathways"),
        ("SP7", "Sp7 Transcription Factor (Osterix)", EntityKind::TranscriptionRegulator, SubcellularLocation::Nucleus, "Master Osteogenic Regulators", "CollecTRI / WikiPathways"),
        ("PPARG", "Peroxisome Proliferator Activated Receptor Gamma", EntityKind::TranscriptionRegulator, SubcellularLocation::Nucleus, "Adipogenic & Anti-Inflammatory NR", "CollecTRI"),
        ("CEBPA", "CCAAT Enhancer Binding Protein Alpha", EntityKind::TranscriptionRegulator, SubcellularLocation::Nucleus, "Myeloid & Metabolic TFs", "CollecTRI"),
        ("KLF2", "Kruppel Like Factor 2", EntityKind::TranscriptionRegulator, SubcellularLocation::Nucleus, "Endothelial Quiescence", "CollecTRI"),
        ("TP53", "Tumor Protein p53", EntityKind::TranscriptionRegulator, SubcellularLocation::Nucleus, "Genome Integrity", "CollecTRI / SIGNOR"),
        ("CCND1", "Cyclin D1", EntityKind::GeneOther, SubcellularLocation::Nucleus, "G1/S Cell Cycle", "Reactome"),
        ("CDKN1A", "Cyclin Dependent Kinase Inhibitor 1A (p21)", EntityKind::GeneOther, SubcellularLocation::Nucleus, "Cell Cycle Checkpoint", "Reactome"),
        ("CDKN1B", "Cyclin Dependent Kinase Inhibitor 1B (p27)", EntityKind::GeneOther, SubcellularLocation::Nucleus, "Cell Cycle Checkpoint", "Reactome"),
        ("PCNA", "Proliferating Cell Nuclear Antigen", EntityKind::GeneOther, SubcellularLocation::Nucleus, "DNA Replication", "Reactome"),
        ("MKI67", "Marker of Proliferation Ki-67", EntityKind::GeneOther, SubcellularLocation::Nucleus, "Mitotic Proliferation", "Reactome"),

        // Downstream Biological Processes (GO) & Disease Phenotypes (Open Targets)
        ("GO:0050900", "Leukocyte Migration & Extravasation", EntityKind::BiologicalProcess, SubcellularLocation::PhenotypeLayer, "Immune Cell Trafficking", "Gene Ontology"),
        ("GO:0006954", "Acute Inflammatory Response", EntityKind::BiologicalProcess, SubcellularLocation::PhenotypeLayer, "Inflammatory & Immune Response", "Gene Ontology"),
        ("GO:0007155", "Endothelial Leukocyte Adhesion", EntityKind::BiologicalProcess, SubcellularLocation::PhenotypeLayer, "Cell-to-Cell Signaling", "Gene Ontology"),
        ("GO:0008283", "Epithelial & Tumor Cell Proliferation", EntityKind::BiologicalProcess, SubcellularLocation::PhenotypeLayer, "Cellular Growth & Proliferation", "Gene Ontology"),
        ("GO:0006915", "Apoptotic Cell Death Execution", EntityKind::BiologicalProcess, SubcellularLocation::PhenotypeLayer, "Cell Death & Survival", "Gene Ontology"),
        ("GO:0001649", "Osteoblast Differentiation & Mineralization", EntityKind::BiologicalProcess, SubcellularLocation::PhenotypeLayer, "Skeletal & Connective Tissue", "Gene Ontology"),
        ("GO:0045444", "Adipocyte Differentiation", EntityKind::BiologicalProcess, SubcellularLocation::PhenotypeLayer, "Lipid & Metabolic Homeostasis", "Gene Ontology"),
        ("EFO:0000305", "ER-Positive Breast Carcinoma Progression", EntityKind::DiseasePhenotype, SubcellularLocation::PhenotypeLayer, "Oncological Conditions", "Open Targets"),
        ("EFO:0003898", "Ankylosing Spondylitis Osteoproliferation", EntityKind::DiseasePhenotype, SubcellularLocation::PhenotypeLayer, "Skeletal & Rheumatic Disorders", "Open Targets"),
        ("EFO:0003767", "Vascular Endothelial Dysfunction", EntityKind::DiseasePhenotype, SubcellularLocation::PhenotypeLayer, "Cardiovascular Pathobiology", "Open Targets"),
    ];

    for &(id, name, kind, loc, cat, src) in reg_nodes {
        kg.upsert_node(EntityNode {
            id: id.to_string(),
            display_name: name.to_string(),
            kind,
            location: loc,
            family_category: cat.to_string(),
            external_id: None,
            provenance_source: src.to_string(),
        });
    }

    let add_edges = |kg: &mut KnowledgeGraph,
                     cls: EdgeClass,
                     items: &[(&str, &str, i8, f64, &str, &str, &str)]| {
        for &(src, dst, s, w, mech, db, pmid) in items {
            kg.add_edge(CausalEdge {
                source: src.to_string(),
                target: dst.to_string(),
                edge_class: cls,
                sign: EdgeSign::from_i8(s),
                weight: w,
                mechanism: mech.to_string(),
                source_db: db.to_string(),
                pmids: vec![pmid.to_string()],
            });
        }
    };

    // 1. Signaling Activation / Post-Translational Cascade Edges (E_a)
    let a_edges: &[(&str, &str, i8, f64, &str, &str, &str)] = &[
        ("TNF", "TNFRSF1A", 1, 0.98, "receptor binding", "SIGNOR", "12612050"),
        ("TNFRSF1A", "TRAF2", 1, 0.95, "complex recruitment", "SIGNOR", "8663075"),
        ("TRAF2", "MAP3K7", 1, 0.90, "polyubiquitination activation", "SIGNOR", "19815554"),
        ("MAP3K7", "IKBKB", 1, 0.95, "phosphorylation", "SIGNOR", "11460167"),
        ("IKBKB", "NFKBIA", -1, 0.98, "phosphorylation-induced degradation", "SIGNOR", "9346484"),
        ("NFKBIA", "NFKB1", -1, 0.96, "cytoplasmic sequestration", "SIGNOR", "8939934"),
        ("NFKBIA", "RELA", -1, 0.96, "cytoplasmic sequestration", "SIGNOR", "8939934"),
        ("IKBKB", "NFKB1", 1, 0.94, "canonical IKK release", "SIGNOR", "11460167"),
        ("IKBKB", "RELA", 1, 0.95, "Ser536 phosphorylation", "SIGNOR", "10597290"),
        ("TNF", "NFKB1", 1, 0.95, "cytokine cascade activation", "OmniPath", "16380449"),
        ("TNF", "RELA", 1, 0.95, "cytokine cascade activation", "OmniPath", "16380449"),
        ("TNF", "STAT1", 1, 0.80, "autocrine interferon loop", "OmniPath", "18329370"),
        ("LPS", "TLR4", 1, 0.98, "PAMP receptor activation", "SIGNOR", "9851930"),
        ("TLR4", "MAP3K7", 1, 0.92, "MyD88-TRAF6 cascade", "SIGNOR", "11460167"),
        ("LPS", "NFKB1", 1, 0.92, "innate TLR4 signaling", "OmniPath", "10430869"),
        ("IFNG", "IFNGR1", 1, 0.98, "ligand binding", "SIGNOR", "8197180"),
        ("IFNGR1", "JAK1", 1, 0.96, "receptor trans-activation", "SIGNOR", "8197180"),
        ("JAK1", "STAT1", 1, 0.97, "Tyr701 phosphorylation", "SIGNOR", "8197180"),
        ("IL6", "JAK1", 1, 0.92, "gp130 activation", "SIGNOR", "7518515"),
        ("JAK1", "STAT3", 1, 0.94, "Tyr705 phosphorylation", "SIGNOR", "7518515"),
        ("ESTRADIOL", "ESR1", 1, 0.99, "nuclear hormone agonist binding", "SIGNOR", "9342352"),
        ("ESR1", "NCOA3", 1, 0.92, "coactivator complex assembly", "SIGNOR", "9267029"),
        ("FULVESTRANT", "ESR1", -1, 0.99, "selective receptor degradation", "ChEMBL", "11507036"),
        ("DEXAMETHASONE", "NR3C1", 1, 0.99, "glucocorticoid agonist", "ChEMBL", "15987703"),
        ("NR3C1", "NFKB1", -1, 0.90, "transrepression", "SIGNOR", "10535948"),
        ("NR3C1", "RELA", -1, 0.90, "transrepression", "SIGNOR", "10535948"),
        ("AREG", "EGFR", 1, 0.92, "EGF-family ligand activation", "SIGNOR", "2491267"),
        ("EGFR", "AKT1", 1, 0.91, "PI3K recruitment", "SIGNOR", "10391249"),
        ("AKT1", "FOXO3", -1, 0.95, "Thr32 phosphorylation nuclear exclusion", "SIGNOR", "10102273"),
        // Multi-hop Wnt / SOST / BMP osteoproliferative signaling (demonstrates K=2 and K=3 CNA!)
        ("SOST", "LRP5", -1, 0.97, "extracellular LRP5/6 antagonism", "SIGNOR", "15878890"),
        ("SOST", "BMPR2", -1, 0.85, "BMP antagonism", "SIGNOR", "14657062"),
        ("SOST", "SMAD1", -1, 0.88, "BMP-SMAD inhibition", "SIGNOR", "14657062"),
        ("DKK1", "LRP5", -1, 0.96, "Wnt co-receptor blockade", "SIGNOR", "11357136"),
        ("WNT3A", "LRP5", 1, 0.96, "Frizzled-LRP5 co-activation", "SIGNOR", "11274176"),
        ("LRP5", "GSK3B", -1, 0.92, "destruction complex inhibition", "SIGNOR", "15260973"),
        ("GSK3B", "CTNNB1", -1, 0.97, "phosphorylation-proteasomal turnover", "SIGNOR", "9545243"),
        ("LRP5", "CTNNB1", 1, 0.94, "canonical Wnt stabilization", "SIGNOR", "15260973"),
        ("BMP2", "BMPR2", 1, 0.96, "receptor kinase activation", "SIGNOR", "8626525"),
        ("BMP6", "BMPR2", 1, 0.95, "receptor kinase activation", "SIGNOR", "8626525"),
        ("BMPR2", "SMAD1", 1, 0.96, "C-terminal phosphorylation", "SIGNOR", "8626525"),
        ("SMAD1", "RUNX2", 1, 0.93, "heteromeric osteogenic complex", "SIGNOR", "10995775"),
        ("CTNNB1", "RUNX2", 1, 0.91, "osteoblast lineage commitment", "WikiPathways", "16199551"),
        ("SMAD1", "EGFR", 1, 0.78, "cross-talk activation", "OmniPath", "22586108"),
    ];
    add_edges(&mut kg, EdgeClass::Activation, a_edges);

    // 2. Transcriptional & Expression Regulon Edges (E_t)
    let t_edges: &[(&str, &str, i8, f64, &str, &str, &str)] = &[
        // TNF / NF-kB / LPS / IFNG inflammatory regulons
        ("TNF", "ICAM1", 1, 0.98, "transcriptional induction", "SIGNOR", "1679050"),
        ("TNF", "VCAM1", 1, 0.98, "transcriptional induction", "SIGNOR", "1382068"),
        ("TNF", "SELE", 1, 0.98, "transcriptional induction", "SIGNOR", "1707827"),
        ("TNF", "CCL2", 1, 0.97, "transcriptional induction", "SIGNOR", "1903794"),
        ("TNF", "CXCL8", 1, 0.98, "transcriptional induction", "SIGNOR", "2663987"),
        ("TNF", "CXCL10", 1, 0.92, "transcriptional induction", "OmniPath", "8617890"),
        ("TNF", "IL6", 1, 0.97, "transcriptional induction", "SIGNOR", "2468689"),
        ("TNF", "IL1B", 1, 0.94, "transcriptional induction", "SIGNOR", "3487394"),
        ("TNF", "PTGS2", 1, 0.95, "transcriptional induction", "SIGNOR", "8144696"),
        ("TNF", "NFKBIA", 1, 0.96, "negative-feedback transcription", "SIGNOR", "8455628"),
        ("TNF", "TNFAIP3", 1, 0.97, "negative-feedback transcription", "SIGNOR", "2118878"),
        ("TNF", "BIRC3", 1, 0.93, "anti-apoptotic transcription", "SIGNOR", "8684456"),
        ("TNF", "SOD2", 1, 0.94, "mitochondrial antioxidant induction", "SIGNOR", "2847949"),
        ("TNF", "MMP9", 1, 0.92, "metalloproteinase induction", "SIGNOR", "8349689"),
        ("TNF", "CSF2", 1, 0.93, "cytokine induction", "SIGNOR", "2425322"),
        ("TNF", "IRF1", 1, 0.90, "TF induction", "OmniPath", "2545613"),
        ("TNF", "NOS3", -1, 0.93, "mRNA destabilization & repression", "SIGNOR", "7504900"),
        ("TNF", "KLF2", -1, 0.94, "endothelial repression", "OmniPath", "15138255"),
        ("TNF", "THBD", -1, 0.92, "endothelial repression", "SIGNOR", "2836506"),
        ("TNF", "PPARG", -1, 0.89, "nuclear receptor down-regulation", "OmniPath", "10593901"),

        ("NFKB1", "ICAM1", 1, 0.99, "direct promoter binding", "CollecTRI", "7512965"),
        ("NFKB1", "VCAM1", 1, 0.99, "direct promoter binding", "CollecTRI", "1382068"),
        ("NFKB1", "SELE", 1, 0.98, "direct promoter binding", "CollecTRI", "1707827"),
        ("NFKB1", "CCL2", 1, 0.97, "direct promoter binding", "CollecTRI", "8163541"),
        ("NFKB1", "CXCL8", 1, 0.99, "direct promoter binding", "CollecTRI", "8083207"),
        ("NFKB1", "IL6", 1, 0.98, "direct promoter binding", "CollecTRI", "2185423"),
        ("NFKB1", "IL1B", 1, 0.96, "direct promoter binding", "CollecTRI", "8264618"),
        ("NFKB1", "PTGS2", 1, 0.96, "direct promoter binding", "CollecTRI", "9285428"),
        ("NFKB1", "NFKBIA", 1, 0.99, "autoregulatory transcription", "CollecTRI", "8455628"),
        ("NFKB1", "TNFAIP3", 1, 0.98, "direct promoter binding", "CollecTRI", "9689070"),
        ("NFKB1", "BIRC3", 1, 0.95, "direct promoter binding", "CollecTRI", "9799235"),
        ("NFKB1", "SOD2", 1, 0.94, "intronic enhancer activation", "CollecTRI", "10330159"),
        ("NFKB1", "MMP9", 1, 0.94, "direct promoter binding", "CollecTRI", "8349689"),
        ("NFKB1", "CSF2", 1, 0.93, "direct promoter binding", "CollecTRI", "8474595"),
        ("NFKB1", "CCND1", 1, 0.88, "direct promoter binding", "CollecTRI", "10490625"),
        ("NFKB1", "MYC", 1, 0.86, "direct promoter binding", "CollecTRI", "2196409"),

        ("RELA", "ICAM1", 1, 0.99, "transactivation", "CollecTRI", "7512965"),
        ("RELA", "VCAM1", 1, 0.99, "transactivation", "CollecTRI", "1382068"),
        ("RELA", "SELE", 1, 0.98, "transactivation", "CollecTRI", "1707827"),
        ("RELA", "CXCL8", 1, 0.99, "transactivation", "CollecTRI", "8083207"),
        ("RELA", "CCL2", 1, 0.97, "transactivation", "CollecTRI", "8163541"),
        ("RELA", "IL6", 1, 0.98, "transactivation", "CollecTRI", "2185423"),
        ("RELA", "NFKBIA", 1, 0.99, "transactivation", "CollecTRI", "8455628"),
        ("RELA", "TNFAIP3", 1, 0.98, "transactivation", "CollecTRI", "9689070"),
        ("RELA", "PTGS2", 1, 0.95, "transactivation", "CollecTRI", "9285428"),
        ("RELA", "NOS3", -1, 0.85, "endothelial repression", "OmniPath", "15138255"),

        ("LPS", "ICAM1", 1, 0.92, "TLR4 transcriptional response", "SIGNOR", "10430869"),
        ("LPS", "VCAM1", 1, 0.92, "TLR4 transcriptional response", "SIGNOR", "10430869"),
        ("LPS", "CCL2", 1, 0.94, "TLR4 transcriptional response", "SIGNOR", "10430869"),
        ("LPS", "CXCL8", 1, 0.95, "TLR4 transcriptional response", "SIGNOR", "10430869"),
        ("LPS", "IL6", 1, 0.96, "TLR4 transcriptional response", "SIGNOR", "10430869"),
        ("LPS", "IL1B", 1, 0.96, "TLR4 transcriptional response", "SIGNOR", "10430869"),
        ("LPS", "PTGS2", 1, 0.95, "TLR4 transcriptional response", "SIGNOR", "10430869"),
        ("LPS", "CXCL10", 1, 0.91, "TRIF-IRF3 induction", "SIGNOR", "11714794"),

        ("IFNG", "CXCL10", 1, 0.99, "ISRE/GAS activation", "SIGNOR", "2111842"),
        ("IFNG", "IRF1", 1, 0.99, "GAS promoter activation", "SIGNOR", "8402893"),
        ("IFNG", "STAT1", 1, 0.96, "positive autoregulation", "SIGNOR", "8197180"),
        ("IFNG", "ICAM1", 1, 0.94, "GAS promoter activation", "SIGNOR", "7512965"),
        ("IFNG", "CCL2", 1, 0.90, "chemokine induction", "SIGNOR", "8163541"),

        ("STAT1", "CXCL10", 1, 0.98, "direct promoter binding", "CollecTRI", "9488479"),
        ("STAT1", "IRF1", 1, 0.99, "direct promoter binding", "CollecTRI", "8402893"),
        ("STAT1", "ICAM1", 1, 0.93, "direct promoter binding", "CollecTRI", "7512965"),
        ("STAT1", "CCL2", 1, 0.89, "direct promoter binding", "CollecTRI", "9488479"),

        ("DEXAMETHASONE", "IL6", -1, 0.95, "glucocorticoid repression", "ChEMBL", "15987703"),
        ("DEXAMETHASONE", "CXCL8", -1, 0.95, "glucocorticoid repression", "ChEMBL", "15987703"),
        ("DEXAMETHASONE", "CCL2", -1, 0.93, "glucocorticoid repression", "ChEMBL", "15987703"),
        ("DEXAMETHASONE", "PTGS2", -1, 0.94, "glucocorticoid repression", "ChEMBL", "15987703"),
        ("DEXAMETHASONE", "SELE", -1, 0.91, "glucocorticoid repression", "ChEMBL", "15987703"),
        ("DEXAMETHASONE", "ICAM1", -1, 0.90, "glucocorticoid repression", "ChEMBL", "15987703"),
        ("DEXAMETHASONE", "MMP9", -1, 0.89, "glucocorticoid repression", "ChEMBL", "15987703"),

        ("KLF2", "NOS3", 1, 0.96, "direct promoter transactivation", "CollecTRI", "15051742"),
        ("KLF2", "THBD", 1, 0.95, "direct promoter transactivation", "CollecTRI", "15713798"),
        ("KLF2", "VCAM1", -1, 0.92, "anti-inflammatory repression", "CollecTRI", "15051742"),
        ("KLF2", "SELE", -1, 0.92, "anti-inflammatory repression", "CollecTRI", "15051742"),

        // ESR1 / Estradiol / Fulvestrant hormone receptor regulons
        ("ESR1", "TFF1", 1, 0.99, "canonical ERE transactivation", "CollecTRI", "17369815"),
        ("ESR1", "GREB1", 1, 0.99, "canonical ERE transactivation", "CollecTRI", "16009731"),
        ("ESR1", "PGR", 1, 0.98, "canonical ERE transactivation", "CollecTRI", "12618429"),
        ("ESR1", "CXCL12", 1, 0.95, "ERE transactivation", "CollecTRI", "14551289"),
        ("ESR1", "MYC", 1, 0.96, "enhancer activation", "CollecTRI", "1984316"),
        ("ESR1", "CCND1", 1, 0.96, "mitogenic induction", "CollecTRI", "9482727"),
        ("ESR1", "AREG", 1, 0.94, "paracrine mitogen induction", "CollecTRI", "7530338"),
        ("ESR1", "BCL2", 1, 0.93, "anti-apoptotic ERE induction", "CollecTRI", "8521391"),
        ("ESR1", "CTSD", 1, 0.94, "canonical ERE transactivation", "CollecTRI", "2566996"),
        ("ESR1", "XBP1", 1, 0.91, "luminal transcriptional program", "CollecTRI", "15665279"),
        ("ESR1", "E2F1", 1, 0.90, "cell-cycle progression", "CollecTRI", "10490633"),
        ("ESR1", "CDKN1A", -1, 0.92, "transcriptional repression", "CollecTRI", "15642340"),
        ("ESR1", "CDKN1B", -1, 0.89, "cell-cycle brake repression", "CollecTRI", "11416145"),
        ("ESR1", "FOXO3", -1, 0.86, "pro-survival repression", "OmniPath", "16954366"),
        ("ESR1", "ERBB2", -1, 0.90, "intronic enhancer repression", "CollecTRI", "18997789"),
        ("ESR1", "IL6", -1, 0.88, "NF-kB transrepression", "CollecTRI", "7544006"),

        ("ESTRADIOL", "TFF1", 1, 0.99, "ER-alpha agonist induction", "SIGNOR", "17369815"),
        ("ESTRADIOL", "GREB1", 1, 0.99, "ER-alpha agonist induction", "SIGNOR", "16009731"),
        ("ESTRADIOL", "PGR", 1, 0.98, "ER-alpha agonist induction", "SIGNOR", "12618429"),
        ("ESTRADIOL", "CXCL12", 1, 0.95, "ER-alpha agonist induction", "SIGNOR", "14551289"),
        ("ESTRADIOL", "MYC", 1, 0.96, "ER-alpha agonist induction", "SIGNOR", "1984316"),
        ("ESTRADIOL", "CCND1", 1, 0.95, "ER-alpha agonist induction", "SIGNOR", "9482727"),
        ("ESTRADIOL", "AREG", 1, 0.93, "ER-alpha agonist induction", "SIGNOR", "7530338"),
        ("ESTRADIOL", "BCL2", 1, 0.92, "ER-alpha agonist induction", "SIGNOR", "8521391"),
        ("ESTRADIOL", "CTSD", 1, 0.93, "ER-alpha agonist induction", "SIGNOR", "2566996"),
        ("ESTRADIOL", "CDKN1A", -1, 0.91, "G1 checkpoint suppression", "SIGNOR", "15642340"),

        ("FULVESTRANT", "TFF1", -1, 0.98, "SERD antagonism", "ChEMBL", "11507036"),
        ("FULVESTRANT", "GREB1", -1, 0.98, "SERD antagonism", "ChEMBL", "11507036"),
        ("FULVESTRANT", "PGR", -1, 0.97, "SERD antagonism", "ChEMBL", "11507036"),
        ("FULVESTRANT", "MYC", -1, 0.94, "SERD antagonism", "ChEMBL", "11507036"),
        ("FULVESTRANT", "CCND1", -1, 0.94, "SERD antagonism", "ChEMBL", "11507036"),
        ("FULVESTRANT", "CDKN1A", 1, 0.90, "cell-cycle arrest derepression", "ChEMBL", "11507036"),

        ("NCOA3", "TFF1", 1, 0.92, "coactivation", "CollecTRI", "9267029"),
        ("NCOA3", "GREB1", 1, 0.92, "coactivation", "CollecTRI", "9267029"),
        ("NCOA3", "CCND1", 1, 0.90, "coactivation", "CollecTRI", "10567357"),
        ("NCOA3", "MYC", 1, 0.89, "coactivation", "CollecTRI", "10567357"),

        ("E2F1", "PCNA", 1, 0.97, "S-phase transactivation", "CollecTRI", "8622684"),
        ("E2F1", "MKI67", 1, 0.95, "mitotic gene induction", "CollecTRI", "15199124"),
        ("E2F1", "CCND1", 1, 0.88, "cell-cycle feedback", "CollecTRI", "9053844"),
        ("TP53", "CDKN1A", 1, 0.99, "direct p53 response element", "CollecTRI", "8242752"),
        ("TP53", "BAX", 1, 0.97, "pro-apoptotic transactivation", "CollecTRI", "7874169"),
        ("TP53", "BCL2", -1, 0.92, "transcriptional repression", "CollecTRI", "8183579"),
        ("TP53", "MKI67", -1, 0.88, "DREAM complex repression", "CollecTRI", "25258670"),

        // Wnt / CTNNB1 / SMAD1 / RUNX2 osteogenic regulons
        ("CTNNB1", "AXIN2", 1, 0.99, "TCF/LEF direct target", "CollecTRI", "11850410"),
        ("CTNNB1", "LEF1", 1, 0.96, "TCF/LEF positive feedback", "CollecTRI", "11564747"),
        ("CTNNB1", "MYC", 1, 0.97, "Wnt enhancer activation", "CollecTRI", "9738008"),
        ("CTNNB1", "CCND1", 1, 0.96, "TCF binding element", "CollecTRI", "10097128"),
        ("CTNNB1", "RUNX2", 1, 0.95, "direct promoter activation", "CollecTRI", "16199551"),
        ("CTNNB1", "BMP2", 1, 0.89, "osteoblast induction", "WikiPathways", "17412763"),
        ("CTNNB1", "PPARG", -1, 0.94, "adipogenesis suppression", "CollecTRI", "12006603"),
        ("CTNNB1", "CEBPA", -1, 0.91, "adipogenesis suppression", "CollecTRI", "12006603"),

        ("RUNX2", "SP7", 1, 0.99, "master osteoblast cascade", "CollecTRI", "11792318"),
        ("RUNX2", "ALPL", 1, 0.98, "OSE2 promoter activation", "CollecTRI", "9179270"),
        ("RUNX2", "BGLAP", 1, 0.99, "osteocalcin OSE2 activation", "CollecTRI", "9179270"),
        ("RUNX2", "COL1A1", 1, 0.96, "bone matrix transcription", "CollecTRI", "10869395"),
        ("RUNX2", "MMP13", 1, 0.96, "hypertrophic remodeling", "CollecTRI", "10446182"),
        ("RUNX2", "BMP2", 1, 0.90, "autocrine BMP loop", "WikiPathways", "15210659"),
        ("RUNX2", "CDKN1A", -1, 0.85, "osteoprogenitor expansion", "CollecTRI", "16537904"),

        ("SMAD1", "RUNX2", 1, 0.95, "BMP-responsive enhancer", "CollecTRI", "10995775"),
        ("SMAD1", "SP7", 1, 0.93, "BMP osteogenic target", "CollecTRI", "12855683"),
        ("SMAD1", "ALPL", 1, 0.92, "BMP osteogenic target", "CollecTRI", "10995775"),
        ("SMAD1", "BMP6", 1, 0.88, "BMP positive feedback", "WikiPathways", "15210659"),
        ("EGFR", "LCN2", 1, 0.87, "epithelial-inflammatory target", "OmniPath", "22586108"),
        ("EGFR", "PTGS2", 1, 0.89, "mitogenic COX-2 induction", "SIGNOR", "10523663"),
    ];
    add_edges(&mut kg, EdgeClass::Transcription, t_edges);

    // 3. Downstream Phenotype & Disease Edges (E_p)
    let p_edges: &[(&str, &str, i8, f64, &str, &str, &str)] = &[
        ("ICAM1", "GO:0050900", 1, 0.96, "leukocyte firm arrest", "Gene Ontology", "1679050"),
        ("VCAM1", "GO:0050900", 1, 0.96, "VLA-4 transmigration", "Gene Ontology", "1382068"),
        ("SELE", "GO:0050900", 1, 0.96, "endothelial tethering & rolling", "Gene Ontology", "1707827"),
        ("CCL2", "GO:0050900", 1, 0.97, "monocyte chemotaxis", "Gene Ontology", "1903794"),
        ("CXCL8", "GO:0050900", 1, 0.98, "neutrophil chemotaxis", "Gene Ontology", "2663987"),
        ("CXCL10", "GO:0050900", 1, 0.95, "Th1 lymphocyte recruitment", "Gene Ontology", "8617890"),
        ("MMP9", "GO:0050900", 1, 0.90, "basement membrane diapedesis", "Gene Ontology", "8349689"),

        ("IL6", "GO:0006954", 1, 0.98, "acute-phase cytokine", "Gene Ontology", "2468689"),
        ("IL1B", "GO:0006954", 1, 0.99, "pyrogenic cytokine", "Gene Ontology", "3487394"),
        ("PTGS2", "GO:0006954", 1, 0.96, "prostaglandin E2 synthesis", "Gene Ontology", "8144696"),
        ("CXCL8", "GO:0006954", 1, 0.95, "inflammatory mediator", "Gene Ontology", "2663987"),
        ("CCL2", "GO:0006954", 1, 0.94, "inflammatory chemokine", "Gene Ontology", "1903794"),
        ("PPARG", "GO:0006954", -1, 0.90, "anti-inflammatory resolution", "Gene Ontology", "9419369"),
        ("KLF2", "GO:0006954", -1, 0.91, "endothelial anti-inflammatory tone", "Gene Ontology", "15051742"),

        ("ICAM1", "GO:0007155", 1, 0.98, "integrin adhesion", "Gene Ontology", "1679050"),
        ("VCAM1", "GO:0007155", 1, 0.98, "integrin adhesion", "Gene Ontology", "1382068"),
        ("SELE", "GO:0007155", 1, 0.97, "selectin adhesion", "Gene Ontology", "1707827"),

        ("NOS3", "EFO:0003767", -1, 0.97, "protective nitric oxide bioavailability", "Open Targets", "9576111"),
        ("KLF2", "EFO:0003767", -1, 0.95, "laminar shear protection", "Open Targets", "15051742"),
        ("THBD", "EFO:0003767", -1, 0.93, "anticoagulant endothelial surface", "Open Targets", "2836506"),
        ("VCAM1", "EFO:0003767", 1, 0.92, "atherogenic activation", "Open Targets", "1382068"),
        ("ICAM1", "EFO:0003767", 1, 0.91, "endothelial activation", "Open Targets", "1679050"),
        ("PTGS2", "EFO:0003767", 1, 0.89, "vascular prostanoid stress", "Open Targets", "8144696"),

        ("MYC", "GO:0008283", 1, 0.98, "ribosome & cell-cycle drive", "Gene Ontology", "1984316"),
        ("CCND1", "GO:0008283", 1, 0.98, "CDK4/6 G1-S progression", "Gene Ontology", "9482727"),
        ("E2F1", "GO:0008283", 1, 0.96, "S-phase DNA synthesis", "Gene Ontology", "8622684"),
        ("PCNA", "GO:0008283", 1, 0.94, "processive replication", "Gene Ontology", "8622684"),
        ("MKI67", "GO:0008283", 1, 0.95, "mitotic chromosome dispersion", "Gene Ontology", "15199124"),
        ("AREG", "GO:0008283", 1, 0.92, "EGFR mitogenic stimulation", "Gene Ontology", "7530338"),
        ("CDKN1A", "GO:0008283", -1, 0.97, "CDK2/4 arrest", "Gene Ontology", "8242752"),
        ("CDKN1B", "GO:0008283", -1, 0.96, "G1 CDK inhibition", "Gene Ontology", "11416145"),
        ("FOXO3", "GO:0008283", -1, 0.90, "antiproliferative program", "Gene Ontology", "10102273"),

        ("TFF1", "EFO:0000305", 1, 0.92, "luminal estrogen growth marker", "Open Targets", "17369815"),
        ("GREB1", "EFO:0000305", 1, 0.96, "estrogen-dependent proliferation", "Open Targets", "16009731"),
        ("PGR", "EFO:0000305", 1, 0.90, "hormone-responsive luminal axis", "Open Targets", "12618429"),
        ("CXCL12", "EFO:0000305", 1, 0.91, "CXCR4 autocrine survival", "Open Targets", "14551289"),
        ("MYC", "EFO:0000305", 1, 0.95, "oncogenic amplification/drive", "Open Targets", "1984316"),
        ("CCND1", "EFO:0000305", 1, 0.95, "luminal cyclin D1 addiction", "Open Targets", "9482727"),
        ("BCL2", "EFO:0000305", 1, 0.91, "apoptosis evasion", "Open Targets", "8521391"),
        ("CDKN1A", "EFO:0000305", -1, 0.89, "cell-cycle arrest", "Open Targets", "15642340"),

        ("BAX", "GO:0006915", 1, 0.98, "MOMP pore formation", "Gene Ontology", "7874169"),
        ("FOXO3", "GO:0006915", 1, 0.90, "BIM/PUMA induction", "Gene Ontology", "10102273"),
        ("BCL2", "GO:0006915", -1, 0.98, "BAX/BAK sequestration", "Gene Ontology", "8521391"),
        ("BIRC3", "GO:0006915", -1, 0.93, "caspase inhibition", "Gene Ontology", "8684456"),

        ("RUNX2", "GO:0001649", 1, 0.99, "osteoblast lineage specification", "Gene Ontology", "9179270"),
        ("SP7", "GO:0001649", 1, 0.99, "mature osteoblast commitment", "Gene Ontology", "11792318"),
        ("ALPL", "GO:0001649", 1, 0.97, "hydroxyapatite matrix mineralization", "Gene Ontology", "9179270"),
        ("BGLAP", "GO:0001649", 1, 0.97, "bone matrix mineralization", "Gene Ontology", "9179270"),
        ("COL1A1", "GO:0001649", 1, 0.95, "osteoid scaffold synthesis", "Gene Ontology", "10869395"),
        ("BMP2", "GO:0001649", 1, 0.96, "osteogenic induction", "Gene Ontology", "8626525"),
        ("BMP6", "GO:0001649", 1, 0.94, "osteogenic induction", "Gene Ontology", "8626525"),

        ("RUNX2", "EFO:0003898", 1, 0.95, "syndesmophyte bone formation", "Open Targets", "23136895"),
        ("SP7", "EFO:0003898", 1, 0.94, "ectopic spinal ossification", "Open Targets", "23136895"),
        ("ALPL", "EFO:0003898", 1, 0.92, "pathologic mineralization", "Open Targets", "23136895"),
        ("BMP2", "EFO:0003898", 1, 0.91, "entheseal bone remodeling", "Open Targets", "23136895"),
        ("AXIN2", "EFO:0003898", 1, 0.89, "hyperactive Wnt signaling", "Open Targets", "23136895"),
        ("LCN2", "EFO:0003898", 1, 0.88, "ankylosing spondylitis biomarker", "Open Targets", "23136895"),
        ("MMP13", "EFO:0003898", 1, 0.88, "chondrocyte hypertrophy", "Open Targets", "23136895"),

        ("PPARG", "GO:0045444", 1, 0.99, "master adipogenic nuclear receptor", "Gene Ontology", "8178145"),
        ("CEBPA", "GO:0045444", 1, 0.98, "adipocyte terminal differentiation", "Gene Ontology", "8178145"),
        ("AXIN2", "GO:0045444", -1, 0.85, "Wnt anti-adipogenic brake", "Gene Ontology", "12006603"),
    ];
    add_edges(&mut kg, EdgeClass::Phenotype, p_edges);

    // 4. Canonical Pathways with Topological Reaction Edges (for MAP simulation & interactive diagram)
    let make_signs = |pairs: &[(&str, i8)]| -> BTreeMap<String, i8> {
        pairs.iter().map(|&(k, v)| (k.to_string(), v)).collect()
    };
    let make_coords = |pts: &[(&str, f32, f32)]| -> BTreeMap<String, (f32, f32)> {
        pts.iter().map(|&(k, x, y)| (k.to_string(), (x, y))).collect()
    };
    let make_topo = |triples: &[(&str, &str, i8, &str)]| -> Vec<PathwayReactionEdge> {
        triples
            .iter()
            .map(|&(a, b, s, lbl)| PathwayReactionEdge {
                from: a.to_string(),
                to: b.to_string(),
                sign: EdgeSign::from_i8(s),
                interaction_label: lbl.to_string(),
            })
            .collect()
    };

    kg.pathways = vec![
        PathwayDef {
            id: "R-HSA-75893".to_string(),
            name: "TNF Signaling & Canonical NF-kB Cascade".to_string(),
            source_db: "Reactome".to_string(),
            category: "Immune System".to_string(),
            url: "https://reactome.org/content/detail/R-HSA-75893".to_string(),
            members: vec![
                "TNF", "TNFRSF1A", "TRAF2", "MAP3K7", "IKBKB", "NFKBIA", "TNFAIP3", "NFKB1",
                "RELA", "ICAM1", "VCAM1", "SELE", "CCL2", "CXCL8", "IL6", "IL1B", "PTGS2",
                "BIRC3", "SOD2", "MMP9", "CSF2",
            ]
            .into_iter()
            .map(String::from)
            .collect(),
            member_expected_signs: make_signs(&[
                ("ICAM1", 1), ("VCAM1", 1), ("SELE", 1), ("CCL2", 1), ("CXCL8", 1),
                ("IL6", 1), ("IL1B", 1), ("PTGS2", 1), ("NFKBIA", 1), ("TNFAIP3", 1),
                ("BIRC3", 1), ("SOD2", 1), ("MMP9", 1), ("CSF2", 1),
            ]),
            node_coords: make_coords(&[
                ("TNF", 0.48, 0.06), ("TNFRSF1A", 0.48, 0.18), ("TRAF2", 0.48, 0.30),
                ("TNFAIP3", 0.24, 0.34), ("MAP3K7", 0.48, 0.42), ("IKBKB", 0.48, 0.53),
                ("NFKBIA", 0.27, 0.58), ("BIRC3", 0.82, 0.50), ("SOD2", 0.82, 0.62),
                ("PTGS2", 0.68, 0.56), ("MMP9", 0.16, 0.07), ("CSF2", 0.72, 0.07),
                ("IL1B", 0.30, 0.07), ("SELE", 0.20, 0.18),
                ("RELA", 0.38, 0.80), ("NFKB1", 0.58, 0.80),
                ("ICAM1", 0.18, 0.85), ("VCAM1", 0.30, 0.85), ("CXCL8", 0.48, 0.85),
                ("CCL2", 0.68, 0.85), ("IL6", 0.84, 0.85),
            ]),
            phospho_nodes: vec!["IKBKB".to_string(), "NFKBIA".to_string(), "RELA".to_string()],
            cell_context_label: "Vascular endothelial cell".to_string(),
            topology_edges: make_topo(&[
                ("TNF", "TNFRSF1A", 1, "binding"),
                ("TNFRSF1A", "TRAF2", 1, "recruitment"),
                ("TRAF2", "MAP3K7", 1, "K63-Ub"),
                ("MAP3K7", "IKBKB", 1, "phospho"),
                ("IKBKB", "NFKBIA", -1, "degradation"),
                ("NFKBIA", "RELA", -1, "inhibition"),
                ("NFKBIA", "NFKB1", -1, "inhibition"),
                ("IKBKB", "RELA", 1, "phospho"),
                ("RELA", "ICAM1", 1, "transcription"),
                ("RELA", "VCAM1", 1, "transcription"),
                ("RELA", "SELE", 1, "transcription"),
                ("RELA", "CXCL8", 1, "transcription"),
                ("NFKB1", "CCL2", 1, "transcription"),
                ("NFKB1", "IL6", 1, "transcription"),
                ("NFKB1", "PTGS2", 1, "transcription"),
                ("NFKB1", "TNFAIP3", 1, "feedback"),
                ("TNFAIP3", "TRAF2", -1, "deubiquitination"),
            ]),
        },
        PathwayDef {
            id: "WP445".to_string(),
            name: "Endothelial Cell Activation & Leukocyte Diapedesis".to_string(),
            source_db: "WikiPathways".to_string(),
            category: "Cardiovascular & Vascular Biology".to_string(),
            url: "https://www.wikipathways.org/pathways/WP445.html".to_string(),
            members: vec![
                "TNF", "IL1B", "ICAM1", "VCAM1", "SELE", "CCL2", "CXCL8", "CXCL10", "MMP9",
                "NOS3", "KLF2", "THBD", "PTGS2",
            ]
            .into_iter()
            .map(String::from)
            .collect(),
            member_expected_signs: make_signs(&[
                ("ICAM1", 1), ("VCAM1", 1), ("SELE", 1), ("CCL2", 1), ("CXCL8", 1),
                ("CXCL10", 1), ("MMP9", 1), ("PTGS2", 1), ("NOS3", -1), ("KLF2", -1), ("THBD", -1),
            ]),
            node_coords: make_coords(&[
                ("IL1B", 0.24, 0.06), ("TNF", 0.54, 0.06), ("MMP9", 0.82, 0.07),
                ("SELE", 0.20, 0.18), ("ICAM1", 0.36, 0.18), ("VCAM1", 0.52, 0.18), ("THBD", 0.76, 0.18),
                ("PTGS2", 0.38, 0.46), ("NOS3", 0.74, 0.50),
                ("KLF2", 0.64, 0.81), ("CXCL8", 0.24, 0.84), ("CCL2", 0.42, 0.84), ("CXCL10", 0.82, 0.84),
            ]),
            phospho_nodes: vec!["NOS3".to_string()],
            cell_context_label: "Microvascular endothelium".to_string(),
            topology_edges: make_topo(&[
                ("TNF", "KLF2", -1, "repression"),
                ("KLF2", "NOS3", 1, "transcription"),
                ("KLF2", "THBD", 1, "transcription"),
                ("KLF2", "VCAM1", -1, "repression"),
                ("IL1B", "SELE", 1, "induction"),
                ("TNF", "ICAM1", 1, "induction"),
                ("TNF", "MMP9", 1, "induction"),
            ]),
        },
        PathwayDef {
            id: "R-HSA-913531".to_string(),
            name: "Interferon-Gamma & JAK-STAT1 Antiviral Signaling".to_string(),
            source_db: "Reactome".to_string(),
            category: "Immune System".to_string(),
            url: "https://reactome.org/content/detail/R-HSA-913531".to_string(),
            members: vec![
                "IFNG", "IFNGR1", "JAK1", "STAT1", "STAT3", "IRF1", "CXCL10", "ICAM1", "CCL2",
            ]
            .into_iter()
            .map(String::from)
            .collect(),
            member_expected_signs: make_signs(&[
                ("STAT1", 1), ("IRF1", 1), ("CXCL10", 1), ("ICAM1", 1), ("CCL2", 1),
            ]),
            node_coords: make_coords(&[
                ("IFNG", 0.46, 0.06), ("CXCL10", 0.76, 0.07), ("CCL2", 0.88, 0.07),
                ("IFNGR1", 0.46, 0.18), ("ICAM1", 0.20, 0.18),
                ("JAK1", 0.46, 0.33), ("STAT1", 0.46, 0.52), ("STAT3", 0.66, 0.50),
                ("IRF1", 0.46, 0.82),
            ]),
            phospho_nodes: vec!["JAK1".to_string(), "STAT1".to_string(), "STAT3".to_string()],
            cell_context_label: "Endothelial / Glial cell".to_string(),
            topology_edges: make_topo(&[
                ("IFNG", "IFNGR1", 1, "binding"),
                ("IFNGR1", "JAK1", 1, "activation"),
                ("JAK1", "STAT1", 1, "Tyr701-P"),
                ("STAT1", "IRF1", 1, "transcription"),
                ("STAT1", "CXCL10", 1, "transcription"),
                ("IRF1", "ICAM1", 1, "co-activation"),
            ]),
        },
        PathwayDef {
            id: "R-HSA-8939211".to_string(),
            name: "ESR1 (Estrogen-Dependent) Nuclear Gene Expression".to_string(),
            source_db: "Reactome".to_string(),
            category: "Endocrine & Nuclear Receptor Signaling".to_string(),
            url: "https://reactome.org/content/detail/R-HSA-8939211".to_string(),
            members: vec![
                "ESTRADIOL", "ESR1", "NCOA3", "TFF1", "GREB1", "PGR", "CXCL12", "MYC",
                "CCND1", "AREG", "BCL2", "CTSD", "XBP1", "E2F1", "PCNA", "MKI67", "CDKN1A",
                "CDKN1B", "FOXO3", "ERBB2",
            ]
            .into_iter()
            .map(String::from)
            .collect(),
            member_expected_signs: make_signs(&[
                ("TFF1", 1), ("GREB1", 1), ("PGR", 1), ("CXCL12", 1), ("MYC", 1),
                ("CCND1", 1), ("AREG", 1), ("BCL2", 1), ("CTSD", 1), ("XBP1", 1),
                ("E2F1", 1), ("PCNA", 1), ("MKI67", 1), ("CDKN1A", -1), ("CDKN1B", -1),
                ("FOXO3", -1), ("ERBB2", -1),
            ]),
            node_coords: make_coords(&[
                ("ESTRADIOL", 0.42, 0.06), ("AREG", 0.18, 0.06), ("CXCL12", 0.76, 0.06),
                ("EGFR", 0.18, 0.18), ("ERBB2", 0.30, 0.18),
                ("GREB1", 0.42, 0.42), ("BCL2", 0.64, 0.46), ("CTSD", 0.80, 0.46), ("FOXO3", 0.20, 0.54), ("XBP1", 0.52, 0.56),
                ("ESR1", 0.40, 0.80), ("NCOA3", 0.52, 0.80), ("TFF1", 0.15, 0.84), ("PGR", 0.27, 0.84),
                ("MYC", 0.65, 0.79), ("CCND1", 0.77, 0.78), ("CDKN1A", 0.89, 0.78), ("CDKN1B", 0.89, 0.86),
                ("E2F1", 0.74, 0.86), ("PCNA", 0.62, 0.89), ("MKI67", 0.82, 0.91),
            ]),
            phospho_nodes: vec!["ESR1".to_string(), "FOXO3".to_string()],
            cell_context_label: "Luminal epithelial cell".to_string(),
            topology_edges: make_topo(&[
                ("ESTRADIOL", "ESR1", 1, "ligand agonist"),
                ("ESR1", "NCOA3", 1, "coactivator"),
                ("ESR1", "TFF1", 1, "ERE"),
                ("ESR1", "GREB1", 1, "ERE"),
                ("ESR1", "PGR", 1, "ERE"),
                ("ESR1", "MYC", 1, "enhancer"),
                ("ESR1", "CCND1", 1, "mitogenic"),
                ("ESR1", "AREG", 1, "ERE"),
                ("AREG", "EGFR", 1, "paracrine"),
                ("CCND1", "E2F1", 1, "Rb release"),
                ("E2F1", "PCNA", 1, "S-phase"),
                ("E2F1", "MKI67", 1, "mitosis"),
                ("ESR1", "CDKN1A", -1, "repression"),
                ("ESR1", "CDKN1B", -1, "repression"),
            ]),
        },
        PathwayDef {
            id: "WP399".to_string(),
            name: "Canonical Wnt / Beta-Catenin & BMP Osteoblastogenesis".to_string(),
            source_db: "WikiPathways".to_string(),
            category: "Developmental & Skeletal Signaling".to_string(),
            url: "https://www.wikipathways.org/pathways/WP399.html".to_string(),
            members: vec![
                "WNT3A", "SOST", "DKK1", "LRP5", "GSK3B", "CTNNB1", "AXIN2", "LEF1", "BMP2",
                "BMP6", "BMPR2", "SMAD1", "RUNX2", "SP7", "ALPL", "BGLAP", "COL1A1", "MMP13",
                "MYC", "CCND1", "PPARG", "CEBPA",
            ]
            .into_iter()
            .map(String::from)
            .collect(),
            member_expected_signs: make_signs(&[
                ("AXIN2", 1), ("LEF1", 1), ("BMP2", 1), ("BMP6", 1), ("RUNX2", 1),
                ("SP7", 1), ("ALPL", 1), ("BGLAP", 1), ("COL1A1", 1), ("MMP13", 1),
                ("MYC", 1), ("CCND1", 1), ("SOST", -1), ("DKK1", -1), ("PPARG", -1), ("CEBPA", -1),
            ]),
            node_coords: make_coords(&[
                ("SOST", 0.18, 0.06), ("WNT3A", 0.34, 0.06), ("DKK1", 0.48, 0.06),
                ("BMP2", 0.72, 0.06), ("BMP6", 0.85, 0.06),
                ("LRP5", 0.34, 0.18), ("BMPR2", 0.76, 0.18),
                ("GSK3B", 0.34, 0.35), ("AXIN2", 0.18, 0.48), ("CTNNB1", 0.34, 0.55),
                ("SMAD1", 0.76, 0.44), ("ALPL", 0.55, 0.55), ("MMP13", 0.62, 0.07),
                ("PPARG", 0.12, 0.82), ("CEBPA", 0.12, 0.89), ("LEF1", 0.26, 0.81),
                ("MYC", 0.38, 0.80), ("CCND1", 0.48, 0.80), ("RUNX2", 0.60, 0.81),
                ("SP7", 0.74, 0.84), ("BGLAP", 0.86, 0.84), ("COL1A1", 0.58, 0.89),
            ]),
            phospho_nodes: vec!["GSK3B".to_string(), "SMAD1".to_string()],
            cell_context_label: "Mesenchymal osteoblast".to_string(),
            topology_edges: make_topo(&[
                ("SOST", "LRP5", -1, "antagonism"),
                ("DKK1", "LRP5", -1, "antagonism"),
                ("WNT3A", "LRP5", 1, "binding"),
                ("LRP5", "GSK3B", -1, "inhibition"),
                ("GSK3B", "CTNNB1", -1, "degradation"),
                ("CTNNB1", "LEF1", 1, "complex"),
                ("CTNNB1", "AXIN2", 1, "transcription"),
                ("CTNNB1", "RUNX2", 1, "transcription"),
                ("CTNNB1", "PPARG", -1, "repression"),
                ("BMP2", "BMPR2", 1, "binding"),
                ("BMPR2", "SMAD1", 1, "phospho"),
                ("SMAD1", "RUNX2", 1, "co-activation"),
                ("RUNX2", "SP7", 1, "transcription"),
                ("RUNX2", "ALPL", 1, "transcription"),
                ("RUNX2", "BGLAP", 1, "transcription"),
                ("RUNX2", "COL1A1", 1, "transcription"),
            ]),
        },
        PathwayDef {
            id: "R-HSA-69278".to_string(),
            name: "G1/S Mitotic Cell Cycle & p53 Checkpoint".to_string(),
            source_db: "Reactome".to_string(),
            category: "Cell Cycle".to_string(),
            url: "https://reactome.org/content/detail/R-HSA-69278".to_string(),
            members: vec![
                "MYC", "CCND1", "E2F1", "PCNA", "MKI67", "TP53", "CDKN1A", "CDKN1B", "BCL2", "BAX",
            ]
            .into_iter()
            .map(String::from)
            .collect(),
            member_expected_signs: make_signs(&[
                ("MYC", 1), ("CCND1", 1), ("E2F1", 1), ("PCNA", 1), ("MKI67", 1),
                ("CDKN1A", -1), ("CDKN1B", -1),
            ]),
            node_coords: make_coords(&[
                ("BCL2", 0.26, 0.45), ("BAX", 0.44, 0.45),
                ("TP53", 0.24, 0.78), ("CDKN1A", 0.38, 0.78), ("CDKN1B", 0.38, 0.86),
                ("MYC", 0.62, 0.76), ("CCND1", 0.76, 0.76), ("E2F1", 0.62, 0.84),
                ("PCNA", 0.50, 0.89), ("MKI67", 0.76, 0.89),
            ]),
            phospho_nodes: vec!["TP53".to_string()],
            cell_context_label: "Proliferating somatic cell".to_string(),
            topology_edges: make_topo(&[
                ("MYC", "CCND1", 1, "transcription"),
                ("CCND1", "E2F1", 1, "CDK4/6"),
                ("E2F1", "PCNA", 1, "transcription"),
                ("E2F1", "MKI67", 1, "transcription"),
                ("TP53", "CDKN1A", 1, "checkpoint"),
                ("CDKN1A", "E2F1", -1, "CDK arrest"),
                ("CDKN1B", "E2F1", -1, "CDK arrest"),
            ]),
        },
    ];

    // 5. Open-Source Perturbation Connectivity Signatures (SigCom LINCS / L1000 consensus profiles)
    kg.signatures = vec![
        PerturbationSignature {
            id: "LINCS-CP-DEX-001".to_string(),
            perturbagen: "Dexamethasone".to_string(),
            kind: PerturbationKind::SmallMoleculeCompound,
            cell_context: "A549 / Endothelial Consensus (10 uM, 24h)".to_string(),
            mechanism_of_action: "Glucocorticoid receptor (NR3C1) agonist; NF-kB transrepressor".to_string(),
            source_db: "SigCom LINCS L1000".to_string(),
            source_url: "https://maayanlab.cloud/sigcom-lincs".to_string(),
            up_genes: vec!["NFKBIA", "KLF2", "THBD", "CDKN1A"]
                .into_iter()
                .map(String::from)
                .collect(),
            down_genes: vec![
                "ICAM1", "VCAM1", "SELE", "CCL2", "CXCL8", "IL6", "IL1B", "PTGS2", "MMP9", "CSF2",
            ]
            .into_iter()
            .map(String::from)
            .collect(),
        },
        PerturbationSignature {
            id: "LINCS-CP-IKK16-002".to_string(),
            perturbagen: "BMS-345541 (Selective IKK-beta Inhibitor)".to_string(),
            kind: PerturbationKind::SmallMoleculeCompound,
            cell_context: "HUVEC / HeLa Consensus (5 uM, 6h)".to_string(),
            mechanism_of_action: "Allosteric IKBKB inhibitor blocking NF-kB nuclear translocation".to_string(),
            source_db: "SigCom LINCS L1000".to_string(),
            source_url: "https://maayanlab.cloud/sigcom-lincs".to_string(),
            up_genes: vec!["NOS3", "KLF2", "PPARG"]
                .into_iter()
                .map(String::from)
                .collect(),
            down_genes: vec![
                "ICAM1", "VCAM1", "SELE", "CXCL8", "CCL2", "CXCL10", "IL6", "BIRC3", "TNFAIP3",
            ]
            .into_iter()
            .map(String::from)
            .collect(),
        },
        PerturbationSignature {
            id: "LINCS-LIG-TNF-003".to_string(),
            perturbagen: "Recombinant Human TNF-alpha".to_string(),
            kind: PerturbationKind::LigandStimulus,
            cell_context: "HUVEC / A375 Consensus (10 ng/mL, 6h)".to_string(),
            mechanism_of_action: "TNFRSF1A agonist activating canonical IKK/NF-kB & JNK".to_string(),
            source_db: "SigCom LINCS L1000".to_string(),
            source_url: "https://maayanlab.cloud/sigcom-lincs".to_string(),
            up_genes: vec![
                "ICAM1", "VCAM1", "SELE", "CCL2", "CXCL8", "IL6", "NFKBIA", "TNFAIP3", "BIRC3", "SOD2",
            ]
            .into_iter()
            .map(String::from)
            .collect(),
            down_genes: vec!["NOS3", "KLF2", "THBD"]
                .into_iter()
                .map(String::from)
                .collect(),
        },
        PerturbationSignature {
            id: "LINCS-CP-FULV-004".to_string(),
            perturbagen: "Fulvestrant (ICI 182,780)".to_string(),
            kind: PerturbationKind::SmallMoleculeCompound,
            cell_context: "MCF-7 Luminal Breast Carcinoma (1 uM, 24h)".to_string(),
            mechanism_of_action: "Selective Estrogen Receptor Downregulator (SERD)".to_string(),
            source_db: "SigCom LINCS L1000".to_string(),
            source_url: "https://maayanlab.cloud/sigcom-lincs".to_string(),
            up_genes: vec!["CDKN1A", "CDKN1B", "FOXO3", "BAX", "ERBB2"]
                .into_iter()
                .map(String::from)
                .collect(),
            down_genes: vec![
                "TFF1", "GREB1", "PGR", "CXCL12", "MYC", "CCND1", "AREG", "BCL2", "CTSD", "E2F1", "MKI67",
            ]
            .into_iter()
            .map(String::from)
            .collect(),
        },
        PerturbationSignature {
            id: "LINCS-CP-TAM-005".to_string(),
            perturbagen: "Tamoxifen / All-Trans Retinoic Acid".to_string(),
            kind: PerturbationKind::SmallMoleculeCompound,
            cell_context: "MCF-7 Consensus (10 uM, 24h)".to_string(),
            mechanism_of_action: "SERM / RAR agonist opposing ER-alpha proliferative transcription".to_string(),
            source_db: "SigCom LINCS L1000".to_string(),
            source_url: "https://maayanlab.cloud/sigcom-lincs".to_string(),
            up_genes: vec!["CDKN1A", "CDKN1B", "FOXO3"]
                .into_iter()
                .map(String::from)
                .collect(),
            down_genes: vec!["TFF1", "GREB1", "MYC", "CCND1", "PGR", "PCNA"]
                .into_iter()
                .map(String::from)
                .collect(),
        },
        PerturbationSignature {
            id: "LINCS-KO-SOST-006".to_string(),
            perturbagen: "Romosozumab-Mimetic / SOST CRISPR Knockout".to_string(),
            kind: PerturbationKind::CrisprKnockout,
            cell_context: "Osteoprogenitor / Mesenchymal Consensus".to_string(),
            mechanism_of_action: "Loss of sclerostin derepresses LRP5/6-Wnt and BMP/SMAD osteogenesis".to_string(),
            source_db: "SigCom LINCS CRISPR".to_string(),
            source_url: "https://maayanlab.cloud/sigcom-lincs".to_string(),
            up_genes: vec![
                "RUNX2", "SP7", "ALPL", "BGLAP", "COL1A1", "AXIN2", "LEF1", "BMP2", "CCND1",
            ]
            .into_iter()
            .map(String::from)
            .collect(),
            down_genes: vec!["PPARG", "CEBPA", "CDKN1A"]
                .into_iter()
                .map(String::from)
                .collect(),
        },
        PerturbationSignature {
            id: "LINCS-CP-LGK974-007".to_string(),
            perturbagen: "LGK-974 (Porcupine Wnt Secretion Inhibitor)".to_string(),
            kind: PerturbationKind::SmallMoleculeCompound,
            cell_context: "HT29 / Osteosarcoma Consensus (1 uM, 24h)".to_string(),
            mechanism_of_action: "PORCN inhibitor blocking Wnt palmitoylation and beta-catenin targets".to_string(),
            source_db: "SigCom LINCS L1000".to_string(),
            source_url: "https://maayanlab.cloud/sigcom-lincs".to_string(),
            up_genes: vec!["CDKN1A", "PPARG", "CEBPA"]
                .into_iter()
                .map(String::from)
                .collect(),
            down_genes: vec![
                "AXIN2", "LEF1", "RUNX2", "SP7", "ALPL", "BGLAP", "MYC", "CCND1", "MMP13",
            ]
            .into_iter()
            .map(String::from)
            .collect(),
        },
    ];

    kg.rebuild_indices();
    kg
}
