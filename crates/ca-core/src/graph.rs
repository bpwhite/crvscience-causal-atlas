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

use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet, HashMap};

/// Biological entity category in the directed causal multigraph.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize, Deserialize)]
pub enum EntityKind {
    Cytokine,
    GrowthFactor,
    TransmembraneReceptor,
    Kinase,
    Phosphatase,
    TranscriptionRegulator,
    Enzyme,
    IonChannel,
    MicroRna,
    Complex,
    ProteinFamily,
    ChemicalDrug,
    GeneOther,
    BiologicalProcess,
    DiseasePhenotype,
}

impl EntityKind {
    pub fn label(self) -> &'static str {
        match self {
            Self::Cytokine => "Cytokine",
            Self::GrowthFactor => "Growth Factor",
            Self::TransmembraneReceptor => "Transmembrane Receptor",
            Self::Kinase => "Kinase",
            Self::Phosphatase => "Phosphatase",
            Self::TranscriptionRegulator => "Transcription Regulator",
            Self::Enzyme => "Enzyme",
            Self::IonChannel => "Ion Channel",
            Self::MicroRna => "microRNA",
            Self::Complex => "Complex",
            Self::ProteinFamily => "Group / Family",
            Self::ChemicalDrug => "Chemical / Drug",
            Self::GeneOther => "Other Gene / Protein",
            Self::BiologicalProcess => "Biological Process",
            Self::DiseasePhenotype => "Disease / Phenotype",
        }
    }

    pub fn is_molecular_regulator(self) -> bool {
        !matches!(self, Self::BiologicalProcess | Self::DiseasePhenotype)
    }

    pub fn is_gene_or_product(self) -> bool {
        matches!(
            self,
            Self::Cytokine
                | Self::GrowthFactor
                | Self::TransmembraneReceptor
                | Self::Kinase
                | Self::Phosphatase
                | Self::TranscriptionRegulator
                | Self::Enzyme
                | Self::IonChannel
                | Self::MicroRna
                | Self::GeneOther
        )
    }
}

/// Subcellular compartment used for layered cellular pathway rendering.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize, Deserialize)]
pub enum SubcellularLocation {
    Extracellular,
    PlasmaMembrane,
    Cytoplasm,
    Nucleus,
    Organelle,
    PhenotypeLayer,
}

impl SubcellularLocation {
    pub fn band_index(self) -> usize {
        match self {
            Self::Extracellular => 0,
            Self::PlasmaMembrane => 1,
            Self::Cytoplasm | Self::Organelle => 2,
            Self::Nucleus => 3,
            Self::PhenotypeLayer => 4,
        }
    }

    pub fn label(self) -> &'static str {
        match self {
            Self::Extracellular => "Extracellular Space",
            Self::PlasmaMembrane => "Plasma Membrane",
            Self::Cytoplasm => "Cytoplasm",
            Self::Nucleus => "Nucleus",
            Self::Organelle => "Mitochondria / Organelle",
            Self::PhenotypeLayer => "Biological Function / Phenotype",
        }
    }
}

/// Causal edge class following Krämer et al. (2014):
/// - `Transcription` (`E_t`): transcriptional/expression regulation of target gene abundance.
/// - `Activation` (`E_a`): post-translational / signaling activation or inhibition of protein activity.
/// - `Phenotype` (`E_p`): causal modulation of downstream biological processes or disease states.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize, Deserialize)]
pub enum EdgeClass {
    Transcription,
    Activation,
    Phenotype,
}

/// Signed direction of causal regulation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize, Deserialize)]
pub enum EdgeSign {
    Activating,
    Inhibiting,
    Ambiguous,
}

impl EdgeSign {
    pub fn as_i8(self) -> i8 {
        match self {
            Self::Activating => 1,
            Self::Inhibiting => -1,
            Self::Ambiguous => 0,
        }
    }

    pub fn as_f64(self) -> f64 {
        f64::from(self.as_i8())
    }

    pub fn from_i8(v: i8) -> Self {
        match v.cmp(&0) {
            std::cmp::Ordering::Greater => Self::Activating,
            std::cmp::Ordering::Less => Self::Inhibiting,
            std::cmp::Ordering::Equal => Self::Ambiguous,
        }
    }
}

/// Node metadata in the open-source knowledge graph.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EntityNode {
    pub id: String,
    pub display_name: String,
    pub kind: EntityKind,
    pub location: SubcellularLocation,
    pub family_category: String,
    pub external_id: Option<String>,
    pub provenance_source: String,
}

/// Directed causal relationship curated from open-source streams (SIGNOR, OmniPath, CollecTRI, GO, Open Targets).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CausalEdge {
    pub source: String,
    pub target: String,
    pub edge_class: EdgeClass,
    pub sign: EdgeSign,
    pub weight: f64,
    pub mechanism: String,
    pub source_db: String,
    pub pmids: Vec<String>,
}

/// Directed reaction/signaling step inside a canonical pathway diagram.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PathwayReactionEdge {
    pub from: String,
    pub to: String,
    pub sign: EdgeSign,
    pub interaction_label: String,
}

/// Canonical pathway definition (from Reactome, WikiPathways GPML, Pathway Commons, or SIGNOR).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PathwayDef {
    pub id: String,
    pub name: String,
    pub source_db: String,
    pub category: String,
    pub url: String,
    pub members: Vec<String>,
    /// Positive (+1) or negative (-1) expected expression direction of key pathway readout members when activated.
    pub member_expected_signs: BTreeMap<String, i8>,
    /// Topological edges used for interactive pathway rendering and Molecule Activity Predictor (MAP).
    pub topology_edges: Vec<PathwayReactionEdge>,
    /// Optional normalized `(x, y)` coordinates in `[0.0, 1.0]` (from GPML CenterX/CenterY or curated vertical lanes).
    #[serde(default)]
    pub node_coords: BTreeMap<String, (f32, f32)>,
    /// Optional subset of phosphorylated signaling nodes rendered with a `(P)` badge.
    #[serde(default)]
    pub phospho_nodes: Vec<String>,
    /// Optional cellular context label displayed below the right end of the plasma membrane bar.
    #[serde(default)]
    pub cell_context_label: String,
}

impl PathwayDef {
    /// Returns spacious `(x, y)` normalized coordinates in `[0.02, 0.98]`, preserving explicit GPML/curated
    /// coordinates when available and computing a deterministic non-overlapping force-directed
    /// subcellular compartment layout otherwise.
    pub fn resolve_layout_coords(&self, kg: &KnowledgeGraph) -> BTreeMap<String, (f32, f32)> {
        crate::layout::resolve_normalized_pathway_coords(self, kg)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum PerturbationKind {
    SmallMoleculeCompound,
    CrisprKnockout,
    LigandStimulus,
    Overexpression,
}

impl PerturbationKind {
    pub fn label(self) -> &'static str {
        match self {
            Self::SmallMoleculeCompound => "Small Molecule",
            Self::CrisprKnockout => "CRISPR KO / shRNA",
            Self::LigandStimulus => "Cytokine / Ligand",
            Self::Overexpression => "cDNA Overexpression",
        }
    }
}

/// Open-source perturbation signature (e.g. SigCom LINCS / L1000 consensus up & down gene sets).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PerturbationSignature {
    pub id: String,
    pub perturbagen: String,
    pub kind: PerturbationKind,
    pub cell_context: String,
    pub mechanism_of_action: String,
    pub source_db: String,
    pub source_url: String,
    pub up_genes: Vec<String>,
    pub down_genes: Vec<String>,
    #[serde(default)]
    pub target_clamps: Vec<(String, i8)>,
}

/// Unified in-memory knowledge graph assembled exclusively from open-source streams.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct KnowledgeGraph {
    pub nodes: BTreeMap<String, EntityNode>,
    pub edges: Vec<CausalEdge>,
    pub pathways: Vec<PathwayDef>,
    pub signatures: Vec<PerturbationSignature>,

    #[serde(skip)]
    pub out_edges_t: HashMap<String, Vec<usize>>,
    #[serde(skip)]
    pub out_edges_at: HashMap<String, Vec<usize>>,
    #[serde(skip)]
    pub in_edges_p: HashMap<String, Vec<usize>>,
    #[serde(skip)]
    pub regulated_gene_universe: BTreeSet<String>,
    #[serde(skip)]
    pub pathway_gene_universe: BTreeSet<String>,
    #[serde(skip)]
    pub in_degree_t: HashMap<String, usize>,
}

impl KnowledgeGraph {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn upsert_node(&mut self, node: EntityNode) {
        self.nodes.entry(node.id.clone()).or_insert(node);
    }

    pub fn ensure_gene_node(&mut self, symbol: &str, source: &str) {
        let key = symbol.trim().to_uppercase();
        if key.is_empty() {
            return;
        }
        self.nodes.entry(key.clone()).or_insert_with(|| EntityNode {
            id: key.clone(),
            display_name: key,
            kind: EntityKind::GeneOther,
            location: SubcellularLocation::Cytoplasm,
            family_category: "General Proteome".to_string(),
            external_id: None,
            provenance_source: source.to_string(),
        });
    }

    pub fn add_edge(&mut self, mut edge: CausalEdge) {
        edge.source = edge.source.trim().to_uppercase();
        edge.target = edge.target.trim().to_uppercase();
        if edge.source.is_empty() || edge.target.is_empty() || edge.source == edge.target {
            return;
        }
        edge.weight = edge.weight.clamp(0.05, 1.0);
        self.edges.push(edge);
    }

    /// Rebuilds fast adjacency indices and gene universes after loading or merging streams.
    pub fn rebuild_indices(&mut self) {
        self.out_edges_t.clear();
        self.out_edges_at.clear();
        self.in_edges_p.clear();
        self.regulated_gene_universe.clear();
        self.pathway_gene_universe.clear();
        self.in_degree_t.clear();

        for (idx, edge) in self.edges.iter().enumerate() {
            match edge.edge_class {
                EdgeClass::Transcription => {
                    self.out_edges_t
                        .entry(edge.source.clone())
                        .or_default()
                        .push(idx);
                    self.out_edges_at
                        .entry(edge.source.clone())
                        .or_default()
                        .push(idx);
                    self.regulated_gene_universe.insert(edge.target.clone());
                    *self.in_degree_t.entry(edge.target.clone()).or_insert(0) += 1;
                }
                EdgeClass::Activation => {
                    self.out_edges_at
                        .entry(edge.source.clone())
                        .or_default()
                        .push(idx);
                }
                EdgeClass::Phenotype => {
                    self.in_edges_p
                        .entry(edge.target.clone())
                        .or_default()
                        .push(idx);
                }
            }
        }

        for pw in &self.pathways {
            for m in &pw.members {
                self.pathway_gene_universe.insert(m.clone());
            }
        }
        for (id, node) in &self.nodes {
            if node.kind.is_gene_or_product() {
                self.pathway_gene_universe.insert(id.clone());
            }
        }
    }

    /// Merges another `KnowledgeGraph` bundle, deduplicating nodes, pathways, and identical causal edges.
    pub fn merge(&mut self, other: KnowledgeGraph) {
        for (id, node) in other.nodes {
            self.nodes.entry(id).or_insert(node);
        }

        let mut existing_edge_keys: HashMap<(String, String, EdgeClass), usize> = self
            .edges
            .iter()
            .enumerate()
            .map(|(idx, e)| ((e.source.clone(), e.target.clone(), e.edge_class), idx))
            .collect();

        for edge in other.edges {
            let key = (edge.source.clone(), edge.target.clone(), edge.edge_class);
            if let Some(&idx) = existing_edge_keys.get(&key) {
                let cur = &mut self.edges[idx];
                if edge.weight > cur.weight {
                    cur.weight = edge.weight;
                }
                if cur.sign == EdgeSign::Ambiguous && edge.sign != EdgeSign::Ambiguous {
                    cur.sign = edge.sign;
                }
                for pmid in edge.pmids {
                    if !cur.pmids.contains(&pmid) {
                        cur.pmids.push(pmid);
                    }
                }
            } else {
                let idx = self.edges.len();
                self.edges.push(edge);
                existing_edge_keys.insert(key, idx);
            }
        }

        let existing_pw_ids: BTreeSet<String> =
            self.pathways.iter().map(|p| p.id.clone()).collect();
        for pw in other.pathways {
            if !existing_pw_ids.contains(&pw.id) {
                self.pathways.push(pw);
            }
        }

        let existing_sig_ids: BTreeSet<String> =
            self.signatures.iter().map(|s| s.id.clone()).collect();
        for sig in other.signatures {
            if !existing_sig_ids.contains(&sig.id) {
                self.signatures.push(sig);
            }
        }

        self.rebuild_indices();
    }
}
