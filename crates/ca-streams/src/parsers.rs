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

use crate::{StreamError, StreamFormat};
use ca_core::{
    CausalEdge, EdgeClass, EdgeSign, EntityKind, EntityNode, KnowledgeGraph, PathwayDef,
    PathwayReactionEdge, PerturbationSignature, SubcellularLocation,
};
use std::collections::BTreeMap;

/// Dispatches raw UTF-8 payload text to the corresponding open-source biological format parser.
pub fn parse_stream_payload(
    format: StreamFormat,
    source_label: &str,
    payload: &str,
) -> Result<(KnowledgeGraph, usize), StreamError> {
    match format {
        StreamFormat::GmtPathways => parse_gmt_pathways(source_label, payload),
        StreamFormat::GpmlXml => parse_gpml_pathway(source_label, payload),
        StreamFormat::SignorTsv => parse_signor_tsv(payload),
        StreamFormat::OmniPathTsv => parse_omnipath_tsv(source_label, payload),
        StreamFormat::StringTsv => parse_string_tsv(payload),
        StreamFormat::GoGafTsv => parse_go_gaf(payload),
        StreamFormat::OpenTargetsGraphQlJson => parse_open_targets_json(payload),
        StreamFormat::LincsSignaturesJson => parse_lincs_signatures_json(payload),
    }
}

/// Parses GMT (`Gene Matrix Transposed`) files from Reactome or WikiPathways.
/// Format per line: `pathway_name_or_id <TAB> description_or_url <TAB> GENE1 <TAB> GENE2 ...`
pub fn parse_gmt_pathways(
    source_db: &str,
    raw: &str,
) -> Result<(KnowledgeGraph, usize), StreamError> {
    let mut kg = KnowledgeGraph::new();
    let mut count = 0usize;

    for (line_idx, line) in raw.lines().enumerate() {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = trimmed.split('\t').collect();
        if cols.len() < 3 {
            continue;
        }
        let raw_title = cols[0].trim();
        let desc_or_url = cols[1].trim();

        // WikiPathways GMT headers often encode `Name%WikiPathways_20261010%WP445%Homo sapiens`
        let (pw_id, pw_name) = if raw_title.contains('%') {
            let parts: Vec<&str> = raw_title.split('%').collect();
            let name = parts[0].trim().to_string();
            let id = parts
                .iter()
                .find(|p| p.starts_with("WP") || p.starts_with("R-HSA-"))
                .copied()
                .unwrap_or(raw_title)
                .to_string();
            (id, name)
        } else {
            let id = if desc_or_url.starts_with("R-HSA-") {
                desc_or_url.to_string()
            } else {
                format!("{}-{:04}", source_db.replace(' ', "_"), line_idx + 1)
            };
            (id, raw_title.to_string())
        };

        let members: Vec<String> = cols[2..]
            .iter()
            .map(|s| s.trim().to_uppercase())
            .filter(|s| !s.is_empty())
            .collect();

        if members.len() < 2 {
            continue;
        }

        for m in &members {
            kg.ensure_gene_node(m, source_db);
        }

        let url = if desc_or_url.starts_with("http://") || desc_or_url.starts_with("https://") {
            desc_or_url.to_string()
        } else if pw_id.starts_with("WP") {
            format!("https://www.wikipathways.org/pathways/{pw_id}.html")
        } else if pw_id.starts_with("R-HSA-") {
            format!("https://reactome.org/content/detail/{pw_id}")
        } else {
            String::new()
        };

        kg.pathways.push(PathwayDef {
            id: pw_id,
            name: pw_name,
            source_db: source_db.to_string(),
            category: "Imported GMT Pathway".to_string(),
            url,
            members,
            member_expected_signs: BTreeMap::new(),
            topology_edges: Vec::new(),
            node_coords: BTreeMap::new(),
            phospho_nodes: Vec::new(),
            cell_context_label: String::new(),
        });
        count += 1;
    }

    kg.rebuild_indices();
    Ok((kg, count))
}

/// Lightweight zero-C-dependency parser for WikiPathways GPML XML files (`<DataNode>` + `<Interaction>`).
pub fn parse_gpml_pathway(
    source_db: &str,
    xml: &str,
) -> Result<(KnowledgeGraph, usize), StreamError> {
    let mut kg = KnowledgeGraph::new();
    let pw_name = extract_xml_attr(xml, "<Pathway", "Name")
        .unwrap_or_else(|| "Imported GPML Pathway".to_string());
    let organism = extract_xml_attr(xml, "<Pathway", "Organism").unwrap_or_default();

    let mut graph_id_to_symbol: BTreeMap<String, String> = BTreeMap::new();
    let mut members = Vec::new();

    for chunk in xml.split("<DataNode").skip(1) {
        let Some(tag_end) = chunk.find('>') else {
            continue;
        };
        let header = format!("<DataNode{}", &chunk[..=tag_end]);
        let Some(label) = extract_xml_attr(&header, "<DataNode", "TextLabel") else {
            continue;
        };
        let sym = label.trim().to_uppercase();
        if sym.is_empty() {
            continue;
        }
        let node_type = extract_xml_attr(&header, "<DataNode", "Type").unwrap_or_default();
        let kind = match node_type.as_str() {
            "Metabolite" => EntityKind::ChemicalDrug,
            "Complex" => EntityKind::Complex,
            "Pathway" => EntityKind::BiologicalProcess,
            _ => EntityKind::GeneOther,
        };

        if let Some(gid) = extract_xml_attr(&header, "<DataNode", "GraphId") {
            graph_id_to_symbol.insert(gid, sym.clone());
        }
        if kind != EntityKind::BiologicalProcess && !members.contains(&sym) {
            members.push(sym.clone());
        }
        kg.upsert_node(EntityNode {
            id: sym.clone(),
            display_name: sym,
            kind,
            location: SubcellularLocation::Cytoplasm,
            family_category: if organism.is_empty() {
                "GPML Pathway".to_string()
            } else {
                organism.clone()
            },
            external_id: None,
            provenance_source: source_db.to_string(),
        });
    }

    let mut topology_edges = Vec::new();
    for ic in xml.split("<Interaction").skip(1) {
        let end_pos = ic.find("</Interaction>").unwrap_or(ic.len());
        let body = &ic[..end_pos];
        let mut refs = Vec::new();
        let mut sign = EdgeSign::Activating;

        for pt in body.split("<Point").skip(1) {
            let Some(pe) = pt.find('>') else {
                continue;
            };
            let ptag = format!("<Point{}", &pt[..=pe]);
            if let Some(gref) = extract_xml_attr(&ptag, "<Point", "GraphRef") {
                refs.push(gref);
            }
            if let Some(arrow) = extract_xml_attr(&ptag, "<Point", "ArrowHead") {
                let al = arrow.to_lowercase();
                if al.contains("inhibition") || al.contains("tbar") {
                    sign = EdgeSign::Inhibiting;
                }
            }
        }
        if refs.len() >= 2 {
            let src_opt = graph_id_to_symbol.get(&refs[0]);
            let dst_opt = refs.last().and_then(|r| graph_id_to_symbol.get(r));
            if let (Some(src), Some(dst)) = (src_opt, dst_opt) {
                if src != dst {
                    topology_edges.push(PathwayReactionEdge {
                        from: src.clone(),
                        to: dst.clone(),
                        sign,
                        interaction_label: "GPML interaction".to_string(),
                    });
                }
            }
        }
    }

    // Extract GPML <Graphics CenterX="..." CenterY="..."/> inside each <DataNode> so author coordinates are preserved!
    let mut raw_xy: Vec<(String, f32, f32)> = Vec::new();
    for chunk in xml.split("<DataNode").skip(1) {
        let end_dn = chunk.find("</DataNode>").unwrap_or(chunk.len());
        let dn_block = &chunk[..end_dn];
        let Some(tag_end) = dn_block.find('>') else { continue; };
        let header = format!("<DataNode{}", &dn_block[..=tag_end]);
        let Some(label) = extract_xml_attr(&header, "<DataNode", "TextLabel") else { continue; };
        let sym = label.trim().to_uppercase();
        if let Some(gfx_idx) = dn_block.find("<Graphics") {
            let gfx_slice = &dn_block[gfx_idx..];
            if let (Some(cx_s), Some(cy_s)) = (
                extract_xml_attr(gfx_slice, "<Graphics", "CenterX"),
                extract_xml_attr(gfx_slice, "<Graphics", "CenterY"),
            ) {
                if let (Ok(cx), Ok(cy)) = (cx_s.parse::<f32>(), cy_s.parse::<f32>()) {
                    raw_xy.push((sym, cx, cy));
                }
            }
        }
    }

    let mut node_coords = BTreeMap::new();
    if raw_xy.len() >= 2 {
        let min_x = raw_xy.iter().map(|p| p.1).fold(f32::INFINITY, f32::min);
        let max_x = raw_xy.iter().map(|p| p.1).fold(f32::NEG_INFINITY, f32::max);
        let min_y = raw_xy.iter().map(|p| p.2).fold(f32::INFINITY, f32::min);
        let max_y = raw_xy.iter().map(|p| p.2).fold(f32::NEG_INFINITY, f32::max);
        let span_x = (max_x - min_x).max(1.0);
        let span_y = (max_y - min_y).max(1.0);
        for (sym, cx, cy) in raw_xy {
            let nx = 0.08 + ((cx - min_x) / span_x) * 0.84;
            let ny = 0.07 + ((cy - min_y) / span_y) * 0.84;
            node_coords.insert(sym, (nx, ny));
        }
    }

    let edge_cnt = topology_edges.len();
    kg.pathways.push(PathwayDef {
        id: format!(
            "GPML-{}",
            pw_name
                .chars()
                .filter(|c| c.is_alphanumeric())
                .take(12)
                .collect::<String>()
        ),
        name: pw_name,
        source_db: source_db.to_string(),
        category: "GPML Pathway".to_string(),
        url: "https://www.wikipathways.org/".to_string(),
        members,
        member_expected_signs: BTreeMap::new(),
        topology_edges,
        node_coords,
        phospho_nodes: Vec::new(),
        cell_context_label: organism,
    });

    kg.rebuild_indices();
    Ok((kg, edge_cnt + 1))
}

fn extract_xml_attr(snippet: &str, _tag: &str, attr: &str) -> Option<String> {
    let needle = format!("{attr}=\"");
    let start = snippet.find(&needle)? + needle.len();
    let rest = &snippet[start..];
    let end = rest.find('"')?;
    Some(rest[..end].to_string())
}

/// Parses SIGNOR 3.0 tab-delimited causal interaction files (`ENTITYA`, `TYPEA`, `ENTITYB`, `TYPEB`, `EFFECT`, `MECHANISM`, `PMID`, `SCORE`).
pub fn parse_signor_tsv(raw: &str) -> Result<(KnowledgeGraph, usize), StreamError> {
    let mut kg = KnowledgeGraph::new();
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .flexible(true)
        .has_headers(true)
        .from_reader(raw.as_bytes());

    let headers = rdr
        .headers()
        .map_err(|e| StreamError::Parse(e.to_string()))?
        .iter()
        .map(|h| h.trim().to_uppercase())
        .collect::<Vec<_>>();

    let idx_of = |names: &[&str]| -> Option<usize> {
        headers
            .iter()
            .position(|h| names.iter().any(|n| h == n))
    };

    let (Some(ia), Some(ib), Some(ieff)) = (
        idx_of(&["ENTITYA", "SOURCE"]),
        idx_of(&["ENTITYB", "TARGET"]),
        idx_of(&["EFFECT"]),
    ) else {
        return Err(StreamError::Parse(
            "SIGNOR TSV missing ENTITYA / ENTITYB / EFFECT columns".to_string(),
        ));
    };

    let itype_b = idx_of(&["TYPEB"]);
    let imech = idx_of(&["MECHANISM"]);
    let ipmid = idx_of(&["PMID"]);
    let iscore = idx_of(&["SCORE"]);

    let mut ingested = 0usize;
    for rec in rdr.records().flatten() {
        let src = rec.get(ia).unwrap_or("").trim().to_uppercase();
        let dst = rec.get(ib).unwrap_or("").trim().to_uppercase();
        if src.is_empty() || dst.is_empty() || src == dst {
            continue;
        }
        let eff = rec.get(ieff).unwrap_or("").to_lowercase();
        let mech = imech
            .and_then(|i| rec.get(i))
            .unwrap_or("causal interaction")
            .trim()
            .to_lowercase();
        let type_b = itype_b
            .and_then(|i| rec.get(i))
            .unwrap_or("protein")
            .to_lowercase();

        let sign = if eff.contains("up-regulates") || eff.contains("activates") {
            EdgeSign::Activating
        } else if eff.contains("down-regulates") || eff.contains("inhibits") {
            EdgeSign::Inhibiting
        } else {
            EdgeSign::Ambiguous
        };

        let edge_class = if type_b.contains("phenotype") || type_b.contains("stimulus") {
            EdgeClass::Phenotype
        } else if eff.contains("quantity")
            || mech.contains("transcriptional")
            || mech.contains("translation")
        {
            EdgeClass::Transcription
        } else {
            EdgeClass::Activation
        };

        let weight = iscore
            .and_then(|i| rec.get(i))
            .and_then(|s| s.trim().parse::<f64>().ok())
            .unwrap_or(0.75)
            .clamp(0.15, 1.0);

        let pmids: Vec<String> = ipmid
            .and_then(|i| rec.get(i))
            .map(|s| {
                s.split([';', ',', '|'])
                    .map(str::trim)
                    .filter(|p| !p.is_empty())
                    .map(String::from)
                    .collect()
            })
            .unwrap_or_default();

        kg.ensure_gene_node(&src, "SIGNOR 3.0");
        if edge_class == EdgeClass::Phenotype {
            kg.upsert_node(EntityNode {
                id: dst.clone(),
                display_name: dst.clone(),
                kind: EntityKind::BiologicalProcess,
                location: SubcellularLocation::PhenotypeLayer,
                family_category: "SIGNOR Phenotype".to_string(),
                external_id: None,
                provenance_source: "SIGNOR 3.0".to_string(),
            });
        } else {
            kg.ensure_gene_node(&dst, "SIGNOR 3.0");
        }

        kg.add_edge(CausalEdge {
            source: src,
            target: dst,
            edge_class,
            sign,
            weight,
            mechanism: mech,
            source_db: "SIGNOR 3.0".to_string(),
            pmids,
        });
        ingested += 1;
    }

    kg.rebuild_indices();
    Ok((kg, ingested))
}

/// Parses OmniPath / CollecTRI / DoRothEA TSV webservice responses:
/// Expected columns: `source_genesymbol` (or `source`), `target_genesymbol` (or `target`),
/// `is_stimulation`, `is_inhibition`, `consensus_stimulation`, `consensus_inhibition`, `curation_effort`, `sources`, `references`.
pub fn parse_omnipath_tsv(
    stream_name: &str,
    raw: &str,
) -> Result<(KnowledgeGraph, usize), StreamError> {
    let mut kg = KnowledgeGraph::new();
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .flexible(true)
        .has_headers(true)
        .from_reader(raw.as_bytes());

    let headers = rdr
        .headers()
        .map_err(|e| StreamError::Parse(e.to_string()))?
        .iter()
        .map(|h| h.trim().to_lowercase())
        .collect::<Vec<_>>();

    let find_col = |candidates: &[&str]| -> Option<usize> {
        candidates
            .iter()
            .find_map(|cand| headers.iter().position(|h| h == cand))
    };

    let (Some(isrc), Some(idst)) = (
        find_col(&["source_genesymbol", "source"]),
        find_col(&["target_genesymbol", "target"]),
    ) else {
        return Err(StreamError::Parse(
            "OmniPath TSV missing source/target columns".to_string(),
        ));
    };

    let istim = find_col(&["consensus_stimulation", "is_stimulation"]);
    let iinh = find_col(&["consensus_inhibition", "is_inhibition"]);
    let ieffort = find_col(&["curation_effort"]);
    let irefs = find_col(&["references"]);

    let is_transcriptional_stream = stream_name.to_lowercase().contains("collectri")
        || stream_name.to_lowercase().contains("dorothea")
        || stream_name.to_lowercase().contains("regulon");

    let mut ingested = 0usize;
    for rec in rdr.records().flatten() {
        let src = rec.get(isrc).unwrap_or("").trim().to_uppercase();
        let dst = rec.get(idst).unwrap_or("").trim().to_uppercase();
        if src.is_empty() || dst.is_empty() || src == dst {
            continue;
        }

        let stim = istim
            .and_then(|i| rec.get(i))
            .map(|v| matches!(v.trim(), "1" | "True" | "true"))
            .unwrap_or(false);
        let inh = iinh
            .and_then(|i| rec.get(i))
            .map(|v| matches!(v.trim(), "1" | "True" | "true"))
            .unwrap_or(false);

        let sign = match (stim, inh) {
            (true, false) => EdgeSign::Activating,
            (false, true) => EdgeSign::Inhibiting,
            _ => EdgeSign::Ambiguous,
        };

        let effort = ieffort
            .and_then(|i| rec.get(i))
            .and_then(|v| v.trim().parse::<f64>().ok())
            .unwrap_or(2.0);
        let weight = (0.55 + (effort.ln_1p() * 0.15)).clamp(0.45, 0.99);

        let pmids: Vec<String> = irefs
            .and_then(|i| rec.get(i))
            .map(|s| {
                s.split(';')
                    .filter_map(|tok| tok.rsplit(':').next())
                    .map(str::trim)
                    .filter(|p| !p.is_empty())
                    .take(5)
                    .map(String::from)
                    .collect()
            })
            .unwrap_or_default();

        kg.upsert_node(EntityNode {
            id: src.clone(),
            display_name: src.clone(),
            kind: if is_transcriptional_stream {
                EntityKind::TranscriptionRegulator
            } else {
                EntityKind::GeneOther
            },
            location: if is_transcriptional_stream {
                SubcellularLocation::Nucleus
            } else {
                SubcellularLocation::Cytoplasm
            },
            family_category: stream_name.to_string(),
            external_id: None,
            provenance_source: stream_name.to_string(),
        });
        kg.ensure_gene_node(&dst, stream_name);

        kg.add_edge(CausalEdge {
            source: src,
            target: dst,
            edge_class: if is_transcriptional_stream {
                EdgeClass::Transcription
            } else {
                EdgeClass::Activation
            },
            sign,
            weight,
            mechanism: if is_transcriptional_stream {
                "transcriptional regulation".to_string()
            } else {
                "post-translational signaling".to_string()
            },
            source_db: stream_name.to_string(),
            pmids,
        });
        ingested += 1;
    }

    kg.rebuild_indices();
    Ok((kg, ingested))
}

/// Parses STRING API TSV interaction responses (`preferredName_A`, `preferredName_B`, `score`).
pub fn parse_string_tsv(raw: &str) -> Result<(KnowledgeGraph, usize), StreamError> {
    let mut kg = KnowledgeGraph::new();
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .flexible(true)
        .has_headers(true)
        .from_reader(raw.as_bytes());

    let headers = rdr
        .headers()
        .map_err(|e| StreamError::Parse(e.to_string()))?
        .iter()
        .map(|h| h.trim().to_lowercase())
        .collect::<Vec<_>>();

    let ia = headers
        .iter()
        .position(|h| h == "preferredname_a" || h == "protein1");
    let ib = headers
        .iter()
        .position(|h| h == "preferredname_b" || h == "protein2");
    let iscore = headers
        .iter()
        .position(|h| h == "score" || h == "combined_score");

    let (Some(ca), Some(cb)) = (ia, ib) else {
        return Err(StreamError::Parse(
            "STRING TSV missing preferredName_A / preferredName_B".to_string(),
        ));
    };

    let mut count = 0usize;
    for rec in rdr.records().flatten() {
        let a = rec.get(ca).unwrap_or("").trim().to_uppercase();
        let b = rec.get(cb).unwrap_or("").trim().to_uppercase();
        if a.is_empty() || b.is_empty() || a == b {
            continue;
        }
        let mut sc = iscore
            .and_then(|i| rec.get(i))
            .and_then(|s| s.trim().parse::<f64>().ok())
            .unwrap_or(0.7);
        if sc > 1.0 {
            sc /= 1000.0;
        }
        if sc < 0.4 {
            continue;
        }

        kg.ensure_gene_node(&a, "STRING v12");
        kg.ensure_gene_node(&b, "STRING v12");
        kg.add_edge(CausalEdge {
            source: a,
            target: b,
            edge_class: EdgeClass::Activation,
            sign: EdgeSign::Ambiguous,
            weight: sc.clamp(0.4, 0.99),
            mechanism: "functional association".to_string(),
            source_db: "STRING v12".to_string(),
            pmids: Vec::new(),
        });
        count += 1;
    }

    kg.rebuild_indices();
    Ok((kg, count))
}

/// Parses Gene Ontology GAF 2.1 annotation lines (`!` comments skipped, Biological Process `P` aspect extracted).
pub fn parse_go_gaf(raw: &str) -> Result<(KnowledgeGraph, usize), StreamError> {
    let mut kg = KnowledgeGraph::new();
    let mut count = 0usize;

    for line in raw.lines() {
        if line.starts_with('!') || line.trim().is_empty() {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 9 {
            continue;
        }
        let symbol = cols[2].trim().to_uppercase();
        let qualifier = cols[3].trim().to_lowercase();
        let go_id = cols[4].trim().to_uppercase();
        let aspect = cols[8].trim();

        if symbol.is_empty() || !go_id.starts_with("GO:") || aspect != "P" || qualifier.contains("not") {
            continue;
        }

        let sign = if qualifier.contains("positively_regulates") {
            EdgeSign::Activating
        } else if qualifier.contains("negatively_regulates") {
            EdgeSign::Inhibiting
        } else {
            EdgeSign::Activating
        };

        kg.ensure_gene_node(&symbol, "GOA GAF");
        kg.upsert_node(EntityNode {
            id: go_id.clone(),
            display_name: go_id.clone(),
            kind: EntityKind::BiologicalProcess,
            location: SubcellularLocation::PhenotypeLayer,
            family_category: "Gene Ontology Biological Process".to_string(),
            external_id: Some(go_id.clone()),
            provenance_source: "GOA GAF 2.1".to_string(),
        });
        kg.add_edge(CausalEdge {
            source: symbol,
            target: go_id,
            edge_class: EdgeClass::Phenotype,
            sign,
            weight: 0.80,
            mechanism: "GO biological process annotation".to_string(),
            source_db: "Gene Ontology".to_string(),
            pmids: Vec::new(),
        });
        count += 1;
    }

    kg.rebuild_indices();
    Ok((kg, count))
}

/// Parses Open Targets Platform v4 GraphQL JSON disease-target association responses.
pub fn parse_open_targets_json(raw: &str) -> Result<(KnowledgeGraph, usize), StreamError> {
    let val: serde_json::Value = serde_json::from_str(raw)?;
    let mut kg = KnowledgeGraph::new();
    let mut count = 0usize;

    if let Some(disease) = val.pointer("/data/disease") {
        let disease_id = disease
            .get("id")
            .and_then(|v| v.as_str())
            .unwrap_or("EFO:0000001")
            .to_uppercase();
        let disease_name = disease
            .get("name")
            .and_then(|v| v.as_str())
            .unwrap_or(&disease_id)
            .to_string();

        kg.upsert_node(EntityNode {
            id: disease_id.clone(),
            display_name: disease_name,
            kind: EntityKind::DiseasePhenotype,
            location: SubcellularLocation::PhenotypeLayer,
            family_category: "Open Targets Disease".to_string(),
            external_id: Some(disease_id.clone()),
            provenance_source: "Open Targets v4".to_string(),
        });

        if let Some(rows) = disease
            .pointer("/associatedTargets/rows")
            .and_then(|v| v.as_array())
        {
            for row in rows {
                let Some(sym) = row
                    .pointer("/target/approvedSymbol")
                    .and_then(|v| v.as_str())
                else {
                    continue;
                };
                let score = row.get("score").and_then(|v| v.as_f64()).unwrap_or(0.6);
                if score < 0.2 {
                    continue;
                }
                let sym_up = sym.trim().to_uppercase();
                kg.ensure_gene_node(&sym_up, "Open Targets v4");
                kg.add_edge(CausalEdge {
                    source: sym_up,
                    target: disease_id.clone(),
                    edge_class: EdgeClass::Phenotype,
                    sign: EdgeSign::Activating,
                    weight: score.clamp(0.2, 0.99),
                    mechanism: "Open Targets integrated disease association".to_string(),
                    source_db: "Open Targets v4".to_string(),
                    pmids: Vec::new(),
                });
                count += 1;
            }
        }
    }

    kg.rebuild_indices();
    Ok((kg, count))
}

/// Parses JSON arrays of `PerturbationSignature` records (SigCom LINCS / open drug profiles).
pub fn parse_lincs_signatures_json(raw: &str) -> Result<(KnowledgeGraph, usize), StreamError> {
    let sigs: Vec<PerturbationSignature> = serde_json::from_str(raw)?;
    let count = sigs.len();
    let mut kg = KnowledgeGraph::new();
    for mut s in sigs {
        s.up_genes = s.up_genes.into_iter().map(|g| g.trim().to_uppercase()).collect();
        s.down_genes = s
            .down_genes
            .into_iter()
            .map(|g| g.trim().to_uppercase())
            .collect();
        kg.signatures.push(s);
    }
    kg.rebuild_indices();
    Ok((kg, count))
}

#[cfg(test)]
mod tests {
    use super::*;
    use ca_core::PerturbationKind;

    #[test]
    fn parses_wikipathways_gmt() {
        let sample = "TNF Signaling%WikiPathways%WP231%Homo sapiens\thttps://www.wikipathways.org/instance/WP231\tTNF\tTNFRSF1A\tNFKB1\tICAM1\n";
        let (kg, cnt) = parse_gmt_pathways("WikiPathways", sample).unwrap();
        assert_eq!(cnt, 1);
        assert_eq!(kg.pathways[0].id, "WP231");
        assert_eq!(kg.pathways[0].members.len(), 4);
    }

    #[test]
    fn parses_signor_and_omnipath_tsvs() {
        let signor = "ENTITYA\tTYPEA\tENTITYB\tTYPEB\tEFFECT\tMECHANISM\tPMID\tSCORE\nTNF\tprotein\tICAM1\tprotein\tup-regulates quantity by expression\ttranscriptional regulation\t1679050\t0.95\n";
        let (kg1, c1) = parse_signor_tsv(signor).unwrap();
        assert_eq!(c1, 1);
        assert_eq!(kg1.edges[0].edge_class, EdgeClass::Transcription);
        assert_eq!(kg1.edges[0].sign, EdgeSign::Activating);

        let omni = "source_genesymbol\ttarget_genesymbol\tconsensus_stimulation\tconsensus_inhibition\tcuration_effort\treferences\nESR1\tCDKN1A\t0\t1\t5\tCollecTRI:15642340\n";
        let (kg2, c2) = parse_omnipath_tsv("CollecTRI", omni).unwrap();
        assert_eq!(c2, 1);
        assert_eq!(kg2.edges[0].sign, EdgeSign::Inhibiting);
    }

    #[test]
    fn parses_perturbation_json() {
        let j = r#"[{"id":"S1","perturbagen":"TestCompound","kind":"SmallMoleculeCompound","cell_context":"MCF7","mechanism_of_action":"HDAC inhibitor","source_db":"SigCom LINCS","source_url":"https://lincsproject.org/","up_genes":["CDKN1A"],"down_genes":["MYC"]}]"#;
        let (kg, c) = parse_lincs_signatures_json(j).unwrap();
        assert_eq!(c, 1);
        assert_eq!(kg.signatures[0].kind, PerturbationKind::SmallMoleculeCompound);
    }
}
