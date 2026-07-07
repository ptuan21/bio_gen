use crate::analysis::search::search;
use crate::sequence::Sequence;

/// A restriction enzyme with an IUPAC recognition site and the cut offset on
/// the forward strand (index into the site where the phosphodiester bond breaks).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Enzyme {
    pub name: &'static str,
    pub site: &'static str,
    pub cut_offset: usize,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SiteHit {
    pub enzyme: &'static str,
    /// 0-based start of the recognition site.
    pub start: usize,
    /// 0-based position of the cut (`start + cut_offset`).
    pub cut: usize,
}

/// A small built-in panel of common enzymes.
pub const ENZYMES: &[Enzyme] = &[
    Enzyme { name: "EcoRI", site: "GAATTC", cut_offset: 1 },
    Enzyme { name: "BamHI", site: "GGATCC", cut_offset: 1 },
    Enzyme { name: "HindIII", site: "AAGCTT", cut_offset: 1 },
    Enzyme { name: "EcoRV", site: "GATATC", cut_offset: 3 },
    Enzyme { name: "PstI", site: "CTGCAG", cut_offset: 5 },
    Enzyme { name: "SmaI", site: "CCCGGG", cut_offset: 3 },
    Enzyme { name: "NotI", site: "GCGGCCGC", cut_offset: 2 },
];

pub fn find_by_name(name: &str) -> Option<&'static Enzyme> {
    ENZYMES.iter().find(|e| e.name.eq_ignore_ascii_case(name))
}

/// All forward-strand occurrences of one enzyme's recognition site.
pub fn find_sites(seq: &Sequence, enzyme: &Enzyme) -> Vec<SiteHit> {
    search(seq, enzyme.site, false)
        .into_iter()
        .map(|m| SiteHit {
            enzyme: enzyme.name,
            start: m.start,
            cut: m.start + enzyme.cut_offset,
        })
        .collect()
}

/// Digest with every enzyme in `enzymes`, sorted by cut position.
pub fn digest(seq: &Sequence, enzymes: &[Enzyme]) -> Vec<SiteHit> {
    let mut hits: Vec<SiteHit> = enzymes
        .iter()
        .flat_map(|e| find_sites(seq, e))
        .collect();
    hits.sort_by_key(|h| (h.cut, h.start));
    hits
}
