# Chapter 1 — Reference Reuse Suggestions

**Source corpus:** 410 annotated references in `outputs/relevance_notes.json` (308 API-generated, 99 hand-written, 3 manual PDF for predecessor theses).
**Already-cited check:** 152 references currently appear in Chapter 1 (`outputs/ch1_citations_in_text.csv`). The `⚲ already cited in Ch1` flag marks those — reusing them for an *additional* uncited claim is "free" citation value.
**Matching method:** keyword/concept retrieval against the corpus, then rationale grounded in each candidate's relevance note (not just titles).

---

## Summary

| Metric | Result |
|---|---:|
| Claims processed in this pass (claims 6-26 in document order) | **21** |
| High-yield subsections (≥3 confident matches across the subsection's claims) | §1.5, §1.7, §1.9, §1.12, §1.15 |
| Low-yield subsections (0-1 confident match across the subsection's claims) | §1.2 (alcohol/IARC), §1.3 (models), §1.6 (rare diseases) |
| New-citation reuses surfaced (in corpus, NOT yet cited in Ch1) | **17** — see consolidated list below |
| Already-cited reuses surfaced (free additional use of existing Ch1 citations) | **21** |
| Gaps flagged (canonical references absent from corpus) | **9** — see consolidated list below |

**Note on §1.13:** the user's brief mentioned a "methyl substituents at positions 4 and 11b" claim about aphidicolin's chemical structure. That claim does not appear in the current `ch1_citation_opportunities.md` (the v2 filter caught it, presumably). If it resurfaces in a future regeneration, it should be flagged as a likely figure-caption / structural-description false positive that needs no citation rather than searched in the bibliography.

### New citations (in corpus, not yet cited in Ch1) — highest value wins

These are references whose relevance notes describe direct fit to a specific Ch1 claim and which are not currently used in the chapter. Adding them is "free" in the sense that no new literature search is needed.

| Surname | Year | Claim it supports |
|---|---|---|
| Choudhury | 2016 | §1.7 — CRISPR/Cas9 in cancer research (BRCA1 demethylation via dCas9-TET1) |
| Hills | 2014 | §1.12 — oncogene-induced replication stress / genomic instability |
| Ekholm-Reed | 2004 | §1.12 — cyclin E disrupts MCM loading + pre-RC assembly |
| Macheret | 2015 | §1.12 — replication stress as a cancer hallmark |
| Srinivasan | 2013 | §1.12 — c-Myc-driven origin firing + fork collapse |
| Schuijers | 2018 | §1.7 — dCas9 epigenome editing at MYC super-enhancer |
| Álvarez-Fernández | 2017 | §1.7 — inducible CRISPR MASTL knockout in breast cancer cells |
| Banerjee | 2020 | §1.9 — multi-region WGS / Darwinian intratumour heterogeneity (CRC) |
| Martelotto | 2014 | §1.9 — review of intratumour heterogeneity in breast cancer |
| McCart | 2019 | §1.2 / §1.9 — breast-cancer genomics + classification |
| Nedeljković | 2019 | §1.9 — clonal heterogeneity in TNBC (background only) |
| Hason | 2019 | §1.3 — zebrafish cancer model review |
| Iñiguez-Muñoz | 2024 | §1.15 — non-coding mutations in regulatory elements |
| Jung | 2021 | §1.15 — pan-cancer WGS of intronic mis-splicing |
| Lee | 2016 | §1.14 / §1.15 — cancer-genome-instability mechanisms via WGS |
| Pagni | 2022 | §1.14 — non-coding regulatory variant interpretation |
| Simpson | 2021 | §1.15 — pan-cancer SV signatures including fragile sites |
| Steele | 2022 | §1.15 — pan-cancer CNV signatures (TCGA 9,873 tumours) |

### Gaps — canonical references absent from the corpus

These are references CC's matching cannot supply because the relevant paper is either not in the bibliography at all, or is in the bibliography but has no abstract (so not in the annotated corpus). Recommend acquiring these manually for the relevant claims.

| Gap | Claim it would support |
|---|---|
| Bagnardi 2015 (or similar alcohol/multi-site cancer meta-analysis) | §1.2 — alcohol as multi-site carcinogen |
| IARC Monographs Volume 44 (1988) or 100E (2012) — alcoholic beverages | §1.2 — IARC 1988 alcohol declaration |
| MacMahon 1970 / Collaborative Group Hormonal Factors 2002 — parity meta-analysis | §1.2 — early childbearing reduces breast-cancer risk (pilot claim 5) |
| Holen 2017 / Park 2018 — animal models of breast cancer review | §1.3 — wide range of animal models |
| Lee 2007 / Goessling 2018 — zebrafish cancer model primary reviews | §1.3 — zebrafish as translational model |
| Ishino 1987 (the original) — *iap* gene paper | §1.5 — Ishino 1987 historical citation (in bibliography but no abstract) |
| Jansen 2002 — original CRISPR-acronym paper | §1.5 — Jansen mobile-element proposal |
| Mojica 2005 / Sander & Joung 2014 — CRISPR-and-rare-disease therapy review | §1.6 — rare diseases policy review |
| PCAWG Consortium 2020 — Nature pan-cancer WGS atlas | §1.15 — pan-cancer non-coding drivers (numbered residue [6]) |

---

## §1.2 Breast Cancer

### Claim 6 — There is much evidence that alcohol is a carcinogen affecting the oesophagus, mouth, breast, colon, and the liver.

**Topic:** meta-analysis of alcohol consumption and multi-site cancer risk

**Suggestions:**

1. **Zhao 2017** [confidence: medium]
   Fit: Mechanistic study of alcohol-induced oxidative DNA damage and p53-dependent DDR in MCF-7 cells. The note flags it as MCF-7 + γH2AX + p53 work, not a multi-site epidemiology meta-analysis. Useful as a *mechanism* citation if §1.2 discusses how alcohol contributes to breast cancer specifically, but does not support the multi-site (oesophagus / mouth / colon / liver) claim.

2. **Yue 2013** [confidence: low] ⚲ already cited in Ch1
   Fit: Reviews ER-dependent and ER-independent mechanisms of breast-cancer carcinogenesis, including genotoxic oestrogen metabolites. Breast-only; does not address the multi-site spectrum.

> **Gap.** The 410-corpus contains no multi-site alcohol-carcinogenicity meta-analysis. A Bagnardi-type meta-analysis (e.g. *Alcohol consumption and site-specific cancer risk: a comprehensive dose-response meta-analysis*, British Journal of Cancer 2015) or the IARC alcohol monograph itself is the proper citation here. Recommend manual addition.

### Claim 7 — Alcohol was declared a carcinogen by the International Agency for Research on Cancer (IARC) in 1988…

**Topic:** IARC Monograph (Volume 44, 1988) classifying alcoholic beverages as carcinogenic

**Suggestions:**

> **No strong matches in the 410-entry corpus.** The IARC 1988 monograph is a primary regulatory document; no corpus paper substitutes for it. Closest candidates considered:
>
> - **Zhao 2017** — alcohol + MCF-7 + DDR; mechanism, not declaration.
> - **Rehm 2020** ⚲ already cited — alcohol-consumption monitoring methodology (the relevance note flags it as having "no relevance to this thesis"; not a useful substitute).
>
> Recommend manual addition of the original IARC Monographs Volume 44 (1988) or the updated Volume 100E (2012) reaffirmation.

---

## §1.3 Animal Models of Breast Cancer

### Claim 8 — There are a wide range of models used for breast cancer research, including mice, rabbits, cats, and rats.

**Topic:** review of animal models used in breast cancer research

**Suggestions:**

1. **Attalla 2021** [confidence: medium] ⚲ already cited in Ch1
   Fit: Reviews the MMTV-PyMT transgenic mouse as the most-used GEMM of breast cancer. Note explicitly recommends it for §1.3 ("animal-models-of-breast-cancer section (1.3) where mouse models are surveyed alongside other organisms"). Mouse-only; doesn't address rabbits, cats, rats.

2. **Hason 2019** [confidence: low]
   Fit: Zebrafish cancer-model review (see claim 10). Could pair with Attalla to broaden the model-organism survey, but the claim specifically lists mammalian models (mice, rabbits, cats, rats) not zebrafish.

> **Gap.** No multi-mammal animal-model review is in the corpus. Recommend adding a Holen 2017 / Park 2018-style breast-cancer animal-models survey.

### Claim 9 — They found that the ketogenic diet induced lower serum insulin and that it could act as an anti-cancer means and that it could increase the effects of rapamycin.

**Topic:** primary study of ketogenic diet plus rapamycin in a mouse breast cancer model

**Suggestions:**

1. **Zou 2020** [confidence: high] ⚲ already cited in Ch1
   Fit: Note explicitly describes "preclinical study examining ketogenic diet and rapamycin/mTOR inhibition in a spontaneous mouse breast-cancer model, demonstrating reduced tumour growth and extended survival." This is precisely the paper the claim points to. (Note's overall relevance assessment to the thesis is "negligible/omit", but for *this specific §1.3 claim* it is the exact target.)

### Claim 10 — Due to the evolutionary conservation of cancer affected programs between zebra fish and Homo Sapiens it is possible to extrapolate research garnered from zebrafish back to humans.

**Topic:** review of zebrafish as a translational model for human cancer biology

**Suggestions:**

1. **Hason 2019** [confidence: high] ⚲ already cited in Ch1
   Fit: Note describes a "review of zebrafish as a cancer model, covering transgenic, transplantation, and xenograft approaches, chemical screening, and in vivo tumour visualisation". Directly on-topic for the §1.3 zebrafish claim.

2. **Lei 2020** [confidence: low] ⚲ already cited in Ch1
   Fit: A specific zebrafish + HUVEC angiogenesis study (note flags it as "no meaningful relevance to this thesis"). Could be cited as an example of zebrafish-derived cancer biology, but Hason 2019 is the cleaner review-level reference.

---

## §1.5 CRISPR and the era of precisely engineered human disease models

### Claim 11 — The term CRISPR has dominated several areas of biomedical and biotechnology in recent times.

**Topic:** broad review of CRISPR's impact on biomedicine and biotechnology

**Suggestions:**

1. **Charpentier 2015** [confidence: medium] ⚲ already cited in Ch1
   Fit: Note describes crRNA biogenesis review across CRISPR-Cas Types I/II/III with explicit framing as the mechanistic precursor to engineered CRISPR/Cas9. Useful as a "how CRISPR became a tool" citation in §1.5.

2. **Ishino 2018** [confidence: medium]
   Fit: Note describes a "concise historical account of CRISPR-Cas systems, tracing the 1987 discovery of the first CRISPR sequence in E. coli through to the development of Cas9-based" — exactly the chronological-impact framing this claim needs.

3. **Choudhury 2016** [confidence: medium]
   Fit: Note describes a CRISPR-dCas9-TET1 engineering application; an example of CRISPR's diversification into epigenome editing supports the "dominated several areas" claim. Not currently cited in Ch1.

### Claim 12 — In 1987, in Japan, the first ever identification of the specific repeats was published based on analysis of the iap gene in Escherichia coli during investigation of the phosphate metabolism.

**Topic:** Ishino 1987 original report of CRISPR repeats in the E. coli iap gene

**Suggestions:**

1. **Ishino 1987** [confidence: high]
   Fit: Note describes it as the "1987 paper characteris[ing] the *E. coli iap* gene and its protein product, notable historically because the sequenced flanking region contained what were later recognised as CRISPR repeat sequences, making it an inadvertent first observation of the CRISPR locus." The exact target paper.
   ⚠ Note this is currently *NOT* among the 152 already-cited refs, but a historical Ishino 1987 citation in Ch1 has already been audited as a known in-text citation. The note text is API-generated from the abstract retrieved in the previous merge; the paper itself is in the bibliography.

2. **Ishino 2018** [confidence: medium]
   Fit: Historical review tracing the 1987 discovery. A complementary "looking back" citation that summarises the original Ishino 1987 contribution; could be cited alongside the primary paper.

### Claim 13 — Ruud Jansen suggested that these structures were frequently found in archaea and bacterial chromosomes in multiple copies and could be mobile elements.

**Topic:** Jansen 2002 paper coining the CRISPR acronym and proposing mobility

**Suggestions:**

1. **Bolotin 2005** [confidence: medium] ⚲ already cited in Ch1
   Fit: Note describes the canonical pre-Cas9 paper showing CRISPR spacers derive from extrachromosomal elements, motivating the adaptive-immunity model. Adjacent to Jansen 2002 in the historical chain ("Cite when tracing the chronology from Ishino 1987 to Jansen 2002 to Bolotin 2005…"). Not a substitute for Jansen 2002 itself but the next link in the chain.

2. **Pourcel 2005** [confidence: medium] ⚲ already cited in Ch1
   Fit: Early CRISPR characterisation in *Yersinia pestis* showing phage origins of spacers. Also adjacent to Jansen 2002.

> **Gap.** The Jansen 2002 paper itself (*Identification of genes that are associated with DNA repeats in prokaryotes*, Mol Microbiol, 43:1565-75) does not appear in the corpus. Recommend manual addition; the chronological CRISPR history hinges on this citation.

---

## §1.6 CRISPR technology and Rare Diseases

### Claim 14 — Rare diseases present a unique problem because their aggregated impact is significant, but due to their individual rarity and low occurrence there is little research and therapeutic development due to lack of resourcing and expertise in any one single disease.

**Topic:** policy review on rare-disease research and therapeutic-development burden

**Suggestions:**

> **No strong matches in the 410-entry corpus.** The 410 references are dominated by cancer biology, replication stress, and CRISPR mechanism. None addresses rare-disease research-and-funding policy. Closest candidates considered:
>
> - **Choudhury 2016** — CRISPR-dCas9 application for BRCA1 demethylation (cancer-specific, not rare-disease).
> - **Solaiman 2021 / Sidali 2023 / Teotia 2024** (predecessor theses) — CRISPR-based investigation of ZFP36L1; not rare-disease policy.
>
> Recommend manual addition of a rare-disease policy review (e.g. Boycott 2017, *Nature Reviews Genetics*; or the EURORDIS / IRDiRC reports).

---

## §1.7 CRISPR and Breast Cancer

### Claim 15 — CRISPR/Cas9 is used widely in cancer research, both for translational and basic research.

**Topic:** review of CRISPR/Cas9 applications in translational cancer research

**Suggestions:**

1. **Choudhury 2016** [confidence: high]
   Fit: Note: "engineers a CRISPR-dCas9-TET1 fusion to drive locus-specific demethylation at the BRCA1 promoter… demonstrates the breadth of CRISPR-Cas9 applications beyond knock-out". Exactly the kind of cancer-research-CRISPR application the claim references. *Not currently cited in Ch1* — high-value reuse.

2. **Teotia 2024** [confidence: high] ⚲ already cited in Ch1 (as predecessor)
   Fit: The MCF-7 B4 clone that this thesis sequences was generated via CRISPR/Cas9 by Teotia. A direct example of CRISPR/Cas9 used translationally in breast-cancer research — the citation is self-referential to the thesis's own experimental design.

3. **Solaiman 2021** [confidence: high] ⚲ already cited in Ch1 (as predecessor)
   Fit: U2OS ZFP36L1-/- KO via CRISPR/Cas9. Equivalent to Teotia for the U2OS arm of the thesis.

4. **Álvarez-Fernández 2017** [confidence: medium]
   Fit: Inducible CRISPR knockout of MASTL in breast-cancer cells — explicit example of CRISPR/Cas9 use in BC research. *Not currently cited in Ch1*.

### Claim 16 — CRISPR/Cas9 has been used to downregulate MYC which is known to have higher expression in high grade BC.

**Topic:** primary study of CRISPR-mediated MYC silencing in breast cancer cells

**Suggestions:**

1. **Schuijers 2018** [confidence: high] ⚲ already cited in Ch1
   Fit: Note: "Mechanistic study of MYC transcriptional dysregulation via cancer-specific super-enhancers… dCas9-DNMT epigenetic editing is used to validate the mechanism." Closest available paper to the claim — uses CRISPR-derived tools to manipulate MYC. The original §1.15 citation was about driver mutations; reuse here for the §1.7 MYC-CRISPR claim is a "free" additional use.

2. **Srinivasan 2013** [confidence: low] ⚲ already cited in Ch1
   Fit: c-Myc-driven replication stress (no CRISPR). Mechanistically about MYC in cancer but doesn't address CRISPR manipulation. Only useful for the "MYC is dysregulated in BC" half of the claim.

### Claim 17 — MASTL kinase activity can be inhibited using CRISPR/Cas9 in human cell lines, which then reduces cell proliferation.

**Topic:** primary study of CRISPR MASTL knockout reducing proliferation in human cells

**Suggestions:**

1. **Álvarez-Fernández 2017** [confidence: high]
   Fit: Note: "This study of MASTL/Greatwall kinase in breast cancer… the inducible CRISPR knockout system in breast cancer cells offering a methodological parallel." Exact target. *Not currently cited in Ch1* — direct new citation reuse.

---

## §1.9 Breast Cancer tumour evolution and clonal variations

### Claim 18 — In a study of 104 triple negative breast cancers (TNBCs) it was shown that there is significant variation in the frequencies of clones.

**Topic:** TNBC clonal-frequency study with an n=104 cohort (likely Shah et al. or similar)

**Suggestions:**

1. **Shah 2012** [confidence: high] ⚲ already cited in Ch1
   Fit: Note: "Landmark WGS-based study demonstrating that primary TNBCs display a continuous spectrum of clonal and mutational evolution at diagnosis, with clonal frequencies of driver mutations… varying widely across cases." This is almost certainly the *exact* paper the n=104 claim references (Shah et al. 2012 sequenced 104 primary TNBCs).

2. **Aftimos 2021** [confidence: medium] ⚲ already cited in Ch1
   Fit: AURORA primary/metastatic BC paired sequencing across 381 patients; documents clonality changes between primary and metastasis. Different cohort/design but adjacent finding.

3. **Banerjee 2020** [confidence: medium]
   Fit: WES of 206 multi-region CRC samples documenting Darwinian intratumour heterogeneity. Disease context (CRC, not TNBC) limits the fit but methodologically parallel. *Not currently cited in Ch1*.

### Claim 19 — The subclones formed due to this somatic heterogeneity have a wide and differing range of biological attributes and capabilities, and subclones with certain advantages will thrive and expand in certain ecosystems…

**Topic:** review of clonal selection and subclonal fitness in tumour evolution

**Suggestions:**

1. **Nowell 1976** [confidence: high] ⚲ already cited in Ch1
   Fit: The conceptual bedrock — "tumours arise clonally and progress through sequential selection of genetically unstable sublines." The ecosystems-and-fitness language in the claim is essentially the Nowell model. Reuse from wherever it's currently cited.

2. **Martelotto 2014** [confidence: high]
   Fit: Note: "Review of intra-tumour heterogeneity in breast cancer, covering clonal evolution, cancer stem cell models, and massively parallel sequencing as a tool for detecting genetic diversity." Directly relevant to §1.9 by the note's own assessment. *Not currently cited in Ch1*.

3. **Banerjee 2020** [confidence: medium]
   Fit: Darwinian intratumour heterogeneity with subclone selection. Methodologically and conceptually on-topic. *Not currently cited in Ch1*.

4. **Shah 2012** [confidence: medium] ⚲ already cited in Ch1
   Fit: TNBC clonal spectrum. Supports the "subclones with differing capabilities" framing.

---

## §1.12 DNA replication and the human genome

### Claim 20 — The fact that genomic instability - a key determining factor of oncogenesis - can be induced by replication stress increases the necessity to understand its underlying mechanisms with more clarity.

**Topic:** review linking replication stress to genomic instability in cancer

**Suggestions:**

1. **Macheret 2015** [confidence: high]
   Fit: Note: "Authoritative review arguing that oncogene-driven DNA replication stress generates genomic instability and selects for apoptosis escape, proposing replication stress as a cancer hallmark in its own right." Direct one-to-one fit. *Not currently cited in Ch1* — major new-citation reuse.

2. **Hills 2014** [confidence: high]
   Fit: Note: "Concise review of how oncogene activation perturbs replication initiation, induces replicative stress and DDR, and drives genomic instability early in tumorigenesis." Directly on-topic. *Not currently cited in Ch1*.

3. **Bartkova 2005** [confidence: high] ⚲ already cited in Ch1
   Fit: Foundational evidence of DDR activation in human precursor lesions implicating replication stress as the engaged stress type. Already a strong §1.8/§1.12 anchor; reuse for this specific instability-from-RS claim is direct.

4. **Bartkova 2006** [confidence: high] ⚲ already cited in Ch1
   Fit: Companion to Bartkova 2005 — oncogene-induced senescence linked to RS-induced DSBs. Same logic.

5. **Briu 2021** [confidence: medium] ⚲ already cited in Ch1
   Fit: Reviews replication-timing / RS / genomic-instability three-way relationship; supports the mechanism-level framing.

### Claim 21 — When levels of cyclin E are elevated, this leads to the obstruction of MCM protein loading, which impedes the formation of pre-replication complexes and also impacts replication initiation in H[uman cells].

**Topic:** primary study of cyclin E overexpression disrupting MCM loading / pre-RC assembly

**Suggestions:**

1. **Ekholm-Reed 2004** [confidence: high]
   Fit: Note: "Mechanistic study showing that cyclin E overexpression in human cells impairs Mcm4 / Mcm7 (and partially Mcm2) loading onto chromatin during telophase / early G1, defectively licensing replication origins." Precisely the paper this claim describes. *Not currently cited in Ch1* — direct new-citation reuse.

2. **Pruitt 2007** [confidence: medium] ⚲ already cited in Ch1
   Fit: MCM2 hypomorph causes chronic RS and tumorigenesis — complementary mechanism (MCM insufficiency rather than cyclin E excess but same downstream consequence). Useful as a paired citation.

3. **Arentson 2002** [confidence: medium] ⚲ already cited in Ch1
   Fit: CDT1 over-licensing oncogene biology — same pre-RC dysregulation theme, different licensing factor.

4. **Seo 2005** [confidence: medium] ⚲ already cited in Ch1
   Fit: CDT1 overexpression drives genomic instability — also a licensing-factor-excess paper analogous in mechanism to cyclin E-driven licensing failure.

### Claim 22 — More recent investigations suggest that the enrichment of yeast heterochromatic DSB marker H2A (H2AX) depends on replication.

**Topic:** yeast study of replication-dependent γH2AX enrichment at heterochromatin

**Suggestions:**

> **No strong matches in the 410-entry corpus.** The corpus is mammalian-cell-heavy. Closest candidates considered:
>
> - **Taymaz-Nikerel 2018** — yeast doxorubicin response with Rad53 activation; addresses yeast DDR but not γH2AX-at-heterochromatin specifically.
> - **Zhu 2017** — yeast rare-variant proxy study for mutational spectra; methodological yeast paper, not γH2AX biology.
> - **Sidali 2023** ⚲ already cited (predecessor) — mammalian γH2AX work, not yeast.
>
> Recommend manual addition of a *Schizosaccharomyces pombe* or *Saccharomyces cerevisiae* γH2AX paper (e.g. Kim 2007, *Nature*; or Chambers & Downs 2012).

---

## §1.14 Whole Genome Sequencing and Genome Instability under Replication Stress

### Claim 23 — Although the technology is aiding science and research in immeasurable ways, one issue with WGS is that the sheer abundance of data can mean that certain discovered variants are poorly understood or not understood at all at the moment.

**Topic:** review of variants-of-unknown-significance interpretation in clinical WGS

**Suggestions:**

1. **Foley 2015** [confidence: high] ⚲ already cited in Ch1
   Fit: Note: "highlighting the challenge that the majority of novel missense variants are classified as variants of unknown significance (VUS) and that large-scale data consolidation is required for accurate interpretation." Direct one-to-one fit to the VUS-burden claim.

2. **Lee 2016** [confidence: medium]
   Fit: Comprehensive review of cancer-genome-instability mechanisms via WGS/WES — provides the broader context the VUS-interpretation problem sits within. *Not currently cited in Ch1*.

3. **Pagni 2022** [confidence: low]
   Fit: Non-coding regulatory variant review — the VUS problem is particularly acute for non-coding variants. Loose fit. *Not currently cited in Ch1*.

4. **Iñiguez-Muñoz 2024** [confidence: low]
   Fit: Somatic + germline non-coding variants in regulatory elements; adjacent to the VUS-interpretation challenge. *Not currently cited in Ch1*.

---

## §1.15 WGS in Cancer Research

### Claim 24 — A major application of WGS in cancer research is the comprehensive assessment of mutational burden across the genome.

**Topic:** review of WGS-derived tumour mutational burden in oncology

**Suggestions:**

1. **Fumet 2020** [confidence: high] ⚲ already cited in Ch1
   Fit: Note explicitly reviews TMB as a biomarker including assessment methods and quantification approaches. Directly on the claim.

2. **Kandoth 2013** [confidence: high] ⚲ already cited in Ch1
   Fit: Pan-Cancer somatic-mutation atlas across 3,281 TCGA tumours — the foundational TMB-quantification resource the claim implicitly references.

3. **Kim 2024** [confidence: medium] ⚲ already cited in Ch1
   Fit: Prospective WGS study across 120 solid-tumour patients reporting TMB and actionability — a clinical-application example complementing the methodology citations.

4. **Lee 2016** [confidence: medium]
   Fit: Comprehensive review of WGS-detected mutational and chromosomal instability — provides the broader conceptual frame for TMB. *Not currently cited in Ch1*.

5. **Aftimos 2021** [confidence: medium] ⚲ already cited in Ch1
   Fit: AURORA's TMB-shorter-relapse finding is a useful BC-specific illustration of the mutational-burden claim.

### Claim 25 — Driver mutations confer a selective advantage and contribute directly to tumor development or progression, whereas passenger mutations accumulate without directly influencing cellular fitness [1].

**Topic:** review distinguishing driver vs passenger mutations in cancer genomes _(numbered residue [1] — convert to Harvard)_

**Suggestions:**

1. **Kandoth 2013** [confidence: high] ⚲ already cited in Ch1
   Fit: Pan-Cancer catalogue of 127 significantly mutated genes (i.e. drivers) across 3,281 tumours. The canonical driver-identification paper of its era.

2. **Beroukhim 2010** [confidence: high]
   Fit: Note: "pan-cancer analysis of focal somatic copy-number alterations across 3,131 specimens identifying 158 recurrently altered regions, of which 122 lack a known cancer-driver gene". Explicit driver-vs-passenger framing for CNV calls. *Not currently cited in Ch1*.

3. **Ellis 2013** [confidence: high] ⚲ already cited in Ch1
   Fit: Cancer Discovery perspective on the breast-cancer omics revolution — synthesises the driver-discovery progress across multiple landmark studies.

4. **Macheret 2015** [confidence: medium]
   Fit: Replication-stress-as-hallmark framing helps explain *why* passenger mutations accumulate (RS-mediated incidental damage). *Not currently cited in Ch1*.

### Claim 26 — Large-scale WGS studies have identified characteristic patterns of driver mutations across cancer types, revealing both frequent alterations in canonical oncogenes and tumor suppressors, as well as unexpected drivers in non-coding regions and regulatory elements that influence cancer development [6].

**Topic:** pan-cancer WGS atlas of driver alterations including non-coding regions (e.g. PCAWG) _(numbered residue [6] — convert to Harvard)_

**Suggestions:**

1. **Kandoth 2013** [confidence: high] ⚲ already cited in Ch1
   Fit: Pan-Cancer TCGA driver atlas — canonical citation for the pan-cancer-WGS-revealed-drivers framing.

2. **Beroukhim 2010** [confidence: high]
   Fit: Pan-cancer CNV atlas across 3,131 specimens including unexplained recurrent regions (i.e. potential non-coding-driver loci). *Not currently cited in Ch1*.

3. **Steele 2022** [confidence: high]
   Fit: Note: "Foundational resource for the WGS and copy-number analysis chapters. This pan-cancer framework of 21 copy number signatures, derived from 9,873 TCGA tumours across 33 cancer types using WGS, whole-exome, and SNP-array data." Recent pan-cancer atlas-level reference. *Not currently cited in Ch1*.

4. **Simpson 2021** [confidence: high]
   Fit: Pan-cancer SV signatures across 38 tumour types from WGS, including fragile-site SVs — particularly relevant given the thesis's CFS theme. *Not currently cited in Ch1*.

5. **Iñiguez-Muñoz 2024** [confidence: medium]
   Fit: Reviews somatic + germline non-coding mutations in regulatory elements — directly supports the "non-coding regions and regulatory elements" half of the claim. *Not currently cited in Ch1*.

6. **Jung 2021** [confidence: medium]
   Fit: Pan-cancer WGS of 1,134 genomes identifying 678 somatic intronic mis-splicing mutations — non-coding-driver evidence. *Not currently cited in Ch1*.

> **Gap.** The actual PCAWG flagship paper (*Pan-cancer analysis of whole genomes*, Nature 578, 2020) is the canonical citation for this claim and is *not* in the corpus. Recommend manual addition.
