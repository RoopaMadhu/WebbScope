"""Curated gene-signature definitions. Each is a candidate gene list; only genes present on the Xenium panel are used (see signatures.py).

The lists are curated from the literature and standard gene sets (MSigDB Hallmark, Schoggins et al. 2011 ISG screens, Interferome, kidney
single-cell atlases); they are NOT derived from this dataset, so review them before relying on a result. `note` states the limits of each.
"""
# Genes removed from EVERY signature and every multi-gene score / analysis (excluded by the user: CCL5 is an unreliable gene here).
EXCLUDED_GENES = frozenset({"CCL5"})
# Older / newer HGNC symbols: the panel may use either, so a candidate is looked up under both.
ALT_SYMBOLS = {"WARS1": "WARS", "STING1": "TMEM173", "SEPTIN4": "SEPT4", "PLAAT4": "RARRES3", "EPRS1": "EPRS", "CGAS": "MB21D1", "TNFSF13B": "BAFF"}
ALT_SYMBOLS.update({v: k for k, v in list(ALT_SYMBOLS.items())})
ALIASES = {"CTGF": "CCN2", "WARS": "WARS1", "HLA-DRA": "HLA-DRA", "NKG2D": "KLRK1", "PD1": "PDCD1", "PDL1": "CD274", "IL28RA": "IFNLR1",
           "SERPINH1": "SERPINH1", "MARCKS": "MARCKS", "TMEM173": "STING1", "LGALS3BP": "LGALS3BP", "IFNB": "IFNB1"}

# --- interferon-stimulated genes ----------------------------------------------------------------------------------------------------------
# Interferon ligands, receptors and JAK kinases are made constitutively (not induced), so they form a separate pathway signature.
IFN_PATHWAY = ["IFNA1", "IFNA2", "IFNB1", "IFNAR1", "IFNAR2", "IFNG", "IFNGR1", "IFNGR2", "IFNL1", "IFNL2", "IFNL3", "IFNL4", "IFNLR1", "IL10RB", "JAK1", "JAK2", "TYK2"]
IFN_III_PATHWAY = ["IFNL1", "IFNL2", "IFNL3", "IFNL4", "IFNLR1", "IL10RB"]
# induced by all three interferon types
IFN_SHARED_CORE = ["STAT1", "STAT2", "IRF9", "IRF1", "IRF7", "ISG15", "IFIT1", "IFIT3", "IFI35", "MX1", "OAS1", "PLSCR1", "IFITM1", "IFITM2", "XAF1", "TAP1",
                   "PSMB9", "SOCS1", "BST2", "SP100", "PARP14", "EPSTI1", "IFI16", "TRIM21", "TRIM22", "USP18", "IFI6", "IFI27", "RSAD2", "HERC6"]
_I = ["ISG15", "IFI6", "IFI27", "IFI44", "IFI44L", "IFIT1", "IFIT2", "IFIT3", "IFIT5", "MX1", "MX2", "OAS1", "OAS2", "OAS3", "OASL", "RSAD2", "HERC5", "HERC6",
      "USP18", "XAF1", "EPSTI1", "LY6E", "SAMD9", "SAMD9L", "DDX58", "IFIH1", "DDX60", "TRIM22", "CMPK2", "SIGLEC1", "IFI35", "BST2", "SLFN5", "SLFN11", "SLFN12",
      "MOV10", "ZBP1", "APOBEC3A", "APOBEC3G", "PARP9", "PARP12", "DTX3L", "SP100", "NMI", "IFITM3", "ISG20", "GMPR", "RTP4",
      "DHX58", "CH25H", "HSH2D", "IFI27L2", "NT5C3A", "RNF213", "UBA7", "UBE2L6", "USP41", "IFIT1B"]
_II = ["GBP1", "GBP2", "GBP3", "GBP4", "GBP5", "GBP6", "GBP7", "IDO1", "CXCL9", "CXCL10", "CXCL11", "IRF1", "SOCS3", "CIITA", "NLRC5", "TAP2", "PSMB8", "PSME1", "PSME2",
       "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1", "CD74", "FCGR1A", "WARS1", "CD274", "PDCD1LG2", "APOL1", "APOL2",
       "APOL6", "LAP3", "ETV7", "IL15", "IL15RA", "UBD", "CASP1", "CASP4", "SERPING1", "TNFSF10", "IFI30", "B2M", "HLA-A", "HLA-B", "HLA-C", "HLA-E",
       "TRIM69", "RARRES3", "PLAAT4", "GCH1", "BATF2", "FGL2", "IRF8", "IRF2", "TAPBP", "ELF1", "LGALS9", "SECTM1", "VAMP5", "VAMP8", "SLC15A3"]
_ISG_EXTRA = ["ADAR", "SAMHD1", "EIF2AK2", "NAMPT", "OGFR", "MVB12A", "CASP7", "CASP8", "DDX41", "STING1", "CGAS", "TLR3", "TLR7", "TLR9", "MYD88", "IRF3", "IRF5", "CD38",
              "CD69", "LGALS3BP", "ICAM1", "VCAM1", "CYBB", "NOS2", "CEACAM1", "TNFSF13B", "SOCS3", "PML", "OPTN", "TRIM25", "LAMP3", "TDRD7", "C1S", "C1R", "CFB", "C2", "C3", "NCOA7", "TNFAIP6", "FAM111A", "HELZ2", "SP110", "SP140", "TRIM5", "TRIM14", "TRIM26", "TRIM38", "PARP10",
              "PNPT1", "RBCK1", "ZNFX1", "GVINP1", "IRGM", "MNDA", "IL7", "IL4I1", "SLFN13", "ISG20L2", "TXNIP", "IFITM3"]
_exclude = set(IFN_PATHWAY)
_shared = set(IFN_SHARED_CORE)
IFN_I = [g for g in _I if g not in _shared and g not in set(_II) and g not in _exclude]      # type I only: not shared, not claimed by type II
IFN_II = [g for g in _II if g not in _shared and g not in set(_I) and g not in _exclude]     # type II only
IFN_MERGED = list(dict.fromkeys(IFN_I + IFN_II + IFN_SHARED_CORE + IFN_III_PATHWAY))
ISG_ALL = list(dict.fromkeys(_I + _II + IFN_SHARED_CORE + _ISG_EXTRA))

# User-provided list (SLE ISG amplification); "4-Sep" is the spreadsheet mangling of SEPTIN4, duplicates are removed on resolve.
ISG_AMP_SLE = ["BATF2", "IFI44L", "XAF1", "OASL", "IFI44", "OAS2", "TRIM6", "HES4", "OTOF", "IFI44L", "DDX60", "IFITM3", "IFIT3", "CXCL10", "SERPING1", "SPATS2L", "RSAD2", "LY6E", "OAS1", "IFIT1", "IFIT3", "OAS3", "OAS1", "OASL", "CMPK2", "HERC5", "OAS1", "MX1", "4-Sep", "HERC6", "IFITM1", "TRIM5", "EIF2AK2", "AIM2", "MT1A", "MOV10", "CCL8", "HELZ2", "ZBP1", "WARS", "LAP3", "TNFSF10", "GBP1", "STAT1", "FBXO6", "PARP10", "OAS2", "TRIM22", "PARP10", "GBP3", "ZNF684", "CARD17", "GALM", "DHX58", "CEACAM1", "UBE2L6", "CEACAM1", "APOL6", "SOCS1", "DDX58", "UBE2L6", "TNFAIP6", "IDO1", "GBP6", "STAT2", "TIMM10", "STAT1", "PARP12", "PARP9", "GBP1", "GBP4", "IFIH1", "IRF7", "PARP14", "IFIT2", "IFI35", "STAT1", "IFIT5", "ZNFX1", "PHF11", "ACTA2", "SP140", "ABCA1", "TCN2", "ZC3HAV1", "HSH2D", "GBP2", "TRIM5", "RHBDF2", "TMEM140", "ADAR", "BTN3A1", "LGALS9", "NBN", "TYMP", "SRBD1", "NCOA7", "UNC93B1", "SP100", "DHRS9", "TRANK1", "TRIM5", "NT5C3A", "ASPRV1", "IRF9", "REC8", "RNF213", "ISG20", "DYNLT1", "TRIM56", "SP140", "TRIM38", "PSMB9", "CPT1B", "BST2", "NMI", "RBCK1", "TRAFD1", "TRIM21", "CARD16", "CHMP5", "TAP2", "SP110", "GADD45B", "TAP1"]
ALIASES["4-Sep"] = "SEPTIN4"

SIGNATURES = [
    # id, group, name, description, genes, note
    dict(id="isg_all", group="Interferon", name="ISG (all interferon-stimulated genes)",
         desc="Every interferon-stimulated gene on the panel: type I, type II, shared-core and further ISGs (receptors / ligands excluded).",
         genes=ISG_ALL, note="The broadest ISG set; includes innate nucleic-acid sensors and antigen-processing genes."),
    dict(id="isg_amp_sle", group="Interferon", name="isg_amp_sle (SLE ISG amplification)",
         desc="Interferon-stimulated genes of the SLE ISG-amplification list supplied by you (duplicates removed).",
         genes=ISG_AMP_SLE, note="User-provided gene list; only the genes on the panel are scored (see the list of missing genes)."),
    dict(id="ifn1", group="Interferon", name="IFN type I (α/β) response",
         desc="ISGs preferentially induced by type I interferons (IFN-α/β), without genes shared by all types or claimed by IFN-γ.",
         genes=IFN_I, note="Many of the strongest type I ISGs (ISG15, IFI6, IFI27, IFIT1/3, MX1, OAS1, RSAD2) are induced by every interferon type, so they sit in the shared core, not here."),
    dict(id="ifn2", group="Interferon", name="IFN type II (γ) response",
         desc="Genes preferentially induced by IFN-γ: GBPs, IDO1, CXCL9-11, antigen processing and presentation, CIITA, NLRC5.",
         genes=IFN_II, note="The panel has no HLA genes, so the antigen-presentation arm of the response is under-represented."),
    dict(id="ifn3", group="Interferon", name="IFN type III (λ) pathway",
         desc="IFN-λ ligands and the IFNLR1 / IL10RB receptor chain.",
         genes=IFN_III_PATHWAY, note="The transcriptional response to IFN-λ is largely the same ISG set as type I, so a response signature specific to type III cannot be defined from gene lists; this is a ligand/receptor (pathway) signature and is weakly informative."),
    dict(id="ifn_core", group="Interferon", name="IFN shared core (induced by all types)",
         desc="ISGs induced by type I, II and III interferons alike.", genes=IFN_SHARED_CORE, note="The common-response reference when contrasting type I with type II."),
    dict(id="ifn_merged", group="Interferon", name="IFN merged (type I + II + III + shared core)",
         desc="Union of the type I, type II, type III and shared-core signatures.", genes=IFN_MERGED, note=""),
    dict(id="ifn_pathway", group="Interferon", name="IFN pathway components (ligands, receptors, JAK)",
         desc="Interferon ligands, receptors and JAK kinases (made constitutively; not a response signature).", genes=IFN_PATHWAY, note=""),
    # --- fibrosis ----------------------------------------------------------------------------------------------------------------------------
    dict(id="fibrosis_ecm", group="Fibrosis", name="Fibrotic ECM deposition",
         desc="Fibrillar and associated extracellular-matrix genes laid down in fibrosis.", genes=["COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL5A2", "COL6A1", "COL6A2", "COL6A3", "COL8A1", "COL11A1", "COL12A1", "COL14A1", "FN1", "TNC", "POSTN", "LUM", "DCN", "BGN", "FBN1", "SPARC", "VCAN", "THBS1", "THBS2", "SERPINH1", "ELN", "MGP", "FBLN1", "FBLN2", "FBLN5", "CTSK", "COMP"], note=""),
    dict(id="myofib", group="Fibrosis", name="Myofibroblast activation",
         desc="Activated fibroblast / myofibroblast programme.", genes=["ACTA2", "TAGLN", "CTHRC1", "FAP", "POSTN", "PDGFRB", "PDGFRA", "COL1A1", "COL1A2", "CCN2", "TGFBI", "LOX", "LOXL2", "SERPINE1", "MYL9", "TPM2", "CNN1", "MYH11", "ADAM12", "GREM1", "WNT5A", "INHBA", "SFRP2"], note="ACTA2/TAGLN/MYH11 are also high in vascular smooth muscle."),
    dict(id="tgfb", group="Fibrosis", name="TGF-β signalling response",
         desc="TGF-β ligands, receptors, SMADs and canonical targets.", genes=["TGFB1", "TGFB2", "TGFB3", "TGFBR1", "TGFBR2", "TGFBR3", "SMAD2", "SMAD3", "SMAD4", "SMAD7", "SERPINE1", "CCN2", "CCN1", "TGFBI", "THBS1", "ITGB6", "ITGAV", "LTBP1", "LTBP2", "LTBP4", "INHBA", "FST", "ID1", "ID2", "ID3", "SKIL", "TGIF1", "PMEPA1", "NOX4", "COL1A1"], note=""),
    dict(id="ecm_remodel", group="Fibrosis", name="ECM remodelling (MMP / TIMP / LOX)",
         desc="Matrix-degrading enzymes, their inhibitors and cross-linking enzymes.", genes=["MMP2", "MMP7", "MMP9", "MMP11", "MMP12", "MMP14", "MMP19", "TIMP1", "TIMP2", "TIMP3", "LOX", "LOXL1", "LOXL2", "LOXL3", "ADAMTS1", "ADAMTS4", "ADAMTS5", "ADAM12", "ADAM17", "CTSB", "CTSD", "CTSK", "CTSL", "PLAU", "PLAUR", "SERPINE1", "SERPINE2", "SERPINA1", "HTRA1"], note=""),
    dict(id="emt", group="Fibrosis", name="Epithelial-mesenchymal transition",
         desc="EMT / partial-EMT genes.", genes=["VIM", "SNAI1", "SNAI2", "ZEB1", "ZEB2", "TWIST1", "TWIST2", "CDH2", "FN1", "ITGA5", "ITGB1", "MMP2", "MMP14", "TGFBI", "SPARC", "COL1A1", "COL3A1", "CCN2", "CCN1", "LGALS1", "S100A4", "TNC", "CD44", "PRRX1"], note=""),
    # --- inflammation / immunity ------------------------------------------------------------------------------------------------------------
    dict(id="complement", group="Immune & inflammation", name="Complement",
         desc="Classical, alternative and terminal complement components, receptors and regulators.", genes=["C1QA", "C1QB", "C1QC", "C1R", "C1S", "C2", "C3", "C4A", "C4B", "C5", "C6", "C7", "C8A", "C9", "CFB", "CFD", "CFH", "CFI", "CFP", "C3AR1", "C5AR1", "C5AR2", "CR1", "CR2", "ITGAM", "ITGAX", "CD55", "CD46", "CD59", "SERPING1", "MASP1", "MASP2", "MBL2", "FCN1", "FCN3"], note=""),
    dict(id="nfkb_inflam", group="Immune & inflammation", name="Inflammatory (NF-κB / TNF) response",
         desc="NF-κB targets and pro-inflammatory cytokines / chemokines.", genes=["TNF", "IL1B", "IL1A", "IL6", "CXCL8", "CXCL1", "CXCL2", "CXCL3", "CCL2", "CCL20", "PTGS2", "NFKB1", "NFKB2", "NFKBIA", "NFKBIZ", "RELB", "TNFAIP3", "TNFAIP2", "IER3", "ICAM1", "VCAM1", "SELE", "BIRC3", "TRAF1", "CD83", "IRAK2", "ZC3H12A", "SOD2", "CSF1", "CSF2", "LIF", "OSM"], note=""),
    dict(id="t_cytotoxic", group="Immune & inflammation", name="Cytotoxic T / NK",
         desc="Cytotoxic effector genes.", genes=["GZMA", "GZMB", "GZMH", "GZMK", "GZMM", "PRF1", "NKG7", "GNLY", "CD8A", "CD8B", "CCL5", "CCL4", "KLRD1", "KLRK1", "KLRB1", "KLRC1", "CST7", "FGFBP2", "CTSW", "FASLG", "EOMES", "TBX21", "CX3CR1"], note="Mostly expressed by T/NK cells, so the score tracks their abundance."),
    dict(id="t_exhaust", group="Immune & inflammation", name="T-cell activation / exhaustion",
         desc="Checkpoint and activation markers.", genes=["PDCD1", "CTLA4", "LAG3", "TIGIT", "HAVCR2", "TOX", "TOX2", "ICOS", "CD27", "TNFRSF9", "TNFRSF18", "TNFRSF4", "CD28", "IL2RA", "CXCL13", "ENTPD1", "BATF", "IRF4", "PRDM1", "CD69"], note=""),
    dict(id="b_plasma", group="Immune & inflammation", name="B cell / plasma cell",
         desc="B-cell identity and plasma-cell secretory programme.", genes=["MS4A1", "CD79A", "CD79B", "CD19", "CD22", "CD37", "BANK1", "CR2", "FCRL1", "FCRL2", "FCRL5", "CD72", "MZB1", "JCHAIN", "SDC1", "XBP1", "PRDM1", "TNFRSF17", "IGHG1", "IGKC", "IGHA1", "IGHM", "CXCL13", "CD27", "CD38", "TNFRSF13B", "TNFRSF13C"], note="Mostly B/plasma cells."),
    dict(id="macrophage", group="Immune & inflammation", name="Macrophage / monocyte",
         desc="Myeloid identity and macrophage activation.", genes=["CD68", "CD163", "MRC1", "MSR1", "CD14", "CSF1R", "C1QA", "C1QB", "C1QC", "APOE", "TREM2", "MARCO", "FCGR3A", "FCGR2A", "LYZ", "AIF1", "LST1", "CTSS", "MS4A7", "SPI1", "SIGLEC1", "VSIG4", "CD86", "ITGAM", "S100A8", "S100A9", "FCER1G", "TYROBP"], note="Mostly myeloid cells."),
    dict(id="mhc_present", group="Immune & inflammation", name="Antigen processing & presentation",
         desc="MHC class I/II machinery and peptide processing.", genes=["B2M", "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1", "CD74", "CIITA", "TAP1", "TAP2", "TAPBP", "PSMB8", "PSMB9", "PSME1", "PSME2", "NLRC5", "CTSS", "IFI30", "CD80", "CD86", "CD83"], note="The panel contains no HLA genes, so this score rests on the non-HLA genes."),
    # --- kidney programmes -------------------------------------------------------------------------------------------------------------------
    dict(id="pt_injury", group="Kidney programmes", name="Proximal tubule injury / dedifferentiation",
         desc="Injury and repair state of the proximal tubule.", genes=["HAVCR1", "LCN2", "VCAM1", "SPP1", "KRT8", "KRT18", "KRT19", "SOX9", "CLU", "ITGB6", "MMP7", "PROM1", "CD24", "TNFRSF12A", "ICAM1", "CCL2", "CXCL1", "CXCL2", "CXCL8", "SERPINA1", "SLPI", "ANXA1", "ANXA2", "S100A6", "VIM", "CDH6", "DCDC2", "GPX3", "LRP2"], note="Contains genes (e.g. LRP2, GPX3) that fall with injury; the score is the mean, so mixed directions dilute it."),
    dict(id="pt_identity", group="Kidney programmes", name="Proximal tubule identity",
         desc="Differentiated proximal tubule transport genes (lost in injury).", genes=["LRP2", "CUBN", "SLC34A1", "SLC22A8", "SLC22A6", "SLC22A2", "SLC13A3", "SLC5A2", "SLC5A12", "SLC4A4", "SLC7A13", "SLC12A1", "PCK1", "FBP1", "GATM", "ALDOB", "HNF4A", "HNF1A", "KL", "MIOX", "CYP27B1", "ACSM2A", "ACSM2B", "PDZK1", "AQP1", "UMOD", "SLC27A2"], note="UMOD/SLC12A1 belong to the thick ascending limb and are dropped if absent from the panel; review the gene list."),
    dict(id="podocyte", group="Kidney programmes", name="Podocyte",
         desc="Podocyte identity and slit-diaphragm genes.", genes=["NPHS1", "NPHS2", "PODXL", "WT1", "SYNPO", "PTPRO", "CLIC5", "MAFB", "CD2AP", "ACTN4", "TJP1", "NES", "PLCE1", "TCF21", "MAGI2", "DDN", "EMX2", "COL4A3", "COL4A4", "COL4A5", "LAMB2", "NEPH1", "KIRREL1", "CDKN1C", "ZNF804A"], note=""),
    dict(id="endo_activ", group="Kidney programmes", name="Endothelial activation",
         desc="Activated / inflamed endothelium.", genes=["ICAM1", "VCAM1", "SELE", "SELP", "ACKR1", "VWF", "PLVAP", "EDN1", "ANGPT2", "CCL2", "CXCL8", "IL6", "PECAM1", "CD34", "ESAM", "EMCN", "KDR", "FLT1", "TEK", "NOS3", "ICAM2", "CD74", "HLA-DRA", "IL33", "CLDN5", "VEGFA", "ENG", "CAV1"], note="Mixed identity (PECAM1, CD34) and activation genes."),
    dict(id="senescence", group="Kidney programmes", name="Senescence / SASP",
         desc="Cell-cycle arrest and senescence-associated secretory phenotype.", genes=["CDKN1A", "CDKN2A", "CDKN2B", "TP53", "GLB1", "SERPINE1", "IL6", "CXCL8", "CXCL1", "CCL2", "MMP3", "MMP1", "IGFBP3", "IGFBP7", "IGFBP5", "GDF15", "LMNB1", "TNFRSF10B", "PLAU", "FN1", "CCN2", "CCN1", "SPP1"], note=""),
    dict(id="proliferation", group="Kidney programmes", name="Proliferation",
         desc="Cell-cycle / proliferation genes.", genes=["MKI67", "TOP2A", "PCNA", "CCNB1", "CCNB2", "CCNA2", "CDK1", "CDC20", "BIRC5", "AURKA", "AURKB", "PLK1", "TYMS", "MCM2", "MCM6", "UBE2C", "CENPF", "NUSAP1", "HMGB2", "STMN1", "RRM2", "CDKN3", "TPX2"], note=""),
    dict(id="hypoxia", group="Kidney programmes", name="Hypoxia / stress",
         desc="HIF targets and cellular stress genes.", genes=["HIF1A", "EPAS1", "VEGFA", "SLC2A1", "CA9", "BNIP3", "PDK1", "LDHA", "ENO1", "PGK1", "ALDOA", "ADM", "NDRG1", "ANGPTL4", "EGLN3", "HMOX1", "FOS", "JUN", "ATF3", "DDIT3", "DDIT4", "HSPA1A", "HSPA1B", "DNAJB1"], note=""),
    dict(id="apoptosis", group="Kidney programmes", name="Apoptosis / cell death",
         desc="Apoptosis and regulated cell-death genes.", genes=["CASP3", "CASP7", "CASP8", "CASP9", "BAX", "BAK1", "BCL2", "BCL2L1", "MCL1", "BID", "BAD", "FAS", "FASLG", "TNFRSF10A", "TNFRSF10B", "TNFSF10", "TP53", "PMAIP1", "BBC3", "GSDMD", "MLKL", "RIPK1", "RIPK3", "XIAP", "BIRC2", "BIRC3"], note=""),
]
