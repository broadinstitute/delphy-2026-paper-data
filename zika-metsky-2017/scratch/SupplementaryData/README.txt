Description of files:
---------------------

Beast input and output/
    Ancestral state reconstruction/
        HKY-relaxed-skyline.xml
            - BEAST XML used as input for ancestral state reconstruction
        HKY-relaxed-skyline.ancestral-state.cds-and-noncoding.state_24725000.fasta
            - Ancestral state constructed from output of the state with the MCC
            tree
    Phylogenetic analyses and model selection/
        SRD06-strict-constant.xml
            - BEAST XML used as input for model selection
            - Strict clock, constant size population
        SRD06-relaxed-constant.xml
            - BEAST XML used as input for model selection
            - Relaxed clock, constant size population
        SRD06-strict-exponential.xml
            - BEAST XML used as input for model selection
            - Strict clock, exponential growth
        SRD06-relaxed-exponential.xml
            - BEAST XML used as input for model selection
            - Relaxed clock, exponential growth
        SRD06-strict-skyline.xml
            - BEAST XML used as input for model selection and for results that
            show concordance with SRD06-relaxed-skyline
            - Strict clock, Skyline tree prior
        SRD06-strict-skyline.log
            - Output of BEAST run with SRD06-strict-skyline.xml
        SRD06-relaxed-skyline.xml
            - BEAST XML used as input for model selection, for MCC tree shown,
            and for other key results that are highlighted (e.g., tMRCA estimates)
            - Relaxed clock, Skyline tree prior
        SRD06-relaxed-skyline.log
            - Output of BEAST run with SRD06-relaxed-skyline.xml
        SRD06-relaxed-skyline.summary.formatted.tree
            - Tree file output by TreeAnnotator after discarding burn-in
    Substitution rates estimation/
        GTR-relaxed-skyline.xml
            - BEAST XML file used as input for estimating A-C, A-G, A-T,
            C-G, C-T, and G-T substitution rates

Root-to-tip/
    ml-tree.GTR-G4.tree
        - Tree file output by PhyML, which was used as the ML tree
    root-to-tip-data.tsv
        - Root-to-tip divergence, dates, and residuals produced by TempEst; the
        first two of these values were used in producing the root-to-tip plot

Sequences and alignments/
    sequences-generated-in-this-study.fasta
        - FASTA (unaligned) of sequences produced in this study that were selected
        for analyses
    alignment-used-for-beast-analyses.no-outgroup.fasta
        - FASTA (aligned) with all sequences used in BEAST analyses
        - In most BEAST analyses only the CDS was used (see Methods for details)
    alignment-used-for-ml-tree.with-outgroup.fasta
        - FASTA (aligned) with all sequences, including outgroup, used to produce
        ML tree
        - BEAST analyses that included an outgroup also came from the CDS of these
        sequences
    unique-viral-contigs-from-mosquito-pools.fasta
        - FASTA of de-duplicated (unique) viral contigs from mosquito pools
