# CBASS_Replicability

We are often in a situation where we want to know how well results from one CBASS run are captured in a second run, e.g. are the top- and bottom-ranking colonies consistent between seasons, treatments, species, etc. .. the scripts here help with that in that they assess 1. linear correlation between 2 sets of ED50 values, 2. number of significant pairwise comparisions between subsets of original (how sensitive to no. of replicates?), and 3. number of falsely assigned top- or bottom-ranking colonies (e.g., to control for robustness in screening for thermally tolerant colonies)


Project folder structure:

```text
.
├── input_example
│   ├── example.csv
│   └── example2.csv
├── old_logic
│   ├── combinations.py
│   ├── Top5_Bottom5_analysis_djb.R
│   └── README.md
├── plots
│   ├── example_replicability_1000reps_all.pdf      <-- From running TopVsBottom_analysis.r
│   ├── example2_replicability_1000reps_all.pdf     <-- From running TopVsBottom_analysis.r
│   ├── subsets_agreement_check.pdf                 <-- From running subsets_agreement_check.r
│   ├── scatter_example2.png                        <-- From running subsets_agreement_check.r
│   ├── scatter_example.png                         <-- From running subsets_agreement_check.r
│   └── README.md
├── TopVsBottom_analysis.r
├── subsets_agreement_check.r
└── README.md
```