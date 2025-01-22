This repository collates code and data to accompany the manuscript **Probabilistic classification of late treatment failure in uncomplicated malaria** by Mehra et al, 2025. This is available as a preprint on medRxiv: https://www.medrxiv.org/content/10.1101/2025.01.21.25320790v1.

The results rely in part on our classification framework *PfRecur*, which has been implemented as an R package, openly available in a GitHub repository: https://github.com/somyamehra/PfRecur. Code and data derived from external sources are highlighted below.

The directory structure is as follows:
* *Angola_TES_data*: this contains data from a recent therapeutic efficacy study conducted in Angola in 2021, published in [Dimbu et al (2024)](https://journals.asm.org/doi/full/10.1128/aac.01525-23); genotypic data (*Dimbu_Angola_TES_genotypes.xlsx*) have been retrieved from an accompanying [GitHub repository](https://github.com/MateuszPlucinski/AngolaTES2021), while clinical metadata (*Dimbu_Angola_TES_metadata.xlsx*) are taken from Supplemental Table S4 of Dimbu et al (2024).
* *Angola_TES_reanalysis*: this contains a re-analysis of the data of Dimbu et al (2024) using our classification framework *PfRecur*.
* *Simulation_model*: this contains code implementing our proposed simulation model for mixtures of recrudescent and newly-inoculated clones, and the classification of simulated data using *PfRecur*.
* *Plucinski_et_al_model*: this contains code and results pertaining to the classification model proposed by [Plucinski et al (2015)](https://journals.asm.org/doi/full/10.1128/aac.00072-15); the file *Plucinski_et_al_code.R* has largely been adapted from a R Markdown file retrieved from an accompanying  [GitHub repository](https://github.com/MateuszPlucinski/AngolaTES2021), with very minor input/output modifications.
* *False_positive_recrudescence*: this contains an analysis of false positive recrudescence rates for *PfRecur* vs the model of Plucinski et al (2015), based on permuted pairs of baseline and recurrent isolates in Dimbu et al (2024).


**References**

Dimbu PR, Labuda S, Ferreira CM, Caquece F, André K, Pembele G, Pode D, João MF, Pelenda VM, Nieto Andrade B, Horton B. Therapeutic response to four artemisinin-based combination therapies in Angola, 2021. Antimicrobial Agents and Chemotherapy. 2024 Apr 3;68(4):e01525-23. Available at: <https://journals.asm.org/doi/full/10.1128/aac.01525-23>

Plucinski MM, Morton L, Bushman M, Dimbu PR, Udhayakumar V. Robust algorithm for systematic classification of malaria late treatment failures as recrudescence or reinfection using microsatellite genotyping. Antimicrobial Agents and Chemotherapy. 2015 Oct;59(10):6096-10. Available at <https://journals.asm.org/doi/full/10.1128/aac.00072-15>

Plucinski M. *AngolaTES2021*. GitHub. 2024. <https://github.com/MateuszPlucinski/AngolaTES2021>
