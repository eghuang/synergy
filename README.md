# synergy

[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)  

Synergy is a R CRAN package prototype for synergy theory methods. Synergy theory is concerned with estimating a baseline dose-effect relationship (DER) aka [dose-response relationships](https://en.wikipedia.org/wiki/Dose%E2%80%93response_relationship) without synergistic or antagonistic effects for a chosen combination of DERs and doses. This package implements incremental effect addivity (IEA), a recently developed synergy theory that offers advantages over those previously used. 

## Installation
Run `devtools::install_github("eghuang/synergy")`.

## Repository contents

#### `DESCRIPTION`
Attribution, licensing, contact information, dependencies, and other information. 

#### `R`
R functions.

#### `man` 
Function documentation.

## Attribution

    # Copyright:     (C) 2017-2019 Sachs Undergraduate Research Apprentice Program
    #                (URAP) class at University of California, Berkeley.
    #                This program and its accompanying materials are distributed
    #                under the terms of the GNU General Public License v3. As
    #                detailed in that license no warranty, explicit or implied,
    #                comes with this suite of R scripts.
    # Package name:  synergy
    # Purpose:       Concerns radiogenic mouse Harderian gland (HG) tumorigenesis but 
    #                can be adapted for other uses. Loads ion and tumor prevalence data
    #                from CSV files and models tumor prevalence using the 
    #                Incremental Effect Additivity model. Additional functions are 
    #                provided to support custom dose-effect relationships and create 
    #                confidence intervals for prevalence estimates.
    # Contact:       Edward Greg Huang
    # Website:       https://github.com/eghuang/synergy
    # Mod history:   24 Dec 2019
    # Attribution:   This R script was developed at UC Berkeley. Authors and contributors
    #                include Edward Greg Huang, Dae Woong Ham, Yimin Lin, Mark Ebert,
    #                Yunzhi Zhang and Ray Sachs 2017-2019.

## Additional background
The algorithms in Synergy were developed and used for modeling cancer risk for astronauts on Mars missions or other extended voyages above low earth orbit. In particular, IEA was introduced in [Siranart et al. 2016](https://www.rrjournal.org/doi/full/10.1667/RR14411.1) and [Ham et al., (2018)](https://www.rrjournal.org/doi/full/10.1667/RR14948.1). More details on synergy theory and the application of IEA in radiation research can be found in [Huang et al., (2019)](https://link.springer.com/article/10.1007%2Fs00411-018-00774-x). Other results from applying IEA synergy theory to experimental results in NASA radiobiology is expected to be published in 2020.

## Relevant references and abbreviations:

    #   ".93Alp" = Alpen et al. "Tumorigenic potential of high-Z, high-LET charged-
    #                           particle radiations." Rad Res 136:382-391 (1993).
    #
    #   ".94Alp" = Alpen et al. "Fluence-based relative biological effectiveness for
    #                           charged particle carcinogenesis in mouse Harderian
    #                           gland." Adv Space Res 14(10): 573-581. (1994).
    #
    #   "16Chang" = Chang et al. "Harderian Gland Tumorigenesis: Low-Dose and LET
    #                            Response." Radiat Res 185(5): 449-460. (2016).
    #
    #   "16Srn" = Siranart et al. "Mixed Beam Murine Harderian Gland Tumorigenesis:
    #                             Predicted Dose-Effect Relationships if neither
    #                             Synergism nor Antagonism Occurs."
    #                             Radiat Res 186(6): 577-591 (2016).
    #
    #   "17Cuc" = Cucinotta & Cacao. "Non-Targeted Effects Models Predict
    #                                Significantly Higher Mars Mission Cancer Risk
    #                                than Targeted Effects Models."
    #                                Sci Rep 7(1): 1832. (2017). PMC5431989
    #
    #   "DER"     = Dose-effect relation(ship)"
    #   "HZE"     = High atomic number Z and high energy
    #   "HG"      = Harderian Gland
    #   "IEA"     = Incremental Effect Additivity
    #   "LET"     = Linear energy transfer; stopping power
    #   "NSRL"    = NASA Space Radiation Laboratory in Brookhaven NY
    #   "NTE"     = Non-targeted effects
    #   "SEA"     = Simple Effect Additivity
    #   "TE"      = Targeted effects
    #   "cGy"     = Centigray
    
