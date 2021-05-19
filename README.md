# 19/05/2021

# Benedetti-et-al. (#NCOMMS-20-37764A)
This repository contains parts of the R scripts developed during my postdoc at ETH Zürich, D-USYS, IBP, UP group.

The present R scripts were developed to produce the results of the manuscript submitted to Nature Communications (#NCOMMS-20-37764A).

- Scripts labelled as 'RSCRIPTBATCH' are those that were run on a local cluster to perform those functions defined within the script in parallel using the 'parallel' R package (R Core Team, 2020) within a mclapply().
- Scripts that were given a number correspond to those R scripts were I formated and mined the main datasets involved. They usually contain organized sequences of code where data are being read, examined, reformatted, analyzed and plotted. These scripts often correspond to the preparation. The number given to a script mainly corresponds to a temporal marker, and not necessarily a direct sequence from one numbered script to another. Considering that multiple studies fall under the OVERSEE project, some of the numbered scripts are deposited in another repository. This is why there might be gaps between the numbers of the present R scripts.
- "OVERSEE" just corresponds to name I gave to the project overarching my main research activities at ETH Zürich (modelling and forecast of global plankton biodiversity). Therefore, other repositories corresponding to other studies/manuscripts also fall within the 'OVERSEE' umbrella.
- The end of the scripts' names usually indicate their main purposes. A more detailed list of the scripts' goals and content is usually given in the beginning of the numbered scripts.
