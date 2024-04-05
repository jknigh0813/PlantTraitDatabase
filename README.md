# PlantTraitDatabase

* TRY_ValidateName_BuildNewick_N.R – Imports raw TRY database trait files for trait N, checks names against a World Flora Online (WFO), and build newick tree using V.Phylomaker2.
* XFT_ValidateName_BuildNewick_N.R – Imports raw XFT database trait files for trait N, checks names against a World Flora Online (WFO), and build newick tree using V.Phylomaker2.
* GlobalTrees_ValidateNames.R – Validates/corrects BCGI Global Tree Search database names against WFO
* PhyloSig_HypTest_TRY.R – Pagel's lambda test on cleaned TRY data.
* PhyloSig_HypTest_TRY_GenusLevel.R – Pagel's lambda test on cleaned TRY data clipped to genus level.
* Solve_PEMs_TRY.R - Solve PEMs from cleaned TRY data.
* Solve_PEMs_Global.R - Solve PEMs for BCGI database. Works in batches to avoid memory bottlenecks.
* Solve_PCoAs_TRY.R - Solve PCoAs from cleaned TRY data.
* RF_TraitPredict_TRY.R – RF model to predict traits from PEMs
* PhyloSig_HypTest_RF_Residuals.R – Pagel's lambda test on RF model residuals.
* RF_TraitPredict_Global.R – Imputes tree species traits in batches from PEMs.
* Graphlan_Export_TRY.R – formats TRY data for fancy phylo plots in Graphlan.
* Graphlan_Export_RF_Residuals.R – formats RF model residuals for fancy phylo plots in Graphlan.
* SEP_Parameter_Estimate.m – Estimates SEP distribution parameters.
* SEPrnd.m – generates SEP random numbers.

References:

Knighton, J., Sanchez-Martinez, P., Anderegg, L. (In Prep). A Globally Comprehensive Database of Tree Hydraulic and Structural Traits Imputed from Phylogenetic Relationships.

Kattge, J., Bönisch, G., Díaz, S., Lavorel, S., Prentice, I. C., Leadley, P., ... & Cuntz, M. (2020). TRY plant trait database–enhanced coverage and open access. Global change biology, 26(1), 119-188.
