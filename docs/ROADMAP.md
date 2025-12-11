# Future Development Roadmap

This document outlines planned features and extensions for the `because` package.

## 1. Occupancy Models (Binomial Measurement Error)

**Concept**: Occupancy models are structurally identical to the Measurement Error (Latent Variable) framework already implemented in `because`.
*   **State Process**: True Presence $z_i \sim \text{Bernoulli}(\psi_i)$ (Latent Variable)
*   **Observation Process**: Detection History $y_{ij} \sim \text{Bernoulli}(z_i \cdot p_{ij})$ (Measurement Error)

**Current Limitation**:
*   The `variability` argument currently hardcodes a **Gaussian** observation process: $X_{obs} \sim N(X_{true}, \sigma)$.

**Proposed Implementation**:
*   Extend `variability` to accept a distribution type.
*   Example API: `variability = c(Sighting = list(type="binomial", trials="reps"))`
*   Update JAGS model generation to support `dbin` or `dbern` detection layers.

**Output**:
*   The model would return estimates for:
    *   **Occupancy ($\psi$)**: The probability a site is occupied (Latent State parameter).
    *   **Detection ($p$)**: The probability of detecting the species if present (Observation Error parameter).
*   Both can be modeled as functions of covariates.

**Killer Feature: Causal Multi-Species Occupancy**
*   Because $\psi$ is a latent variable, it can be used as a **predictor** for other species.
*   Example: `Occupancy_SpeciesB ~ Occupancy_SpeciesA + Habitat`
*   This enables **Causal Mediation Analysis** on true occupancy states (e.g., *Habitat -> Species A -> Species B*), correcting for imperfect detection on both A and B.

## 2. Non-Gaussian Optimization
*   Extend the optimized linear algebra (matrix inversion) approach to Binomial and Poisson response models using PQL (Penalized Quasi-Likelihood) or similar approximations, or improved sampling schemes.


## 3. Phylogeography (Spatial-Phylogenetic Interactions)
*   **Concept**: Allow the strength of the phylogenetic signal to vary across space, or spatial autocorrelation to vary by clade.
*   **Implementation**: A tensor product of the Phylogenetic and Spatial precision matrices ($K_{phy} \otimes K_{spatial}$).
*   **Application**: Testing if evolutionary rates differ between regions (e.g., Tropics vs Temperate), effectively bridging Macroecology and Phylogeography.

## 4. Integrated Population Models (IPMs)
*   **Concept**: IPMs are the logical endpoint of this framework. They combine multiple observation datasets (Counts, Mark-Recapture, Productivity) to estimate shared latent population parameters.
*   **Architecture**:
    *   **Latent Process**: Population dynamics ($N_{t+1} = N_t \cdot S + R$).
    *   **Observation Models**: $C_t \sim \text{Poisson}(N_t)$, $y_{ij} \sim \text{Bernoulli}(S)$.
*   **Potential**: Since `because` already explicitly separates Latent States from Observations, extending it to support time-series dynamics ($t \to t+1$) would allow it to fit full IPMs.

**Killer Feature: Causal Demographic Hypothesis Testing**
*   Because demographic rates (Survival, Recruitment) and states (Population Size $N_t$) are latent variables, they can be used in causal tests.
*   Example: Test if **Density Dependence** is real ($N_t \to S_t$) or if **Climate** drives population via Recruitment ($Climate \to R_t \to N_{t+1}$).
*   This moves IPMs from "estimation tools" to "hypothesis testing engines".

## 5. Automatic DAG Visualization
*   **Goal**: Replace text summaries with publication-ready path diagrams.
*   **Implementation**: 
    *   S3 method `plot(fit)`.
    *   Uses `igraph` or `DiagrammeR` to draw the causal graph.
    *   Edges colored by significance (Red/Blue); Width scaled by coefficient magnitude.

## 6. Zero-Inflation Support
*   **Problem**: Ecological count data often has excess zeros (true absence vs. undetected).
*   **Solution**: Support `family = "zip"` (Zero-Inflated Poisson) and `zinb` (Zero-Inflated Negative Binomial).
*   **Mechanism**: A mixture model with a binary "structural zero" process and a count process.

## 7. Alternative Backends (HPC)
*   **Goal**: Handle massive datasets ($N > 10,000$).
*   **Plan**: functionality to export the generated DAG/Model to **Stan** or **Nimble** code for compilation and execution on high-performance clusters.

## 8. Causal Biodiversity Modeling (Hill Numbers)
*   **Concept**: Instead of just estimating diversity, *test causes* of diversity patterns using Hill numbers ($^{q}D$) as response variables.
*   **Integration**:
    *   **Response**: Site-level diversity (Taxonomic, Phylogenetic, or Functional Hill numbers).
    *   **Predictors**: Environmental drivers (Temperature, Productivity).
    *   **Controls**: Spatial autocorrelation ($u_{spatial}$) and Biogeographic history.
*   **Causal Tests**: Test if e.g., *Productivity -> Functional Diversity -> Ecosystem Stability*. 
*   **Innovation**: Most diversity studies are correlative. `because` allows formally testing causal pathways to biodiversity while accounting for spatial-phylogenetic non-independence.

## 9. Causal Utilities (New Functions)

These would be separate helper functions to strengthen the causal claims made by the models.

*   **`because_do(fit, intervention)`**: Counterfactual Simulation.
    *   **Goal**: Answer "What if?" questions.
    *   **Usage**: `because_do(fit, set(Temperature = +2))` 
    *   **Output**: The predicted distribution of the response under the intervention (cutting the incoming arrows to Temperature).

*   **`because_sensitivity(fit)`**: Unmeasured Confounding Analysis.
    *   **Goal**: Assess robustness. "How strong would an unmeasured confounder have to be to explain away this result?" (similar to E-values).
    *   **Implementation**: Post-hoc analysis of the residuals.

*   **`because_iv_strength(fit)`**: Instrumental Variable Diagnostics.
    *   **Goal**: If using an IV (e.g., Rainfall as an IV for Productivity), test its strength (F-statistic equivalent) and validity constraints where possible.
