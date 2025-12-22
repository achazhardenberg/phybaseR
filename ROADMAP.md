# Future Development Roadmap

This document outlines planned features and extensions for the `because` package.

## 1. Custom Priors

-   **Request**: Users need the ability to override default informative/weakly-informative priors, **crucial for mechanistic models**.
-   **Implementation**: Add a `priors` argument to `because()`.
-   **Mechanism**:
    -   Accept a named list: `list(alpha_Birth = "dnorm(0, 0.001)", sigma_Storks = "dunif(0, 10)")`.
    -   During JAGS model generation, replace the default prior string for that specific parameter with the user-supplied string.
-   **Validation**: Check that the provided string is valid JAGS syntax (or catch errors from JAGS).

## 2. Occupancy Models (Binomial Measurement Error)

**Concept**: Occupancy models are structurally identical to the Measurement Error (Latent Variable) framework already implemented in `because`. \* **State Process**: True Presence $z_i \sim \text{Bernoulli}(\psi_i)$ (Latent Variable) \* **Observation Process**: Detection History $y_{ij} \sim \text{Bernoulli}(z_i \cdot p_{ij})$ (Measurement Error)

**Current Limitation**: \* The `variability` argument currently hardcodes a **Gaussian** observation process: $X_{obs} \sim N(X_{true}, \sigma)$.

**Proposed Implementation**: \* Extend `variability` to accept a distribution type. \* Example API: `variability = c(Sighting = list(type="binomial", trials="reps"))` \* Update JAGS model generation to support `dbin` or `dbern` detection layers.

**Output**: \* The model would return estimates for: \* **Occupancy (**$\psi$): The probability a site is occupied (Latent State parameter). \* **Detection (**$p$): The probability of detecting the species if present (Observation Error parameter). \* Both can be modeled as functions of covariates.

**Killer Feature: Causal Multi-Species Occupancy** \* Because $\psi$ (Occupancy) and $p$ (Detection) are typically treated as nuisance parameters in other software, `because` can treat them as **endogenous response variables**. \* **Causal Links on Detection**: Test behavioral hypotheses (e.g., *Predator Presence -\> Prey Detection Probability*). Prey hide when predators are around ($p \downarrow$). \* **Causal Links on Occupancy**: Test demographic hypotheses (e.g., *Predator Presence -\> Prey Occupancy*). Prey avoid predator-rich sites ($z \downarrow$). \* **Confounding Control**: Use latent variables (MAG approach) to control for unmeasured environmental confounders that affect both species, allowing robust D-separation tests on the latent occupancy states. \* This effectively creates a **Causal Occupancy Engine**, moving beyond simple correlation. \* **Backend Strategy**: Adopt the **Rota et al. (2016)** approach (Multivariate Bernoulli). Unlike Latent MVN (Pollock), this models interaction probabilities directly (Log-linear), allowing for covariate-dependent interactions (e.g., Competition varies with Temperature).

## 3. Spatial Modeling (CAR Models)

-   [ ] **Native CAR Model Support**: Implement `dcar_normal()` in JAGS models to support Conditional Autoregressive (CAR) spatial models natively. - Proper estimation of the spatial dependence parameter ($\rho$). - Better alignment with standard ecological spatial modelling practices (e.g. `spdep`).

## 4. Non-Gaussian Optimization

-   Extend the optimized linear algebra (matrix inversion) approach to Binomial and Poisson response models using PQL (Penalized Quasi-Likelihood) or similar approximations, or improved sampling schemes.

## 5. Phylogeography (Spatial-Phylogenetic Interactions)

-   **Concept**: Allow the strength of the phylogenetic signal to vary across space, or spatial autocorrelation to vary by clade.
-   **Implementation**: A tensor product of the Phylogenetic and Spatial precision matrices ($K_{phy} \otimes K_{spatial}$).
-   **Application**: Testing if evolutionary rates differ between regions (e.g., Tropics vs Temperate), effectively bridging Macroecology and Phylogeography.

## 6. Integrated Population Models (IPMs)

-   **Concept**: IPMs are the logical endpoint of this framework. They combine multiple observation datasets (Counts, Mark-Recapture, Productivity) to estimate shared latent population parameters.
-   **Architecture**:
    -   **Latent Process**: Population dynamics ($N_{t+1} = N_t \cdot S + R$).
    -   **Observation Models**: $C_t \sim \text{Poisson}(N_t)$, $y_{ij} \sim \text{Bernoulli}(S)$.
-   **Potential**: Since `because` already explicitly separates Latent States from Observations, extending it to support time-series dynamics ($t \to t+1$) would allow it to fit full IPMs.

**Killer Feature: Causal Demographic Hypothesis Testing** \* Because demographic rates (Survival, Recruitment) and states (Population Size $N_t$) are latent variables, they can be used in causal tests. \* Example: Test if **Density Dependence** is real ($N_t \to S_t$) or if **Climate** drives population via Recruitment ($Climate \to R_t \to N_{t+1}$). \* This moves IPMs from "estimation tools" to "hypothesis testing engines".

## 7. Automatic DAG Visualization (COMPLETE)

-   **Goal**: Replace text summaries with publication-ready path diagrams.
-   **Implementation**:
    -   S3 method `plot(fit)`.
    -   Uses `igraph` or `DiagrammeR` to draw the causal graph.
    -   Edges colored by significance (Red/Blue); Width scaled by coefficient magnitude.

## 8. Zero-Inflation Support (COMPLETE)

-   **Problem**: Ecological count data often has excess zeros (true absence vs. undetected).
-   **Solution**: Support `family = "zip"` (Zero-Inflated Poisson) and `zinb` (Zero-Inflated Negative Binomial).
-   **Mechanism**: A mixture model with a binary "structural zero" process and a count process.

## 9. Alternative Backends (NIMBLE & INLA)

-   **Goal**: Enable high-performance computing for massive datasets ($N > 10,000$) or complex mixtures.
-   **Strategy**: **Multi-Backend API**.
    -   Add a `backend` argument to `because(..., backend = "jags")`.
    -   **"jags"** (Default): Fast startup, no compiler needed. Best for prototyping and standard usage.
    -   **"nimble"**: Compiles to C++. User can customize samplers (e.g. HMC + Gibbs). **Crucial for ODEs and complex mechanistic models**.
    -   **"inla"**: Approximate Bayesian Inference.
        -   **Use Case**: Massive spatial/temporal datasets where MCMC is too slow.
        -   **Role**: The engine for `because_discover` (Structure Learning) and large-scale Spatiotemporal mapping.
        -   **Trade-off**: Extremely fast, but less flexible for very custom non-standard distributions.
-   **Implementation**: Refactor `because.R` to adopt a plugin architecture for model execution.
-   **Note**: Stan support is deferred indefinitely due to the need for a separate code generator.

## 10. Causal Biodiversity Modeling (Hill Numbers)

-   **Concept**: Instead of just estimating diversity, *test causes* of diversity patterns using Hill numbers ($^{q}D$) as response variables.
-   **Integration**:
    -   **Response**: Site-level diversity (Taxonomic, Phylogenetic, or Functional Hill numbers).
    -   **Predictors**: Environmental drivers (Temperature, Productivity).
    -   **Controls**: Spatial autocorrelation ($u_{spatial}$) and Biogeographic history.
-   **Causal Tests**: Test if e.g., *Productivity -\> Functional Diversity -\> Ecosystem Stability*.
-   **Innovation**: Most diversity studies are correlative. `because` allows formally testing causal pathways to biodiversity while accounting for spatial-phylogenetic non-independence.

## 11. Causal Utilities (Interventions & Robustness)

-   **`because_do(fit, ...)`**: Counterfactual Simulation.
    -   **Concept**: Implements Pearl's **do-operator** ($\text{do}(X=x)$) for structural interventions.
    -   **Mechanism**:
        1.  **Surgery**: Removes incoming edges to the intervened variable in the DAG.
        2.  **Propagation**: Propagates the fixed value $x$ to all downstream descendants (children) using the fitted structural equations.
        3.  **Bayesian**: Returns the full posterior distribution of the counterfactual outcome.
    -   **Application**: "What would the population size be if we had reduced poaching by 50%?" (as opposed to just comparing areas with low vs high poaching).
-   **`because_sensitivity(fit)`**: Unmeasured Confounding Analysis.
    -   **Goal**: Assess robustness. "How strong would an unmeasured confounder have to be to explain away this result?" (similar to E-values).
    -   **Implementation**: Post-hoc analysis of the residuals.

------------------------------------------------------------------------

# JAGS Optimization Implementation Progress

*(Archived from previous technical roadmap)*

## Completed (2025-12-03)

-   ✅ **Benchmarking**: Confirmed random effects formulation is \~4.6x - 15.5x faster.
-   ✅ **Code Generation**: Implemented additive random effects (`u_phylo + u_spatial`).
-   ✅ **Variance Partitioning**: Now estimates separate `sigma` components instead of aggregate lambda.
-   ✅ **Multiple Covariance**: Validated with `tests/test_multiple_covariance.R`.
-   ✅ **Non-Gaussian Families**: Validated support for `binomial`, `poisson`, `negbinomial`, `zip`, `zinb` with full multi-tree phylogenetic signal.

## 12. Deterministic Nodes (Formula interface) (COMPLETE)

-   **Status**: **Implemented**.
-   **Feature**: `R/deterministic_nodes.R` automatically parses formulas with interactions (`:`), logic, and functions (e.g., `I(x^2)`).
-   **Mechanism**: Converts them into internal deterministic nodes (`det_...[i] <- ...`) in JAGS, ensuring correct propagation of interventions.
-   **Benefit**: Allows "Mechanism-First" modeling (e.g. `growth ~ Vmax * Light / (Km + Light)`) without manual data pre-processing.

## 13. State Space Models (Time Series)

-   **Concept**: With the addition of **Deterministic Nodes** (Item 12) and **Latent Variables** (Item 1), `because` is structurally capable of fitting State Space Models (SSM).
-   **Missing Piece**: **Time Indexing** ($t, t-1$).
-   **Proposed Implementation**:
    -   Support autoregressive notation in formulas: `N ~ N[t-1] * S + R`.
    -   Internally map this to a JAGS loop: `N[i] <- N[i-1] * ...`.
    -   Combine with `variability` (measurement error) to handle the "Partially Observed" component.
    -   Combine with `variability` (measurement error) to handle the "Partially Observed" component.
    -   **Result**: A general-purpose Bayesian SSM engine that can *also* handle phylogeny and spatial autocorrelation (e.g., Spatiotemporal SSMs). "Process-based" causal inference.
    -   **Strategic Note**: This is the **primary mechanistic engine** for `because`, as most ecological data is discrete (year-to-year), making this far more widely useful than continuous ODEs.
    -   **Temporal D-Separation**:
        -   DAGs are "unrolled" over time ($t-1 \to t$).
        -   Independence tests check for **lagged effects** (e.g., does $X_{t-1}$ predict $Y_t$ given $X_t$?).
        -   This allows testing Granger Causality-like hypotheses within the Bayesian structural framework.

## 14. Structure Learning (Discovery Engine)

-   **Problem**: Users often have many variables and unknown structure.
-   **Solution**: **`because_discover(data)`**
-   **Mechanism**: A "Phylogenetic PC Algorithm".
    -   Starts with a fully connected graph.
    -   Iteratively removes edges based on conditional independence tests (D-separation).
    -   **Innovation**: Uses the `because` backend to perform these tests *while accounting for phylogeny and spatial structure*, avoiding the biases of standard structure learning algorithms.
-   **Hybrid Discovery-Confirmation Workflow**:
    1.  **Exploratory Phase**: Use a fast backend (INLA/Horseshoe JAGS) to "learn" the structure of environmental variables where theory is weak.
    2.  **Theoretical Phase**: Integrate the learned environmental DAG as a fixed subgraph into a more complex biological model (e.g., Occupancy) where user provides strong prior theoretical knowledge.
    3.  **Refinement**: Use d-separation tests on the combined model to identify remaining structural misfits.

## 15. Automated Mediation Analysis (COMPLETE)

-   **Goal**: Decompose effects into "Direct" and "Indirect" pathways.
-   **Solution**: **`because_mediation(fit, exposure, outcome)`**
-   **Mechanism**:
    -   Automatically traces all paths from $X$ to $Y$ in the fitted DAG using `igraph`.
    -   Calculates Indirect Effects (e.g., $X \to M \to Y$) by multiplying coefficients ($\beta_{XM} \cdot \beta_{MY}$) chain-by-chain.
    -   Calculates Total Effects (Sum of all paths).
    -   Answers: *"How much of the effect is mediated by M?"* with full Bayesian uncertainty intervals.

## 16. The Front-Door Criterion

-   **Problem**: Unmeasured confounding between $X$ and $Y$ where back-door adjustment is impossible.
-   **Solution**: Detect and exploit "Front-Door" structures ($X \to M \to Y$, with $M$ isolated from confounders).
-   **Mechanism**: A topological check effectively allowing causal identification even when standard controls are missing.

## 17. Transportability (External Validity)

-   **Problem**: Can a model learned in one region (e.g., Europe) predict in another (e.g., Amazon)?
-   **Solution**: **`because_transport(fit, new_context)`**
-   **Mechanism**: Uses **Pearl's Selection Diagrams**. Users flag which structural equations differ in the new environment (Selection Nodes $S$). The tool computes whether the causal effect can be validly transported or needs recalibration.

## 18. Data Fusion (Grand DAGs)

-   **Problem**: Disparate datasets (e.g., Dataset A has $\{X, Y\}$, Dataset B has $\{Y, Z\}$).
-   **Solution**: **`because_merge(list(data1, data2))`**
-   **Mechanism**: Because `because` treats missing data as latent parameters, it can combine datasets into a single "Grand DAG", imputing the missing variables for each dataset based on the shared structure.

## 19. Differential Equations (ODE Solver Interface)

-   **Concept**: The "Gold Standard" of mechanistic modeling in ecology (e.g., Lotka-Volterra, SIR models).
-   **Problem**: Currently, `because` assumes discrete-time updates or static relationships. Many ecological theories are defined in continuous time (`dN/dt`).
-   **Solution**: Investigate JAGS extensions (like passing C++ modules) or consider alternative backends (Stan/NIMBLE) for this specific feature, as native JAGS lacks a built-in `ode()` solver.
-   **Workaround**: Use discrete approximation (Euler method) within JAGS loops for now.
-   **Priority**: **Low/Niche**. As noted by user, most ecological data is discrete (e.g. annual), so Item 13 (SSMs) covers 90% of use cases. This remains a long-term goal for specific sub-fields (e.g., epidemiology).

## 20. Bayesian Model Averaging (BMA) & Variable Selection (Experimental)

-   **Status**: **Prototype Implemented** (See `drafts/demo_bma_rhino.R`).
-   **Concept**: Use **Gibbs Variable Selection** (Indicator Variables) to calculate Posterior Inclusion Probabilities (PIPs) for edges, effectively performing Bayesian D-separation tests.
-   **Mechanism**:
    -   User specifies a "Supermodel".
    -   `because` automatically injects $z \sim dbern(0.5)$ for optional paths.
    -   Returns PIPs for structural validation.
-   **Plan**: Move prototype logic into `because()` to allow `bma=TRUE`.

## 21. Tidyverse DSL (The "TidySEM" Wrapper)

-   **Concept**: An optional, pipe-based interface inspired by `TidySEM` but powered by `because`.
-   **Goal**: "Steal the Syntax, Keep the Backend". Provide the grammar users love without the frequentist limitations.
-   **Core Principle**: The traditional `because()` function remains the stable engine. This DSL is a lightweight wrapper.
-   **Proposed Syntax** (Matching TidySEM): `r     model <- tidy_because(data) |>       add_paths(Risk ~ Population + Climate) |>       set_family(Risk = "binomial") |>  # 'because' specific extension       estimate_jags()`
-   **Key Features**:
    -   **`add_paths()`**: Incrementally build equations.
    -   **`set_family()`**: Explicitly handle families (Bernoulli, Poisson, etc).
    -   **`set_variability()`**: Explicitly handle measurement error.
    -   **Base Pipe (`|>`)**: Uses the native R pipe (requires R \>= 4.1.0) to avoid `magrittr` dependency.
    -   **`broom` Compatibility**: `tidy(fit)` and `glance(fit)` for seamless plotting.
-   **Priority**: **Medium** (Excellent for complex ecological datasets with many species).

------------------------------------------------------------------------

# Appendix: Technical Implementation Strategies

## A. Strategy for Structure Learning with Missing Data

The "Purist" approach (treating missing data as parameters in structure learning) is computationally infeasible. We will adopt a **"Pragmatic Hybrid" (Scout & Surveyor)** approach:

1.  **Pre-Imputation (The Scout)**:
    -   Use fast tools (e.g., `Rphylopars`, `mice`) to generate $K$ (e.g., 5) completed datasets *before* starting discovery.
    -   This bypasses the missing-predictor limitation of fast engines like INLA.
2.  **Parallel Discovery**:
    -   Run `because_discover(backend="inla")` on each of the $K$ imputed datasets.
    -   This exploits INLA's speed (seconds vs hours) to search the massive graph space.
3.  **Consensus Pooling**:
    -   Identify edges that appear in the majority of the $K$ graphs (robustness check).
    -   This yields the **Candidate Graph**.
4.  **Full Estimation (The Surveyor)**:
    -   Pass the Candidate Graph to `because(backend="jags")`.
    -   JAGS treats the original missing data as parameters (discarding the pre-imputation) to properly estimate uncertainty and parameters in the final model.

## B. The INLA vs. JAGS Role

-   **INLA (The Gaussian Machine)**:
    -   **Role**: *Structure Learning, Spatial Mapping*.
    -   **Constraint**: Cannot handle non-linear logic (polynomials of latent vars) or complex missing data in predictors naturally.
    -   **Advantage**: Speed (1000x faster).
-   **NIMBLE (The Custom Workshop)**:
    -   **Role**: *ODEs, Complex Mechanistic Logic, Custom Samplers*.
    -   **Advantage**: Compiles to C++, allowing for differential equation solvers and arbitrary non-linear functions that are too slow or impossible in JAGS.
    -   **Strategic Value**: The bridge to full "Process-Based" modeling.
-   **JAGS (The General Engine)**:
    -   **Role**: *Final Estimation, `because_do` (Interventions), Complex Logic*.
    -   **Advantage**: Handles everything (interaction of latent vars, mechanism-based equations).
    -   **Constraint**: Slow.
