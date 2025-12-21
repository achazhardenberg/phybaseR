# Tech Note: Bayesian Model Averaging (BMA) in JAGS

## The Goal: Model Probabilities without RJMCMC

You want to know: **"What is the probability that** $X$ causes $Y$?" ($P(M_1 | Data)$).

Traditionally, this requires **Reversible Jump MCMC (RJMCMC)**, where the sampler jumps between models of different dimensions (e.g., Model A has 2 params, Model B has 3). JAGS cannot do this because it requires a fixed number of nodes.

## The Solution: Gibbs Variable Selection (Indicator Variables)

Instead of deleting a parameter, we include it but multiply it by a binary switch $z$.

### 1. The Math

For a predictor $X$ with coefficient $\beta$, we define the **Effective Coefficient** $\beta^*$ as:

$$ \beta^* = z \times \beta $$

Where: \* $\beta \sim \mathcal{N}(0, \sigma^2)$ (The "slab": the effect size *if* it exists) \* $z \sim \text{Bernoulli}(\pi)$ (The "spike": the probability of inclusion, usually $\pi=0.5$)

### 2. The Mechanics (JAGS Code)

The JAGS model looks like this:

``` jags
model {
  # Priors
  beta_raw ~ dnorm(0, 0.001)   # The size of the effect
  z        ~ dbern(0.5)        # The existence of the effect (Coin flip)

  # Effective Coefficient
  beta_final <- z * beta_raw

  # Likelihood
  for (i in 1:N) {
    mu[i] <- alpha + beta_final * X[i]
    Y[i] ~ dnorm(mu[i], tau)
  }
}
```

### 3. How it Works (The "Switch")

-   **Iteration 1**: $z=1$. The model is $Y = \alpha + \beta X$. The sampler estimates $\beta$ based on the data.
-   **Iteration 2**: $z=0$. The model is $Y = \alpha + 0$. The $\beta$ parameter is "detached" from the likelihood (it floats freely or samples from the prior), and the model effectively becomes the Null Model.

The MCMC sampler will "visit" the $z=1$ state more often if including $X$ improves the model fit (Likelihood) enough to justify the complexity cost (Prior).

### 4. Interpretation (The "Gold Standard")

After running 10,000 iterations, you calculate the mean of $z$:

$$ P(X \to Y) = \frac{1}{N} \sum z_i $$

-   **Result = 0.85**: "There is an 85% probability that X causes Y." (Posterior Inclusion Probability).
-   **Result = 0.05**: "There is only a 5% chance this link exists."

This is **structurally identical** to obtaining Bayes Factors regarding the inclusion of the variable, but purely within a single JAGS run.

## Why this is better than "Significance"

-   **P-value**: "Assuming the effect is zero, how weird is this data?" (Backward logic).
-   **BMA (**$z$): "What is the probability the effect is non-zero?" (Forward logic).

## 5. Magnitude vs. Probability (The "Two Dimensions")

You currently have **Standardized Path Coefficients**. These tell you: *"If this arrow exists, how strong is it?"*

With BMA, every edge gets **two** numbers:

1.  **Strength (**$\beta$): The standardized path coefficient (Magnitude).
2.  **Certainty (PIP)**: The probability the arrow exists at all (Existence).

### Visualizing the Result

In a BMA-weighted DAG: \* **Edge Width** = Strength ($\beta$). \* **Edge Transparency** = Certainty (PIP).

### The "Shrinkage" Effect

If you averaged the coefficients across *all* iterations (including the zeros where $z=0$), you get the **BMA Estimate**:

$$ \beta_{BMA} \approx \beta_{strength} \times P(Inclusion) $$

-   A strong effect ($\beta=0.8$) that is uncertain ($P=0.5$) becomes a moderate estimate ($0.4$).
-   A moderate effect ($\beta=0.4$) that is certain ($P=1.0$) stays moderate ($0.4$).

## 6. Computational Cost

**Yes, it is more expensive.**

While the math per iteration is cheap (just one extra multiplication), the **convergence** is slower.

1.  **The Mixing Problem**: The sampler has to "jump" between two distinct peaks ($z=0$ vs $z=1$). If the data strongly supports one model, it gets "stuck" there and rarely visits the other, making it hard to estimate the exact probability (e.g., distinguishing 99% from 99.9%).
2.  **Iteration Count**: A standard model might converge in 10,000 iterations. A BMA model might need **100,000 to 500,000 iterations** to properly explore all the combinations of On/Off switches.

**The Trade-off**: \* **Fixed Model**: fast, but tells you nothing about structural uncertainty. \* **BMA**: 10x slower, but answers "Is this model even right?".

## 7. BMA vs. Structure Learning (The "Scout vs. Judge")

You asked: *"How is this different from Structure Learning?"*

They are two tools for the same problem (finding the DAG), but they work at different scales.

### 1. Structure Learning (The Scout)

-   **Goal**: Find a plausible DAG from scratch (Nothing $\to$ Something).
-   **Method**: Algorithms like PC (Constraint-based) or Greedy Search (Score-based).
-   **Scale**: Can handle 100 variables and billions of possible edges.
-   **Output**: One "Best Guess" DAG (or a few top candidates).
-   **Weakness**: Often ignores parameter uncertainty, assumes linearity, can be biased.

### 2. BMA (The Judge)

-   **Goal**: Rigorously test a specific set of hypotheses (Something $\to$ Truth).
-   **Method**: Full Bayesian MCMC with Indicator Variables.
-   **Scale**: Can handle \~5-10 uncertain arrows. You can't run BMA on 100 variables.
-   **Output**: Exact Probabilities ($P=0.83$) and Parameter Estimates averaged over uncertainty.
-   **Strength**: Handles all the complexity (Phylogeny, Spatial, Non-linear) correctly.

### The "Because" Workflow

1.  **Run `because_discover` (Structure Learning)**: Use the fast "Scout" to scan the data and suggest a candidate DAG.
2.  **Refine with BMA**: Identify the 3-4 edges the Scout was unsure about (or that you theoretically doubt).

## 8. The Ultimate D-separation Test (Direct Probability of Independence)

You asked: *"If we used BMA on a d-separation test, would we get the probability they are independent?"*

**YES.** And this is a massive improvement over standard tests.

### Standard Test (Frequentist-ish)

-   We add the arrow $X \to Y$.
-   We check if the 95% CI includes zero.
-   **Result**: "We failed to reject the null hypothesis." (Ambiguous: Is it truly independent, or just weak data?)

### BMA Test (Bayesian)

-   We add the arrow $X \to Y$ with an indicator $z$.
-   We calculate the PIP (Probability of Inclusion).
    -   $P(Dependent) = \text{Mean}(z)$
    -   $P(Independent) = 1 - \text{Mean}(z)$

### Example

If `z_RS_BM` has a mean of **0.12**: \* There is a 12% chance the variables are connected. \* There is an **88% probability they are conditionally independent**.

This gives you a direct, interpretable number quantifying **support for the Null Hypothesis**, which is impossible in frequentist statistics.
