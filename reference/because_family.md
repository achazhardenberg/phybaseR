# Create a Custom Family (Distribution) for use with because

This helper function allows users to easily create custom distribution
families by providing the JAGS likelihood code. The function handles S3
method registration.

## Usage

``` r
because_family(
  name,
  jags_likelihood,
  link = "identity",
  extra_priors = NULL,
  description = NULL
)
```

## Arguments

- name:

  Character string; the name of your family (e.g., "student_t").

- jags_likelihood:

  A character string with the JAGS likelihood code. Use placeholders:
  `{response}` for variable name, `{mu}` for mean, `{i}` for loop index,
  `{tau}` for precision parameter.

- link:

  Link function name (default: "identity"). Options: "identity", "log",
  "logit".

- extra_priors:

  Character vector of additional JAGS prior statements. Use `{response}`
  placeholder for variable-specific naming.

- description:

  Optional description of the family.

## Value

A constructor function for the family that can be passed to because(...,
family = ...).

## Examples

``` r
if (FALSE) { # \dontrun{
# Create a Student-t family for robust regression
student_t <- because_family(
  name = "student_t",
  jags_likelihood = "{response}[{i}] ~ dt({mu}[{i}], {tau}, df_{response})",
  extra_priors = c("df_{response} ~ dunif(2, 100)")
)

# Use it
fit <- because(Y ~ X, data = data, family = student_t())
} # }
```
