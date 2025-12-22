# Check Storks Estimates for Raw vs Scaled Data

devtools::load_all(".")
data(storks)

# 1. Scaled Estimates (Reference)
print("--- SCALED DATA (Reference) ---")
storks_std <- storks
storks_std[2:5] <- scale(storks[2:5])

fit_std <- because(
    equations = list(
        Storks ~ Area,
        Birth ~ Area,
        Humans ~ Birth
    ),
    data = storks_std,
    n.iter = 1000,
    quiet = TRUE
)
print("Standardized Estimates:")
print(summary(fit_std))


# 2. Raw Estimates (The "Bug" Case)
print("--- RAW DATA (The User Case) ---")
fit_raw <- because(
    equations = list(
        Storks ~ Area,
        Birth ~ Area,
        Humans ~ Birth
    ),
    data = storks,
    n.iter = 1000,
    quiet = TRUE
)
# We expect weird large values for raw data
print("Raw Estimates:")
print(summary(fit_raw))

# Verify if beta_Birth_Area ~ 0.002 is correct for raw data
# LM check
lm_raw <- lm(Birth ~ Area, data = storks)
print("LM Estimate for Birth ~ Area:")
print(coef(lm_raw))
