import sys
import re

file_path = 'R/because.R'
with open(file_path, 'r') as f:
    content = f.read()

# Replace `names(data)` with `safe_names(data)` where safe_names is defined or inline
# Inline: `(if (is.null(names(data))) character(0) else names(data))`

# Replacing globally is risky?
# Only inside the list block?
# The list block starts with `if (is.data.frame(data) || ...)`

# I'll replace `intersect(eq_vars, names(data))` 
# and `setdiff(eq_vars, names(data))`

content = content.replace('intersect(eq_vars, names(data))', 'intersect(eq_vars, (if (is.null(names(data))) character(0) else names(data)))')
content = content.replace('setdiff(eq_vars, names(data))', 'setdiff(eq_vars, (if (is.null(names(data))) character(0) else names(data)))')

with open(file_path, 'w') as f:
    f.write(content)

print("Robust names usage applied.")
