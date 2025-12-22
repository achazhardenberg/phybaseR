import sys
import re

file_path = 'R/because.R'
with open(file_path, 'r') as f:
    content = f.read()

# Fix the merged line: `occ_vars <- ... if (length`
# Regex to find: assignment followed by spaces and if
pattern_err = r'(occ_vars <- names\(family\)\[family == "occupancy"\])\s*(if \(length\(occ_vars\) > 0\) \{)'

# Replace with newline
replacement = r'\1\n        \2'

new_content = re.sub(pattern_err, replacement, content)

with open(file_path, 'w') as f:
    f.write(new_content)

print("Syntax error repaired")
