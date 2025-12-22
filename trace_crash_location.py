import sys
import re

file_path = 'R/because.R'
with open(file_path, 'r') as f:
    lines = f.readlines()

# I want to insert debug messages inside the block starting at line 697.
# I'll search for specific markers inside that block.

markers = [
    ('eq_vars <- unique', 'message("DEBUG: Trace 1 - entering block")'),
    ('available_vars <- intersect', 'message("DEBUG: Trace 2 - vars calculated")'),
    ('if (length(missing_vars) > 0 && length(available_vars) == 0)', 'message("DEBUG: Trace 3 - missing check done")'),
    ('row_ids <- NULL', 'message("DEBUG: Trace 4 - row_ids init")'),
    ('# Also check for variability-related columns', 'message("DEBUG: Trace 5 - before se_cols")'),
    ('dummy_vars <- character(0)', 'message("DEBUG: Trace 6 - dummy vars init")'),
    ('data_list <- list()', 'message("DEBUG: Trace 7 - data_list init")'),
    ('data <- data_list', 'message("DEBUG: Trace 8 - data assigned")')
]

for marker_text, msg in markers:
    found = False
    for i, line in enumerate(lines):
        if marker_text in line:
            # Insert message before the line
            lines.insert(i, f'    if (!quiet) {msg}\n')
            found = True
            break # only first occurrence? careful.
    if not found:
        print(f"Warning: Marker '{marker_text}' not found")

with open(file_path, 'w') as f:
    f.writelines(lines)
    
print("Crash traces added.")
