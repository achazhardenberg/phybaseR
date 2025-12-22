import sys
import re

file_path = 'R/because.R'
with open(file_path, 'r') as f:
    lines = f.readlines()

# Clean up previous traces first? 
# They confuse the line numbers if I don't.
# I'll enable "overwrite" effectively by just finding the lines again.

# Target block: lines 708-716 of original.
# Search for 'available_vars <- intersect'

found_idx = -1
for i, line in enumerate(lines):
    if 'available_vars <- intersect' in line:
        found_idx = i
        break

if found_idx != -1:
    # Insert debugs after found_idx (Trace 2 is already before it)
    
    # After line 708:
    # Check if we already have Trace 2 before it.
    
    # I want to verify step by step
    debug_code = [
       '    if (!quiet) message("DEBUG: Detail A - intersect started")\n',
       '    if (!quiet) message("DEBUG: eq_vars keys: ", paste(head(eq_vars), collapse=","))\n',
       '    if (!quiet) message("DEBUG: data names: ", paste(head(names(data)), collapse=","))\n',
    ]
    # Insert BEFORE 708
    lines[found_idx:found_idx] = debug_code
    
    # Find setdiff line (missing_vars)
    # Re-search because indices shifted
    for j in range(len(lines)):
        if 'missing_vars <- setdiff' in lines[j] and 'latent' not in lines[j]:
             lines.insert(j, '    if (!quiet) message("DEBUG: Detail B - setdiff missing_vars")\n')
             # After setdiff? No, loop updates indices.
             # I'll just do it safely.
             break

    for k in range(len(lines)):
        if 'if (!is.null(latent))' in lines[k]:
             lines.insert(k, '    if (!quiet) message("DEBUG: Detail C - latent check")\n')
             break

    for l in range(len(lines)):
         if 'missing_vars <- setdiff' in lines[l] and 'latent' in lines[l]:
             lines.insert(l, '    if (!quiet) message("DEBUG: Detail D - setdiff latent")\n')
             break
             
with open(file_path, 'w') as f:
    f.writelines(lines)

print("Detailed traces added.")
