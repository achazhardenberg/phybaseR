import sys

file_path = 'R/because.R'
with open(file_path, 'r') as f:
    lines = f.readlines()

# Find the block lines 700-720
# Search for `available_vars <- intersect`
idx = -1
for i, line in enumerate(lines):
    if 'available_vars <- intersect' in line:
        idx = i
        break

if idx != -1:
    # Insert before intersect
    lines.insert(idx, '    if (!quiet) message("DEBUG: Step 1 - before intersect")\n')
    idx += 2 # +1 for insert, +1 for intersect line
    
    # Insert after intersect
    lines.insert(idx, '    if (!quiet) message("DEBUG: Step 2 - after intersect / before setdiff")\n')
    idx += 2
    
    # Insert after setdiff
    lines.insert(idx, '    if (!quiet) message("DEBUG: Step 3 - after setdiff / before latent check")\n')
    # find latent check? assume next few lines
    
    # Search for if (!is.null(latent))
    found = False
    for j in range(idx, idx+10):
        if 'if (!is.null(latent))' in lines[j]:
             lines.insert(j, '    if (!quiet) message("DEBUG: Step 4 - latent check line next")\n')
             idx = j + 2
             found = True
             break
    
    if found:
        # After latent check?
        pass

    with open(file_path, 'w') as f:
        f.writelines(lines)
    print("Manual granular traces added.")
else:
    print("Block not found.")
