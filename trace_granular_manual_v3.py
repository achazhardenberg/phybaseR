import sys

file_path = 'R/because.R'
with open(file_path, 'r') as f:
    lines = f.readlines()

# Locate robust names lines to find the block
# `setdiff(eq_vars, (if ...`

idx = -1
for i, line in enumerate(lines):
    if 'setdiff(eq_vars' in line:
        idx = i
        break

if idx != -1:
    # Adding traces after setdiff
    idx += 1
    lines.insert(idx, '    if (!quiet) message("DEBUG: Trace A - After primary setdiff")\n'); idx += 1
    
    # Locate latent check
    # search forward
    for j in range(idx, idx+20):
        if 'if (!is.null(latent))' in lines[j]:
             lines.insert(j, '    if (!quiet) message("DEBUG: Trace B - Before latent check")\n'); idx = j + 2
             lines.insert(j+2, '    if (!quiet) message("DEBUG: Trace C - Inside latent block")\n'); idx = j + 3
             break
             
    # Locate missing vars check
    for k in range(idx, idx+20):
        if 'if (length(missing_vars) > 0 && length(available_vars)' in lines[k]:
             lines.insert(k, '    if (!quiet) message("DEBUG: Trace D - Before fatal missing check")\n'); idx = k+2
             break
             
    # Locate missing warning check
    for l in range(idx, idx+20):
        if 'if (length(missing_vars) > 0 && !quiet)' in lines[l]:
             lines.insert(l, '    if (!quiet) message("DEBUG: Trace E - Before warning missing check")\n'); idx = l+2
             break
             
    # Locate row_ids check (is.data.frame wrapper)
    for m in range(idx, idx+30):
        if 'if (is.data.frame(data)) {' in lines[m]:
             lines.insert(m, '    if (!quiet) message("DEBUG: Trace F - Before row_ids wrapper")\n'); idx = m+2
             break

    with open(file_path, 'w') as f:
        f.writelines(lines)
    print("V3 traces added.")
else:
    print("Could not find robust names line.")
