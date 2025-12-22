import sys
import re

file_path = 'R/because.R'
with open(file_path, 'r') as f:
    lines = f.readlines()

# Locate the row_ids block
# row_ids <- NULL
# if (!is.null(id_col)) { ... } else { ... }

# I want to wrap the IF/ELSE part in if (is.data.frame(data)) { ... }
# row_ids <- NULL is fine to stay outside.

start_marker = '    # Handle id_col for matching'
# lines around 731.

found_idx = -1
for i, line in enumerate(lines):
    if start_marker in line:
        found_idx = i
        break

if found_idx != -1:
    # row_ids <- NULL is usually the next line or so
    
    # I want to insert if (is.data.frame(data)) { after row_ids <- NULL?
    # Or just wrap the whole logic block?
    
    # line 732: row_ids <- NULL
    # line 733: if (!is.null(id_col)) {
    
    # We want:
    # row_ids <- NULL
    # if (is.data.frame(data)) {
    #   if (!is.null(id_col)) { ... } else { ... rn ... }
    # }
    
    # So insert `if (is.data.frame(data)) {` before line 733.
    # And `}` after line 746 (closing brace of else).
    
    # Find line 733
    idx_if = -1
    for j in range(found_idx, len(lines)):
        if 'if (!is.null(id_col)) {' in lines[j]:
            idx_if = j
            break
            
    if idx_if != -1:
        lines.insert(idx_if, '    if (is.data.frame(data)) {\n')
        
        # Now find the end.
        # It's an if/else block.
        # if (...) { ... } else { ... }
        
        # Scan for closing brace of else block.
        # It ends before "# Also check for variability-related columns"
        
        idx_end = -1
        for k in range(idx_if, len(lines)):
            if '# Also check for variability-related columns' in lines[k]:
                idx_end = k
                break
        
        if idx_end != -1:
            # Check if previous line is closing brace?
            # It usually is.
            # Insert } before comments
            lines.insert(idx_end, '    }\n')
            
            with open(file_path, 'w') as f:
                f.writelines(lines)
            print("Row ID logic wrapped.")
        else:
             print("Could not find end of row_ids block")
    else:
        print("Could not find start of row_ids if-block")
else:
    print("Could not find row_ids marker")
