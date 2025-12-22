import sys
import re

file_path = 'R/because.R'
with open(file_path, 'r') as f:
    lines = f.readlines()

# Target block: lines 700-750.
# I want to insert message("DEBUG: Line X") before every line that isn't empty/comment

start_marker = 'if (is.data.frame(data) || (is.list(data)'
found = False
count = 0

new_lines = []
for line in lines:
    if start_marker in line:
        found = True
    
    if found and count < 60: # Limit scope
        if line.strip() and not line.strip().startswith('#') and 'if (!quiet)' not in line:
             line_stripped = line.strip().replace('"', "'")
             new_lines.append(f'    if (!quiet) message("DEBUG: Executing: {line_stripped}")\n')
        count += 1
    
    new_lines.append(line)

with open(file_path, 'w') as f:
    f.writelines(new_lines)

print("Granular traces added.")
