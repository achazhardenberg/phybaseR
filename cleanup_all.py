import sys
import re

file_path = 'R/because.R'
with open(file_path, 'r') as f:
    lines = f.readlines()

new_lines = []
for line in lines:
    # Remove lines containing our debug tokens
    if 'DEBUG: ' in line or 'if (!quiet) message("DEBUG:' in line:
        continue
    # Remove clean versions of traces if I removed the if(!quiet) part but left message
    if 'message("DEBUG:' in line:
        continue
    new_lines.append(line)

# Also ensure no double blank lines or messed up braces from manual deletions?
# Just preserving is safer.

with open(file_path, 'w') as f:
    f.writelines(new_lines)

print("All debug traces cleaned.")
