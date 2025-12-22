import sys
import re

file_path = 'R/because.R'
with open(file_path, 'r') as f:
    lines = f.readlines()

# I know lines 802 to 814 (1-indexed) are the mess.
# Index 0 is line 1.
# So index 801 to 814.

# Let's verify the content of line 802 to 814
start_idx = 801
end_idx = 814

# Replace with clean code
clean_code = [
    '    if (length(random_data_updates) > 0) {\n',
    '      for (nm in names(random_data_updates)) {\n',
    '        if (!is.null(nm) && nm != "") {\n',
    '          data[[nm]] <- random_data_updates[[nm]]\n',
    '        }\n',
    '      }\n',
    '    }\n'
]

# Check if lines match roughly what I expect to avoid destroying wrong code
# Line 807 was the orphan else
if 'else' in lines[806] or 'else' in lines[807]: # roughly
    # Apply replacement
    new_lines = lines[:start_idx] + clean_code + lines[end_idx:]
    
    with open(file_path, 'w') as f:
        f.writelines(new_lines)
    print("Loop repaired successfully.")
else:
    print("Safety check failed: Did not find 'else' in declared range. Printing lines:")
    for i in range(start_idx, end_idx):
        print(f"{i+1}: {lines[i].strip()}")
    sys.exit(1)
