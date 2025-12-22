import sys
import re

file_path = 'R/because.R'
with open(file_path, 'r') as f:
    content = f.read()

# 1. content is messy because of my previous insertion inside the loop.
# I will locate the specific block I inserted previously and remove it.
# It started with `if (!quiet) { message("DEBUG: Final data keys ...`
# end ended with `}` inside the loop.

pattern_to_remove = r'if \(!quiet\) \{\s*message\("DEBUG: Final data keys.*?\n\s*\}'
content = re.sub(pattern_to_remove, '', content, flags=re.DOTALL)

# 2. Now insert new debug before the loop.
# Look for `if (length(random_data_updates) > 0)`

marker = 'if (length(random_data_updates) > 0)'
idx = content.find(marker)

if idx == -1:
    print("Could not find loop marker")
    sys.exit(1)

# Backtrack to newline
idx_insert = content.rfind('\n', 0, idx) + 1

debug_code = r'''
    if (!quiet) {
       message("DEBUG: Check random_data_updates. Length: ", length(random_data_updates))
       if (length(random_data_updates) > 0) {
          message("DEBUG: Keys in updates: ", paste(names(random_data_updates), collapse=", "))
       } else {
          message("DEBUG: random_data_updates is EMPTY.")
       }
    }
'''

new_content = content[:idx_insert] + debug_code + content[idx_insert:]

with open(file_path, 'w') as f:
    f.write(new_content)

print("Debug traces relocated")
