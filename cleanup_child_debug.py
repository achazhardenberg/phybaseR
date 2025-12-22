import sys
import re

file_path = 'R/because.R'
with open(file_path, 'r') as f:
    content = f.read()

# Remove the debug block I added:
# if (!quiet) {
#    message("DEBUG: Child because call.")
#    ...
# }

# Regex to match the block starting with explicit string and ending with }
pattern = r'\s*if \(!quiet\) \{\s*message\("DEBUG: Child because call\."\).*?\n\s*\}'
new_content = re.sub(pattern, '', content, flags=re.DOTALL)

with open(file_path, 'w') as f:
    f.write(new_content)

print("Child debug traces removed")
