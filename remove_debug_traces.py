import sys

file_path = 'R/because.R'
with open(file_path, 'r') as f:
    content = f.read()

# Remove DEBUG messages I added.
# They look like: `if(!quiet) { message("DEBUG: ...`

# I'll just use the previous backup if I had one? No.
# I'll read the file and strip the specific blocks I added.

# Block 1: Initial potential_latents debug
start_1 = '      if(!quiet) {\n         message("DEBUG: Initial potential_latents:'
end_1 = '      }\n'

# Block 2: occ_vars found
# It was a single line if statement in my python script?
# trace_2 = r'''        if(!quiet) message("DEBUG: occ_vars found: ", paste(occ_vars, collapse=", "))'''
marker_2 = 'message("DEBUG: occ_vars found:'

# Block 3: Calling because_dsep
start_3 = '    if(!quiet) {\n      message("DEBUG: Calling because_dsep")'
end_3 = '    }\n' # The next `}` after start_3

import re

# Remove Block 1
# Regex pattern: if(!quiet) {\n\s*message\("DEBUG: Initial.*?}\n
pattern1 = r'\s*if\(!quiet\) \{\s*message\("DEBUG: Initial.*?\n\s*\}'
content = re.sub(pattern1, '', content, flags=re.DOTALL)

# Remove Block 2 - it's a single line usually
# pattern2 = r'\s*if\(!quiet\) message\("DEBUG: occ_vars found:.*?\n'
pattern2 = r'\s*if\(!quiet\) message\("DEBUG: occ_vars found:.*?\n'
content = re.sub(pattern2, '', content)

# Remove Block 3
pattern3 = r'\s*if\(!quiet\) \{\s*message\("DEBUG: Calling because_dsep.*?\n\s*\}'
content = re.sub(pattern3, '', content, flags=re.DOTALL)

with open(file_path, 'w') as f:
    f.write(content)

print("Debug traces removed")
