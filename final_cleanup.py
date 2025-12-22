import sys
import re

file_path = 'R/because.R'
with open(file_path, 'r') as f:
    content = f.read()

# Define the exact block to remove (based on what I saw in view_file)
# I'll likely use string replacement or regex with DOTALL

block_start = '        # Check if SiteID is in data'
# The block ends with the closing brace of the ELSE block, THEN the closing brace of the orphan IF
# }
# }

marker_end = '        }\n      }'

# Let's find start
idx_start = content.find(block_start)
if idx_start == -1:
    print("Could not find block start")
    sys.exit(1)

# Find the end brace sequence after start
# It looks like:
# ...
#         }
#       }
#       fit <- because(

idx_fit = content.find('fit <- because(', idx_start)

# The content to remove is everything from idx_start to just before idx_fit (minus indentation/newline maybe)
# Actually, let's just splice it out.

new_content = content[:idx_start] + content[idx_fit:]

with open(file_path, 'w') as f:
    f.write(new_content)

print("Removed corrupted debug block successfully.")
