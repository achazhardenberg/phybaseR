import sys

file_path = 'R/because.R'
with open(file_path, 'r') as f:
    content = f.read()

# Locate the data list combination
# data <- data_list
marker = 'data <- data_list'
idx = content.find(marker)

if idx == -1:
    # Try finding the loop
    marker = 'if (length(random_data_updates) > 0) {'
    idx = content.find(marker)

if idx == -1:
    print("Could not find data merge block")
    sys.exit(1)

# I want to insert AFTER the loop closes?
# The loop:
#     if (length(random_data_updates) > 0) {
#       for (nm in names(random_data_updates)) {
#         data[[nm]] <- random_data_updates[[nm]]
#       }
#     }

# Let's verify what I'm targeting.
# I'll search for the loop specifically.

block = r'''if (length(random_data_updates) > 0) {
      for (nm in names(random_data_updates)) {
        data[[nm]] <- random_data_updates[[nm]]
      }
    }'''

idx_block = content.find(block.replace('\n', '').replace(' ', ''))
# Regex search for block because whitespace varies
import re
pattern = r'if \(length\(random_data_updates\) > 0\) \{.*?\}'
match = re.search(pattern, content, re.DOTALL)

if match:
    end_pos = match.end()
    debug_code = r'''
    if (!quiet) {
       message("DEBUG: Final data keys before JAGS: ", paste(names(data), collapse=", "))
       if ("N_SiteID" %in% names(data)) {
          message("DEBUG: N_SiteID value: ", data[["N_SiteID"]])
       } else {
          message("DEBUG: N_SiteID MISSING from data list!")
          message("DEBUG: random_data_updates keys: ", paste(names(random_data_updates), collapse=", "))
       }
    }
'''
    new_content = content[:end_pos] + debug_code + content[end_pos:]
    
    with open(file_path, 'w') as f:
        f.write(new_content)
    print("Deep debug traces added")
else:
    print("Could not find random_data_updates block via regex")
    sys.exit(1)
