import sys
import re

file_path = 'R/because.R'
with open(file_path, 'r') as f:
    content = f.read()

# Locate `fit <- because(` inside `run_single_dsep_test`
marker_fit = 'fit <- because('
idx_fit = content.find(marker_fit)

if idx_fit == -1:
    print("Error: Could not find child because call")
    sys.exit(1)

# Debug block
debug_block = r'''
      if (!quiet) {
        message("DEBUG: Child because call.")
        message("DEBUG: test_eq: ", deparse(test_eq))
        message("DEBUG: Random arg present? ", !is.null(random))
        if (!is.null(random)) message("DEBUG: Random: ", deparse(random))
        
        # Check if SiteID is in data
        # Note: test_data is the dataframe passed
        if ("SiteID" %in% names(test_data)) {
           message("DEBUG: SiteID FOUND in test_data.")
        } else {
           message("DEBUG: SiteID NOT FOUND in test_data. Names: ", paste(names(test_data), collapse=", "))
        }
      }
'''

# Insert before `fit <- because`
new_content = content[:idx_fit] + debug_block + content[idx_fit:]

with open(file_path, 'w') as f:
    f.write(new_content)

print("Debug traces added to run_single_dsep_test")
