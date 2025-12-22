import sys

file_path = 'R/because.R'
with open(file_path, 'r') as f:
    content = f.read()

# Locate beginning of function body after inputs
marker = '# --- Input Validation & Setup ---'
idx = content.find(marker)
if idx == -1:
    print("Could not find marker")
    sys.exit(1)

insert_idx = content.find('\n', idx) + 1

debug_code = r'''
  if (!quiet) {
     # message("DEBUG: because() called.")
     # message("DEBUG: equations length: ", length(equations))
     if (!is.null(random)) {
        # message("DEBUG: random arg: ", deparse(random))
     }
  }
'''

# Better location: After random_terms are parsed.
# Line 398 in view 916.
marker_rt = 'random_terms <- random_terms[!duplicated(keys)]'
idx_rt = content.find(marker_rt)
if idx_rt == -1:
    # Try finding the block end
    marker_rt = '# --- Polynomial Term Extraction ---'
    idx_rt = content.find(marker_rt)

insert_idx_2 = idx_rt

debug_code_2 = r'''
  if (!quiet) {
     message("DEBUG: Parsed random_terms count: ", length(random_terms))
     if (length(random_terms) > 0) {
        msg <- sapply(random_terms, function(x) paste(x$response, x$group, sep="|"))
        message("DEBUG: random_terms: ", paste(msg, collapse=", "))
     } else {
        message("DEBUG: random_terms is EMPTY.")
        if(!is.null(random)) message("DEBUG: global random was: ", deparse(random))
     }
  }
'''

new_content = content[:insert_idx_2] + debug_code_2 + content[insert_idx_2:]

with open(file_path, 'w') as f:
    f.write(new_content)

print("Debug traces added to because()")
