import sys

file_path = 'R/because.R'
with open(file_path, 'r') as f:
    content = f.read()

# Insert at top of function body
marker = '# --- Input Validation & Setup ---'
idx = content.find(marker)
insert_idx = content.find('\n', idx) + 1

debug_code = r'''
  if (!quiet) {
     message("DEBUG: class(data): ", class(data)[1])
     is_list_debug <- is.list(data) && !is.data.frame(data)
     message("DEBUG: is_list_data logic check: ", is_list_debug)
  }
'''

new_content = content[:insert_idx] + debug_code + content[insert_idx:]

with open(file_path, 'w') as f:
    f.write(new_content)

print("Data type debug added")
