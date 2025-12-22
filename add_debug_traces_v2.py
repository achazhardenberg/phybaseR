import sys

file_path = 'R/because.R'
with open(file_path, 'r') as f:
    content = f.read()

# 1. Trace inside latent detection logic
# Look for: `if (!is.null(family)) {` (inside potential_latents block)
# Step 897, line 1482.
marker_latent = 'if (!is.null(family)) {'
# There might be multiple. We want the one inside the auto-detect block.
# Context: `potential_latents <- setdiff`
marker_context = 'potential_latents <- setdiff(vars_in_equations, vars_in_data)'
idx_context = content.find(marker_context)

if idx_context == -1:
    print("Error: Could not find potential_latents logic")
    # Try finding `potential_latents` assignment
    # View shows:
    # 1477:       potential_latents <- setdiff(vars_in_equations, vars_in_data)
    sys.exit(1)

# Find insertion point after `potential_latents` init
idx_insert_1 = content.find('\n', idx_context) + 1
trace_1 = r'''
      if(!quiet) {
         message("DEBUG: Initial potential_latents: ", paste(potential_latents, collapse=", "))
         message("DEBUG: vars_in_equations: ", paste(vars_in_equations, collapse=", "))
         message("DEBUG: vars_in_data: ", paste(vars_in_data, collapse=", "))
      }
'''

# Find the family check inside
check_marker = 'if (!is.null(family)) {'
idx_check = content.find(check_marker, idx_insert_1)
# Find insertion point inside that block
idx_insert_2 = content.find('\n', idx_check) + 1

trace_2 = r'''
        occ_vars <- names(family)[family == "occupancy"]
        if(!quiet) message("DEBUG: occ_vars for exclusion: ", paste(occ_vars, collapse=", "))
'''
# Note: occ_vars is computed inside originally, I should just trace after it.
# View 897 line 1483: `occ_vars <- ...`
# So find `occ_vars <- `
idx_occ = content.find('occ_vars <- names(family)', idx_check)
if idx_occ != -1:
    idx_insert_2 = content.find('\n', idx_occ) + 1

# 2. Trace before because_dsep call
marker_dsep = 'dsep_result <- because_dsep('
idx_dsep = content.find(marker_dsep)
if idx_dsep == -1:
    print("Error: Could not find because_dsep call")
    sys.exit(1)

trace_3 = r'''
    if(!quiet) {
      message("DEBUG: Calling because_dsep...")
      message("DEBUG: Equations[1]: ", deparse(equations[[1]]))
      message("DEBUG: Latent: ", paste(latent, collapse=", "))
    }
'''

# Apply traces
new_content = content[:idx_insert_1] + trace_1 + content[idx_insert_1:idx_dsep] + trace_3 + content[idx_dsep:]
# Oops, idx_dsep shifted.
# Better to build sequentially.

part1 = content[:idx_insert_1]
part2 = trace_1
part2 += content[idx_insert_1:idx_occ] # Wait, idx_occ is absolute index in original
# Let's slice carefully.

# Re-find in steps.
c1 = content
idx_1 = c1.find('potential_latents <- setdiff(vars_in_equations, vars_in_data)') + len('potential_latents <- setdiff(vars_in_equations, vars_in_data)')
trace_1 = r'''
      if(!quiet) {
         message("DEBUG: Initial potential_latents: ", paste(potential_latents, collapse=", "))
         message("DEBUG: vars_in_equations: ", paste(vars_in_equations, collapse=", "))
      }
'''

idx_occ = c1.find('occ_vars <- names(family)[family == "occupancy"]', idx_1)
if idx_occ == -1:
    print("Error finding occ_vars")
    sys.exit(1)
idx_occ_end = c1.find('\n', idx_occ) + 1
trace_2 = r'''
        if(!quiet) message("DEBUG: occ_vars found: ", paste(occ_vars, collapse=", "))
'''

idx_dsep = c1.find('dsep_result <- because_dsep(', idx_occ_end)
trace_3 = r'''
    if(!quiet) {
      message("DEBUG: Calling because_dsep")
      message("DEBUG: Equations passed: ", paste(sapply(equations, function(e) deparse(e)), collapse="; "))
      message("DEBUG: Latent passed: ", paste(latent, collapse=", "))
    }
'''

output = c1[:idx_1] + trace_1 + c1[idx_1:idx_occ_end] + trace_2 + c1[idx_occ_end:idx_dsep] + trace_3 + c1[idx_dsep:]

with open(file_path, 'w') as f:
    f.write(output)

print("Debug traces added")
