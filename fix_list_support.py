import sys
import re

file_path = 'R/because.R'
with open(file_path, 'r') as f:
    content = f.read()

# 1. Modify the IF condition to support lists
# if (is.data.frame(data)) {
# to
# if (is.data.frame(data) || is_list_data) {

pattern_if = r'if \(is\.data\.frame\(data\)\) \{'
repl_if = 'if (is.data.frame(data) || is_list_data) {'

if pattern_if not in content and 'is_list_data' in content:
   # Maybe whitespace differs
   pattern_if = r'if\s*\(is\.data\.frame\(data\)\)\s*\{'

content = re.sub(pattern_if, repl_if, content)

# 2. Replace colnames(data) with names(data) within the function
# (Global replace is probably safe enough in this R file context, but let's be slightly careful)
# Actually, strict global replace of `colnames(data)` -> `names(data)` is safe for dataframes too.

content = content.replace('colnames(data)', 'names(data)')


# 3. Cleanup ALL debug traces I added.
# - add_main_debug.py traces
# - check_data_type.py traces
# - relocate_debug.py traces

# Regex for all if (!quiet) { message("DEBUG: ...") ... } blocks?
# I'll just remove specific strings to be safe.

debug_patterns = [
    r'\s*if \(!quiet\) \{\s*message\("DEBUG: class\(data\):".*?\n\s*\}', # check_data_type
    r'\s*if \(!quiet\) \{\s*# message\("DEBUG: because\(\) called\."\).*?\}', # add_main_debug (commented msg)
    r'\s*if \(!quiet\) \{\s*message\("DEBUG: Parsed random_terms count:".*?\n\s*\}', # add_main_debug part 2
    r'\s*if \(!quiet\) \{\s*message\("DEBUG: Check random_data_updates\..*?\n\s*\}' # relocate_debug
]

for pat in debug_patterns:
    content = re.sub(pat, '', content, flags=re.DOTALL)

with open(file_path, 'w') as f:
    f.write(content)

print("List support fixed and debug traces removed.")
