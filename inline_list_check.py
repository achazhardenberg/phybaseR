import sys
import re

file_path = 'R/because.R'
with open(file_path, 'r') as f:
    content = f.read()

# Replace:
# if (is.data.frame(data) || is_list_data) {
# with
# if (is.data.frame(data) || (is.list(data) && !is.data.frame(data))) {

target = r'if \(is\.data\.frame\(data\) \|\| is_list_data\) \{'
# Handle potential whitespace differences
pattern = r'if\s*\(is\.data\.frame\(data\)\s*\|\|\s*is_list_data\)\s*\{'

replacement = 'if (is.data.frame(data) || (is.list(data) && !is.data.frame(data))) {'

new_content = re.sub(pattern, replacement, content)

if new_content == content:
    print("Could not find pattern to replace.")
    # Check context in file
    idx = content.find('is_list_data')
    if idx != -1:
        print("Found 'is_list_data' at context:")
        print(content[idx-50:idx+50])
    
    # Try simpler replacement if manual edit messed up regex
    target_simple = 'if (is.data.frame(data) || is_list_data)'
    if target_simple in content:
        new_content = content.replace(target_simple, 'if (is.data.frame(data) || (is.list(data) && !is.data.frame(data)))')
    else:
        sys.exit(1)

with open(file_path, 'w') as f:
    f.write(new_content)

print("Inlined list check successfully.")
