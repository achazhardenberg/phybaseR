import sys

file_path = 'R/because.R'
with open(file_path, 'r') as f:
    lines = f.readlines()

# 1-based index 798 corresponds to list index 797
start_idx = 797
end_idx = 801 # inclusive delete range

# Check if line 797 contains 'else' (safety check)
if 'else' in lines[start_idx]:
    # Replace lines 797-800 with clean line
    # Lines involved are 798, 799, 800, 801
    
    # Wait, in view 1263:
    # 798: data <- data_list else {
    # 799:       message(...)
    # 800:    }
    # 801: }
    
    # We want line 798 to be '    data <- data_list\n'
    # And remove 799, 800, 801.
    
    lines[start_idx] = '    data <- data_list\n'
    # Remove lines at idx 798, 799, 800?
    
    # Slicing: replace lines[798:801] with nothing?
    # Index 798 corresponds to line 799.
    # Index 799 corresponds to line 800.
    # Index 800 corresponds to line 801.
    
    del lines[start_idx+1 : start_idx+4] # delete 3 lines following
    
    with open(file_path, 'w') as f:
        f.writelines(lines)
    print("Syntax else repaired.")
else:
    print("Safety check failed: 'else' not found on target line.")
    print(lines[start_idx])
    sys.exit(1)
