# Debug sort(list()) behavior
empty_list <- list()
print("Sorting empty list:")
tryCatch(
    {
        sort(empty_list)
    },
    error = function(e) {
        print(paste("FAILURE:", e$message))
    }
)

print("Sorting numeric(0):")
print(sort(numeric(0)))
