# Setup file for testthat
# This file is automatically sourced before tests run

# Transform function to strip trailing whitespace from snapshot output
# This ensures snapshots are consistent across platforms (macOS, Linux, Windows)
strip_trailing_ws <- function(x) {
  gsub("[[:space:]]+$", "", x)
}
