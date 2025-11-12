# RDatabase: Handle to a loaded VDJdb

`RDatabase` is an object provided by the Rust backend (via extendr) that
holds an in-memory VDJdb table and exposes methods for filtering.

Methods include:

- `new_from_file(path)`

- `new_from_vdjdb(use_fat_db = FALSE)`

- `len()`

- `filter(species = NULL, gene = NULL, min_vdjdb_score = 0)`

- `filter_by_epitope_size(min_size = 1)`

See also:
[match_tcr_df](https://furlan-lab.github.io/vdjmatchR/reference/match_tcr_df.md),
[match_tcr_many_df](https://furlan-lab.github.io/vdjmatchR/reference/match_tcr_many_df.md),
[vdjdb_download](https://furlan-lab.github.io/vdjmatchR/reference/vdjdb_download.md).
