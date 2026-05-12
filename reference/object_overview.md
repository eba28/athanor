# Print an overview of a Seurat object.

Prints a structured summary of a Seurat object, covering its assays,
reductions, neighbor graphs, graphs, WNN configuration (if present), and
Misc slot contents.

## Usage

``` r
object_overview(seurat_obj)
```

## Arguments

- seurat_obj:

  The Seurat object.

## Value

A overview of the Seurat object, and invisibly returns the input object
for piping if desired.

## Details

Each section is only printed if the corresponding slot is non-empty, so
the output degrades gracefully on minimal objects. For RNA, the number
of variable features is shown. For ADT, all marker names are listed. The
WNN section is printed only when a `w.nn` neighbor graph is detected and
the `FindMultiModalNeighbors` command is recorded in the object. Misc
slot values are printed inline for short atomic vectors and summarized
by type for larger or complex objects.
