[sceExplorer:]{.underline}

[Raw data]{.underline}: [ ]{.underline}

-   Input: Count matrix/ .f5 file/ 10X files/ .rds object and annotation
    file

-   Output: .rds object (SingleCellExperiment class)

[QC filter (]{.underline}SingleCellExperiment):

-   Input: .rds object

-   Output: filtered data

[Normalization]{.underline}:

-   Input: QC-filtered data

<!-- -->

-   Output: Normalized counts

[Highly Variable Gene Selection]{.underline}:

-   Input: Normalized counts

-   Output: VGcounts

[Dimension Reduction]{.underline}:

i.  [Linear]{.underline}:

-   Input: VGcounts

-   Output: PCA

ii. [Non-Linear]{.underline}:

-   Input: PCA

-   Output: Umap/tSNE

[Clustering]{.underline}:

-   Input: PCA/Umap/tSNE

-   Output: Clusters

[Batch Correction]{.underline}:

i.  Correction on Expression Matrix

-   Input: VGcounts

-   Output: BEVGcounts

    -   [Dimension Reduction]{.underline}:

        -   [Linear]{.underline}:

            -   Input: BEVGcounts

            -   Output: BEPCA

ii. Correction on Dimension Reduction

-   Input: counts

-   Output: BEPCA

[Non-Linear]{.underline}:

-   Input: BEPCA

-   Output: Umap/tSNE

[Clustering]{.underline}:

-   Input: BEPCA/Umap/tSNE

-   Output: Clusters

[DEG Analysis]{.underline}:

-   Input: Raw counts and comparison vector

-   Output: DEG table

[Cell Type Annotation]{.underline}:

-   Input: Count Matrix and clusters

-   Output: Cell type

[Cell Development]{.underline}:

-   Input: Cell Type

-   Output: Cell Development
