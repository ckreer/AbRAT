### Repertoire Characteristics

#### Gene Segment Usage

Collapsing requires a unique clone identifier in the dataset and reduces gene segment counts to only unique occurrences  
within each clone. Note that if clustering settings are less stringent, a single clone may contain multiple gene segments  
(e.g., multiple light chain V gene segments if not restricted during clustering).

#### Hydrophobicity Values

The **GRAVY (Grand Average of Hydropathy)** score is calculated as follows:

$$
\text{GRAVY} = \frac{\sum_{i=1}^{n} h(a_i)}{n}
$$

where $h(a_i)$ is the hydrophobicity value (using either the **Kyte-Doolittle** or **Eisenberg** scale) for the $i$-th amino acid  
and $n$ is the total number of amino acids in the sequence.
- **Positive GRAVY values** indicate a generally hydrophobic protein.
- **Negative GRAVY values** indicate a more hydrophilic protein.

#### Net Charge Values

The **CDR3 net charge** at pH 7.4 is computed using the Python **peptides** package (via `Peptide.charge(pH=7.4)`), which  
applies the Henderson–Hasselbalch equation with standard pKa values. A positive net charge indicates a basic sequence, whereas  
a negative net charge indicates an acidic sequence.

#### CDR3 Diversity

AbRAT implements two widely used diversity indices to quantify the heterogeneity of CDR3 sequences.  
*Note: The basic implementation of both indices considers only identical sequences, not similar ones, when estimating diversity.*

##### Shannon Index

The **Shannon Index** is calculated as:

$$
H = -\sum_{i=1}^{S} p_i \ln(p_i)
$$

where $p_i$ is the relative frequency of the $i$-th unique CDR3 sequence and $S$ is the total number of unique sequences.
- **Interpretation:** A higher Shannon Index indicates a more diverse repertoire with a more even distribution.

##### Inverse Simpson Index

The **Inverse Simpson Index** is computed as:

$$
D = \frac{1}{\sum_{i=1}^{S} p_i^2}
$$

- **Interpretation:** A higher Inverse Simpson Index indicates a more diverse and evenly distributed repertoire, whereas a lower  
  value suggests that the repertoire is dominated by a few sequences.

##### Subsampling for Fair Comparisons

To compare diversity across different repertoires, each dataset is randomly subsampled to match the size of the smallest dataset.  
The mean diversity index over 20 iterations is reported along with the standard deviation. For the smallest dataset, the index is  
calculated on the complete data (thus, the standard deviation is zero).

#### Additional Repertoire Characteristics

- **Clonal Assignment:** The filtered and clustered repertoire data can be further analyzed to assess clonal relationships and diversity.
- **Visualization:** Various metrics (gene segment usage, hydrophobicity, net charge, diversity indices) are presented on interactive  
  dashboards for in-depth exploration.