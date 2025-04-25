### Data Format

In {{AbRAT}}, the names of _ab1-files_ serve as a **unique identifier** for a specific sequence. Therefore, 
_ab1-files_ must follow a strict naming convention with **13 positions** separated by underscores ('_').

```
COHORT_SUBJECT_TIMEPOINT_SAMPLE_SUBSET_PLATE_WELL_CHAIN_PRIMERSET_SOURCE-SUBSOURCE_SEQPRIMER_SEQCOMPANY_SEQ-REPEAT.ab1
```
    
The file name is divided into 4 blocks that contain information on:
1. **Sample:** Details about the study cohort, subject ID, timepoint, sample material, and subset.
2. **Cell Position:** Information about the physical location (plate and well) of the cell.
3. **PCR Details:** Information on the PCR performed, including chain (must be HC, KC, or LC), primer set, and source. 
   The subsource is derived from the string after the dash '-' and is optional but recommended.
4. **Sequencing:** Information on the sequencing process, including the sequencing primer, company, and sequence run number.

**Sample Block:**
    ```
    COHORT_SUBJECT_TIMEPOINT_SAMPLE_SUBSET_
    ```
- Includes study cohort, subject ID, timepoint (required for longitudinal samples), sample material, and cell subset.
- The sample and subset fields are flexible; for example, the sample field might indicate the tissue type (e.g., PBMCs, lymph node) 
  or a preselection (e.g., CD20-IgG for gating during FACS analysis). The subset might include the bait protein used for antigen-specific sorting.

**Cell Position Block:**
    ```
    PLATE_WELL
    ```
- Contains the cell's physical storage information: the plate number (sequentially numbered) and the well (formatted as A1–H12).
- Together with the sample block, this information forms the unique *B_CELL_ID*.

**PCR Block:**
    ```
    CHAIN_PRIMERSET_SOURCE-SUBSOURCE
    ```
- Encodes information about the PCR, including the chain (HC, KC, or LC), primer set, and source.
- **Important:** _CHAIN_ must be one of **HC**, **KC**, or **LC** to ensure proper assembly of the B cell receptor. 
  The subsource, if present, is the first string after the dash '-' and aids in identifying repeated PCR events or plasmid batches.

**Sequencing Block:**
    ```
    SEQPRIMER_SEQCOMPANY_SEQ-REPEAT.ab1
    ```
- Contains details of the sequencing process: the sequencing primer, the sequencing company or instrument, and the sequence run number.

###### Example
```
SARS2_IDC10_t1_PBMCs_IgG-S488_1_B8_HC_oPR_2ND-1_IgInt_EF_SEQ-1.ab1
SARS2_IDC10_t1_PBMCs_IgG-S488_1_B8_KC_oPR_MIDI-435_Ckrev_EF_SEQ-2.ab1
```
These two files come from the same B cell sampled from Subject IDC10 in a SARS-CoV-2 study (SARS2) at the first timepoint (t1). 
The cell, obtained from PBMCs, was sorted from an IgG and S-Protein-Alexa488 positive gate on plate 1 in well B8. 
The first file represents the heavy chain (HC) sequence from the first nested PCR attempt (2ND-1) using the openPrimeR set, 
sequenced with the heavy chain reverse primer 'IgInt' at Eurofins (EF) for the first time (SEQ-1). 
The second file represents the kappa chain (KC) sequence from the same cell, derived from a Midi-prep plasmid (plasmid 435) 
with the kappa constant region reverse primer 'Ckrev' at Eurofins (EF). However, the first sequence of this plasmid had quality issues and was re-sequenced (SEQ-2).
    