# GRAB-ACh Analysis Pipeline

MATLAB pipeline accompanying:

**Goodwin et al. (2026)**

This repository contains the analysis pipeline used to process 2-photon GRAB-ACh imaging acquired during virtual reality behaviour.

---

## Requirements

- MATLAB R2021b or newer (tested)
- Image Processing Toolbox
- Signal Processing Toolbox

External dependency:

- ScanImageTiffReader (Vidrio Technologies)


## Input data

The pipeline requires

- ScanImage TIFF imaging stack
- Virmen behavioural log files
- Matching imaging and behavioural session information

Specify the imaging and behaviour locations at the top of the main script:

```matlab
pat = '...';
behaviourRoot = '...';
```

---

## Running the analysis

Run

```matlab
main/grabACh_fluorescenceVsSpeed.m
```

The script will

1. Load imaging data
2. Motion-correct the imaging stack
3. Define the GRAB-ACh ROI
4. Perform background subtraction
5. Calculate ΔF/F
6. Load behavioural data
7. Align imaging and behaviour
8. Generate summary figures
9. Save processed variables to `data.mat`

---

## Output

The pipeline generates

- ΔF/F traces
- Behaviour-aligned variables
- Summary figures (.svg)
- `data.mat` containing processed imaging and behavioural variables
