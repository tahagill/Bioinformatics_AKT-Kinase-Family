# AKT Kinase Family Bioinformatics Analysis

### Update: Most recent update integrates Docker setup and auto PDF analysis report generation

## Overview
This project focuses on the AKT kinase family, specifically AKT1, AKT2, and AKT3, which are critical regulators of cell survival, growth, and metabolism. The project performs a comprehensive bioinformatics analysis of these genes, including sequence alignment, phylogenetic tree construction, mutation analysis, 3D structure visualization, and pathway analysis. The goal is to understand the structural, functional, and evolutionary relationships between these genes and their roles in diseases such as cancer and metabolic disorders.

## Genes Involved
- **AKT1**: Plays a key role in cell survival and growth. Mutations in AKT1 are associated with cancers like breast cancer and Proteus syndrome.
- **AKT2**: Involved in insulin signaling and glucose metabolism. Mutations in AKT2 are linked to severe insulin resistance and type 2 diabetes.
- **AKT3**: Primarily expressed in the brain and associated with brain development disorders like megalencephaly.

## Bioinformatics Tools Used
- **ClustalOmega**: For multiple sequence alignment (MSA) of AKT protein and DNA sequences.
- **Phylogenetic Tree Construction**: Using Biopython's DistanceTreeConstructor to build evolutionary trees.
- **3D Structure Visualization**: Using py3Dmol to visualize protein structures and map mutations.
- **Sequence Logo Generation**: Using logomaker to visualize conserved regions in aligned sequences.
- **Heatmaps and Plots**: Using matplotlib, seaborn, and plotly for data visualization.

## Databases Used
- **NCBI (National Center for Biotechnology Information)**: For fetching DNA and protein sequences, ClinVar mutation data, and PDB structures.
- **UniProt**: For domain architecture and functional annotation of AKT proteins.
- **KEGG (Kyoto Encyclopedia of Genes and Genomes)**: For pathway analysis of the PI3K-AKT-mTOR signaling pathway.
- **GTEx (Genotype-Tissue Expression)**: For gene expression analysis across human tissues.

## Tech Stack
- **Python**: Primary programming language for bioinformatics analysis.
- **Biopython**: For sequence manipulation, alignment, and phylogenetic tree construction.
- **Pandas**: For data manipulation and analysis.
- **Matplotlib/Seaborn/Plotly**: For data visualization.
- **py3Dmol**: For interactive 3D protein structure visualization.
- **Docker**: Containerization for reproducible execution environment
- **GitHub Actions**: CI/CD pipeline for automated testing and deployment
- **Google Colab**: For running the project in a cloud environment.

## Project Structure
The project is organized into the following directories:
```
Bioinformatics_AKT-Kinase-Family/
├── Dockerfile # Docker configuration
├── .dockerignore # Docker exclusion rules
├── codebase.py # Functions pipeline for the analysis
├── main.py # Main script to execute the pipeline
├── requirements.txt # Required Python packages
├── README.md # This file
├── .gitignore # Git ignore file
├── temp_input_protein.fasta # Temporary input file for protein sequences
└── kinase_project/ # Created after running the code
├── data/
│ ├── sequences/ # DNA and protein sequences for AKT1, AKT2, AKT3
│ ├── mutations/ # ClinVar mutation data
│ ├── structures/ # PDB files for AKT1, AKT2, AKT3
└── results/ # Output files (plots, alignments, trees, etc.)
```

## Questions Answered by the Codebase

### Sequence Analysis:
- What are the conserved regions in the AKT family proteins?
- How do the DNA and protein sequences of AKT1, AKT2, and AKT3 compare?

### Phylogenetic Analysis:
- What is the evolutionary relationship between AKT1, AKT2, and AKT3?

### Mutation Analysis:
- What are the clinically significant mutations in the AKT genes?
- How do these mutations affect protein structure and function?

### 3D Structure Visualization:
- Where are the mutations located in the 3D structure of AKT proteins?
- What is the impact of mutations on protein domains (e.g., PH domain, kinase domain)?

### Domain Architecture:
- What are the functional domains in AKT1, AKT2, and AKT3?
- How do these domains contribute to protein function?

### Expression Analysis:
- How are AKT1, AKT2, and AKT3 expressed across different human tissues?
- Which tissues show the highest expression of each AKT gene?

### Pathway Analysis:
- How do AKT genes fit into the PI3K-AKT-mTOR signaling pathway?
- What are the key interactions and regulatory mechanisms in this pathway?


# AKT-Kinase Family Analysis

This repository contains tools for analyzing the AKT kinase protein family.

## How to Run the Project

### Option 1: Google Colab (No Dependency Issues but might be outdated or maintained poorly)

Direct Link to Google Colab Notebook (Preffereed as less dependancy issues): https://colab.research.google.com/drive/1jPHj7jZNLtjZW6wW1zwzeK8uDN3d8QRS?usp=sharing]

### Option 2: Docker (Recommended)

#### Using Pre-built Image from Docker Hub (you need docker cli / hub installed)

To pull the image from DockerHub:

```bash
docker pull tahagill/akt-analysis:latest
```

```bash
docker run -it --rm \
  -v "$(pwd)/results:/app/kinase_project/results" \
  -e ENTREZ_EMAIL="your@email.com" \
  tahagill/akt-analysis
```

#### If you want to build from source !!  (you need docker cli / hub installed)

```bash
git clone https://github.com/tahagill/Bioinformatics_AKT-Kinase-Family.git
cd Bioinformatics_AKT-Kinase-Family
docker build -t akt-analysis .
docker run -it --rm \
  -v "$(pwd)/results:/app/kinase_project/results" \
  -e ENTREZ_EMAIL="your@email.com" \
  akt-analysis
```


### Option 3: Local Execution (old school)
```bash
git clone https://github.com/tahagill/Bioinformatics_AKT-Kinase-Family.git
cd Bioinformatics_AKT-Kinase-Family
pip install -r requirements.txt
python main.py
```

## Docker Image Verification
The pre-built Docker image has been successfully pushed to Docker Hub:
```bash
docker push tahagill/akt-analysis:latest
```

Output:
```
The push refers to repository [docker.io/tahagill/akt-analysis]
0824de02b755: Pushed
ee4ab3e14ed5: Pushed
ff1399ac0930: Mounted from library/python
4b785e93aa71: Mounted from library/python
73b62020d887: Pushed
dac1d3453b30: Mounted from library/python
1fc0a4ca81a8: Pushed
7cf63256a31a: Mounted from library/python
c947297ba602: Pushed
latest: digest: sha256:bb1d64b8a14ff5f2fafed57eeb1f0e5dda3041e39fa9409ed529c9463d7ca641 size: 2272
```

## Requirements
* Docker (for Option 2 )
* Python 3.8+ (for Option 3)
* See `requirements.txt` for Python dependencies

4. **View Results**:
   - All outputs (plots, alignments, trees, etc.) will be saved in the `results/` directory.
   - Interactive 3D visualizations will be displayed in the notebook or saved as HTML files.
   - All data will be saved in `data/` directory 

## Key Features of the Codebase

### Sequence Fetching and Alignment:
- Fetches DNA and protein sequences from NCBI.
- Performs multiple sequence alignment (MSA) using ClustalOmega.

### Phylogenetic Tree Construction:
- Builds a phylogenetic tree to visualize evolutionary relationships.

### Multiple Sequence Alignment & Sequence Alignment:
-  Uses ClustalOmega 

### Mutation Analysis:
- Fetches ClinVar mutation data and maps mutations to protein structures.
- Highlights mutations in 3D protein structures.

### Domain Architecture:
- Fetches domain information from UniProt and visualizes domain architecture.

### Expression Analysis:
- Analyzes gene expression data from GTEx and Human Protein Atlas.

### Pathway Analysis:
- Visualizes the PI3K-AKT-mTOR pathway using KEGG data.

## Future Enhancements
- Integrate more mutation databases (e.g., COSMIC, dbSNP).
- Add machine learning models to predict the impact of mutations.
- Expand pathway analysis to include other signaling pathways.


## Acknowledgments
- NCBI, UniProt, KEGG, and GTEx for providing the data used in this project.
- Biopython, Pandas, and other open-source libraries for enabling this analysis.
