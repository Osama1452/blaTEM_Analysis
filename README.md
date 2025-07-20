TEM-1 β-lactamase Mutation and Structural Analysis


This repository contains an in-depth bioinformatics analysis of the TEM-1 β-lactamase gene, focusing on mutation detection, structural visualization, nucleotide diversity, and sequence diversity metrics. The project leverages sequence data retrieved from NCBI, performs multiple sequence alignments, calculates nucleotide diversity, and visualizes structural mutations using 3D modeling techniques. The included scripts and visualizations aim to provide a thorough understanding of the genetic and structural variations in TEM-1 β-lactamase sequences.

Overview

The project includes Python scripts and Jupyter notebooks that:

Retrieve TEM-1 sequences from the NCBI database.
Conduct multiple sequence alignments using Clustal Omega.
Calculate nucleotide diversity (π) and identify single nucleotide polymorphisms (SNPs).
Generate a variety of visualizations, including heatmaps, 3D plots, sliding window analyses, and bubble charts.
Visualize the 3D structure of the TEM-1 protein with highlighted mutation sites using NGL Viewer or similar tools.

The accompanying PDF document (3d_bar_nucleotide-merged.pdf) contains the results of these analyses, including diversity plots, SNP distributions, nucleotide frequency bar graphs, and 3D protein structure reconstructions.

Installation

Dependencies

To run the scripts and reproduce the results, install the following dependencies:

1- Python 3.x

2- Biopython: For sequence handling and NCBI queries (pip install biopython)

3- NumPy: For numerical computations (pip install numpy)

4- Matplotlib: For plotting and visualizations (pip install matplotlib)

5- Pandas: For data manipulation and CSV handling (pip install pandas)

6- Scikit-learn: For PCA analysis (pip install scikit-learn)

7- UMAP-learn: For dimensionality reduction (pip install umap-learn)

8- Seaborn: For enhanced visualizations (pip install seaborn)

9- NGLView: For 3D protein structure visualization (pip install nglview)

10- IPywidgets: For interactive Jupyter notebook support (pip install ipywidgets)

11- Clustal Omega: For multiple sequence alignment (download and install from Clustal Omega)

Setup

Clone the repository:
git clone <repository-url>
cd mutation_structural_analysis_ncbi


Install the required Python packages:
pip install -r requirements.txt


Ensure Clustal Omega is installed and the path to clustalo.exe is correctly set in the script (e.g., D:\Osama\biotec\Bioinformatics\mutation_structural_analysis_ncbi\clustal-omega-1.2.2-win64\clustalo.exe).





Usage

Run the Jupyter Notebook or individual Python scripts in the mutation_structural_analysis_ncbi directory.

The script will:

Download TEM-1 sequences from NCBI.
Perform alignment and diversity calculations.
Generate plots and save them as PNG files or the compiled PDF (3d_bar_nucleotide-merged.pdf).
Create a 3D structure visualization and save it as an HTML file or image.


View the generated files:

Sequence data: tem1_sequences.fasta, filtered_sequences.fasta
Alignment: aligned_sequences.aln
Plots: Various .png files (e.g., sliding_window_pi.png, 3d_nucleotide_diversity.png) and 3d_bar_nucleotide-merged.pdf
3D Viewer: mutations_structure.html



Results
The analysis of TEM-1 β-lactamase sequences has yielded a rich set of results, detailed in the 3d_bar_nucleotide-merged.pdf file. Below is an exhaustive breakdown of the findings:
1. 3D Plot of Nucleotide Frequency (Page 1)

This plot provides a three-dimensional representation of nucleotide frequency across the TEM-1 sequences. The visualization highlights the distribution of nucleotides (A, T, C, G) in a spatial context, with peaks indicating regions of high frequency. The presence of "A" annotations suggests a focus on adenine-rich regions, which may correlate with functional or structural significance.

2. Bubble Chart of Sequence Diversity and SNPs (Page 3)

A bubble chart illustrates the diversity of sequences and the distribution of single nucleotide polymorphisms (SNPs). The size and position of bubbles reflect the magnitude of diversity and the frequency of SNPs at specific positions. This chart is crucial for identifying hotspots of genetic variation within the TEM-1 gene.

3. 3D Protein Structure Reconstruction (Page 5)

The PDF includes a 3D reconstruction of the TEM-1 protein structure, likely derived from cryo-electron microscopy (cryo-EM) data. The irregular, folded appearance with a dense structure suggests a complex protein conformation. This visualization is instrumental for studying the protein's function, its interactions with other molecules (e.g., β-lactam antibiotics), and its spatial arrangement at near-atomic resolution. The repetition of descriptive text emphasizes the reliability of cryo-EM as a technique for such analyses.

4. Distribution of Nucleotides in blaTEM Sequences (Page 9)

This section presents a detailed bar graph of nucleotide frequencies across the blaTEM sequences. The graph quantifies the proportion of A, T, C, and G at each position, providing insights into sequence conservation and variability. The extensive numerical data (e.g., "0.0" repetitions) indicates a baseline or low-variance region, with occasional deviations highlighting mutation-prone sites.

5. Bar Graph of Nucleotide Frequencies (Page 11)

A dedicated bar graph further elaborates on nucleotide frequencies, offering a clear visual comparison of A, T, C, and G distributions. This complements the 3D plot and helps in pinpointing regions with significant nucleotide bias, which could be linked to evolutionary pressures or functional domains.

6. Distribution of SNPs in Water Sequences (Page 18)

This analysis focuses on the distribution of SNPs, possibly in a water-related context (e.g., aqueous environment effects on sequence stability). The repeated "0%" values suggest minimal SNP presence in certain regions, while deviations indicate areas of potential interest for further investigation into environmental influences on mutation rates.

7. Nucleotide Helix Representation (Page 21)

The "Nucleotide Helix" label suggests a helical projection or model of the nucleotide sequence, possibly integrating structural and sequence data. This could represent a secondary structure prediction or a visual aid for understanding the 3D arrangement of the DNA or protein.

Additional Observations

Pages with truncated or sparse data (e.g., Pages 2, 4, 6-8, 10, 12-17, 19-20, 22-26) contain partial visualizations or metadata (e.g., sequence indices, "The" annotations, "# 1" markers), which may indicate incomplete rendering or additional context not fully captured in the OCR output. These pages likely serve as supporting material or placeholders for further detailed analysis.

The results collectively provide a multi-dimensional view of TEM-1 β-lactamase, combining genetic diversity metrics with structural insights, making this dataset valuable for researchers studying antibiotic resistance mechanisms.


Contributing

Contributions are welcome! Please fork the repository and submit pull requests with improvements or bug fixes. Ensure to update the documentation and tests accordingly.

License

This project is licensed under the MIT License. See the LICENSE file for details.
