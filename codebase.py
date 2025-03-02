# !! Install required packages if running on Google Colab !! (remove the #)

#!pip install biopython requests pandas matplotlib seaborn plotly py3Dmol -q
#!apt-get install clustalo
#!pip install logomaker

# If not running on google colab just do pip install -r requirements.txt on ur terminal

import os
import re
import requests
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.graph_objects as go
import logomaker as lm
import py3Dmol

from Bio import Entrez, SeqIO, Phylo, AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.PDB.PDBParser import PDBParser
from IPython.display import HTML, display, Image
from tqdm.auto import tqdm  
from matplotlib.lines import Line2D

Entrez.email = "tahagill99@gmail.com"  ## you can add ur email too.

'''from google.colab import drive
drive.mount('/content/drive')   !!! IF YOU WANT TO STORE DATA ON GOOGLE DRIVE UNCOMMENT THIS !!!
BASE_DIR = "/content/drive/My Drive/kinase_project" '''
BASE_DIR = os.path.join(os.getcwd(), "kinase_project")  # keep this if you want data to be stored on local machine, if you want to work on cloud see above



# Define global directory paths
sequence_dir = os.path.join(BASE_DIR, "data/sequences")
mutation_dir = os.path.join(BASE_DIR, "data/mutations")
structure_dir = os.path.join(BASE_DIR, "data/structures")
results_dir = os.path.join(BASE_DIR, "results")

def create_project_directories():
    """Create project directories if they don't exist."""
    directories = [sequence_dir, mutation_dir, structure_dir, results_dir]
    for directory in directories:
        os.makedirs(directory, exist_ok=True)
        print(f"Directory created or exists: {directory}")

GENES = {
    'AKT1': {
        'dna_refseq': 'NM_001014431',
        'protein_refseq': 'NP_001014431',
        'pdb': '3O96',
        'uniprot': 'P31749'
    },
    'AKT2': {
        'dna_refseq': 'NM_001626',
        'protein_refseq': 'NP_001626',
        'pdb': '2JDR',
        'uniprot': 'P31751'
    },
    'AKT3': {
        'dna_refseq': 'NM_005465',
        'protein_refseq': 'NP_005465',
        'pdb': '4EKK',
        'uniprot': 'Q9Y243'
    }
}
def fetch_sequence(gene, seq_type='dna'):
    """Fetch DNA/protein sequences with validation"""
    seq_dir = os.path.join(sequence_dir, seq_type)
    os.makedirs(seq_dir, exist_ok=True)

    refseq_id = GENES[gene][f'{seq_type}_refseq']
    filename = os.path.normpath(os.path.join(seq_dir, f"{gene}_{seq_type}.fasta"))

    if os.path.exists(filename):
        print(f"Existing {seq_type.upper()} sequence for {gene}")
        return SeqIO.read(filename, "fasta")

    try:
        db = 'nucleotide' if seq_type == 'dna' else 'protein'
        handle = Entrez.efetch(db=db, id=refseq_id, rettype="fasta", retmode="text")
        record = SeqIO.read(handle, "fasta")

        # Sequence validation
        valid_chars = set("ATCGN" if seq_type == 'dna' else "ACDEFGHIKLMNPQRSTVWY")
        if not all(c in valid_chars for c in record.seq):
            raise ValueError(f"Invalid characters in {seq_type.upper()} sequence")

        with open(filename, "w") as f:
            SeqIO.write(record, f, "fasta")
        print(f"{seq_type.upper()} sequence saved for {gene}: {filename}")
        return record

    except Exception as e:
        print(f"Error fetching {seq_type.upper()} sequence for {gene}: {str(e)}")
        return None

def fetch_all_sequences():
    """Test this in a separate cell"""


    for gene in GENES:
        print(f"\n=== Processing {gene} ===")
        dna_seq = fetch_sequence(gene, 'dna')
        protein_seq = fetch_sequence(gene, 'protein')

def fetch_genomic_data():
    """Fetch and validate nucleoti,de sequences from NCBI"""
    all_exists = True
    for gene, config in GENES.items():
        fasta_path = os.path.join(sequence_dir, f"{gene}_sequence.fasta")
        if not os.path.exists(fasta_path):
            all_exists = False
            try:
                handle = Entrez.efetch(db="nucleotide", id=config['dna_refseq'], rettype="fasta", retmode="text")
                record = SeqIO.read(handle, "fasta")
                with open(fasta_path, "w") as f:
                    f.write(f">{gene}\n{str(record.seq)}\n")
                print(f"{gene} nucleotide sequence stored in {fasta_path}")
            except Exception as e:
                print(f" {gene} sequence error: {str(e)}")
        else:
            print(f"{gene} sequence already exists at {fasta_path}")

    if all_exists:
        print("All genomic sequences already exist. Skipping download.")

def fetch_pdb():
    """Download PDB structures"""
    all_exists = True
    for gene, config in GENES.items():
        pdb_id = config['pdb']
        pdb_path = os.path.join(structure_dir, f"{gene}_structure.pdb")
        if not os.path.exists(pdb_path):
            all_exists = False
            try:
                url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
                response = requests.get(url)
                response.raise_for_status()
                with open(pdb_path, "w") as f:
                    f.write(response.text)
                print(f"{gene} structure ({pdb_id}) saved to {pdb_path}")
            except Exception as e:
                print(f"Error fetching {pdb_id}: {str(e)}")
        else:
            print(f"{gene} structure already exists at {pdb_path}")

    if all_exists:
        print(" All PDB structures already exist. Skipping download.")

import os
import subprocess

def align_sequences():
    combined_fasta = os.path.normpath(os.path.join(sequence_dir, "AKT_family.fasta"))
    aligned_fasta = os.path.normpath(os.path.join(results_dir, "AKT_family_aligned.fasta"))

    if not os.path.exists(combined_fasta):
        print("Creating combined FASTA file...")
        sequences = []
        for gene in GENES:
            gene_file = os.path.normpath(os.path.join(sequence_dir, f"{gene}_sequence.fasta"))
            if os.path.exists(gene_file):
                sequences.append(SeqIO.read(gene_file, "fasta"))
            else:
                print(f"Missing sequence file for {gene}: {gene_file}")
                return

        # Write combined file
        with open(combined_fasta, "w") as f:
            SeqIO.write(sequences, f, "fasta")
        print(f"Created combined FASTA at {combined_fasta}")

    # Run alignment
    cmd = f'clustalo -i "{combined_fasta}" -o "{aligned_fasta}" --auto --force'
    try:
        subprocess.run(cmd, check=True, shell=True)
        print(f"✅ MSA completed! Aligned sequences saved to {aligned_fasta}")
    except Exception as e:
        print(f"ClustalOmega error: {str(e)}")
        if not os.path.exists(combined_fasta):
            print("Combined FASTA file was not created successfully")
#  MSA

def perform_msa(genes, seq_type='protein'):
    """Perform multiple sequence alignment using ClustalOmega"""
    import subprocess
    
    # Collect sequences
    sequences = []
    for gene in genes:
        seq = fetch_sequence(gene, seq_type)  
        if seq:
            sequences.append(seq)
    
    if not sequences:
        print("No sequences found for MSA")
        return None

    # Write temporary input file
    input_file = os.path.join(sequence_dir, f"temp_{seq_type}_input.fasta")
    SeqIO.write(sequences, input_file, "fasta")
    
    # Run ClustalOmega
    output_file = os.path.join(results_dir, f"AKT_{seq_type}_aligned.fasta")
    try:
        subprocess.run(
            ["clustalo", "-i", input_file, "-o", output_file, "--auto", "--force"],
            check=True,
            stderr=subprocess.PIPE
        )
        return AlignIO.read(output_file, "fasta")
    except subprocess.CalledProcessError as e:
        print(f"ClustalOmega error: {e.stderr.decode().strip()}")
        return None
    except Exception as e:
        print(f"MSA failed: {str(e)}")
        return None
# Phylogenetic Tree
def build_phylogenetic_tree():
    tree_path = os.path.join(results_dir, "AKT_phylogenetic_tree.newick")
    tree_img_path = os.path.join(results_dir, "phylogenetic_tree.png")

    # Force rebuilding the tree to fix visualization issues
    aligned_file = os.path.join(results_dir, "AKT_family_aligned.fasta")

    # Check if alignment file exists
    if not os.path.exists(aligned_file):
        print(f"Aligned file not found at {aligned_file}. Run align_sequences() first.")
        return

    alignment = AlignIO.read(aligned_file, "fasta")

    calculator = DistanceCalculator('identity')
    dm = calculator.get_distance(alignment)

    constructor = DistanceTreeConstructor(calculator, 'upgma')
    tree = constructor.build_tree(alignment)

    # Save Newick tree
    Phylo.write([tree], tree_path, "newick") # type: ignore

    #labels and branches
    plt.figure(figsize=(10, 6))
    ax = plt.subplot(111)

    # Draw tree with customizations
    Phylo.draw(tree, axes=ax, show_confidence=False, do_show=False) # type: ignore

    plt.title("AKT Family Phylogenetic Tree")

    # Save to file
    plt.savefig(tree_img_path, dpi=300, bbox_inches='tight')
    display(plt.gcf())
    plt.close()

    print("Phylogenetic tree saved as PNG and Newick file!")

    # Display tree data in text form
    print("\nTree in Newick format:")
    with open(tree_path, 'r') as f:
        newick = f.read().strip()
        print(newick)



def map_mutations_to_structure():
    """Map mutations onto 3D structures for AKT1, AKT2, and AKT3."""
    structure_view_path = os.path.join(results_dir, "akt_mutations_3d.html")
    structure_img_path = os.path.join(results_dir, "akt_mutations_3d.png")
    mutation_data_path = os.path.join(mutation_dir, "akt_mutations.csv")

    # Load mutation data
    if os.path.exists(mutation_data_path):
        mutations_df = pd.read_csv(mutation_data_path)
    else:
        print("Mutation data not found. Run fetch_clinvar_data() first.")
        return

    # Check if mutation data has required columns
    if mutations_df.empty or 'ProteinChange' not in mutations_df.columns:
        print("No mutation data available for 3D visualization")
        # Just display the structures without mutations
        for gene in GENES:
            pdb_path = os.path.join(structure_dir, f"{gene}_structure.pdb")
            if not os.path.exists(pdb_path):
                print(f"Missing structure for {gene}")
                continue

            print(f"Displaying structure for {gene} (without mutations)")
            viewer = py3Dmol.view()
            viewer.addModel(open(pdb_path).read(), 'pdb')
            viewer.setStyle({'cartoon': {'color': 'spectrum'}})
            viewer.zoomTo()
            display(viewer.show())
        return

    # If we have valid mutation data
    for gene in GENES:
        pdb_path = os.path.join(structure_dir, f"{gene}_structure.pdb")
        if not os.path.exists(pdb_path):
            print(f"Missing structure for {gene}")
            continue

        viewer = py3Dmol.view()
        viewer.addModel(open(pdb_path).read(), 'pdb')
        viewer.setStyle({'cartoon': {'color': 'spectrum'}})

        # Highlight mutations
        if 'GeneSymbol' in mutations_df.columns:
            mutations = mutations_df[mutations_df['GeneSymbol'] == gene]
        else:
            print(f"'GeneSymbol' column not found. Assuming all mutations belong to {gene}.")
            mutations = mutations_df  # Use all mutations for this gene

        if mutations.empty:
            print(f"No mutations found for {gene}")
        else:
            print(f"Highlighting {len(mutations)} mutations for {gene}")
            for _, row in mutations.iterrows():
                try:
                    match = re.search(r'\d+', row['ProteinChange'])
                    if match:
                        pos = int(match.group())
                        viewer.addStyle({'resi': str(pos)},
                                      {'sphere': {'radius': 1.5, 'color': 'red'}})
                except Exception as e:
                    print(f"Error highlighting mutation {row['ProteinChange']}: {str(e)}")

        viewer.zoomTo()
        display(viewer.show())
        print(f"{gene} structure rendered")

    # Save HTML representation
    html_content = viewer._make_html()
    with open(structure_view_path, 'w') as f:
        f.write(html_content)

    print(f"Interactive 3D structure with mutations saved to {structure_view_path}")

# Mutation Analysis
def mutation_analysis():
    logo_path = os.path.join(results_dir, 'sequence_logo.png')
    structure_img_path = os.path.join(results_dir, 'AKT1_3D_structure.png')
    structure_jpg_path = os.path.join(results_dir, 'AKT1_3D_structure.jpg')

    # We'll rebuild the visualizations to ensure they display properly
    aligned_file = os.path.join(results_dir, "AKT_family_aligned.fasta")

    # Check if alignment file exists
    if not os.path.exists(aligned_file):
        print(f"Aligned file not found at {aligned_file}. Run align_sequences() first.")
        return

    alignment = AlignIO.read(aligned_file, "fasta")

    # Generate and display sequence logo
    plt.figure(figsize=(15, 4))
    matrix = lm.alignment_to_matrix([str(rec.seq) for rec in alignment])
    logo = lm.Logo(matrix, color_scheme='chemistry')
    plt.title("AKT Family Sequence Logo")

    # Save and display
    plt.savefig(logo_path, bbox_inches='tight', dpi=300)
    display(plt.gcf())
    plt.close()
    print(f"Sequence logo saved to {logo_path}")

    # Use Plotly for better 3D structure visualization (more compatible and interactive)
    pdb_path = os.path.join(structure_dir, "AKT1_structure.pdb")

    # Check if PDB file exists
    if not os.path.exists(pdb_path):
        print(f"PDB file not found at {pdb_path}. Run fetch_pdb() first.")
        return

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("AKT1", pdb_path)
    if structure is None:
        print("Failed to parse PDB structure")
        return

    # Extract coordinates and atom types
    coords = []
    atom_types = []
    for atom in structure.get_atoms():
        coords.append(atom.get_coord())
        atom_types.append(atom.element)

    coords = np.array(coords)

    # Create color mapping for different atom types
    color_map = {
        'C': '#333333',  # Carbon
        'N': '#3050F8',  # Nitrogen
        'O': '#FF0D0D',  # Oxygen
        'S': '#FFFF30',  # Sulfur
        'P': '#FF8000',  # Phosphorus
        'H': '#FFFFFF',  # Hydrogen
        ' ': '#654321'   # Default
    }

    colors = [color_map.get(atom, color_map[' ']) for atom in atom_types]

    # Create interactive 3D plot with Plotly
    fig = go.Figure(data=[go.Scatter3d(
        x=coords[:, 0],
        y=coords[:, 1],
        z=coords[:, 2],
        mode='markers',
        marker=dict(
            size=2,
            color=colors,
            opacity=0.8
        ),
        text=atom_types,
        hoverinfo='text'
    )])

    fig.update_layout(
        title="AKT1 Protein Structure (Interactive)",
        scene=dict(
            xaxis_title='X',
            yaxis_title='Y',
            zaxis_title='Z'
        ),
        width=800,
        height=800,
        margin=dict(l=0, r=0, b=0, t=30)
    )
    fig.show()

    fig_static = plt.figure(figsize=(10, 10))
    ax = fig_static.add_subplot(111, projection='3d')

    # custom colors based on atom type
    for i, (x, y, z) in enumerate(coords):
        ax.scatter(x, y, z, color=colors[i], s=2, alpha=0.7) # type: ignore

    ax.set_title("AKT1 Structure")
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    # In mutation_analysis() function:

    ax.set_zlabel('Z')  # type: ignore # Now works on 3D axis

    # Save in multiple formats
    plt.savefig(structure_img_path, dpi=300, bbox_inches='tight')
    plt.savefig(structure_jpg_path, dpi=300, bbox_inches='tight')

    # Display the static version
    display(plt.gcf())
    plt.close()

    print(f"3D structure visualizations saved to {structure_img_path} and {structure_jpg_path}")
    print("Interactive 3D structure also displayed in notebook")

# CLINVAR INTEGRATION
def fetch_clinvar_data():
    """Fetch ClinVar data only if not already downloaded"""
    filename = os.path.join(mutation_dir, "akt_mutations.csv")
    
    # Check if file exists and is valid (prevents redundant download)
    if os.path.exists(filename):
        try:
            existing_df = pd.read_csv(filename)
            
            # Basic validation of required columns
            required_cols = ['GeneSymbol', 'ProteinChange', 'ClinicalSignificance']
            if all(col in existing_df.columns for col in required_cols):
                print(f"✓ Using existing ClinVar data from {filename}")
                return existing_df
            else:
                print("Existing file missing required columns, redownloading...")
        except Exception as e:
            print(f"Corrupted existing file: {str(e)}, redownloading...")

    # File doesn't exist or is invalid - download fresh
    try:
        print("Downloading latest ClinVar data...")
        clinvar_url = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
        df = pd.read_csv(clinvar_url, sep='\t', compression='gzip', dtype={'Name': str})

        akt_mutations = df[
            (df['GeneSymbol'].isin(['AKT1', 'AKT2', 'AKT3'])) &
            (df['Type'] == 'single nucleotide variant') &
            (df['Assembly'] == 'GRCh38')
        ].copy()

        akt_mutations['ProteinChange'] = akt_mutations['Name'].str.extract(r'(p\.[A-Za-z]+\d+[A-Za-z]+)')
        akt_mutations.dropna(subset=['ProteinChange'], inplace=True)
        
        # Save validated data
        akt_mutations.to_csv(filename, index=False)
        print(f"Saved new ClinVar data to {filename}")
        return akt_mutations

    except Exception as e:
        print(f"ClinVar error: {str(e)}")
        return pd.DataFrame(columns=['GeneSymbol', 'ProteinChange', 'ClinicalSignificance'])

def validate_clinvar_data():
    """Test this in separate cell"""
    filename = os.path.join(mutation_dir, "akt_mutations.csv")

    try:
        df = pd.read_csv(filename)
        if 'ProteinChange' not in df.columns and 'Name' in df.columns:
            df['ProteinChange'] = df['Name'].str.extract(r'(p\.[A-Za-z]+\d+[A-Za-z]+)')
            df.dropna(subset=['ProteinChange'], inplace=True)
            df.to_csv(filename, index=False)
        return df
    except:
        return pd.DataFrame()
#   Conservation Analysis
def calculate_conservation():
    conservation_plot_path = os.path.join(results_dir, 'conservation_plot.png')
    aligned_file = os.path.join(results_dir, "AKT_family_aligned.fasta")

    # Check if alignment file exists
    if not os.path.exists(aligned_file):
        print(f"Aligned file not found at {aligned_file}. Run align_sequences() first.")
        return

    alignment = AlignIO.read(aligned_file, "fasta")

    conservation = []
    for i in range(len(alignment[0])):
        column = [rec.seq[i] for rec in alignment]
        conservation.append(sum(1 for aa in column if aa == column[0]) / len(column))

    plt.figure(figsize=(12, 4))
    sns.lineplot(x=range(len(conservation)), y=conservation)
    plt.title("Position-wise Conservation Score")
    plt.xlabel("Alignment Position")
    plt.ylabel("Conservation (%)")

    # Add annotations for highly conserved regions
    high_cons_threshold = 0.9
    for i, cons in enumerate(conservation):
        if cons >= high_cons_threshold:
            plt.plot(i, cons, 'ro', markersize=4)

    # legend
    custom_lines = [Line2D([0], [0], color='red', marker='o', linestyle='None')]
    plt.legend(custom_lines, ['Highly conserved (≥90%)'])

    plt.tight_layout()
    plt.savefig(conservation_plot_path, dpi=300)

    display(plt.gcf())
    plt.close()

    print(f"Conservation plot saved to {conservation_plot_path}")

    # stats
    print(f"\nConservation Statistics:")
    print(f"Average conservation: {np.mean(conservation):.2f}")
    print(f"Number of highly conserved positions (≥90%): {sum(1 for c in conservation if c >= 0.9)}")
    print(f"Number of weakly conserved positions (<50%): {sum(1 for c in conservation if c < 0.5)}")

#  Domain Architecture
def fetch_domain_architecture():
    for gene, config in GENES.items():
        domain_plot_path = os.path.join(results_dir, f"{gene}_domains.png")

        try:
            url = f"https://rest.uniprot.org/uniprotkb/{config['uniprot']}.json"
            res = requests.get(url)
            res.raise_for_status()
            res_json = res.json()

            # Get all feature types for richer visualization
            domains = [ft for ft in res_json['features'] if ft['type'] == 'Domain']
            motifs = [ft for ft in res_json['features'] if ft['type'] == 'Motif']
            regions = [ft for ft in res_json['features'] if ft['type'] == 'Region']
            active_sites = [ft for ft in res_json['features'] if ft['type'] == 'Active site']
            binding_sites = [ft for ft in res_json['features'] if ft['type'] == 'Binding site']

            plt.figure(figsize=(12, 4))

            # Base protein as gray bar
            protein_length = res_json.get('sequence', {}).get('length', 500)
            plt.hlines(y=1, xmin=0, xmax=protein_length, linewidth=10, color='gray', alpha=0.3)

            # Counter for vertical positioning of labels
            y_pos = 1

            # Plot domains with distinct colors
            domain_colors = ['teal', 'blue', 'purple', 'green']
            for i, dom in enumerate(domains):
                start = dom['location']['start']['value']
                end = dom['location']['end']['value']
                color = domain_colors[i % len(domain_colors)]
                plt.hlines(y=y_pos, xmin=start, xmax=end, linewidth=20, color=color, alpha=0.8)
                plt.text((start+end)/2, y_pos+0.3, dom.get('description', 'Domain'),
                        ha='center', va='bottom', fontsize=9, fontweight='bold')

            # Plot motifs
            for motif in motifs:
                start = motif['location']['start']['value']
                end = motif['location']['end']['value']
                plt.hlines(y=y_pos-0.2, xmin=start, xmax=end, linewidth=10, color='orange', alpha=0.7)
                plt.text((start+end)/2, y_pos-0.4, motif.get('description', 'Motif'),
                        ha='center', va='top', fontsize=8)

            # Mark active sites with red diamonds
            for site in active_sites:
                pos = site['location']['start']['value']
                plt.scatter(pos, y_pos, marker='D', s=80, color='red')
                plt.text(pos, y_pos+0.25, site.get('description', 'Active site'),
                        ha='center', va='bottom', fontsize=8, rotation=45)

            # Mark binding sites with yellow stars
            for site in binding_sites:
                pos = site['location']['start']['value']
                plt.scatter(pos, y_pos, marker='*', s=100, color='gold')
                plt.text(pos, y_pos-0.25, site.get('description', 'Binding site'),
                        ha='center', va='top', fontsize=8, rotation=45)

            plt.yticks([])
            plt.title(f"{gene} Domain Architecture")
            plt.xlabel("Amino Acid Position")
            plt.xlim(0, protein_length + 10)
            plt.ylim(0, 2)

            # Add a legend
            from matplotlib.patches import Patch
            from matplotlib.lines import Line2D
            legend_elements = [
                Patch(facecolor='teal', alpha=0.8, label='Domain'),
                Patch(facecolor='orange', alpha=0.7, label='Motif'),
                Line2D([0], [0], marker='D', color='w', markerfacecolor='red', markersize=10, label='Active site'),
                Line2D([0], [0], marker='*', color='w', markerfacecolor='gold', markersize=10, label='Binding site')
            ]
            plt.legend(handles=legend_elements, loc='upper right')

            plt.tight_layout()
            plt.savefig(domain_plot_path, dpi=300, bbox_inches='tight')
            display(plt.gcf())
            plt.close()

            print(f"Enhanced domain plot for {gene} saved to {domain_plot_path}")

            # Print domain information as text for reference
            print(f"\n{gene} Domain Information:")
            for i, dom in enumerate(domains):
                print(f"  • {dom.get('description', 'Domain')}: Position {dom['location']['start']['value']}-{dom['location']['end']['value']}")

        except Exception as e:
            print(f"Domain error for {gene}: {str(e)}")

    print("All domain plots processed!")
def analyze_disease_associations():
    """Associate AKT mutations with disease phenotypes"""

    disease_plot_path = os.path.join(results_dir, "akt_disease_associations.png")
    mutation_heatmap_path = os.path.join(results_dir, "akt_mutation_heatmap.png")
    disease_data_path = os.path.join(results_dir, "akt_disease_data.csv")

    # synthetic data based on known associations !! [will change in future commits]

    # Known disease associations for AKT genes !!!  [will change in future commits]
    disease_data = [
        # AKT1
        {"Gene": "AKT1", "Disease": "Proteus syndrome", "Mutation": "E17K", "Significance": "Pathogenic", "Count": 27},
        {"Gene": "AKT1", "Disease": "Breast cancer", "Mutation": "E17K", "Significance": "Likely pathogenic", "Count": 18},
        {"Gene": "AKT1", "Disease": "Colorectal cancer", "Mutation": "E17K", "Significance": "Likely pathogenic", "Count": 15},
        {"Gene": "AKT1", "Disease": "Ovarian cancer", "Mutation": "Various", "Significance": "Uncertain", "Count": 12},
        {"Gene": "AKT1", "Disease": "Endometrial cancer", "Mutation": "Various", "Significance": "Likely pathogenic", "Count": 9},

        # AKT2
        {"Gene": "AKT2", "Disease": "Hypoglycemia", "Mutation": "E17K", "Significance": "Pathogenic", "Count": 8},
        {"Gene": "AKT2", "Disease": "Severe insulin resistance", "Mutation": "Various", "Significance": "Pathogenic", "Count": 13},
        {"Gene": "AKT2", "Disease": "Type 2 Diabetes", "Mutation": "Various", "Significance": "Risk factor", "Count": 22},
        {"Gene": "AKT2", "Disease": "Colorectal cancer", "Mutation": "Various", "Significance": "Uncertain", "Count": 7},

        # AKT3
        {"Gene": "AKT3", "Disease": "Megalencephaly", "Mutation": "E17K", "Significance": "Pathogenic", "Count": 14},
        {"Gene": "AKT3", "Disease": "Hemimegalencephaly", "Mutation": "Various", "Significance": "Pathogenic", "Count": 11},
        {"Gene": "AKT3", "Disease": "Brain development disorders", "Mutation": "Various", "Significance": "Likely pathogenic", "Count": 19},
        {"Gene": "AKT3", "Disease": "Neurological disorders", "Mutation": "Various", "Significance": "Uncertain", "Count": 16}
    ]

    # Convert to DataFrame
    disease_df = pd.DataFrame(disease_data)

    # Save to CSV
    disease_df.to_csv(disease_data_path, index=False)
    print(f"Disease association data saved to {disease_data_path}")

    plt.figure(figsize=(14, 8))

    # Define a custom palette for significance
    significance_colors = {
        "Pathogenic": "#FF5555",
        "Likely pathogenic": "#FFAA55",
        "Uncertain": "#FFFF55",
        "Risk factor": "#55AAFF"
    }

    # Create grouped bar chart
    g = sns.catplot(
        data=disease_df,
        kind="bar",
        x="Gene",
        y="Count",
        hue="Significance",
        palette=significance_colors,
        height=6,
        aspect=1.5
    )

    g.despine(left=True)
    g.set_axis_labels("AKT Gene", "Number of Reported Cases")
    if g.legend:
        g.legend.set_title("Clinical Significance")


    plt.title("AKT Gene Family Disease Associations", fontsize=16, pad=20)
    plt.tight_layout()
    plt.savefig(disease_plot_path, dpi=300, bbox_inches='tight')

    display(g.fig)
    plt.close()

    print(f"Disease association plot saved to {disease_plot_path}")

    # Create a heatmap of mutation types across genes
    plt.figure(figsize=(12, 8))

    # Added some synthetic mutation data for a more complete heatmap
    mutation_types = ['Missense', 'Frameshift', 'Nonsense', 'Splice-site', 'Indel']
    mutation_data = []

    np.random.seed(42)  # For reproducibility !!!!!!

    for gene in ["AKT1", "AKT2", "AKT3"]:
        for m_type in mutation_types:
            # Different genes have different mutation profiles
            if gene == "AKT1" and m_type == "Missense":
                base = 35
            elif gene == "AKT2" and m_type == "Frameshift":
                base = 28
            elif gene == "AKT3" and m_type == "Splice-site":
                base = 22
            else:
                base = 10

            count = max(0, int(base + np.random.normal(0, 5)))

            mutation_data.append({
                "Gene": gene,
                "Mutation_Type": m_type,
                "Count": count
            })

    # Convert to DataFrame
    mutation_df = pd.DataFrame(mutation_data)
    mutation_pivot = mutation_df.pivot(index="Gene", columns="Mutation_Type", values="Count")

    # Create heatmap
    sns.heatmap(
        mutation_pivot,
        cmap="YlOrRd",
        annot=True,
        fmt="d",
        linewidths=0.5
    )

    plt.title("Mutation Type Distribution in AKT Genes", fontsize=14)
    plt.tight_layout()
    plt.savefig(mutation_heatmap_path, dpi=300, bbox_inches='tight')

    # Display in notebook
    display(plt.gcf())
    plt.close()

    print(f"Mutation type heatmap saved to {mutation_heatmap_path}")

    # Create a sunburst chart for AKT1 mutations and diseases
    plt.figure(figsize=(10, 10))

    # Filter for AKT1 data
    akt1_data = disease_df[disease_df["Gene"] == "AKT1"]

    # Create a simple pie chart as an alternative to sunburst
    plt.pie(
        akt1_data["Count"],
        labels=akt1_data["Disease"], # type: ignore
        autopct='%1.1f%%',
        startangle=90,
        shadow=True,
        explode=[0.1 if "cancer" in disease else 0 for disease in akt1_data["Disease"]],
        colors=sns.color_palette("Set3", len(akt1_data))
    )

    plt.axis('equal')
    plt.title("AKT1 Disease Distribution", fontsize=16)

    akt1_disease_path = os.path.join(results_dir, "akt1_disease_distribution.png")
    plt.tight_layout()
    plt.savefig(akt1_disease_path, dpi=300, bbox_inches='tight')

    # Display in notebook
    display(plt.gcf())
    plt.close()

    print(f"AKT1 disease distribution plot saved to {akt1_disease_path}")

    return disease_df

def analyze_expression_data():
    """Analyze gene expression levels of AKT genes across tissues using real GTEx data"""

    expression_plot_path = os.path.join(results_dir, "akt_expression_heatmap.png")
    expression_boxplot_path = os.path.join(results_dir, "akt_expression_boxplot.png")
    tissue_expr_path = os.path.join(results_dir, "akt_tissue_expression.png")
    expression_data_path = os.path.join(results_dir, "akt_expression_data.csv")

    if os.path.exists(expression_data_path):
        print(f"Loading existing expression data from {expression_data_path}")
        expr_df = pd.read_csv(expression_data_path)
    else:
        print("Fetching real GTEx data for AKT genes...")

        # GTEx data portal API endpoint for querying gene expression
        base_url = "https://gtexportal.org/rest/v1/expression/medianGeneExpression"

        # Gene symbols we want to query
        genes = ["AKT1", "AKT2", "AKT3"]

        # tissue types to query
        tissues = [
            "Adipose Tissue", "Adrenal Gland", "Bladder", "Blood", "Blood Vessel",
            "Brain", "Breast", "Cervix Uteri", "Colon", "Esophagus", "Fallopian Tube",
            "Heart", "Kidney", "Liver", "Lung", "Muscle", "Nerve", "Ovary", "Pancreas",
            "Pituitary", "Prostate", "Skin", "Small Intestine", "Spleen", "Stomach",
            "Testis", "Thyroid", "Uterus", "Vagina"
        ]

        #  empty list to store our expression data
        expression_data = []

        # fetch expression data across tissues for each gene
        for gene in tqdm(genes, desc="Fetching gene data"):
            # Construct API call parameters
            params = {
                "format": "json",
                "geneId": gene
            }

            try:
                # Make API request
                response = requests.get(base_url, params=params)
                response.raise_for_status()  # Raise an exception for HTTP errors

                # Parse the response
                data = response.json()

                # Process the data - extract median TPM for each tissue
                for tissue_data in data["medianGeneExpression"]:
                    tissue_name = tissue_data["tissueSiteDetailId"]

                    # Extract the main tissue type
                    main_tissue = next((t for t in tissues if t in tissue_name), tissue_name)

                    # Get the expression value (median TPM)
                    expression = tissue_data["median"]

                    # Append to our data list
                    expression_data.append({
                        "Gene": gene,
                        "Tissue": main_tissue,
                        "Tissue_Detail": tissue_name,
                        "Expression": expression
                    })

                print(f"Successfully fetched data for {gene}")

            except requests.exceptions.RequestException as e:
                print(f"Error fetching data for {gene}: {e}")

                # If API fails, using backup data from literature !!!
                # These are approximate values based on published studies !!!
                backup_data = []

                if gene == "AKT1":
                    # AKT1 is broadly expressed with higher levels in blood and skin
                    for tissue in tissues:
                        base_expr = 20  # TPM

                        if tissue in ["Blood", "Skin"]:
                            expr = base_expr * 1.5
                        elif tissue in ["Brain", "Heart"]:
                            expr = base_expr * 0.8
                        else:
                            expr = base_expr

                        backup_data.append({
                            "Gene": gene,
                            "Tissue": tissue,
                            "Tissue_Detail": tissue,
                            "Expression": expr
                        })

                elif gene == "AKT2":
                    # AKT2 has higher expression in insulin-responsive tissues
                    for tissue in tissues:
                        base_expr = 15  # TPM

                        if tissue in ["Muscle", "Liver", "Adipose Tissue"]:
                            expr = base_expr * 2.0
                        elif tissue == "Pancreas":
                            expr = base_expr * 1.7
                        else:
                            expr = base_expr

                        backup_data.append({
                            "Gene": gene,
                            "Tissue": tissue,
                            "Tissue_Detail": tissue,
                            "Expression": expr
                        })

                elif gene == "AKT3":
                    # AKT3 has higher expression in neural tissues
                    for tissue in tissues:
                        base_expr = 10  # TPM

                        if tissue == "Brain":
                            expr = base_expr * 3.0
                        elif tissue in ["Heart", "Lung"]:
                            expr = base_expr * 1.5
                        else:
                            expr = base_expr

                        backup_data.append({
                            "Gene": gene,
                            "Tissue": tissue,
                            "Tissue_Detail": tissue,
                            "Expression": expr
                        })

                print(f"Using backup data for {gene}")
                expression_data.extend(backup_data)

        # Alternative data source if GTEx API failed
        if not expression_data:
            print("GTEx API failed, using alternative data source: Human Protein Atlas")

            # Human Protein Atlas data URLs for each gene
            hpa_urls = {
                "AKT1": "https://www.proteinatlas.org/ENSG00000142208-AKT1/tissue",
                "AKT2": "https://www.proteinatlas.org/ENSG00000105221-AKT2/tissue",
                "AKT3": "https://www.proteinatlas.org/ENSG00000117020-AKT3/tissue"
            }

            # Since direct API access to HPA might be restricted i useed literature-based values
            # These are normalized expression values (NX) based on HPA data

            hpa_data = {
                "AKT1": {
                    "Adipose Tissue": 21.8, "Adrenal Gland": 35.5, "Brain": 17.4, "Blood": 43.7,
                    "Blood Vessel": 22.4, "Bone Marrow": 44.2, "Breast": 23.5, "Colon": 21.8,
                    "Esophagus": 15.9, "Heart": 13.2, "Kidney": 20.6, "Liver": 17.9,
                    "Lung": 24.1, "Muscle": 14.7, "Nerve": 18.3, "Ovary": 25.2,
                    "Pancreas": 11.7, "Prostate": 21.5, "Skin": 31.9, "Small Intestine": 22.8,
                    "Spleen": 40.5, "Stomach": 14.3, "Testis": 22.8, "Thyroid": 28.2,
                    "Uterus": 19.7
                },
                "AKT2": {
                    "Adipose Tissue": 38.9, "Adrenal Gland": 22.6, "Brain": 5.4, "Blood": 21.3,
                    "Blood Vessel": 35.2, "Bone Marrow": 18.7, "Breast": 28.4, "Colon": 19.8,
                    "Esophagus": 18.9, "Heart": 22.7, "Kidney": 30.8, "Liver": 52.7,
                    "Lung": 16.2, "Muscle": 48.3, "Nerve": 9.4, "Ovary": 22.3,
                    "Pancreas": 35.7, "Prostate": 22.4, "Skin": 18.3, "Small Intestine": 21.4,
                    "Spleen": 28.4, "Stomach": 13.9, "Testis": 6.8, "Thyroid": 33.1,
                    "Uterus": 31.5
                },
                "AKT3": {
                    "Adipose Tissue": 9.4, "Adrenal Gland": 18.7, "Brain": 64.3, "Blood": 5.7,
                    "Blood Vessel": 12.4, "Bone Marrow": 4.6, "Breast": 15.2, "Colon": 6.9,
                    "Esophagus": 6.2, "Heart": 12.7, "Kidney": 7.2, "Liver": 3.8,
                    "Lung": 13.5, "Muscle": 7.8, "Nerve": 22.6, "Ovary": 18.9,
                    "Pancreas": 4.8, "Prostate": 12.3, "Skin": 18.6, "Small Intestine": 5.2,
                    "Spleen": 9.1, "Stomach": 5.7, "Testis": 21.8, "Thyroid": 17.3,
                    "Uterus": 15.8
                }
            }

            # Convert HPA data to our format
            for gene, tissue_data in hpa_data.items():
                for tissue, expression in tissue_data.items():
                    expression_data.append({
                        "Gene": gene,
                        "Tissue": tissue,
                        "Tissue_Detail": tissue,
                        "Expression": expression
                    })

            print("Using Human Protein Atlas data as fallback")

        # Convert to DataFrame
        expr_df = pd.DataFrame(expression_data)

        # If tissue-detail granularity is not needed, aggregate by main tissue and gene
        expr_df_aggregated = expr_df.groupby(['Gene', 'Tissue'])['Expression'].mean().reset_index()

        # Sort tissues by average expression for better visualization
        tissue_order = expr_df_aggregated.groupby('Tissue')['Expression'].mean().sort_values(ascending=False).index.tolist()
        expr_df_aggregated['Tissue'] = pd.Categorical(expr_df_aggregated['Tissue'], categories=tissue_order, ordered=True)

        # Sort and save to CSV
        expr_df_aggregated = expr_df_aggregated.sort_values(['Tissue', 'Gene'])
        expr_df_aggregated.to_csv(expression_data_path, index=False)
        print(f"Expression data saved to {expression_data_path}")

        # Use the aggregated dataframe for visualization
        expr_df = expr_df_aggregated

    # Create heatmap
    plt.figure(figsize=(16, 8))

    # Pivot the data for the heatmap
    heatmap_data = expr_df.pivot(index="Gene", columns="Tissue", values="Expression")

    # Create heatmap with custom colormap
    ax = sns.heatmap(heatmap_data, cmap="YlOrRd", annot=True, fmt=".1f",
                linewidths=0.5, cbar_kws={'label': 'Expression Level (TPM)'})

    # Improve readability of x-axis labels
    plt.xticks(rotation=45, ha='right')

    plt.title("AKT Gene Family Expression Across Human Tissues (GTEx Data)", fontsize=14)
    plt.tight_layout()
    plt.savefig(expression_plot_path, dpi=300, bbox_inches='tight')

    display(plt.gcf())
    plt.close()

    print(f"Expression heatmap saved to {expression_plot_path}")

    # Create boxplot to show expression distribution
    plt.figure(figsize=(10, 6))

    # Create boxplot
    ax = sns.boxplot(x="Gene", y="Expression", data=expr_df, palette="Set2")
    sns.swarmplot(x="Gene", y="Expression", data=expr_df, color="0.25", alpha=0.5)

    # Add statistics
    gene_stats = expr_df.groupby('Gene')['Expression'].agg(['mean', 'median', 'std']).reset_index()

    for i, row in gene_stats.iterrows():
        plt.text(i, expr_df['Expression'].max() * 0.9, # type: ignore
                f"Mean: {row['mean']:.1f}\nMedian: {row['median']:.1f}",
                ha='center', va='top', fontsize=9)

    plt.title("Distribution of AKT Gene Expression Across Human Tissues", fontsize=14)
    plt.ylabel("Expression Level (TPM)")
    plt.grid(axis='y', linestyle='--', alpha=0.7)

    plt.savefig(expression_boxplot_path, dpi=300, bbox_inches='tight')

    # Display in notebook
    display(plt.gcf())
    plt.close()

    print(f"Expression boxplot saved to {expression_boxplot_path}")

    # tissue-specific comparison plot for top tissues
    plt.figure(figsize=(16, 10))

    # Get top tissues by expression
    top_tissues = expr_df.groupby('Tissue')['Expression'].mean().nlargest(15).index.tolist()
    tissue_subset = expr_df[expr_df['Tissue'].isin(top_tissues)]

    # grouped bar chart
    ax = sns.barplot(x="Tissue", y="Expression", hue="Gene", data=tissue_subset, palette="Set1")

    plt.title("AKT Gene Expression in Top 15 Tissues", fontsize=14)
    plt.ylabel("Expression Level (TPM)")
    plt.xlabel("Tissue")
    plt.xticks(rotation=45, ha='right')
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.legend(title="Gene")

    plt.tight_layout()
    plt.savefig(tissue_expr_path, dpi=300, bbox_inches='tight')

    display(plt.gcf())
    plt.close()

    print(f"Tissue-specific expression plot saved to {tissue_expr_path}")

    # tissue-specific comparison line plot
    plt.figure(figsize=(16, 8))

    # Sort tissues by AKT1 expression for consistent ordering
    tissue_order = expr_df[expr_df['Gene'] == 'AKT1'].sort_values('Expression', ascending=False)['Tissue'].tolist()
    expr_df['Tissue'] = pd.Categorical(expr_df['Tissue'], categories=tissue_order, ordered=True)
    expr_df = expr_df.sort_values('Tissue')

    # Plot expression profiles as lines
    sns.lineplot(data=expr_df, x='Tissue', y='Expression', hue='Gene',
                 markers=True, dashes=False, palette='Set1', linewidth=2.5, markersize=8)

    plt.title("AKT Gene Expression Profiles Across Human Tissues", fontsize=14)
    plt.ylabel("Expression Level (TPM)")
    plt.xlabel("Tissue")
    plt.xticks(rotation=45, ha='right')
    plt.grid(linestyle='--', alpha=0.7)
    plt.legend(title="Gene", loc='best')

    profile_plot_path = os.path.join(results_dir, "akt_expression_profile.png")
    plt.tight_layout()
    plt.savefig(profile_plot_path, dpi=300, bbox_inches='tight')

    display(plt.gcf())
    plt.close()

    print(f"Expression profile plot saved to {profile_plot_path}")

    # Add summary statistics
    print("\nSummary Statistics for AKT Gene Expression:")
    summary = expr_df.groupby('Gene')['Expression'].agg(['mean', 'median', 'min', 'max', 'std'])
    print(summary)

    # Find tissues with highest expression for each gene
    highest_expr = pd.DataFrame(columns=['Gene', 'Highest_Tissue', 'Expression'])

    for gene in expr_df['Gene'].unique():
        gene_data = expr_df[expr_df['Gene'] == gene]
        max_tissue = gene_data.loc[gene_data['Expression'].idxmax()]
        highest_expr = pd.concat([highest_expr, pd.DataFrame({
            'Gene': [gene],
            'Highest_Tissue': [max_tissue['Tissue']],
            'Expression': [max_tissue['Expression']]
        })], ignore_index=True)

    print("\nTissues with highest expression for each AKT gene:")
    print(highest_expr)

    # Save these insights to a text file
    insights_path = os.path.join(results_dir, "akt_expression_insights.txt")

    with open(insights_path, 'w') as f:
        f.write("AKT GENE FAMILY EXPRESSION ANALYSIS\n")
        f.write("===================================\n\n")
        f.write("Summary Statistics:\n")
        f.write(str(summary) + "\n\n")
        f.write("Tissues with highest expression for each AKT gene:\n")
        f.write(str(highest_expr) + "\n\n")

        f.write("Gene-specific insights:\n")

        # AKT1 insights
        akt1_data = expr_df[expr_df['Gene'] == 'AKT1']
        akt1_high = akt1_data.nlargest(3, 'Expression')
        akt1_low = akt1_data.nsmallest(3, 'Expression')

        f.write("AKT1:\n")
        f.write(f"- Highest expression in: {', '.join(akt1_high['Tissue'].astype(str))}\n")
        f.write(f"- Lowest expression in: {', '.join(akt1_low['Tissue'].astype(str))}\n")
        f.write(f"- Mean expression: {akt1_data['Expression'].mean():.2f} TPM\n\n")

        # AKT2 insights
        akt2_data = expr_df[expr_df['Gene'] == 'AKT2']
        akt2_high = akt2_data.nlargest(3, 'Expression')
        akt2_low = akt2_data.nsmallest(3, 'Expression')

        f.write("AKT2:\n")
        f.write(f"- Highest expression in: {', '.join(akt2_high['Tissue'].astype(str))}\n")
        f.write(f"- Lowest expression in: {', '.join(akt2_low['Tissue'].astype(str))}\n")
        f.write(f"- Mean expression: {akt2_data['Expression'].mean():.2f} TPM\n\n")

        # AKT3 insights
        akt3_data = expr_df[expr_df['Gene'] == 'AKT3']
        akt3_high = akt3_data.nlargest(3, 'Expression')
        akt3_low = akt3_data.nsmallest(3, 'Expression')

        f.write("AKT3:\n")
        f.write(f"- Highest expression in: {', '.join(akt3_high['Tissue'].astype(str))}\n")
        f.write(f"- Lowest expression in: {', '.join(akt3_low['Tissue'].astype(str))}\n")
        f.write(f"- Mean expression: {akt3_data['Expression'].mean():.2f} TPM\n\n")

        # Comparative insights
        f.write("Comparative Insights:\n")

        # Identify tissues where AKT1 is dominant
        akt1_dominant = []
        for tissue in expr_df['Tissue'].unique():
            tissue_data = expr_df[expr_df['Tissue'] == tissue]
            if tissue_data[tissue_data['Gene'] == 'AKT1']['Expression'].values[0] > \
               max(tissue_data[tissue_data['Gene'] == 'AKT2']['Expression'].values[0],
                   tissue_data[tissue_data['Gene'] == 'AKT3']['Expression'].values[0]):
                akt1_dominant.append(tissue)

        f.write(f"- AKT1 is the dominant isoform in: {', '.join(akt1_dominant[:5])}{'...' if len(akt1_dominant) > 5 else ''}\n")

        # Identify tissues where AKT2 is dominant
        akt2_dominant = []
        for tissue in expr_df['Tissue'].unique():
            tissue_data = expr_df[expr_df['Tissue'] == tissue]
            if tissue_data[tissue_data['Gene'] == 'AKT2']['Expression'].values[0] > \
               max(tissue_data[tissue_data['Gene'] == 'AKT1']['Expression'].values[0],
                   tissue_data[tissue_data['Gene'] == 'AKT3']['Expression'].values[0]):
                akt2_dominant.append(tissue)

        f.write(f"- AKT2 is the dominant isoform in: {', '.join(akt2_dominant[:5])}{'...' if len(akt2_dominant) > 5 else ''}\n")

        # Identify tissues where AKT3 is dominant
        akt3_dominant = []
        for tissue in expr_df['Tissue'].unique():
            tissue_data = expr_df[expr_df['Tissue'] == tissue]
            if tissue_data[tissue_data['Gene'] == 'AKT3']['Expression'].values[0] > \
               max(tissue_data[tissue_data['Gene'] == 'AKT1']['Expression'].values[0],
                   tissue_data[tissue_data['Gene'] == 'AKT2']['Expression'].values[0]):
                akt3_dominant.append(tissue)

        f.write(f"- AKT3 is the dominant isoform in: {', '.join(akt3_dominant[:5])}{'...' if len(akt3_dominant) > 5 else ''}\n")

    print(f"Expression analysis insights saved to {insights_path}")
    return expr_df

def analyze_promoter_regions():
    """Analyze promoter regions and regulatory elements of AKT genes."""


    promoter_plot_path = os.path.join(results_dir, "akt_promoter_elements.png")
    promoter_data_path = os.path.join(results_dir, "akt_promoter_elements.csv")

    # Define upstream region (e.g., 1000 bp upstream of transcription start site)
    UPSTREAM_LENGTH = 1000

    # Fetch promoter regions from Ensembl
    promoter_data = []
    for gene, config in GENES.items():
        try:
            # Fetch sequence for upstream region
            handle = Entrez.efetch(db="nucleotide", id=config['dna_refseq'], rettype="fasta", retmode="text", seq_start=1, seq_stop=UPSTREAM_LENGTH)
            record = SeqIO.read(handle, "fasta")
            promoter_seq = str(record.seq)

            # Count GC content (a simple metric for promoter analysis)
            gc_content = (promoter_seq.count("G") + promoter_seq.count("C")) / len(promoter_seq)

            promoter_data.append({
                "Gene": gene,
                "Promoter_Sequence": promoter_seq,
                "GC_Content": gc_content
            })
            print(f"Fetched promoter region for {gene}")
        except Exception as e:
            print(f"Error fetching promoter region for {gene}: {str(e)}")

    # Convert to DataFrame
    promoter_df = pd.DataFrame(promoter_data)

    # Save to CSV
    promoter_df.to_csv(promoter_data_path, index=False)
    print(f"Promoter data saved to {promoter_data_path}")

    # Create visualization , GC content bar plot
    plt.figure(figsize=(8, 6))
    sns.barplot(x="Gene", y="GC_Content", data=promoter_df, palette="Set3")
    plt.title("GC Content in AKT Promoter Regions", fontsize=14)
    plt.xlabel("Gene")
    plt.ylabel("GC Content (%)")
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.savefig(promoter_plot_path, dpi=300, bbox_inches='tight')
    display(plt.gcf())
    plt.close()

    print(f"Promoter analysis plot saved to {promoter_plot_path}")

    return promoter_df
def analyze_splice_variants():
    """Analyze splice variants of AKT genes using Ensembl data."""

    splice_plot_path = os.path.join(results_dir, "akt_splice_variants.png")
    splice_data_path = os.path.join(results_dir, "akt_splice_variants.csv")

    # Fetch splice variant data from GenBank records
    splice_data = []
    for gene, config in GENES.items():
        try:
            handle = Entrez.efetch(db="nucleotide", id=config['dna_refseq'], rettype="gb", retmode="text")
            record = SeqIO.read(handle, "genbank")

            # Count exons in mRNA features
            for feature in record.features:
                if feature.type == "mRNA":
                    exons = [part for part in feature.location.parts]
                    splice_data.append({
                        "Gene": gene,
                        "Transcript": feature.qualifiers.get("product", ["Unknown"])[0],
                        "Exon_Count": len(exons),
                        "Transcript_Length": len(feature)
                    })
            print(f"Fetched splice variants for {gene}")
        except Exception as e:
            print(f"Error fetching splice variants for {gene}: {str(e)}")

    # Create DataFrame and handle empty case
    splice_df = pd.DataFrame(splice_data)

    if splice_df.empty:
        print("No splice variant data found. Skipping visualization.")
        return splice_df

    # Save to CSV
    splice_df.to_csv(splice_data_path, index=False)
    print(f"Splice variant data saved to {splice_data_path}")

    # Create visualization
    plt.figure(figsize=(10, 6))
    sns.boxplot(x="Gene", y="Exon_Count", data=splice_df, palette="Set2")
    plt.title("Exon Count Distribution in AKT Splice Variants")
    plt.xlabel("Gene")
    plt.ylabel("Number of Exons")
    plt.tight_layout()
    plt.savefig(splice_plot_path, dpi=300)
    plt.close()

    print(f"Splice variant plot saved to {splice_plot_path}")
    return splice_df



def analyze_pathways():
    """Analyze AKT genes in the PI3K-AKT-mTOR pathway.  [just for showoff or if you want to understand the pathway of the PI3K-AKT-mTOR pathway]"""

    # Fetch pathway data from KEGG
    pathway_id = "hsa04151"  # PI3K-AKT-mTOR pathway
    pathway_img_url = "https://www.kegg.jp/kegg/pathway/hsa/hsa04151.png"
    pathway_img_path = os.path.join(results_dir, "pi3k_akt_mtor_pathway.png")

    try:
        response = requests.get(pathway_img_url)
        response.raise_for_status()
        # Save image
        with open(pathway_img_path, "wb") as f:
            f.write(response.content)
        # Display
        display(Image(filename=pathway_img_path))
        print(f"Pathway visualization saved to {pathway_img_path}")

    except Exception as e:
        print(f"Pathway error: {str(e)}")

    # Build document properly

def run_pipeline():
    """Run the entire Project pipeline for AKT gene analysis."""
    try:
        print("🚀 Starting AKT gene analysis pipeline...")
        create_project_directories()

        # Step 1: Fetch genomic data
        print("\n1️⃣ Fetching nucleotide sequences...")
        fetch_genomic_data()

        # Step 2: Fetch PDB structures
        print("\n2️⃣ Fetching PDB structures...")
        fetch_pdb()

        # Step 3: Fetch ClinVar data (even if synthetic)
        print("\n3️⃣ Fetching ClinVar mutations...")
        clinvar_df = fetch_clinvar_data()
        print("\nClinVar Data Preview:")
        display(clinvar_df.head())  # Show first 5 entries

        # Step 4: Perform old alignment (ClustalO)
        print("\n4️⃣ Performing old sequence alignment (ClustalO)...")
        align_sequences()
        aligned_file_old = os.path.join(results_dir, "AKT_family_aligned.fasta")
        if os.path.exists(aligned_file_old):
            print("\nOld Alignment Preview:")
            alignment_old = AlignIO.read(aligned_file_old, "fasta")
            print(alignment_old)

        # Step 5: Perform improved MSA (ClustalOmega)
        print("\n5️⃣ Performing improved MSA (ClustalOmega)...")
        protein_alignment = perform_msa(GENES.keys(), seq_type='protein')
        if protein_alignment:
            print("\nImproved MSA Preview:")
            print(protein_alignment)

        # Step 6: Build phylogenetic tree
        print("\n6️⃣ Building phylogenetic tree...")
        build_phylogenetic_tree()

        # Step 7: Mutation analysis
        print("\n7️⃣ Analyzing mutations...")
        mutation_analysis()

        # Step 8: Map mutations to 3D structure
        print("\n8️⃣ Mapping mutations to 3D structure...")
        map_mutations_to_structure()

        # Step 9: Conservation analysis
        print("\n9️⃣ Calculating conservation...")
        calculate_conservation()

        # Step 10: Domain architecture
        print("\n🔟 Plotting domain architecture...")
        fetch_domain_architecture()

        # Step 11: Splice variant analysis
        print("\n1️⃣1️⃣ Analyzing splice variants...")
        analyze_splice_variants()

        # Step 12: Promoter region analysis
        print("\n1️⃣2️⃣ Analyzing promoter regions...")
        analyze_promoter_regions()

        # Step 13: Disease associations
        print("\n1️⃣3️⃣ Analyzing disease associations...")
        analyze_disease_associations()

        # Step 14: Expression data analysis
        print("\n1️⃣4️⃣ Analyzing expression data...")
        analyze_expression_data()

        # Step 15: Analyze pathways
        print("\n1️⃣5️⃣ Analyzing pathways...")
        analyze_pathways()

        print("\n✅ Pipeline completed! All visualizations have been saved.")
        print(f"📁 Results directory: {results_dir}")

        # Generate PDF report using ReportGenerator
        from report_generator import ReportGenerator  # Import ReportGenerator
        ReportGenerator().generate_report()
        print(" PDF Report Created")

    except Exception as e:
        print(f"❌ Pipeline failed with error: {str(e)}")

