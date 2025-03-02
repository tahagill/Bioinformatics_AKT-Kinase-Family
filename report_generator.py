import os
from typing import List
from reportlab.lib.pagesizes import letter
from reportlab.platypus import (
    BaseDocTemplate, Paragraph, Spacer, Image, 
    Table, TableStyle, PageTemplate, Frame, Flowable
)
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle, StyleSheet1
from reportlab.lib.units import inch
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib import colors

# Define directory paths
BASE_DIR = os.path.join(os.getcwd(), "kinase_project")
sequence_dir = os.path.join(BASE_DIR, "data/sequences")
mutation_dir = os.path.join(BASE_DIR, "data/mutations")
results_dir = os.path.join(BASE_DIR, "results")

# Register professional fonts
pdfmetrics.registerFont(TTFont("Roboto", "Roboto-Regular.ttf"))
pdfmetrics.registerFont(TTFont("Roboto-Bold", "Roboto-Bold.ttf"))

class ReportGenerator:
    def __init__(self, results_dir: str = results_dir):
        self.results_dir = results_dir
        self.styles = self._create_styles()
        self.elements: List[Flowable] = []

    def _create_styles(self) -> StyleSheet1:
        """Create custom styles without redefining existing ones"""
        styles = getSampleStyleSheet()
        
        # Custom title style
        styles.add(ParagraphStyle(
            name="BioTitle",
            fontName="Roboto-Bold",
            fontSize=22,
            leading=26,
            spaceAfter=18,
            alignment=1,
            textColor=colors.darkblue
        ))
        
        # Custom section header style
        styles.add(ParagraphStyle(
            name="BioSectionHeader",
            fontName="Roboto-Bold",
            fontSize=16,
            leading=20,
            spaceAfter=12,
            textColor=colors.darkcyan
        ))
        
        # Modify existing BodyText style
        body_style = styles["BodyText"]
        body_style.fontName = "Roboto"
        body_style.fontSize = 12
        body_style.leading = 14
        body_style.spaceAfter = 6
        
        return styles

    def _add_title(self, title: str) -> None:
        self.elements.append(Paragraph(title, self.styles["BioTitle"]))
        self.elements.append(Spacer(1, 12))

    def _add_section_header(self, header: str) -> None:
        self.elements.append(Paragraph(header, self.styles["BioSectionHeader"]))
        self.elements.append(Spacer(1, 12))

    def _add_image(self, image_path: str, width: float = 6*inch) -> None:
        if os.path.exists(image_path):
            try:
                img = Image(image_path, width=width, height=4*inch)
                self.elements.append(img)
                self.elements.append(Spacer(1, 12))
            except Exception as e:
                print(f"⚠️ Could not load image {image_path}: {str(e)}")
        else:
            print(f"⚠️ Image not found: {image_path}")

    def _add_text(self, text: str) -> None:
        self.elements.append(Paragraph(text, self.styles["BodyText"]))
        self.elements.append(Spacer(1, 6))

    def _add_table(self, data: List[List[str]], headers: List[str]) -> None:
        """Add table with type-safe data"""
        try:
            # Convert all data to strings
            str_data = [[str(item) for item in row] for row in data]
            str_headers = [str(h) for h in headers]
            
            table = Table([str_headers] + str_data)
            table.setStyle(TableStyle([
                ('BACKGROUND', (0,0), (-1,0), colors.darkcyan),
                ('TEXTCOLOR', (0,0), (-1,0), colors.whitesmoke),
                ('ALIGN', (0,0), (-1,-1), 'CENTER'),
                ('FONTNAME', (0,0), (-1,0), 'Roboto-Bold'),
                ('FONTSIZE', (0,0), (-1,0), 12),
                ('BOTTOMPADDING', (0,0), (-1,0), 12),
                ('BACKGROUND', (0,1), (-1,-1), colors.beige),
                ('GRID', (0,0), (-1,-1), 1, colors.black)
            ]))
            self.elements.append(table)
            self.elements.append(Spacer(1, 12))
        except Exception as e:
            print(f"⚠️ Failed to create table: {str(e)}")

    def _add_fasta_previews(self):
        """Add DNA and protein sequence previews"""
        # DNA Sequences
        self._add_section_header("1. DNA Sequence Previews")
        for gene in ["AKT1", "AKT2", "AKT3"]:
            dna_path = os.path.join(sequence_dir, f"{gene}_sequence.fasta")
            if os.path.exists(dna_path):
                try:
                    with open(dna_path, "r") as f:
                        preview = "".join([next(f) for _ in range(5)])
                    self._add_text(f"<b>{gene} DNA Sequence Preview:</b>")
                    self._add_text(preview)
                except Exception as e:
                    print(f"⚠️ Could not read DNA sequence for {gene}: {str(e)}")
            else:
                print(f"⚠️ DNA sequence file not found: {dna_path}")
        
        # Protein Sequences
        self._add_section_header("2. Protein Sequence Previews")
        for gene in ["AKT1", "AKT2", "AKT3"]:
            protein_path = os.path.join(sequence_dir, "protein", f"{gene}_protein.fasta")
            if os.path.exists(protein_path):
                try:
                    with open(protein_path, "r") as f:
                        preview = "".join([next(f) for _ in range(5)])
                    self._add_text(f"<b>{gene} Protein Sequence Preview:</b>")
                    self._add_text(preview)
                except Exception as e:
                    print(f"⚠️ Could not read protein sequence for {gene}: {str(e)}")
            else:
                print(f"⚠️ Protein sequence file not found: {protein_path}")

    def _add_alignment_text(self) -> None:
        """Add sequence alignment text from file"""
        alignment_path = os.path.join(results_dir, "AKT_family_aligned.fasta")
        if os.path.exists(alignment_path):
            try:
                with open(alignment_path, "r") as f:
                    alignment = f.read()
                self._add_section_header("3. Multiple Sequence Alignment")
                self._add_text(alignment)
            except Exception as e:
                print(f"⚠️ Could not read alignment file: {str(e)}")
        else:
            print(f"⚠️ Alignment file not found: {alignment_path}")



    def _add_full_analysis(self):
        """Add all analysis components"""
        # Existing analysis
        self._add_section_header("5. Phylogenetic Analysis")
        self._add_image(os.path.join(results_dir, "phylogenetic_tree.png"))
        
        # Disease analysis
        self._add_section_header("6. Disease Analysis")
        self._add_image(os.path.join(results_dir, "akt1_disease_distribution.png"))
        self._add_image(os.path.join(results_dir, "akt_disease_associations.png"))
        self._add_image(os.path.join(results_dir, "akt_mutation_heatmap.png"))
        
        # Expression analysis
        self._add_section_header("7. Expression Analysis")
        self._add_image(os.path.join(results_dir, "akt_expression_heatmap.png"))
        self._add_image(os.path.join(results_dir, "akt_expression_boxplot.png"))
        self._add_image(os.path.join(results_dir, "akt_tissue_expression.png"))
        self._add_image(os.path.join(results_dir, "akt_expression_profile.png"))
        
        # Structural features
        self._add_section_header("8. Structural Features")
        for gene in ["AKT1", "AKT2", "AKT3"]:
            self._add_image(os.path.join(results_dir, f"{gene}_domains.png"))
        self._add_image(os.path.join(results_dir, "akt_promoter_elements.png"))

    def generate_report(self, output_file: str = "AKT_Report.pdf") -> None:
        """Generate full PDF report"""
        doc = BaseDocTemplate(
            output_file,
            pagesize=letter,
            leftMargin=0.5*inch,
            rightMargin=0.5*inch,
            topMargin=0.5*inch,
            bottomMargin=0.5*inch
        )
        
        frame = Frame(
            doc.leftMargin, doc.bottomMargin,
            doc.width, doc.height,
            id='normal'
        )
        doc.addPageTemplates([PageTemplate(id='AllPages', frames=[frame])])

        try:
            self._add_title("AKT Kinase Family Analysis Report \n Author: Taha Ahmad")
            self._add_fasta_previews()
            self._add_alignment_text()

            self._add_full_analysis()
            
            doc.build(self.elements)
            print(f"✅ Successfully generated report: {output_file}")
            
        except PermissionError:
            alt_path = os.path.join(os.path.expanduser("~"), "Desktop", output_file)
            doc = BaseDocTemplate(alt_path, pagesize=letter)
            doc.build(self.elements)
            print(f"✅ Report saved to desktop: {alt_path}")
        except Exception as e:
            print(f"❌ Failed to generate report: {str(e)}")
            print("Check the error messages above for missing files or other issues.")

