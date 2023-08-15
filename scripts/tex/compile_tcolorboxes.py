import regex
import argparse
import os
import subprocess

# Define and parse command line arguments
parser = argparse.ArgumentParser(description="Extract tcolorboxes from a LaTeX source file.")
parser.add_argument("source_file", help="LaTeX source file to extract from.")
parser.add_argument("output_dir", help="Directory to write output files to.")
args = parser.parse_args()

# Read the content of the LaTeX source file
with open(args.source_file, "r") as source_file:
    source_content = source_file.read()

# Specify the preamble content
preamble = r"""
\documentclass{standalone}
\usepackage{amsmath, amsthm, amssymb, amsfonts}
\usepackage{tcolorbox}

\newtheorem{theorem}{Theorem}
\newtheorem{conjecture}{Conjecture}
\newtheorem{lemma}{Lemma}
\newtheorem{corollary}{Corollary}
\newtheorem{proposition}{Proposition}
\newtheorem{observation}{Observation}
\newtheorem{definition}{Definition}
\newtheorem{problem}{Problem}
\newtheorem{example}{Example}
\begin{document}
\begin{tcolorbox}
"""

# Specify the postamble content
postamble = "\\end{tcolorbox}\n\\end{document}"

# Use a regular expression to find each labeled tcolorbox
tcolorbox_regex = regex.compile(r"\\exportabletcolorbox\{(.*?)\}\{((?:\{(?:[^\{\}]*|(?2))*\}|[^\{\}])*)\}", regex.DOTALL)


print(list(tcolorbox_regex.finditer(source_content)))
# For each labeled tcolorbox...
for i, match in enumerate(tcolorbox_regex.finditer(source_content)):
    label, content = match.groups()

    # Create the content of a standalone document
    standalone_content = preamble + content + postamble

    # Write the standalone document to a .tex file in the output directory
    with open(os.path.join(args.output_dir, f"{label}.tex"), "w") as standalone_file:
        standalone_file.write(standalone_content)

    # Compile the standalone document to a PDF
    # This requires pdflatex to be installed and in your PATH
    subprocess.run(["pdflatex", "-output-directory=" + args.output_dir, os.path.join(args.output_dir, f"{label}.tex")])
