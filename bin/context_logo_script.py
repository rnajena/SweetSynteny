# (C) 2024, Tom Eulenfeld & Maria Schreiber, MIT license
import argparse
from collections import Counter
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

GAP_REPLACEMENT_COLOR = '#FFFFFF'

def setup_parser():
    """Configure command-line argument parser."""
    parser = argparse.ArgumentParser(description='Process color sequences and generate motif plot.')
    parser.add_argument('-i', '--input', required=True, help='Path to the input TSV file with color sequences')
    parser.add_argument('-o', '--output', required=True, help='Directory where output files (e.g. plots) will be saved')
    return parser

def _read_color_seqs(tsv):
    """Reads a Tab-Separated Values (TSV) file expected to contain ID and color sequence data and returns color_dic."""
    id_to_colors_dict = {}
    with open(tsv) as f:
        for line in f:
            if '#' not in line:
                continue
            # e.g. line: 'id\t#color1_details#color2_details#...'
            id_, color_str = line.split()
            color_str_lst = color_str.split("#")
            # Slice string into color codes every 7 characters (each hex color)\
            colors = ['#' + c.split("_")[0] for c in color_str_lst]

            id_to_colors_dict[id_] = colors[1:]
    # Keys are the IDs (str) and values are lists of hex color codes (str), e.g., {'sample_id': ['#RRGGBB', '#AABBCC', ...]}
    return id_to_colors_dict

def _color_mapping(colors, gapcolor='-'):
    """Creates two complementary dictionaries for mapping hex color codes to single-character symbols."""
    import string
    charfrom = {} # charfrom: Maps color code -> character symbol.
    colorfrom = {} # colorfrom: Maps character symbol -> color code.
    #for color, char in zip([gapcolor] + list(colors), '-' + string.ascii_letters[:len(colors)], strict=True):
    for color, char in zip([gapcolor] + list(colors), '-' + string.ascii_letters[:len(colors)]):
        charfrom[color] = char
        colorfrom[char] = color    
    return charfrom, colorfrom

def _write_mafftinput(id_to_colors_dict, out, charfrom):
    """Writes sequences from a dictionary of encoded color sequences to a file in FASTA format."""
    with open(out, 'w', encoding='latin1') as fout:
        for id_, colors in id_to_colors_dict.items():
            chars = ''.join(charfrom[color] for color in colors)
            fout.write(f'>{id_}\n{chars}\n')

def _read_mafftoutput(fname, colorfrom):
    """Reads a sequence alignment file in MAFFT output and converts the sequence characters into a list of corresponding color codes."""
    id_ = None
    ali = {}
    with open(fname) as f:
        for line in f:
            if line.startswith('>'):
                id_ = line.lstrip('>').strip()
            elif line.strip():
                assert id_ is not None
                colors = [colorfrom[char] for char in line.strip()]
                ali[id_] = colors
                id_ = None
    return ali

def textmafft(input_file):
    """Performs a Multiple Sequence Alignment on color sequences using the MAFFT program via os."""
    import os
    seqs = _read_color_seqs(input_file)
    all_colors = sorted(set(c for seq in seqs.values() for c in seq))
    charfrom, colorfrom = _color_mapping(all_colors, gapcolor='-' * len(all_colors[0]))
    _write_mafftinput(seqs, 'mafftin.fasta', charfrom)
    os.system('mafft --auto --text --op 1.53 mafftin.fasta > mafftout.fasta')
    ali = _read_mafftoutput('mafftout.fasta', colorfrom)
    return ali

def create_fig(args, all_colors, num_positions, pos_color_bitscores):
    """Generates and saves a stacked bar plot visualizing the "bit score" (information 
    content) of different color-coded features at each position of an aligned sequence."""
    # Prepare the bottom array to stack bar segments
    bottom = np.zeros(num_positions)

    fig, ax = plt.subplots(figsize=(6, 3))

    for color in all_colors:
        heights = [pos_color_bitscores[i].get(color, 0) for i in range(num_positions)]
        ax.bar(range(num_positions), heights, bottom=bottom, color=color, edgecolor='black', width=0.8)
        bottom += heights

    # Add gene_of_interest text above the middle bar
    middle_pos = num_positions // 2
    height = bottom[middle_pos]
    ax.text(middle_pos, height / 2, 'gene_of_interest', rotation=90, va='center', ha='center', fontsize=10, color='white')

    ax.set_xticks(range(num_positions))
    ax.set_xticklabels(range(num_positions))
    ax.set_xlabel('Position in Color List')
    ax.set_ylabel('Bit Score (Information Content)')
    ax.set_title('Stacked Bar Plot of Bit Scores by Position and Color')
    
    handles = [plt.Rectangle((0,0),1,1, color=color) for color in all_colors]
    ax.legend(handles, all_colors, title="Color/Feature", bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout(rect=[0, 0, 0.8, 1])

    plt.savefig(args.output + '/' + args.input.split('/')[-1].replace('.tsv','') + 'cluster_motif.png')
    plt.savefig(args.output + '/' + args.input.split('/')[-1].replace('.tsv','') + 'cluster_motif.svg')

def calculate_bit_scores(color_lists):
    """Calculates the information content (bit score) for each unique color at every 
    position across a list of aligned color sequences.    
    Assuming bit score per color at each position:
    1. Calculate the frequency p_i of each color i at a given position.
    2. Use the formula for information content per color:
        bit_score_i = p_i x log_2 p_i/q_i
    Where:
        p_i = frequency of color i at that position (number of times the color appears divided by total sequences)
        q_i = background frequency (assume uniform background = 1 / number_of_colors)
    The height of each stacked bar segment is the bit score of that color at that position."""
    num_positions = len(color_lists[0])
    num_seqs = len(color_lists)

    # Prepare frequency counts per position
    pos_color_counts = []
    for i in range(num_positions):
        colors_at_pos = [color_list[i] for color_list in color_lists]
        counts = Counter(colors_at_pos)
        pos_color_counts.append(counts)

    # Collect all unique colors appearing anywhere (for consistent color assignment)
    all_colors = sorted({color for counts in pos_color_counts for color in counts.keys()})

    # Assume uniform background frequency for each color
    bg_freq = 1.0 / len(all_colors)
    
    # Calculate bit scores for each position and color
    pos_color_bitscores = []
    for counts in pos_color_counts:
        total = sum(counts.values())
        bitscores = {}
        for color in all_colors:
            p = counts.get(color, 0) / total if total > 0 else 0
            if p > 0:
                bitscores[color] = p * np.log2(p / bg_freq)
            else:
                bitscores[color] = 0
        pos_color_bitscores.append(bitscores)

    return all_colors, num_positions, pos_color_bitscores

def main():
    parser = setup_parser()
    args = parser.parse_args()
    
    data = textmafft(args.input)
    
    # Replace '-------' with '#FFFFFF' in the lists inside the dict
    for key, color_list in data.items():
        data[key] = [GAP_REPLACEMENT_COLOR if color == '-------' else color for color in color_list]

    color_lists = list(data.values())
    
    all_colors, num_positions, pos_color_bitscores = calculate_bit_scores(color_lists)

    create_fig(args, all_colors, num_positions, pos_color_bitscores)

if __name__ == '__main__':
    main()
