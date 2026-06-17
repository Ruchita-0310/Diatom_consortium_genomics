#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import re

# 1. Load Data
# Note: Using your specific filename 'metaeuk_output_polyp_taxonomy_tax_per_pred.tsv'
df_meta = pd.read_csv('metaeuk_output_polyp_taxonomy_tax_per_pred.tsv', sep='\t')
df_bins = pd.read_csv('contig_to_bin.txt', sep='\t', names=['bin_name', 'contig ID'])

# 2. Updated Organelle-Aware Labeling Logic
def get_organelle_label(classification):
    if pd.isna(classification) or str(classification).strip() == "": 
        return 'Unclassified'
    
    clean_class = str(classification).strip()
    
    # Eukaryotic check
    if 'd_Eukaryota' in clean_class:
        return 'Eukaryota'
    
    # Mitochondria check
    if 'o_Rickettsiales' in clean_class or 'o__Rickettsiales' in clean_class:
        return 'Mitochondria-derived'
    
    # Chloroplast check (Cyanobacteria proxy)
    if 'p_Cyanobacteria' in clean_class:
        return 'Chloroplast-derived'
    
    # Standard Root check
    if clean_class == "_cellular organisms":
        return 'Ambiguous (Cellular Org)'
    
    # Standard Bacteria check
    if 'd_Bacteria' in clean_class:
        return 'Bacteria'
    
    # Archaea/Viruses
    return 'Other'

df_meta['label'] = df_meta['Classification'].apply(get_organelle_label)

# 3. Aggregate hits per contig
contig_stats = df_meta.groupby('Contig_ID')['label'].value_counts().unstack(fill_value=0)

# 4. Final Decision Logic (Priority Hierarchy)
def final_decision_hierarchy(row):
    euk = row.get('Eukaryota', 0)
    mito = row.get('Mitochondria-derived', 0)
    chloro = row.get('Chloroplast-derived', 0)
    bac = row.get('Bacteria', 0)
    other = row.get('Other', 0)
    ambig = row.get('Ambiguous (Cellular Org)', 0)
    
    total_bio = euk + mito + chloro + bac + other + ambig
    
    if total_bio == 0:
        return 'Unclassified'

    # STEP A: The 30% Rule (HIGHEST PRIORITY)
    # If this is met, we don't care if Bacteria has more hits.
    euk_score = (euk + mito + chloro) / total_bio
    if euk_score > 0.30:
        return 'Eukaryota'
    
    # STEP B: Bacterial Majority
    # Must be strictly greater than Euk-group and Other.
    if bac > (euk + mito + chloro) and bac > other:
        return 'Bacteria'
    
    # STEP C: Other Majority
    if other > (euk + mito + chloro) and other > bac:
        return 'Other'
    
    # STEP D: TIE-BREAK / NO MAJORITY
    # If we reached here, it failed 30% Euk and no domain is a clear winner.
    return 'Ambiguous (Cellular Org)'

contig_stats['Final_Label'] = contig_stats.apply(final_decision_hierarchy, axis=1)

# 5. Merge and Normalize
final_df = pd.merge(df_bins, contig_stats.reset_index(), left_on='contig ID', right_on='Contig_ID', how='left')
final_df['Final_Label'] = final_df['Final_Label'].fillna('Unclassified')

# 6. Sorting and Export
def bin_sort_key(name):
    match = re.search(r'bin\.(\d+)', str(name))
    return (0, int(match.group(1))) if match else (1, str(name))

bin_summary = final_df.groupby(['bin_name', 'Final_Label']).size().unstack(fill_value=0)
bin_summary_norm = bin_summary.div(bin_summary.sum(axis=1), axis=0)
bin_summary_norm = bin_summary_norm.reindex(sorted(bin_summary_norm.index, key=bin_sort_key))

final_df.to_csv('contig_classification_final_priority.csv', index=False)

# 7. Plotting
color_map = {
    'Bacteria': "#2ac92a", 'Eukaryota': "#3197fd", 
    'Ambiguous (Cellular Org)': "#fa9734", 'Unclassified': '#999999', 'Other': "#db7a7a"
}
cat_order = ['Bacteria', 'Eukaryota', 'Ambiguous (Cellular Org)', 'Unclassified', 'Other']

for cat in cat_order:
    if cat not in bin_summary_norm.columns:
        bin_summary_norm[cat] = 0.0
bin_summary_norm = bin_summary_norm[cat_order]

plt.figure(figsize=(18, 8))
bin_summary_norm.plot(kind='bar', stacked=True, figsize=(18, 8), 
                      color=[color_map[c] for c in cat_order], 
                      edgecolor='white', linewidth=0.1)

plt.ylim(0, 1.0)
plt.title('Normalized Contig Distribution (Priority Hierarchy)')
plt.legend(title='Classification', bbox_to_anchor=(1.02, 1), loc='upper left')
plt.tight_layout()
plt.savefig('normalized_histogram_final.png')

print("Task complete. Created 'contig_classification_final_priority.csv'.")
