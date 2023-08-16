import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import argparse

def load_data(csv_path):
    df = pd.read_csv(csv_path)
    return df

def set_alpha(ax, alpha):
    for patch in ax.patches:
        r, g, b, _ = patch.get_facecolor()
        patch.set_facecolor((r, g, b, alpha))

def main():
    parser = argparse.ArgumentParser(description='Plot false positive rate from a CSV file using boxplot.')
    parser.add_argument('csv_path', type=str, help='Path to the CSV file containing the DataFrame.')
    parser.add_argument('--output', type=str, help='Output prefix.')
    args = parser.parse_args()

    df = load_data(args.csv_path)

    fig, axes = plt.subplots(figsize=(8, 4), nrows=1, ncols=2)

    sns.boxplot(data=df, x='clones', y='false_positive_rate', hue='algorithm', fliersize=0, ax=axes[0])
    set_alpha(axes[0], 0.5)

    sns.stripplot(data=df, x='clones', y='false_positive_rate', hue='algorithm', dodge=True, jitter=True, marker='o', alpha=0.5, legend=False, ax=axes[0])
    axes[0].set_xlabel('Number of Clones')
    axes[0].set_ylabel('False Positive Rate')
    axes[0].get_legend().remove()

    sns.boxplot(data=df, x='clones', y='false_negative_rate', hue='algorithm', fliersize=0, ax=axes[1])
    set_alpha(axes[1], 0.5)
    sns.stripplot(data=df, x='clones', y='false_negative_rate', hue='algorithm', dodge=True, jitter=True, marker='o', alpha=0.5, legend=False, ax=axes[1])
    axes[1].set_xlabel('Number of Clones')
    axes[1].set_ylabel('False Negative Rate')
    axes[1].legend(title='Algorithm')

    handles, labels = axes[1].get_legend_handles_labels()
    labels = list(map(lambda x: {'pairtree': 'Pairtree', 'allele_minima': 'AlleleMinima'}[x], labels))
    axes[1].legend(handles, labels, title='Algorithm', loc='upper right')
    fig.suptitle('FPR and FNR Versus Number of Clones')
    fig.tight_layout()

    fig.savefig(args.output + '_fprfnr_vs_clones.pdf')

    fig, axes = plt.subplots(figsize=(8, 4), nrows=1, ncols=2)

    sns.boxplot(data=df[df['clones'] == 50], x='samples', y='false_positive_rate', hue='algorithm', fliersize=0, ax=axes[0])
    set_alpha(axes[0], 0.5)
    sns.stripplot(data=df[df['clones'] == 50], x='samples', y='false_positive_rate', hue='algorithm', dodge=True, jitter=True, marker='o', alpha=0.5, legend=False, ax=axes[0])

    axes[0].set_xlabel('Number of Samples')
    axes[0].set_ylabel('False Positive Rate')
    axes[0].get_legend().remove()

    sns.boxplot(data=df[df['clones'] == 50], x='samples', y='false_negative_rate', hue='algorithm', fliersize=0, ax=axes[1])
    set_alpha(axes[1], 0.5)
    sns.stripplot(data=df[df['clones'] == 50], x='samples', y='false_negative_rate', hue='algorithm', dodge=True, jitter=True, marker='o', alpha=0.5, legend=False, ax=axes[1])

    axes[1].set_xlabel('Number of Samples')
    axes[1].set_ylabel('False Negative Rate')
    axes[1].legend(title='Algorithm')

    handles, labels = axes[1].get_legend_handles_labels()
    labels = list(map(lambda x: {'pairtree': 'Pairtree', 'allele_minima': 'AlleleMinima'}[x], labels))
    axes[1].legend(handles, labels, title='Algorithm', loc='upper right')
    fig.suptitle('FPR and FNR Versus Number of Samples w/ 50 Clones')
    
    fig.tight_layout()
    fig.savefig(args.output + '_50clones_fprfnr_vs_samples.pdf')

    fig, axes = plt.subplots(figsize=(10, 5), nrows=1, ncols=2)
    
    df['clone_and_sample'] = '(' + df['clones'].astype(str) + ', ' + df['samples'].astype(str) + ')'
    df = df.sort_values(by='clone_and_sample', key=lambda x: x.map(eval))

    sns.boxplot(data=df[df['clones'] >= 30], x='clone_and_sample', y='false_positive_rate', hue='algorithm', fliersize=0, ax=axes[0])
    set_alpha(axes[0], 0.5)
    sns.stripplot(
        data=df[df['clones'] >= 30], x='clone_and_sample', y='false_positive_rate', hue='algorithm', 
        dodge=True, jitter=True, marker='o', alpha=0.5, legend=False, ax=axes[0]
    )

    for tick in axes[0].get_xticklabels():
        tick.set_rotation(60)

    axes[0].set_xlabel('Number of Clones and Samples')
    axes[0].set_ylabel('False Positive Rate')
    axes[0].get_legend().remove()

    sns.boxplot(data=df[df['clones'] >= 30], x='clone_and_sample', y='false_negative_rate', hue='algorithm', fliersize=0, ax=axes[1])
    set_alpha(axes[1], 0.5)
    sns.stripplot(
        data=df[df['clones'] >= 30], x='clone_and_sample', y='false_negative_rate', hue='algorithm', 
        dodge=True, jitter=True, marker='o', alpha=0.5, legend=False, ax=axes[1]
    )

    for tick in axes[1].get_xticklabels():
        tick.set_rotation(60)

    axes[1].set_xlabel('Number of Clones and Samples')
    axes[1].set_ylabel('False Negative Rate')
    axes[1].legend(title='Algorithm')

    handles, labels = axes[1].get_legend_handles_labels()
    labels = list(map(lambda x: {'pairtree': 'Pairtree', 'allele_minima': 'AlleleMinima'}[x], labels))
    axes[1].legend(handles, labels, title='Algorithm', loc='upper right')
    fig.suptitle('FPR and FNR Versus Number of Clones and Samples w/ >= 30 Clones')

    fig.tight_layout()
    fig.savefig(args.output + '_fprfnr_vs_clonesandsamples.pdf')

    df['U_error'] = df['U_error'].astype(float) / (df['clones'] * df['samples'])
    df['F_error'] = df['F_error'].astype(float) / (df['clones'] * df['samples'])

    fig, axes = plt.subplots(figsize=(8, 4), nrows=1, ncols=2)

    sns.boxplot(data=df[df['clones'] == 50], x='samples', y='U_error', hue='algorithm', fliersize=0, ax=axes[0])
    set_alpha(axes[0], 0.5)
    sns.stripplot(data=df[df['clones'] == 50], x='samples', y='U_error', hue='algorithm', dodge=True, jitter=True, marker='o', alpha=0.5, legend=False, ax=axes[0])

    axes[0].set_xlabel('Number of Samples')
    axes[0].set_ylabel(r'$\frac{1}{mn}\|U - \hat{U}\|_1$')
    axes[0].get_legend().remove()

    sns.boxplot(data=df[df['clones'] == 50], x='samples', y='F_error', hue='algorithm', fliersize=0, ax=axes[1])
    set_alpha(axes[1], 0.5)
    sns.stripplot(data=df[df['clones'] == 50], x='samples', y='F_error', hue='algorithm', dodge=True, jitter=True, marker='o', alpha=0.5, legend=False, ax=axes[1])

    axes[1].set_xlabel('Number of Samples')
    axes[1].set_ylabel(r'$\frac{1}{mn}\|F - \hat{F}\|_1$')
    axes[1].legend(title='Algorithm')

    handles, labels = axes[1].get_legend_handles_labels()
    labels = list(map(lambda x: {'pairtree': 'Pairtree', 'allele_minima': 'AlleleMinima'}[x], labels))
    axes[1].legend(handles, labels, title='Algorithm', loc='upper right')
    fig.suptitle('Matrix Error Versus Number of Samples w/ 50 Clones')

    fig.tight_layout()
    fig.savefig(args.output + '_50clones_matrix_error_vs_samples.pdf')

    fig, axes = plt.subplots(figsize=(8, 4), nrows=1, ncols=2)
    sns.boxplot(data=df, x='clones', y='U_error', hue='algorithm', fliersize=0, ax=axes[0])
    set_alpha(axes[0], 0.5)
    sns.stripplot(data=df, x='clones', y='U_error', hue='algorithm', dodge=True, jitter=True, marker='o', alpha=0.5, legend=False, ax=axes[0])

    axes[0].set_xlabel('Number of Clones')
    axes[0].set_ylabel(r'$\frac{1}{mn}\|U - \hat{U}\|_1$')
    axes[0].get_legend().remove()

    sns.boxplot(data=df, x='clones', y='F_error', hue='algorithm', fliersize=0, ax=axes[1])
    set_alpha(axes[1], 0.5)
    sns.stripplot(data=df, x='clones', y='F_error', hue='algorithm', dodge=True, jitter=True, marker='o', alpha=0.5, legend=False, ax=axes[1])

    axes[1].set_xlabel('Number of Clones')
    axes[1].set_ylabel(r'$\frac{1}{mn}\|F - \hat{F}\|_1$')
    axes[1].legend(title='Algorithm')

    handles, labels = axes[1].get_legend_handles_labels()
    labels = list(map(lambda x: {'pairtree': 'Pairtree', 'allele_minima': 'AlleleMinima'}[x], labels))
    axes[1].legend(handles, labels, title='Algorithm', loc='upper right')
    fig.suptitle('Matrix Error Versus Number of Clones')
    
    fig.tight_layout()
    fig.savefig(args.output + '_matrix_error_vs_clones.pdf')

    fig, axes = plt.subplots(figsize=(10, 5), nrows=1, ncols=2)

    sns.boxplot(data=df[df['clones'] >= 30], x='clone_and_sample', y='U_error', hue='algorithm', fliersize=0, ax=axes[0])
    set_alpha(axes[0], 0.5)
    sns.stripplot(
        data=df[df['clones'] >= 30], x='clone_and_sample', y='U_error', hue='algorithm', 
        dodge=True, jitter=True, marker='o', alpha=0.5, legend=False, ax=axes[0]
    )

    for tick in axes[0].get_xticklabels():
        tick.set_rotation(60)

    axes[0].set_xlabel('Number of Clones and Samples')
    axes[0].set_ylabel(r'$\frac{1}{mn}\|U - \hat{U}\|_1$')
    axes[0].get_legend().remove()

    sns.boxplot(data=df[df['clones'] >= 30], x='clone_and_sample', y='F_error', hue='algorithm', fliersize=0, ax=axes[1])
    set_alpha(axes[1], 0.5)
    sns.stripplot(
        data=df[df['clones'] >= 30], x='clone_and_sample', y='F_error', hue='algorithm', 
        dodge=True, jitter=True, marker='o', alpha=0.5, legend=False, ax=axes[1]
    )

    for tick in axes[1].get_xticklabels():
        tick.set_rotation(60)

    axes[1].set_xlabel('Number of Clones and Samples')
    axes[1].set_ylabel(r'$\frac{1}{mn}\|F - \hat{F}\|_1$')
    axes[1].legend(title='Algorithm')

    handles, labels = axes[1].get_legend_handles_labels()
    labels = list(map(lambda x: {'pairtree': 'Pairtree', 'allele_minima': 'AlleleMinima'}[x], labels))
    axes[1].legend(handles, labels, title='Algorithm', loc='upper right')
    fig.suptitle('Matrix Error Versus Number of Clones and Samples w/ >= 30 Clones')

    fig.tight_layout()
    fig.savefig(args.output + '_matrix_error_vs_clonesandsamples.pdf')

    plt.show()

if __name__ == '__main__':
    main()

