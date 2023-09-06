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

algorithm_name_map = {
    'pairtree': 'Pairtree',
    'allele_minima': 'AlleleMinima',
    'calder': 'CALDER',
    'citup': 'CITUP'
}

def plot_matrix_error(df, output=None, xaxis='samples', xaxislabel='Number of Samples', rotate=False):
    fig, axes = plt.subplots(figsize=(8, 4), nrows=1, ncols=2)

    sns.boxplot(data=df, x=xaxis, y='U_error', hue='algorithm', fliersize=0, ax=axes[0])
    set_alpha(axes[0], 0.5)
    sns.stripplot(data=df, x=xaxis, y='U_error', hue='algorithm', dodge=True, jitter=True, marker='o', alpha=0.5, legend=False, ax=axes[0])

    axes[0].set_xlabel(xaxislabel)
    axes[0].set_ylabel(r'$\frac{1}{mn}\|U - \hat{U}\|_1$')
    axes[0].get_legend().remove()
    
    sns.boxplot(data=df, x=xaxis, y='F_error', hue='algorithm', fliersize=0, ax=axes[1])
    set_alpha(axes[1], 0.5)
    sns.stripplot(data=df, x=xaxis, y='F_error', hue='algorithm', dodge=True, jitter=True, marker='o', alpha=0.5, legend=False, ax=axes[1])

    axes[1].set_xlabel(xaxislabel)
    axes[1].set_ylabel(r'$\frac{1}{mn}\|F - \hat{F}\|_1$')
    axes[1].legend(title='Algorithm')

    handles, labels = axes[1].get_legend_handles_labels()
    labels = list(map(lambda x: algorithm_name_map[x], labels))
    axes[1].legend(handles, labels, title='Algorithm', loc='upper right')

    if rotate:
        for tick in axes[0].get_xticklabels():
            tick.set_rotation(60)

        for tick in axes[1].get_xticklabels():
            tick.set_rotation(60)

    fig.tight_layout()
    if output:
        fig.savefig(output)

def plot_fpr_fnr(df, output=None, xaxis='clones', xaxislabel='Number of Clones', rotate=False):
    fig, axes = plt.subplots(figsize=(8, 4), nrows=1, ncols=2)

    sns.boxplot(data=df, x=xaxis, y='false_positive_rate', hue='algorithm', fliersize=0, ax=axes[0])
    set_alpha(axes[0], 0.5)

    sns.stripplot(data=df, x=xaxis, y='false_positive_rate', hue='algorithm', dodge=True, jitter=True, marker='o', alpha=0.5, legend=False, ax=axes[0])
    axes[0].set_xlabel(xaxislabel)
    axes[0].set_ylabel('False Positive Rate')
    axes[0].get_legend().remove()

    sns.boxplot(data=df, x=xaxis, y='false_negative_rate', hue='algorithm', fliersize=0, ax=axes[1])
    set_alpha(axes[1], 0.5)
    sns.stripplot(data=df, x=xaxis, y='false_negative_rate', hue='algorithm', dodge=True, jitter=True, marker='o', alpha=0.5, legend=False, ax=axes[1])
    axes[1].set_xlabel(xaxislabel)
    axes[1].set_ylabel('False Negative Rate')
    axes[1].legend(title='Algorithm')

    handles, labels = axes[1].get_legend_handles_labels()
    labels = list(map(lambda x: algorithm_name_map[x], labels))
    axes[1].legend(handles, labels, title='Algorithm', loc='upper right')

    if rotate:
        for tick in axes[0].get_xticklabels():
            tick.set_rotation(60)

        for tick in axes[1].get_xticklabels():
            tick.set_rotation(60)

    fig.tight_layout()
    if output:
        fig.savefig(output)

def main():
    parser = argparse.ArgumentParser(description='Plot false positive rate from a CSV file using boxplot.')
    parser.add_argument('csv_path', type=str, help='Path to the CSV file containing the DataFrame.')
    parser.add_argument('--output', type=str, help='Output prefix.')
    args = parser.parse_args()

    df = load_data(args.csv_path)

    df['U_error'] = df['U_error'].astype(float) / (df['clones'] * df['samples'])
    df['F_error'] = df['F_error'].astype(float) / (df['clones'] * df['samples'])

    df['clone_and_sample'] = '(' + df['clones'].astype(str) + ', ' + df['samples'].astype(str) + ')'
    df = df.sort_values(by='clone_and_sample', key=lambda x: x.map(eval))
    
    plot_fpr_fnr(
        df[~df['algorithm'].isin(['calder', 'citup'])], 
        output=args.output + 'large_fpr_fnr.pdf'
    )

    plot_fpr_fnr(
        df[df['clones'] <= 10], 
        output=args.output + '_10clones_fpr_fnr.pdf'
    )

    plot_matrix_error(
        df[~df['algorithm'].isin(['calder', 'citup'])], 
        output=args.output + 'large_matrix_error.pdf', 
        xaxis='clone_and_sample', 
        xaxislabel='(Number of Clones, Number of Samples)', 
        rotate=True
    )

    plot_matrix_error(
        df[df['clones'] <= 10], 
        output=args.output + '_10clones_matrix_error.pdf', 
        xaxis='clone_and_sample', 
        xaxislabel='(Number of Clones, Number of Samples)', 
        rotate=True
    )

    plot_matrix_error(
        df[df['clones'] == 50], 
        output=args.output + '_50samples_matrix_error.pdf', 
        xaxis='samples', 
        xaxislabel='Number of Samples'
    )

    plt.show()

if __name__ == '__main__':
    main()

