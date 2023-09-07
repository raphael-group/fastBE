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

# use seaborn color pallete
algorithm_color_map = {
    'pairtree': sns.color_palette()[0],
    'allele_minima': sns.color_palette()[1],
    'calder': sns.color_palette()[2],
    'citup': sns.color_palette()[3]
}

def plot_matrix_error(df, output=None, xaxis='samples', xaxislabel='Number of Samples', rotate=False):
    fig, axes = plt.subplots(figsize=(8, 4), nrows=1, ncols=2)

    sns.boxplot(data=df, x=xaxis, y='U_error', hue='algorithm', fliersize=0, ax=axes[0], palette=algorithm_color_map)
    set_alpha(axes[0], 0.5)
    sns.stripplot(
        data=df, x=xaxis, y='U_error', hue='algorithm', dodge=True, jitter=True, marker='o', alpha=0.5, 
        legend=False, ax=axes[0], palette=algorithm_color_map
    )

    axes[0].set_xlabel(xaxislabel)
    axes[0].set_ylabel(r'$\frac{1}{mn}\|U - \hat{U}\|_1$')
    axes[0].get_legend().remove()
    
    sns.boxplot(data=df, x=xaxis, y='F_error', hue='algorithm', fliersize=0, ax=axes[1], palette=algorithm_color_map)
    set_alpha(axes[1], 0.5)
    sns.stripplot(
        data=df, x=xaxis, y='F_error', hue='algorithm', dodge=True, jitter=True, marker='o', alpha=0.5, legend=False, ax=axes[1],
        palette=algorithm_color_map
    )

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

    print(df.groupby(['algorithm', xaxis])[['U_error', 'F_error']].mean())

    fig.tight_layout()
    if output:
        fig.savefig(output, transparent=True)

def plot_fpr_fnr(df, output=None, xaxis='clones', xaxislabel='Number of Clones', rotate=False):
    fig, axes = plt.subplots(figsize=(8, 4), nrows=1, ncols=2)

    sns.boxplot(data=df, x=xaxis, y='false_positive_rate', hue='algorithm', fliersize=0, ax=axes[0], palette=algorithm_color_map)
    set_alpha(axes[0], 0.5)

    sns.stripplot(data=df, x=xaxis, y='false_positive_rate', hue='algorithm', dodge=True, jitter=True, marker='o', alpha=0.5, legend=False, ax=axes[0], palette=algorithm_color_map)
    axes[0].set_xlabel(xaxislabel)
    axes[0].set_ylabel('False Positive Rate')
    axes[0].get_legend().remove()

    sns.boxplot(data=df, x=xaxis, y='false_negative_rate', hue='algorithm', fliersize=0, ax=axes[1], palette=algorithm_color_map)
    set_alpha(axes[1], 0.5)
    sns.stripplot(data=df, x=xaxis, y='false_negative_rate', hue='algorithm', dodge=True, jitter=True, marker='o', alpha=0.5, legend=False, ax=axes[1], palette=algorithm_color_map)
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

    print(df.groupby(['algorithm', xaxis])[['false_positive_rate', 'false_negative_rate']].mean())

    fig.tight_layout()
    if output:
        fig.savefig(output, transparent=True)

def main():
    parser = argparse.ArgumentParser(description='Plot false positive rate from a CSV file using boxplot.')
    parser.add_argument('csv_path', type=str, help='Path to the CSV file containing the DataFrame.')
    parser.add_argument('--output', type=str, help='Output prefix.')
    args = parser.parse_args()

    nmutations = [20, 30, 40, 50]
    nclones    = [3, 5, 10]
    nsamples   = [5, 10, 25]
    seeds      = [0, 1, 2, 3, 4, 5]
    coverage   = [100]

    df = load_data(args.csv_path)
    df = df[
        (df['coverage'].isin(coverage) & 
        df['mutations'].isin(nmutations) & 
        df['clones'].isin(nclones) & 
        df['samples'].isin(nsamples) & 
        df['seed'].isin(seeds)) | 
        (df['clones'] > 10)
    ]

    df['U_error'] = df['U_error'].astype(float) / (df['clones'] * df['samples'])
    df['F_error'] = df['F_error'].astype(float) / (df['clones'] * df['samples'])
    df['true_positives'] = df['positives'] - df['false_negatives']
    df['true_negatives'] = df['negatives'] - df['false_positives']

    df['precision'] = df['true_positives'] / (df['true_positives'] + df['false_positives'])
    df['recall'] = df['true_positives'] / (df['true_positives'] + df['false_negatives'])
    df['f1_score'] = (2 * df['precision'] * df['recall']) / (df['precision'] + df['recall'])

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

    # plot_matrix_error(
        # df[~df['algorithm'].isin(['calder', 'citup'])], 
        # output=args.output + 'large_matrix_error.pdf', 
        # xaxis='clone_and_sample', 
        # xaxislabel='(Number of Clones, Number of Samples)', 
        # rotate=True
    # )

    plot_matrix_error(
        df[df['clones'] <= 10], 
        output=args.output + '_10clones_matrix_error.pdf', 
        xaxis='clones',
        xaxislabel='Number of Clones'
    )

    plot_matrix_error(
        df[~df['algorithm'].isin(['calder', 'citup'])],
        output=args.output + 'large_matrix_error.pdf',
        xaxis='clones',
        xaxislabel='Number of Clones'
    )
    # plot_matrix_error(
        # df[df['clones'] == 50], 
        # output=args.output + '_50samples_matrix_error.pdf', 
        # xaxis='samples', 
        # xaxislabel='Number of Samples'
    # )

    fig, axes = plt.subplots(figsize=(8, 3.5), nrows=1, ncols=2)
    axes = axes[::-1]
    sns.boxplot(data=df[df['algorithm'].isin(['allele_minima', 'pairtree']) & (df['clones'] > 10)], x='clones', y='f1_score', hue='algorithm', fliersize=0, ax=axes[0], palette=algorithm_color_map)
    set_alpha(axes[0], 0.5)

    sns.stripplot(data=df[df['algorithm'].isin(['allele_minima', 'pairtree']) & (df['clones'] > 10)], x='clones', y='f1_score', hue='algorithm', dodge=True, jitter=True, marker='o', alpha=0.5, legend=False, ax=axes[0], palette=algorithm_color_map)
    axes[0].set_xlabel('Number of Clones')
    axes[0].set_ylabel('F1 Score')

    sns.boxplot(data=df[df['clones'] <= 10], x='clones', y='f1_score', hue='algorithm', fliersize=0, ax=axes[1], palette=algorithm_color_map)
    set_alpha(axes[1], 0.5)
    sns.stripplot(data=df[df['clones'] <= 10], x='clones', y='f1_score', hue='algorithm', dodge=True, jitter=True, marker='o', alpha=0.5, legend=False, ax=axes[1], palette=algorithm_color_map)
    axes[1].set_xlabel('Number of Clones')
    axes[1].set_ylabel('F1 Score')
    axes[1].legend(title='Algorithm')
    axes[1].get_legend().remove()


    handles, labels = axes[1].get_legend_handles_labels()
    labels = list(map(lambda x: algorithm_name_map[x], labels))
    axes[0].legend(handles, labels, title='Algorithm', loc='upper right')

    fig.tight_layout()
    fig.savefig(args.output + '_f1_score.pdf', transparent=True)

    print(df[df['clones'] <= 10].groupby(['algorithm', 'clones'])['f1_score'].mean())
    print(df[df['clones'] <= 10].algorithm.value_counts())
    print(df[~df['algorithm'].isin(['calder', 'citup'])].algorithm.value_counts())

    plt.show()

if __name__ == '__main__':
    main()

