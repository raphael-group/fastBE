import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import argparse

def main():
    parser = argparse.ArgumentParser(description='Plot false positive rate from a CSV file using boxplot.')
    parser.add_argument('csv_path', type=str, help='Path to the CSV file containing the DataFrame.')
    parser.add_argument('--output', type=str, help='Output prefix.')
    args = parser.parse_args()

    df = pd.read_csv(args.csv_path)
    df['clone_and_sample'] = '(' + df['clones'].astype(str) + ', ' + df['samples'].astype(str) + ')'
    df = df.sort_values(by='clone_and_sample', key=lambda x: x.map(eval))

    fig, axes = plt.subplots(ncols=2, figsize=(8, 5))
    sns.boxplot(x='clone_and_sample', y='elapsed_time', hue='iteration', data=df, ax=axes[0])
    sns.boxplot(x='clone_and_sample', y='objective_value', hue='iteration', data=df, ax=axes[1])

    for ax in axes:
        ax.set_xlabel('(Clone, Sample)')
        ax.set(yscale='log')
        for tick in ax.get_xticklabels():
            tick.set_rotation(60)

    axes[0].set_ylabel('Elapsed Time (s)')
    axes[1].set_ylabel('Objective Value')

    fig.tight_layout()
    plt.show()

if __name__ == '__main__':
    main()
