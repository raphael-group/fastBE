import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plots timing results.')
    parser.add_argument('cpp_results', type=str, help='The path to the cpp timing results CSV.')
    parser.add_argument('python_results', type=str, help='The path to the python timing results CSV.')
    parser.add_argument('--output', type=str, help='Plot output prefix.')
    
    args = parser.parse_args()

    cpp_df = pd.read_csv(args.cpp_results)
    cpp_df = cpp_df.rename(columns={"time (s)": "cpp_dp_time"})

    python_df = pd.read_csv(args.python_results)
    python_timing_df = python_df[[
        'm', 'n', 's', 'c', 'r', 'lp_time_without_building', 'lp_time_with_building', 
        'dual_time_without_building', 'dual_time_with_building', 'dp_time',
    ]]

    df = pd.merge(cpp_df, python_timing_df, on=['m', 'n', 's', 'c', 'r'])
    df = df.melt(
        id_vars=['m', 'n', 's', 'c', 'r'], 
        value_vars=['cpp_dp_time', 'lp_time_without_building', 'lp_time_with_building', 'dual_time_without_building', 'dual_time_with_building', 'dp_time'], 
        var_name='algorithm',  value_name='time (s)'
    )

    df = df[df['algorithm'].isin(['lp_time_without_building', 'cpp_dp_time'])]

    df['algorithm'] = df['algorithm'].replace({
        'cpp_dp_time': 'Handcrafted DP',
        'dual_time_without_building': 'Dual LP (Gurobi)',
        'lp_time_without_building': 'LP (Gurobi)',
    })

    fig, ax = plt.subplots(figsize=(7, 5))

    df['(n, s)'] = df.apply(lambda row: (row["n"], row["s"]), axis=1)

    col_order = df.apply(lambda row: (row["n"], row["s"]), axis=1).sort_values().unique()

    sns.boxplot(data=df, x='(n, s)', y='time (s)', hue='algorithm', ax=ax, order=col_order)

    ax.set_xlabel('(# Clones, # Samples)')
    ax.set_ylabel('Time (s)')
    ax.legend(title='Algorithm', loc='upper right')

    for tick in ax.get_xticklabels():
        tick.set_rotation(60)

    plt.tight_layout()

    if args.output is not None:
        plt.savefig(f'{args.output}_vs_parameters.pdf')

    fig, axes = plt.subplots(ncols=2, figsize=(10, 5))

    speedup_df = python_timing_df.merge(cpp_df, on=['m', 'n', 's', 'c', 'r'])
    speedup_df['speedup_no_build'] = speedup_df['lp_time_without_building'] / speedup_df['cpp_dp_time'] 
    speedup_df['speedup_with_build'] = speedup_df['lp_time_with_building'] / speedup_df['cpp_dp_time']
    speedup_df['(n, s)'] = speedup_df.apply(lambda row: f'({row["n"]}, {row["s"]})', axis=1)

    sns.histplot(data=speedup_df['speedup_with_build'], ax=axes[0])
    sns.histplot(data=speedup_df['speedup_no_build'], ax=axes[1])
    
    axes[0].set_title('Relative Speedup With Model Building Time')
    axes[1].set_title('Relative Speedup Without Model Building Time')

    axes[0].set_xlabel('Relative Speedup (Time of LP (Gurobi) / Handcrafted DP)')
    axes[1].set_xlabel('Relative Speedup (Time of LP (Gurobi) / Handcrafted DP)')

    if args.output is not None:
        plt.savefig(f'{args.output}_relative_speedup.pdf')

    plt.tight_layout()

    cpp_df['cpp_obj'] = cpp_df['objective_value']
    df = pd.merge(cpp_df, python_df, on=['m', 'n', 's', 'c', 'r'])
    df = df[['m', 'n', 's', 'c', 'r', 'cpp_obj', 'lp_obj']]

    fig, ax = plt.subplots(figsize=(7, 5))

    sns.scatterplot(data=df, x='cpp_obj', y='lp_obj', ax=ax)
    ax.set_xlabel('Handcrafted DP Objective Value')
    ax.set_ylabel('LP (Gurobi) Objective Value')

    plt.tight_layout()

    if args.output is not None:
        plt.savefig(f'{args.output}_objective_comparison.pdf')

    # df['error'] = df['lp_obj'] - df['cpp_obj']

    # plot error as histogram
    # fig, ax = plt.subplots(figsize=(7, 5))

    # sns.histplot(data=df['error'], ax=ax)
    # ax.set_xlabel('LP Objective Value - Handcrafted DP Objective Value')
    # ax.set_ylabel('Count')
# 
    # plt.tight_layout()
    # if args.output is not None:
        # plt.savefig(f'{args.output}_objective_error.pdf')

    plt.show()
