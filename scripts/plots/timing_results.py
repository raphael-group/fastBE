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
        'm', 'n', 's', 'c', 'r', 'gurobi_lp_time_without_building', 'gurobi_lp_time_with_building', 
        'cplex_lp_time_without_building', 'cplex_lp_time_with_building',
        'gurobi_dual_time_without_building', 'gurobi_dual_time_with_building', 'dp_time',
    ]]

    df = pd.merge(cpp_df, python_timing_df, on=['m', 'n', 's', 'c', 'r'])
    df = df.melt(
        id_vars=['m', 'n', 's', 'c', 'r'], 
        value_vars=['cpp_dp_time', 
                    'gurobi_lp_time_without_building', 
                    'gurobi_lp_time_with_building', 
                    'gurobi_dual_time_without_building', 
                    'gurobi_dual_time_with_building', 
                    'cplex_lp_time_without_building',
                    'cplex_lp_time_with_building',
                    'dp_time'], 
        var_name='algorithm',  value_name='time (s)'
    )

    df = df[df['algorithm'].isin(['gurobi_lp_time_without_building', 'cplex_lp_time_without_building', 'cpp_dp_time'])]

    df['algorithm'] = df['algorithm'].replace({
        'cpp_dp_time': 'Handcrafted DP',
        'dual_time_without_building': 'Dual LP (Gurobi)',
        'gurobi_lp_time_without_building': 'LP (Gurobi)',
        'cplex_lp_time_without_building': 'LP (CPLEX)',
    })

    fig, ax = plt.subplots(figsize=(7, 5))

    df['(n, s)'] = df.apply(lambda row: (row["n"], row["s"]), axis=1)

    col_order = df.apply(lambda row: (row["n"], row["s"]), axis=1).sort_values().unique()

    sns.boxplot(data=df, x='(n, s)', y='time (s)', hue='algorithm', ax=ax, order=col_order)

    ax.set_xlabel('(# Clones, # Samples)')
    ax.set_ylabel('Time (s)')
    ax.set_yscale('log')
    ax.legend(title='Algorithm', loc='upper right')

    for tick in ax.get_xticklabels():
        tick.set_rotation(60)

    plt.tight_layout()

    if args.output is not None:
        plt.savefig(f'{args.output}_vs_parameters.pdf')

    fig, axes = plt.subplots(ncols=2, figsize=(10, 5))

    speedup_df = python_timing_df.merge(cpp_df, on=['m', 'n', 's', 'c', 'r'])
    speedup_df['gurobi_speedup_no_build'] = speedup_df['gurobi_lp_time_without_building'] / speedup_df['cpp_dp_time'] 
    speedup_df['gurobi_speedup_with_build'] = speedup_df['gurobi_lp_time_with_building'] / speedup_df['cpp_dp_time']
    speedup_df['cplex_speedup_no_build'] = speedup_df['cplex_lp_time_without_building'] / speedup_df['cpp_dp_time'] 
    speedup_df['cplex_speedup_with_build'] = speedup_df['cplex_lp_time_with_building'] / speedup_df['cpp_dp_time']
    speedup_df['(n, s)'] = speedup_df.apply(lambda row: f'({row["n"]}, {row["s"]})', axis=1)

    speedup_df = speedup_df.melt(
        id_vars=['m', 'n', 's', 'c', 'r'], 
        value_vars=['gurobi_speedup_no_build', 'gurobi_speedup_with_build', 'cplex_speedup_no_build', 'cplex_speedup_with_build'],
        var_name='algorithm', value_name='speedup'
    )

    sns.histplot(
        data=speedup_df[speedup_df['algorithm'].str.contains('with_build')].replace({'algorithm': {
            'gurobi_speedup_with_build': 'LP (Gurobi)', 'cplex_speedup_with_build': 'LP (CPLEX)'}
        }),
        x='speedup', hue='algorithm', ax=axes[0]
    )

    sns.histplot(
        data=speedup_df[speedup_df['algorithm'].str.contains('no_build')].replace({'algorithm': {
            'gurobi_speedup_no_build': 'LP (Gurobi)', 'cplex_speedup_no_build': 'LP (CPLEX)'}
        }),
        x='speedup', hue='algorithm', ax=axes[1]
    )

    # axes[0].set_title('Relative Runtime With Model Building Time')
    # axes[1].set_title('Relative Runtime Without Model Building Time')

    axes[0].set_xlabel('Relative Runtime (Time of LP / Handcrafted DP)')
    axes[1].set_xlabel('Relative Runtime (Time of LP / Handcrafted DP)')

    axes[0].get_legend().set_title('Algorithm')
    axes[1].get_legend().set_title('Algorithm')

    if args.output is not None:
        plt.savefig(f'{args.output}_relative_speedup.pdf')

    plt.tight_layout()

    cpp_df['cpp_obj'] = cpp_df['objective_value']
    df = pd.merge(cpp_df, python_df, on=['m', 'n', 's', 'c', 'r'])
    df = df[['m', 'n', 's', 'c', 'r', 'cpp_obj', 'cplex_lp_obj']]

    fig, ax = plt.subplots(figsize=(7, 5))

    sns.scatterplot(data=df, x='cpp_obj', y='cplex_lp_obj', ax=ax)
    ax.set_xlabel('Handcrafted DP Objective Value')
    ax.set_ylabel('LP (CPLEX) Objective Value')

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
