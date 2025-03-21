import numpy as np
import pandas as pd
import re
import sqlite3

from decimal import Decimal

def fexp(number):
    (sign, digits, exponent) = Decimal(number).as_tuple()
    return len(digits) + exponent - 1

def fman(number):
    return Decimal(number).scaleb(-fexp(number)).normalize()

def format_as_latex(df, header, sidebar, footer):

    def exp_formatter(x):
        return '${:.1f}\\times10^{{{}}}$'.format(fman(x), fexp(x))

    df_string = df.to_string(header=None, index=None, float_format=exp_formatter)

    df_string = re.sub(' +', ' & ', df_string)        # insert & as separators
    df_string = re.sub('\n', '\\\\\\\n', df_string) # insert \\ 
    df_string = re.sub('\Z', '\\\\\\\n', df_string)   # insert \\ and a newline in the end

    df_lines = re.split('\n', df_string)
    sidebar_lines = re.split('\n', sidebar)

    mid_part = ''

    for sb_line, df_line in zip(sidebar_lines, df_lines):
        mid_part += sb_line + df_line + '\n'

    return header + '\n' + mid_part + '\n' +  footer

def get_df_w_z_results():
    db_name = '../data/ms_data.db'

    point_id = 105 # id of an example point
    with sqlite3.connect(db_name) as db:
        query = """SELECT inf_theory, inf_num FROM RESULTS
                                WHERE point_id == {}
                                AND basis == 'z'
                                AND carrier_opt == {}
                                AND ham_type == {}
                                AND init_spin_2 == 1"""

        df_noopt_ld = pd.read_sql_query(query.format(point_id, "FALSE", "'LD'"),   db)
        df_opt_ld   = pd.read_sql_query(query.format(point_id, "TRUE",  "'LD'"),   db)
        df_noopt_full = pd.read_sql_query(query.format(point_id, "FALSE", "'full'"), db)
        df_opt_full = pd.read_sql_query(query.format(point_id, "TRUE",  "'full'"), db)

    df = pd.DataFrame(
                {'opt' :   [  df_opt_ld['inf_theory'][0],   df_opt_ld['inf_num'][0],   df_opt_full['inf_num'][0]],
                 'noopt' : [df_noopt_ld['inf_theory'][0], df_noopt_ld['inf_num'][0], df_noopt_full['inf_num'][0]]}
            )


    return df


if __name__ == '__main__':

    header =\
r"""
\begin{tabular}{c|c|c}
    \toprule
    & Optimized & Non-optimized\\
    \midrule"""

    sidebar = \
r"""    $1-F_0$ &
    Num. (LD)   &
    Num. (full) &
"""

    footer = """\\bottomrule \n\end{tabular}"""

    df = get_df_w_z_results()


    print(format_as_latex(df, header, sidebar, footer))
