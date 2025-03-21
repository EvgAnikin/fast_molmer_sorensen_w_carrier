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


def exp_formatter(x):
    return '${:.1f}\\times10^{{{}}}$'.format(fman(x), fexp(x))


def format_as_latex(df, header, sidebar, footer):

    df_string = df.to_string(header=None, index=None, float_format=exp_formatter)

    df_string = re.sub(' +', ' & ', df_string)        # insert & as separators
    df_string = re.sub('\n', '\\\\\\\n\n', df_string) # insert \\ and an additional newline in the end of each line
    df_string = re.sub('\Z', '\\\\\\\n', df_string)   # insert \\ and a newline in the end

    df_lines = re.split('\n', df_string)
    sidebar_lines = re.split('\n', sidebar)

    mid_part = ''

    for sb_line, df_line in zip(sidebar_lines, df_lines):
        mid_part += sb_line + df_line + '\n'

    return header + '\n' + mid_part + '\n' +  footer


def get_df_w_x_results():
    db_name = '../data/ms_data.db'

    point_id = 105 # id of an example point

    with sqlite3.connect(db_name) as db:
        query = """SELECT ph_ex_theory, ph_ex_num, spin_flip_theory, spin_flip_num, inf_theory, inf_num FROM RESULTS
                                WHERE point_id == {}
                                AND basis == 'x'
                                AND carrier_opt == {}
                                AND ham_type == {}
                                AND init_spin_2 == {}"""

        df_noopt_ld_p   = pd.read_sql_query(query.format(point_id, "FALSE", "'LD'",   1), db)
        df_opt_ld_p     = pd.read_sql_query(query.format(point_id, "TRUE",  "'LD'",   1), db)
        df_noopt_full_p = pd.read_sql_query(query.format(point_id, "FALSE", "'full'", 1), db)
        df_opt_full_p   = pd.read_sql_query(query.format(point_id, "TRUE",  "'full'", 1), db)

        df_noopt_ld_m   = pd.read_sql_query(query.format(point_id, "FALSE", "'LD'",  -1), db)
        df_opt_ld_m     = pd.read_sql_query(query.format(point_id, "TRUE",  "'LD'",  -1), db)
        df_noopt_full_m = pd.read_sql_query(query.format(point_id, "FALSE", "'full'",-1), db)
        df_opt_full_m   = pd.read_sql_query(query.format(point_id, "TRUE",  "'full'",-1), db)


    opt_theory = [df_opt_ld_p['ph_ex_theory'][0], df_opt_ld_p['spin_flip_theory'][0], df_opt_ld_p['inf_theory'][0],
                  df_opt_ld_m['ph_ex_theory'][0], df_opt_ld_m['spin_flip_theory'][0], df_opt_ld_m['inf_theory'][0]]

    opt_num_ld = [df_opt_ld_p['ph_ex_num'][0], df_opt_ld_p['spin_flip_num'][0], df_opt_ld_p['inf_num'][0],
                  df_opt_ld_m['ph_ex_num'][0], df_opt_ld_m['spin_flip_num'][0], df_opt_ld_m['inf_num'][0]]

    opt_num_full = [df_opt_full_p['ph_ex_num'][0], df_opt_full_p['spin_flip_num'][0], df_opt_full_p['inf_num'][0],
                    df_opt_full_m['ph_ex_num'][0], df_opt_full_m['spin_flip_num'][0], df_opt_full_m['inf_num'][0]]

    noopt_theory = [df_noopt_ld_p['ph_ex_theory'][0], df_noopt_ld_p['spin_flip_theory'][0], df_noopt_ld_p['inf_theory'][0],
                    df_noopt_ld_m['ph_ex_theory'][0], df_noopt_ld_m['spin_flip_theory'][0], df_noopt_ld_m['inf_theory'][0]]

    noopt_num_ld = [df_noopt_ld_p['ph_ex_num'][0], df_noopt_ld_p['spin_flip_num'][0], df_noopt_ld_p['inf_num'][0],
                    df_noopt_ld_m['ph_ex_num'][0], df_noopt_ld_m['spin_flip_num'][0], df_noopt_ld_m['inf_num'][0]]

    noopt_num_full = [df_noopt_full_p['ph_ex_num'][0], df_noopt_full_p['spin_flip_num'][0], df_noopt_full_p['inf_num'][0],
                      df_noopt_full_m['ph_ex_num'][0], df_noopt_full_m['spin_flip_num'][0], df_noopt_full_m['inf_num'][0]]


    df = pd.DataFrame(
                {
                    'opt.theory'     : opt_theory,
                    'opt.num.LD'     : opt_num_ld,
                    'opt.num.full'   : opt_num_full,
                    'noopt.theory'   : noopt_theory,
                    'noopt.num.LD'   : noopt_num_ld,
                    'noopt.num.full' : noopt_num_full, 
                }
            )
    return df


if __name__ == '__main__':

    header = r"""\begin{tabular}{c|c|c|c|c|c|c|c}
    \toprule
    \multicolumn{2}{c|}{} & \multicolumn{3}{|c|}{Optimized} & \multicolumn{3}{|c}{Non-optimized}\\
    \midrule
    \multicolumn{2}{c|}{} & Theory & Num. (LD) & Num. (full) & Theory & Num. (LD) & Num. (full)\\
    \midrule"""

    sidebar = \
r"""\multirow{3}{*}{$|1,1\rangle_x$}  & $P_\mathrm{ph}^s$ &
\cline{2-8}
                                  & $P_\mathrm{flip}^s$   &
\cline{2-8}
                                  &  $1 - F_\mathrm{tot}$ &
\cline{1-8}
\multirow{3}{*}{$|1,-1\rangle_x$}    & $P_\mathrm{ph}^s$ &
\cline{2-8}                                               
                                 & $P_\mathrm{flip}^s$   &
\cline{2-8}                                               
                                 &  $1 - F_\mathrm{tot}$ &"""

    footer = """\\bottomrule \n\end{tabular}"""


    df = get_df_w_x_results()

    print(format_as_latex(df, header, sidebar, footer))
