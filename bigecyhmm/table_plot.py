import os
from typing import Optional
import pandas as pd
import matplotlib.pyplot as plt


def plot_table(display_df: pd.DataFrame, output_path: str = os.path.join('plots', 'group_stats_table.png')) -> None:

    data_frame = display_df.copy()

    #Assume data_stats.py provides exact column names. Fail if not.
    data_frame = data_frame.reset_index(drop=True)
    data_frame.insert(0, 'ID', range(1, len(data_frame) + 1))
    # median columns are produced by data_stats and start with 'Median '
    median_cols = [c for c in data_frame.columns if c.startswith('Median ')]
    desired = ['ID', 'Metabolic Function', 'significant', 'p_bh'] + median_cols
    df = data_frame[desired]

    # create table
    # use the same figure height as the donut plot (20 inches) so combined figures align vertically
    fig, ax = plt.subplots(figsize=(8, 20))
    ax.axis('off')

    table = ax.table(
        cellText=df.values.tolist(),
        colLabels=df.columns.tolist(),
        loc='center',
        cellLoc='center',
        colLoc='center',
    )

    table.auto_set_font_size(False)
    table.set_fontsize(8)
    try:
        table.auto_set_column_width(col=list(range(len(df.columns))))
    except Exception:
        pass

    #draw horizontal bottom edges for every row, and a thicker one for the top row
    try:
        for (row, col), cell in table.get_celld().items():
            cell.set_edgecolor('gray')
            cell.set_linewidth(0.5)
            cell.set_antialiased(False)
            cell.visible_edges = 'B'
            if row == 0:
                cell.set_edgecolor('black')
                cell.set_linewidth(1.0)
                cell.get_text().set_fontweight('bold')
    except Exception:
        # styling failures are non-fatal; leave table as-is
        pass

    plt.subplots_adjust(left=0.05, right=0.95, top=0.98, bottom=0.02)
    
    # draw to populate the renderer
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    try:
        renderer = fig.canvas.get_renderer()
        bbox = table.get_window_extent(renderer)
        bbox_inches = bbox.transformed(fig.dpi_scale_trans.inverted())
        fig.savefig(output_path, dpi=300, bbox_inches=bbox_inches, pad_inches=0.01)
    except Exception:
        fig.savefig(output_path, dpi=300, bbox_inches='tight', pad_inches=0.02)
    plt.close(fig)
