import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from typing import Dict, List, Tuple
from typing import Optional

from PIL import Image


def get_group_col_names(df: pd.DataFrame, groups_tsv: str) -> Tuple[List[str], List[List[str]], dict]:
    # create list of column names per "analysis-group" as defined in the groups TSV. Also return group names and dict for later use.

    mapping_df = pd.read_csv(groups_tsv, sep='\t', dtype=str)
    # trim whitespace to avoid accidental mismatches
    if 'sample' not in mapping_df.columns or 'group' not in mapping_df.columns:
        raise ValueError("Group mapping TSV must contain columns 'sample' and 'group'")

    mapping_df['sample'] = mapping_df['sample'].astype(str).str.strip()
    mapping_df['group'] = mapping_df['group'].astype(str).str.strip()

    groups_dict = mapping_df.groupby('group', sort=False)['sample'].apply(list).to_dict()
    cols = list(df.columns)

    group_names = list(groups_dict.keys())
    resolved_groups: List[List[str]] = []
    for samples in groups_dict.values():
        samples_clean = [str(s).strip() for s in samples if isinstance(s, str)]
        if '' in samples_clean:
            matched = cols.copy()
        else:
            matched = [s for s in samples_clean if s in df.columns]
        # preserve order and remove duplicates
        uniq = list(dict.fromkeys(matched))
        resolved_groups.append(uniq)

    return group_names, resolved_groups, groups_dict


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


def plot_donut(
    df: pd.DataFrame,
    groups: dict,
    metabolic_labels: List[str],
    output_path: str,
    background_path: str,
    background_scale: float = 1.0,
    background_offset: Tuple[float, float] = (0.0, 0.0),
    group_col_names: List[List[str]] | None = None,
):
    #create figure
    fig = plt.figure(figsize=(20, 20))

    #prepare angular ticks
    num_segments = len(metabolic_labels)
    angle_ticks = np.linspace(0, 2 * np.pi, num_segments + 1)

    #optional background image positioned relative to the full figure
    if background_path and os.path.exists(background_path):
        scale = float(background_scale) if background_scale > 0 else 1.0
        width = 1.0 * scale
        height = 1.0 * scale

        ox, oy = background_offset
        ox = float(ox)
        oy = float(oy)

        #center the background in figure coordinates and apply offsets (offsets are in figure fractions)
        left = (1.0 - width) / 2.0 + ox
        bottom = (1.0 - height) / 2.0 + oy

        bg_ax = fig.add_axes([left, bottom, width, height], zorder=0)
        bg_img = plt.imread(background_path)

        #preserve image aspect and center the image
        img_h, img_w = bg_img.shape[:2]
        img_aspect = float(img_w) / float(img_h)
        cell_aspect = float(width) / float(height) if height != 0 else 1.0

        if img_aspect > cell_aspect:
            new_width = width
            new_height = width / img_aspect
            new_left = left
            new_bottom = bottom + (height - new_height) / 2
        else:
            new_height = height
            new_width = height * img_aspect
            new_bottom = bottom
            new_left = left + (width - new_width) / 2

        bg_ax.set_position([new_left, new_bottom, new_width, new_height])
        bg_ax.imshow(bg_img, aspect='auto', origin='upper', alpha=1)
        bg_ax.axis('off')

    #polar axes for donut plot 
    ax = fig.add_subplot(1, 1, 1, projection='polar', zorder=1)
    ax.set_facecolor('none')

    #start at top and go clockwise
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)

    ax.set_rmax(1.05)       #outer rim beyond 1 to see plot lines
    ax.set_rmin(-0.05)      #inner rim inferior 0 to see plot lines
    ax.set_rorigin(-5)      #donut hole size
    ax.set_thetamin(0)
    ax.set_thetamax(360)
    ax.set_rticks([0, 0.25, 0.5, 0.75, 1])  # % markers
    ax.set_rlabel_position(-22.5)
    ax.set_yticklabels([''] * len(ax.get_yticks()))
    ax.set_xticks(angle_ticks)
    ax.set_xticklabels([''] * len(angle_ticks))

    #use equal box aspect for uniform radial scaling
    ax.set_aspect('equal')

    mid_angles = (angle_ticks[:-1] + angle_ticks[1:]) / 2
    inner_labels = [str(i) for i in range(1, num_segments + 1)]
    inner_r = ax.get_rmax() * -0.18
    outer_r = ax.get_rmax() * 1.12

    def split_label_near_middle(label: str) -> str:
        s = label.strip()
        if not s:
            return s
        spaces = [i for i, ch in enumerate(s) if ch.isspace()]
        if not spaces:
            return s
        mid = len(s) // 2
        best = min(spaces, key=lambda i: abs(i - mid))
        return s[:best].rstrip() + '\n' + s[best + 1:].lstrip()

    #get outer label rotation right
    for angle, outer_label, inner_label in zip(mid_angles, metabolic_labels, inner_labels):
        rotation = np.rad2deg(-angle)
        if (angle > (np.pi / 2.0)) and (angle < (3.0 * np.pi / 2.0)):
            rotation = (rotation + 180) % 360
            if rotation > 180:
                rotation -= 360

        wrapped_outer = split_label_near_middle(outer_label)
        #outer labels
        ax.text(angle, 
                outer_r, 
                wrapped_outer, 
                ha='center', 
                va='center', 
                fontsize=9, 
                rotation=rotation, 
                rotation_mode='anchor'
                )
        #inner labels
        ax.text(angle, 
                inner_r, 
                inner_label, 
                ha='center', 
                va='center', 
                fontsize=9, 
                fontweight='bold', 
                rotation=0, 
                rotation_mode='anchor'
                )

    #prepare radial series and draw filled bands + median lines
    angles = np.linspace(0, 2 * np.pi, 1000)
    segment_angles = np.linspace(0, 2 * np.pi, num_segments + 1)
    colors = ["#2169A7", 
              "#B5704B",
              "#9C4486",
              "#2F6B3E", 
              "#F0CA20"
              ]
    n_steps = 30 #steps for min/max shading from median origin

    #If caller supplied precomputed column lists use them; otherwise resolve from `groups` dict by exact sample names
    if group_col_names is None:
        group_col_names = []
        for gname, samples in groups.items():
            matched = []
            if samples is None:
                group_col_names.append([])
                continue
            for s in samples:
                if isinstance(s, str) and s in df.columns:
                    matched.append(s)
            #preserve order and remove duplicates
            seen = set()
            uniq = [c for c in matched if not (c in seen or seen.add(c))]
            group_col_names.append(uniq)

    group_names = list(groups.keys())
    for grp_idx, (gname, cols) in enumerate(zip(group_names, group_col_names)):
        median_series = df[cols].median(axis=1).to_numpy(dtype=float)
        min_series = df[cols].min(axis=1).to_numpy(dtype=float)
        max_series = df[cols].max(axis=1).to_numpy(dtype=float)

        #radial arrays sampled densely across angles
        r_median = np.zeros_like(angles)
        r_min_arr = np.zeros_like(angles)
        r_max_arr = np.zeros_like(angles)
        for i in range(len(segment_angles) - 1):
            mask = (angles >= segment_angles[i]) & (angles < segment_angles[i + 1])
            r_median[mask] = median_series[i]
            r_min_arr[mask] = min_series[i]
            r_max_arr[mask] = max_series[i]
        r_median[angles >= segment_angles[-1]] = median_series[-1]
        r_min_arr[angles >= segment_angles[-1]] = min_series[-1]
        r_max_arr[angles >= segment_angles[-1]] = max_series[-1]

        #upper/lower minmax fill bands
        for step in range(n_steps):
            frac = (step + 1) / n_steps
            alpha = 0.03 * (1 - step / n_steps)
            r_upper = r_median + (r_max_arr - r_median) * frac
            ax.fill_between(angles, r_median, r_upper, color=colors[grp_idx % len(colors)], alpha=alpha)
        for step in range(n_steps):
            frac = (step + 1) / n_steps
            alpha = 0.03 * (1 - step / n_steps)
            r_lower = r_median - (r_median - r_min_arr) * frac
            ax.fill_between(angles, r_lower, r_median, color=colors[grp_idx % len(colors)], alpha=alpha)

        ax.plot(angles, r_median, color=colors[grp_idx % len(colors)], label=gname)

    ax.grid(True, color='lightgrey', linewidth=0.5)
    ax.legend(prop={'size': 12})

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close(fig)


def combine_images_side_by_side(
    left_path: str,
    right_path: str,
    output_path: str,
    padding: int = 20,
    match_height: str = 'left',  # 'left' | 'right' | 'max' | 'min'
    bg_color: Tuple[int, int, int] = (255, 255, 255),
) -> str:
# this script combines the two output images of the table and the donut side-by-side for easier viewing.
  
    os.makedirs(os.path.dirname(output_path) or '.', exist_ok=True)

    left = Image.open(left_path).convert('RGB')
    right = Image.open(right_path).convert('RGB')

    h_left, w_left = left.size[1], left.size[0]
    h_right, w_right = right.size[1], right.size[0]

    if match_height == 'left':
        target_h = h_left
    elif match_height == 'right':
        target_h = h_right
    elif match_height == 'max':
        target_h = max(h_left, h_right)
    else:  # 'min' or fallback
        target_h = min(h_left, h_right)

    def scale_to_height(img: Image.Image, target_h: int) -> Image.Image:
        w, h = img.size
        if h == target_h:
            return img
        new_w = int(round(w * (target_h / h)))
        return img.resize((new_w, target_h), Image.LANCZOS)

    left = scale_to_height(left, target_h)
    right = scale_to_height(right, target_h)

    total_w = left.size[0] + padding + right.size[0]
    out = Image.new('RGB', (total_w, target_h), color=bg_color)
    out.paste(left, (0, 0))
    out.paste(right, (left.size[0] + padding, 0))

    out.save(output_path, dpi=(300, 300))
    return output_path