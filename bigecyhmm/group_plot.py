import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from typing import List, Tuple

from PIL import Image


def plot_table(display_df: pd.DataFrame, output_path: str = os.path.join('plots', 'group_stats_table.png')) -> None:
    """ Plot table showcasing statistical analysis for the different comparison.

    Args:
        display_df (pd.DataFrame): DataFrame containing statistical results as columns and metabolic pathway as rows, modified to display
        output_path (str): path to output table png
    """

    data_frame = display_df.copy()

    # Assume display_df from group_analysis.py provides exact column names. Fail if not.
    data_frame = data_frame.reset_index(drop=True)
    data_frame.insert(0, 'ID', range(1, len(data_frame) + 1))
    # Median columns start with 'Median '.
    median_cols = [c for c in data_frame.columns if c.startswith('Median ')]
    desired_cols = ['ID', 'Metabolic Function', 'p', 'p_bh', 'significant'] + median_cols
    df = data_frame[desired_cols]

    # Use the same figure height as the donut plot (20 inches) so combined figures align vertically.
    fig, ax = plt.subplots(figsize=(8, 20))
    ax.axis('off')

    # Generate the dataframe as a table on an ax object.
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

    # Draw horizontal bottom edges for every row, and a thicker one for the top row.
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
    background_offset: Tuple[float, float] = (0.0, 0.0),
    background_scale: float = 1.0,
    group_col_names: List[List[str]] | None = None):
    """ Draw donut plot from DataFrame containing metabolic abundance in sample, groups and background image.

    Args:
        df (pd.DataFrame): DataFrame containg metabolic function as row, samples as columns and showing abundance of metabolic pathway in sample
        groups (dict): dictionary linking group (keys) to sample (values)
        metabolic_labels (list): list of metabolic pathway in input file
        output_path (str): path to output donut plot png
        background_path (str): path to background image for donut plot
        background_offset (tuple): adjust background image position (right, up)
        background_scale (float): adjust background image scale relative to donut
        group_col_names (list): list of lists, each list is linked to a group and contains the samples associated with the group
    """
    #create figure
    fig = plt.figure(figsize=(20, 20))

    #prepare angular ticks
    num_segments = len(metabolic_labels)
    if num_segments != len(df.index):
        raise ValueError(
            f"Segment count mismatch: {num_segments} labels but {len(df.index)} rows in df"
        )
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
    ax.set_rticks([0.25, 0.5, 0.75])  # % markers
    ax.set_yticklabels([''] * len(ax.get_yticks()))
    ax.set_xticks(angle_ticks)
    ax.set_xticklabels([''] * len(angle_ticks))

    #use equal box aspect for uniform radial scaling
    ax.set_aspect('equal')

    mid_angles = (angle_ticks[:-1] + angle_ticks[1:]) / 2
    # Place custom % labels on the boundary between the last and first segments.
    axis_label_angle = angle_ticks[0]
    for r_val, pct_label in zip([0.25, 0.5, 0.75], ['25%', '50%', '75%']):
        ax.text(
            axis_label_angle,
            r_val,
            pct_label,
            ha='center',
            va='center',
            fontsize=8,
            bbox=dict(boxstyle='round,pad=0.25', facecolor='white', alpha=0.8, edgecolor='none'),
            zorder=5,
        )
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
              "#5FA571", 
              "#F0CA20"
              ]
    n_steps = 30 #steps for IQR shading from median outward

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
        if not cols:
            continue

        median_series = df[cols].median(axis=1).to_numpy(dtype=float)
        q25_series = df[cols].quantile(0.25, axis=1).to_numpy(dtype=float)
        q75_series = df[cols].quantile(0.75, axis=1).to_numpy(dtype=float)

        #Build per-segment arrays and separate segments with NaN to avoid boundary bridging artifacts.
        nan_pad = np.array([np.nan])
        seg_a_parts, seg_med_parts, seg_q25_parts, seg_q75_parts = [], [], [], []
        for i in range(num_segments):
            sa_start = segment_angles[i]
            sa_end = segment_angles[i + 1]
            seg_a = angles[(angles > sa_start) & (angles < sa_end)]
            seg_a = np.concatenate([[sa_start], seg_a, [sa_end]])
            seg_med = np.full_like(seg_a, median_series[i])
            seg_q25 = np.full_like(seg_a, q25_series[i])
            seg_q75 = np.full_like(seg_a, q75_series[i])
            seg_a_parts.append(np.concatenate([seg_a, nan_pad]))
            seg_med_parts.append(np.concatenate([seg_med, nan_pad]))
            seg_q25_parts.append(np.concatenate([seg_q25, nan_pad]))
            seg_q75_parts.append(np.concatenate([seg_q75, nan_pad]))

        r_angles = np.concatenate(seg_a_parts)
        r_median = np.concatenate(seg_med_parts)
        r_q25 = np.concatenate(seg_q25_parts)
        r_q75 = np.concatenate(seg_q75_parts)

        #upper/lower IQR fill bands with exponential alpha decay from the median
        for step in range(n_steps):
            frac = (step + 1) / n_steps
            alpha = 0.20 * np.exp(-4 * step / n_steps)
            r_upper = r_median + (r_q75 - r_median) * frac
            ax.fill_between(r_angles, r_median, r_upper, color=colors[grp_idx % len(colors)], alpha=alpha)
        for step in range(n_steps):
            frac = (step + 1) / n_steps
            alpha = 0.20 * np.exp(-4 * step / n_steps)
            r_lower = r_median - (r_median - r_q25) * frac
            ax.fill_between(r_angles, r_lower, r_median, color=colors[grp_idx % len(colors)], alpha=alpha)

        ax.plot(r_angles, r_median, color=colors[grp_idx % len(colors)], label=gname, zorder=3)

    # Per-sample scatter points with group/sample angular separation in each segment.
    n_groups = len(group_col_names)
    if n_groups > 0:
        segment_width = (2 * np.pi) / num_segments
        inter_group_span = segment_width * 0.70
        group_band_width = inter_group_span / n_groups
        intra_sample_span = group_band_width * 0.65

        for seg_idx in range(num_segments):
            for grp_idx, cols in enumerate(group_col_names):
                if not cols:
                    continue
                vals = pd.to_numeric(df.iloc[seg_idx][cols], errors='coerce').dropna().values.astype(float)
                n_vals = len(vals)
                if n_vals == 0:
                    continue

                group_offset = (grp_idx - (n_groups - 1) / 2.0) * group_band_width
                group_center = mid_angles[seg_idx] + group_offset

                if n_vals == 1:
                    thetas = np.array([group_center])
                else:
                    thetas = group_center + np.linspace(-intra_sample_span / 2, intra_sample_span / 2, n_vals)

                # This adds for each sample a dot linked to a function showing the abundance of this function.
                ax.scatter(
                    thetas,
                    vals,
                    s=10,
                    color=colors[grp_idx % len(colors)],
                    alpha=0.75,
                    edgecolors='none',
                    zorder=2,
                )

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
    bg_color: Tuple[int, int, int] = (255, 255, 255)) -> str:
    """ This script combines the two output images of the table and the donut side-by-side for easier viewing.

    Args:
        left_path (str): path to donut plot png
        right_path (str): path to table png file
        output_path (str): path to output combined png
        padding (int): padding to compute total weight
        match_height (str): how to match height
        bg_color (tuple): colour for background

    Returns:
        output_path (str): path to output combined png
    """
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