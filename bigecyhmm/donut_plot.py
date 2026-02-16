import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from data_stats import get_group_col_names
from typing import Dict, List, Tuple


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

    ax.grid(True)
    ax.legend()

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
