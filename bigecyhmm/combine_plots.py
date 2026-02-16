from PIL import Image
import os
from typing import Tuple

# this script combines the two output images of the table and the donut side-by-side for easier viewing.

def combine_images_side_by_side(
    left_path: str,
    right_path: str,
    output_path: str,
    padding: int = 20,
    match_height: str = 'left',  # 'left' | 'right' | 'max' | 'min'
    bg_color: Tuple[int, int, int] = (255, 255, 255),
) -> str:
  
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