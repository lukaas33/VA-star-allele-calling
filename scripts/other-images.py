import os
from pathlib import Path
import sys
path_root = Path(__file__).parents[1]
sys.path.append(str(path_root))

from modules.assets import generate_images

for name in generate_images.image_configs:
    if os.path.exists(f"C:\\Users\\Lucas-PC\\Downloads\\{name}.svg"):
        continue
    print(name)
    os.system(f"python main.py --example {name}")