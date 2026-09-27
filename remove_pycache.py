import os
import shutil

def remove_pycache(root_dir: str = ".") -> None:
    """Recursively remove all __pycache__ directories."""
    for dirpath, dirnames, _ in os.walk(root_dir):
        if "__pycache__" in dirnames:
            pycache_path = os.path.join(dirpath, "__pycache__")
            shutil.rmtree(pycache_path)
            print(f"Removed: {pycache_path}")

if __name__ == "__main__":
    remove_pycache()