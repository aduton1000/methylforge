"""Enable `python -m methylforge ...` alongside the `methylforge` console script."""
from .cli import main

if __name__ == "__main__":
    import sys
    sys.exit(main())
