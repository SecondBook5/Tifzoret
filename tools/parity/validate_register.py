import sys
from tools.parity.register_schema import load_and_validate
if __name__ == "__main__":
    for p in sys.argv[1:]:
        load_and_validate(p); print(f"OK: {p}")
