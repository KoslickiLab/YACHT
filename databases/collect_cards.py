"""
Copy the data card (README.md) of every packaged database to OUT_DIR/<prefix>.md and write OUT_DIR/README.md, an
index of the cards (docs/database_cards in the repository).

Usage: collect_cards.py PACKAGE_DIR OUT_DIR
"""
import collections
import os
import shutil
import sys


def main():
    package_dir, out_dir = sys.argv[1:3]
    os.makedirs(out_dir, exist_ok=True)
    cards = collections.defaultdict(dict)
    titles = {}
    for summary in sorted(f for f in os.listdir(package_dir) if f.endswith("_pretrained.summary.json")):
        prefix = summary[: -len("_pretrained.summary.json")]
        name, ani = prefix.rsplit("_", 1)
        card = os.path.join(package_dir, f"{prefix}_pretrained", "README.md")
        shutil.copyfile(card, os.path.join(out_dir, f"{prefix}.md"))
        cards[name][ani] = f"{prefix}.md"
        with open(card) as f:
            paragraphs = [p for p in f.read().split("\n\n")[1:] if p.strip() and not p.startswith("**For YACHT")]
        titles[name] = " ".join(paragraphs[0].split()) if paragraphs else ""
    with open(os.path.join(out_dir, "README.md"), "w") as f:
        f.write(
            "# Data cards of the YACHT pretrained databases\n\n"
            "One data card per database and ANI threshold, as in the database's zip and on its Zenodo record "
            "([YACHT community](https://zenodo.org/communities/yacht)). See [../databases.md](../databases.md) for "
            "an overview.\n\n| database | | ANI thresholds |\n|---|---|---|\n"
        )
        for name in sorted(cards):
            links = " · ".join(f"[{ani}]({cards[name][ani]})" for ani in sorted(cards[name], key=float))
            f.write(f"| `{name}` | {titles[name]} | {links} |\n")
    print(f"{sum(len(v) for v in cards.values())} cards of {len(cards)} databases in {out_dir}")


if __name__ == "__main__":
    main()
