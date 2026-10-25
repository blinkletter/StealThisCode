# add_colab.py — insert an "Open in Colab" badge into every notebook
import pathlib, nbformat

USER, REPO, BRANCH = "blinkletter", "StealThisCode", "main"
MARKER = "colab-badge.svg"

for nb_path in pathlib.Path("chapters").rglob("*.ipynb"):
    if ".ipynb_checkpoints" in nb_path.parts:
        continue
    nb = nbformat.read(nb_path, as_version=4)

    # skip notebooks that already have a badge
    if any(MARKER in c.source for c in nb.cells if c.cell_type == "markdown"):
        continue

    url = (f"https://colab.research.google.com/github/"
           f"{USER}/{REPO}/blob/{BRANCH}/{nb_path.as_posix()}")
    badge = nbformat.v4.new_markdown_cell(
        f"[![Open In Colab](https://colab.research.google.com/assets/{MARKER})]({url})"
    )

    # keep a frontmatter cell first, if there is one
    first = nb.cells[0] if nb.cells else None
    pos = 1 if first and first.cell_type == "markdown" and first.source.lstrip().startswith("---") else 0
    nb.cells.insert(pos, badge)

    nbformat.write(nb, nb_path)
    print("added badge:", nb_path)