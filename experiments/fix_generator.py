import json

nb_path = "cic3_seeding_strategies.ipynb"
with open(nb_path, "r") as f:
    nb = json.load(f)

for cell in nb["cells"]:
    if cell["cell_type"] == "code":
        source = "".join(cell["source"])
        if "links.append([])" in source and "triangles.append([])" not in source:
            source = source.replace("links.append([])", "links.append([])\n            triangles.append([])")
            cell["source"] = source.splitlines(True)

with open(nb_path, "w") as f:
    json.dump(nb, f, indent=1)

print("Generator fixed!")
