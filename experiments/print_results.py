import json

nb_path = "cic3_seeding_strategies.ipynb"
with open(nb_path, "r") as f:
    nb = json.load(f)

for cell in nb["cells"]:
    if cell["cell_type"] == "code":
        for out in cell.get("outputs", []):
            if out.get("name") == "stdout":
                text = "".join(out.get("text", []))
                if "Topology    Strategy" in text:
                    print(text)
