input_file = "martini_v3.0.0.itp"  # Replace with your actual filename


start_line = 316064
end_line = 318589

with open(input_file, "r") as f:
    lines = f.readlines()

with open("1.01modified.itp", "w") as f:
    for i, line in enumerate(lines):
        parts = line.split()
        if start_line <= i + 1 <= end_line and len(parts) == 5 and parts[1] !='TW' and parts[1] != 'SW' and parts[1] != 'W' and parts[1] != 'TQ5':
            try:
                parts[-1] = f"{float(parts[-1]) * 1.01:.6e}"
            except ValueError:
                pass  # Skip lines that don't match expected format
        f.write(" ".join(parts) + "\n")

print(f"Modified file saved as 1.01modified.itp")



with open(input_file, "r") as f:
    lines = f.readlines()

with open("1.04modified.itp", "w") as f:
    for i, line in enumerate(lines):
        parts = line.split()
        if start_line <= i + 1 <= end_line and len(parts) == 5 and parts[1] !='TW' and parts[1] != 'SW' and parts[1] != 'W' and parts[1] != 'TQ5':
            try:
                parts[-1] = f"{float(parts[-1]) * 1.04:.6e}"
            except ValueError:
                pass  # Skip lines that don't match expected format
        f.write(" ".join(parts) + "\n")

print(f"Modified file saved as 1.04modified.itp")

with open(input_file, "r") as f:
    lines = f.readlines()

with open("1.08modified.itp", "w") as f:
    for i, line in enumerate(lines):
        parts = line.split()
        if start_line <= i + 1 <= end_line and len(parts) == 5 and parts[1] !='TW' and parts[1] != 'SW' and parts[1] != 'W' and parts[1] != 'TQ5':
            try:
                parts[-1] = f"{float(parts[-1]) * 1.08:.6e}"
            except ValueError:
                pass  # Skip lines that don't match expected format
        f.write(" ".join(parts) + "\n")

print(f"Modified file saved as 1.08modified.itp")

####################
with open(input_file, "r") as f:
    lines = f.readlines()

with open("1.03modified.itp", "w") as f:
    for i, line in enumerate(lines):
        parts = line.split()
        if start_line <= i + 1 <= end_line and len(parts) == 5 and parts[1] !='TW' and parts[1] != 'SW' and parts[1] != 'W' and parts[1] != 'TQ5':
            try:
                parts[-1] = f"{float(parts[-1]) * 1.03:.6e}"
            except ValueError:
                pass  # Skip lines that don't match expected format
        f.write(" ".join(parts) + "\n")

print(f"Modified file saved as 1.03modified.itp")


with open(input_file, "r") as f:
    lines = f.readlines()

with open("1.06modified.itp", "w") as f:
    for i, line in enumerate(lines):
        parts = line.split()
        if start_line <= i + 1 <= end_line and len(parts) == 5 and parts[1] !='TW' and parts[1] != 'SW' and parts[1] != 'W' and parts[1] != 'TQ5':
            try:
                parts[-1] = f"{float(parts[-1]) * 1.06:.6e}"
            except ValueError:
                pass  # Skip lines that don't match expected format
        f.write(" ".join(parts) + "\n")

print(f"Modified file saved as 1.06modified.itp")



with open("0.995modified.itp", "w") as f:
    for i, line in enumerate(lines):
        parts = line.split()
        if start_line <= i + 1 <= end_line and len(parts) == 5 and parts[1] !='TW' and parts[1] != 'SW' and parts[1] != 'W' and parts[1] != 'TQ5':
            try:
                parts[-1] = f"{float(parts[-1]) * 0.995:.6e}"
            except ValueError:
                pass  # Skip lines that don't match expected format
        f.write(" ".join(parts) + "\n")

print(f"Modified file saved as 0.995modified.itp")



