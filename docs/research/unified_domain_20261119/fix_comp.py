from pathlib import Path
p=Path('stiffgwpy_fast/exact_background.py')
s=p.read_text(encoding='utf8').replace('F_nodes[0] = 0.0\n    for i in range(n - 1):\n        F_nodes[i + 1] = F_nodes[i] + integral[i]', '''F_nodes[0] = 0.0
    compensation = 0.0
    for i in range(n - 1):
        if relative_scale:
            increment = integral[i] - compensation
            F_nodes[i + 1] = F_nodes[i] + increment
            compensation = (F_nodes[i + 1] - F_nodes[i]) - increment
        else:
            F_nodes[i + 1] = F_nodes[i] + integral[i]''')
p.write_text(s,encoding='utf8')

