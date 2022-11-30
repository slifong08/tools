import os, sys
from pathlib import Path

mystery_fruit = sys.argv[1]
mystery_fruit_2 = sys.argv[2]

out = os.path.join(os.getcwd(), "fruit_list.txt")
Path(out).touch()

items = ["apples", "pears", "peaches", mystery_fruit, mystery_fruit_2]

print("noodles are not fruit!")

with open(out, "w") as write_file:
    for i in items:
        line=f"{i}\n"
        write_file.write(line)

    write_file.close()
