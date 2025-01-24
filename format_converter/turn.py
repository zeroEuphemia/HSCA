import os, sys

def format_translate(from_path, to_path):
    output_file = open(to_path, "w")
    output_file.write("Model apache\n\nParameters:\n")

    with open(from_path, "r") as input_file:
        lines = input_file.readlines()
    
    state, var_index, value_nums = "", 0, []
    for line in lines:
        line = line.replace("\n", "")
        if line == "[Parameter]":
            state = "Parameter"
            continue
        elif line == "[Constraint]":
            state = "Constraint"
            output_file.write("\nConstraints:\n")
            continue
        if state == "Parameter":
            num = line.count(",") + 1
            if num <= 1:
                continue
            value_nums.append(num)
            if num == 2:
                output_file.write(f"p{var_index}: Boolean\n")
            else:
                output_file.write(f"p{var_index}: {'{'}")
                for i in range(num):
                    output_file.write(f" v{i}")
                output_file.write(" }\n")
            var_index += 1

        elif state == "Constraint":
            arr, first = line.split(), True
            output_file.write("# ")
            for item in arr:
                if item == "||":
                    continue
                if '!' in item:
                    item_1 = item.split("!")[0]
                else:
                    item_1 = item.split("=")[0]
                item_2 = item.split("=")[1]
                idx = int(item_1[1:])
                # print(item_1, idx, item_2, end = " ")

                if first:
                    first = False
                else:
                    output_file.write(" or ")

                if value_nums[idx] == 2:
                    if item_2 == "1":
                        output_file.write(f"!{item_1}")
                    else:
                        output_file.write(item_1)
                else:
                    output_file.write(f"{item_1}!=v{item_2}")
            output_file.write(" #\n")
            # print()
    # print(len(value_nums), value_nums)
    output_file.close()
    return

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 turn.py from_path to_path")
        exit(0)
    from_path = sys.argv[1]
    to_path = sys.argv[2]
    format_translate(from_path, to_path)
