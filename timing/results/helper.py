import os

if __name__ == "__main__":
    current_dir = os.getcwd()

    # Iterate over the files in the current directory
    for filename in sorted(os.listdir(current_dir)):
        if os.path.isfile(os.path.join(current_dir, filename)):
            label = None
            if "-O1" in filename:
                label = "Low"
            elif "tree" in filename:
                label = "Mid"
            else:
                label = "High"
            output = "(\"" + filename + "\", " + label + "),"
            print(output)
