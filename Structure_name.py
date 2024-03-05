# %%
import os

def modify_mae_files(directory):
    """
    For each .mae file in the specified directory, write the file's name into line 38.

    Args:
    - directory (str): Path to the directory containing the .mae files.
    """

    for filename in os.listdir(directory):
        if filename.endswith(".mae"):
            file_path = os.path.join(directory, filename)
            
            with open(file_path, 'r') as file:
                content = file.readlines()

            # Modify line 38 with the filename (without the extension)
            for i, line in enumerate(content):
                if i > 1 and content[i-2].strip() == "i_m_ct_format" and content[i-1].strip() == ":::":
                    # Modify the line with the filename (without the extension)
                    base_name = os.path.splitext(filename)[0]
                    content[i] = f' {base_name}\n'
                    break


            with open(file_path, 'w') as file:
                file.writelines(content)

current_directory = os.getcwd()
modify_mae_files(current_directory)


