import os
import json

# Define the folder prefixes to search for and the settings template file path
def find_folders_with_dune_module(root_folder):
    """Find directories in the root folder that contain the given file."""
    matching_folders = []

    # Iterate through the directories in the root folder
    for entry in os.scandir(root_folder):
        if entry.is_dir():
            # Check if the specific file exists in the current directory
            file_path = os.path.join(entry.path, "dune.module")
            if os.path.exists(file_path):
                matching_folders.append(entry.path)

    if matching_folders:
        print(f"Found directories containing 'dune.module':")
        for folder in matching_folders:
            print(os.path.basename(folder))

    return matching_folders

def add_folder_to_list(root_folder,folder_list):
    """Add a folder to the list if it is not already present."""
    while True:
        folder = input("Enter the folder name (enter 'q' to quit): ")
        folder_path = os.path.join(root_folder,folder)
        if folder.lower() == 'q':
            break
        if folder not in folder_list and os.path.exists(folder_path):
            folder_list.append(folder_path)
        else:
            print("Folder already included or is not valid.")
    return folder_list

def generate_settings_dict(matching_folders):
    config_dict = {
    "cmake.environment": {
        "CMAKE_COLOR_DIAGNOSTICS": "ON"
    },
    "cmake.configureSettings": {
        f"{os.path.basename(folder)}_DIR": f"../{os.path.basename(folder)}/build-cmake"
        for folder in matching_folders
    }
    }

    return config_dict

def generate_settings_file(vscode_folder, matching_folders):
    config_dict = generate_settings_dict(matching_folders)
    setting_file_path = os.path.join(vscode_folder, "settings.json")
    with open(setting_file_path, "w") as file:
        json.dump(config_dict, file, indent=4)

def create_vscode_folder(destination_folder, matching_folders):
    """Create the .vscode folder and generate the settings file within it."""
    vscode_folder = os.path.join(destination_folder, ".vscode")
    os.makedirs(vscode_folder, exist_ok=True)
    generate_settings_file(vscode_folder, matching_folders)

def generate_workspace_setting(root_folder, matching_folders):
    """Generate a workspace settings file for VS Code based on matching folders."""
    folders_list = [
        {"name": os.path.basename(folder).replace("dune-", ""), "path": f"./{os.path.relpath(folder, root_folder)}"}
        for folder in matching_folders
    ]

    workspace_data = {
        "folders": folders_list,
        "settings": {
            "cmake.configureOnOpen": False
        }
    }

    workspace_file_path = os.path.join(root_folder, "vs.code-workspace")
    try:
        with open(workspace_file_path, "w") as file:
            json.dump(workspace_data, file, indent=4)
    except Exception as e:
        print(f"Error writing workspace settings: {e}")

if __name__ == "__main__":
    root_folder = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../"))
    matching_folders = find_folders_with_dune_module(root_folder)

    to_add_folders = input("Do you want to add more folders? (y/n): ").strip().lower()
    if to_add_folders == "y":
        add_folder_to_list(root_folder,matching_folders)

        print("\n\n\nUpdated directories:")
        for folder in matching_folders:
            print(os.path.basename(folder))

    for folder in matching_folders:
        print(f"Creating .vscode folder in {folder}")
        create_vscode_folder(folder, matching_folders)

    generate_workspace_setting(root_folder, matching_folders)
    print("VS Code workspace setup complete!")
    print(f"Workspace file created at: {os.path.join(root_folder, 'vs.code-workspace')}, use this file to open the workspace in VS Code.")
