import os


class Incar:

    def __init__(self, tags=None):
        if tags is None:
            self._tags = {}
        else:
            self._tags = tags

    @property
    def tags(self):
        return self._tags

    @tags.setter
    def tags(self, new_tags):
        if isinstance(new_tags, dict):
            self._tags = new_tags
        else:
            print("New tags variable for INCAR file is not a dictionary")

    def add_tag(self, key, value):
        if key in self._tags:
            print(f"Tag {key} already exists, value: {self._tags[key]}")
        else:
            self._tags[key] = value

    def remove_tag(self, key):
        if key in self._tags:
            del self._tags[key]
        else:
            print(f"Can't remove tag. Tag {key} does not exist")

    def update_tag(self, key, value):
        if key in self._tags:
            self._tags[key] = value
        else:
            print(f"Can't update tag. Tag {key} does not exist")

    def write(self, path='.'):
        initial_directory = os.getcwd()
        os.chdir(path)
        f = open('INCAR', 'w')
        for tag in self._tags:
            f.write(f"{tag} = {self._tags[tag]}\n")
        f.close()
        os.chdir(initial_directory)
