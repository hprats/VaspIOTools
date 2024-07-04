import os


class Incar:
    """A class that represents an INCAR file.

    Attributes:
        tags (dict): a dictionary of INCAR tags.

    Examples:
    >>> tags = {'IBRION': 2, 'EDIFF': 1E-05, 'EDIFFG': -0.01}
    >>> my_incar = Incar(tags)
    """

    def __init__(self, tags=None):
        if tags is None:
            self._tags = {}
        else:
            self._tags = tags

    @classmethod
    def from_dict(cls, dct):
        tags = dct['_tags']

        incar = cls(tags=tags)
        return incar

    @property
    def tags(self):
        return self._tags

    @tags.setter
    def tags(self, new_tags):
        if isinstance(new_tags, dict):
            self._tags = new_tags
        else:
            print("new_tags must be a dictionary")

    def add_tag(self, key, value):
        """Add a new tag."""
        if key in self._tags:
            print(f"Tag {key} already exists, value: {self._tags[key]}")
            print(f"Tag not added, use update_tag to modify an existing tag")
        else:
            self._tags[key] = value

    def remove_tag(self, key):
        """Remove an existing tag."""
        if key in self._tags:
            del self._tags[key]
        else:
            print(f"Can't remove tag. Tag {key} does not exist")

    def update_tag(self, key, value):
        """Modify the value of an existing tag.
        If the tag does not exist, it is added."""
        self._tags[key] = value

    def write(self, path):
        """Write the INCAR file."""
        with open(f"{path}/INCAR", 'w') as infile:
            for tag in self._tags:
                infile.write(f"{tag} = {self._tags[tag]}\n")

    @classmethod
    def from_file(cls, path):
        """Read an existing INCAR file and store it as a INCAR object"""
        with open(f"{path}/INCAR", 'r') as infile:
            lines = infile.readlines()
        tags = {}
        for line in lines:
            if '=' in line:
                tags[line.split('=')[0].strip()] = line.split('=')[1].strip()
        incar = cls(tags=tags)
        return incar
