import os


# class with phenox directories
class PhenoXPaths:
    def __init__(self, base_dir="~/git/phenox/"):
        self.base_dir = base_dir
        self.data_dir = os.path.join(base_dir, 'data')
        self.src_dir = os.path.join(base_dir, 'phenox')
        self.test_dir = os.path.join(base_dir, 'tests')