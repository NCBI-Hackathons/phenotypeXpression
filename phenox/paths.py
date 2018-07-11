import os


# class with phenox directories
class PhenoXPaths:
    def __init__(self, base_dir="~/git/phenotypeXpression/"):
        self.base_dir = os.path.expanduser(base_dir)
        self.data_dir = os.path.join(self.base_dir, 'data')
        self.src_dir = os.path.join(self.base_dir, 'phenox')
        self.test_dir = os.path.join(self.base_dir, 'tests')