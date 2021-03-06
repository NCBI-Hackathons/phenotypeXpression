import os


# class with phenox directories
class PhenoXPaths:
    def __init__(self, outprefix: str, base_dir=os.getcwd()):
        assert os.path.exists(os.path.join(base_dir, 'PhenoX.png'))
        self.base_dir = os.path.expanduser(base_dir)
        self.data_dir = os.path.join(self.base_dir, 'data')
        self.src_dir = os.path.join(self.base_dir, 'phenox')
        self.test_dir = os.path.join(self.base_dir, 'tests')
        self.output_dir = os.path.join(self.base_dir, 'output')
        self.outprefix = outprefix