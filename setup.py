from setuptools import setup

# Version
version = None
with open("./fastq_preprocessor/__init__.py", "r") as f:
    for line in f.readlines():
        line = line.strip()
        if line.startswith("__version__"):
            version = line.split("=")[-1].strip().strip('"')
assert version is not None, "Check version in fastq_preprocessor/__init__.py"

setup(name='fastq_preprocessor',
      version=version,
      description='Fast short read fastq preprocessor with optional contamination removal',
      url='https://github.com/jolespin/fastq_preprocessor',
      author='Josh L. Espinoza',
      author_email='jespinoz@jcvi.org',
      license='GNU AFFERO GENERAL PUBLIC LICENSE',
      packages=["fastq_preprocessor"],
      install_requires=[
      "genopype >=2023.5.15",
      "soothsayer_utils >=2022.6.24",
      "pandas >=0.24",
      "numpy",
      "tqdm",
      "scandir",
      ],
    include_package_data=True,
     scripts=[
         "fastq_preprocessor/fastq_preprocessor",
         "fastq_preprocessor/fastq_preprocessor.py",
         "fastq_preprocessor/fastq_preprocessor_short.py",
         "fastq_preprocessor/fastq_preprocessor_long.py",
         ],

)
