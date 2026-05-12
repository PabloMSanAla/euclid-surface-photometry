from setuptools import setup, find_packages

setup(
    name='euclid-surface-photometry',
    version='0.1.0',
    description='Euclid surface photometry pipeline',
    author='Your Name',
    packages=find_packages(),
    install_requires=[
        'numpy', 'scipy', 'astropy', 'photutils', 'matplotlib', 'pandas', 'pyyaml', 'opencv-python', 'tqdm', 'fitz',
    ],
    entry_points={
        'console_scripts': [
            'euclid-pipeline=euclid_surface_photometry.pipeline:main',
        ],
    },
    include_package_data=True,
    zip_safe=False,
)
