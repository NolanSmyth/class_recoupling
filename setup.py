from skbuild import setup
import os

# define absolute paths
ROOT_FOLDER = os.path.join(os.path.dirname(os.path.abspath(__file__)), ".")
INCLUDE_FOLDER = os.path.join(ROOT_FOLDER, "include")

def get_class_version(include_folder):
    # Recover the CLASS version
    with open(os.path.join(include_folder, 'common.h'), 'r') as v_file:
        version = ""
        for line in v_file:
            if line.find("_VERSION_") != -1:
                # get rid of the " and the v
                version = line.split()[-1][2:-1]
                break
        return version
        
setup(
    name='classy',
    version=get_class_version(INCLUDE_FOLDER),
    description='Python interface to the Cosmological Boltzmann code CLASS',
    url='http://www.class-code.net',
    include_package_data=True,
    packages=["classy"],
)
