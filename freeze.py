import sys
from cx_Freeze import setup, Executable

# Dependencies are automatically detected, but it might need fine tuning.
build_exe_options = {"includes":  ["scipy.special._ufuncs_cxx", "scipy.linalg"]}

# GUI applications require a different base on Windows (the default is for a
# console application).
base = None
# if sys.platform == "win32":
#     base = "Win32GUI"

setup(  name = "GenoWAPv1.2: Post-GWAS Tissue Prioritize",
        version = "1.0.1",
        description = "Prioritize with Tissue-Specific Functional Annotation",
        options = {"build_exe": build_exe_options},
        executables = [Executable("GenoWAP.py", base=base)] )
