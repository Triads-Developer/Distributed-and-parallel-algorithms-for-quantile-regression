import ctypes
import pathlib

if __name__ == "__main__":
    libname = pathlib.Path().absolute() / "libcmult.so"
    c_lib = ctypes.CDLL(libname)
