import os

class pythonAPI:

    def __init__(self, namespace="mmg"):
        self.ns = namespace
        self.modules = []

    def writeAPI(self):
        with open(self.ns + ".py", "w") as f:
            f.write(python_header)


python_header = """import ctypes
import os

libdir  = os.path.dirname(os.path.realpath(__file__))
libpath = os.path.join(libdir,"libtest.dylib")

lib = ctypes.CDLL(libpath)
"""

api = pythonAPI()

api.writeAPI()