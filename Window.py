import struct
from sys import getsizeof


class Window:

    def __init__(self, start, end, value):
        if start >= end:
            print("Size is less than 1")
            exit(-1)

        #self.start = struct.pack('I', start)
        #self.end = struct.pack('I', end)
        #self.value = struct.pack('f', value)
        #self.size = struct.pack('I', end - start)

        self.start = start
        self.end = end
        self.value = value
        self.size = end - start

