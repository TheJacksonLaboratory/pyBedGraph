from Window import Window

START_INDEX = 1
END_INDEX = 2
VALUE_INDEX = 3


class Chromosome:

    def __init__(self, name):
        self.name = name
        self.data_index = 0

        self.last_index = -1

        self.windows = []

    def add_data(self, data):
        window = Window(self.name, int(data[START_INDEX]), int(data[END_INDEX]),
                        float(data[VALUE_INDEX].strip()))

        self.last_index = window.end

        self.windows.append(window)
