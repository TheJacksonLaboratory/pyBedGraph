class Window:

    def __init__(self, chromosome_name, start, end, value):
        self.chromosome_name = chromosome_name
        self.start = start
        self.end = end
        self.value = value

        self.size = end - start

        if self.size < 1:
            print("Size is less than 1")
            exit(-1)
