

class ResultIterator:
    def __init__(self, data, parser, blocks, rows, columns, header=False):
        """

        :type parser: KeyParser
        :type data: BlockReader
        :param data:
        :param parser:
        :param blocks:
        :param rows:
        :param columns:
        """
        self.data = data
        self.parser = parser
        self.blocks = sorted(list(blocks))
        self.rows = set(rows)
        self.columns = set(columns)
        self.header = header

        # Iterator status parameters
        self.__currentBlock = 0
        self.__currentLine = 0
        self.__currentBlockLines = []

    def __iter__(self):
        return self

    def __next__(self):
        """2.6-3.x version"""
        return self.next()

    def next(self):
        """2.5 version"""

        if self.header:
            self.header = False
            return self.data.header()

        if self.__currentLine >= len(self.__currentBlockLines):

            while True:
                if self.__currentBlock >= len(self.blocks):
                    # No more blocks to read
                    raise StopIteration

                # Read next block
                block = self.blocks[self.__currentBlock]
                lines = self.data.read(block)

                # Filter the lines
                self.__currentBlockLines = list(filter(self.filter, lines))
                self.__currentLine = 0

                # Point to the next block
                self.__currentBlock += 1

                # Exit block reading while if this block has lines
                # to return
                if len(self.__currentBlockLines) > 0:
                    break

        # Current line
        line = self.__currentBlockLines[self.__currentLine]

        # Iterate to the next line
        self.__currentLine += 1

        return line

    def filter(self, line):
        row, col = self.parser.keys(line)
        return row in self.rows and col in self.columns

    def read(self, n=0):
        try:
            return next(self) + '\n'
        except StopIteration:
            return ''