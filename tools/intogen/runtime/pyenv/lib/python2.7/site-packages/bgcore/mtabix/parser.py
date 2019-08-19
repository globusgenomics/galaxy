

class KeyParser:
    def __init__(self, separator='\t', row=1, column=0):
        """

        :param separator: Field separator to split one line in several values
        :param row: The position in a line of the row identifier (starting at zero)
        :param column: The position in a line of the column identifier (starting at zero)
        """
        self.separator = separator
        self.row = row
        self.column = column

    def keys(self, line):
        """

        Return the row and column identifiers of this line

        :param line:
        :return:
        """
        values = line.split(self.separator)
        return values[self.row], values[self.column]