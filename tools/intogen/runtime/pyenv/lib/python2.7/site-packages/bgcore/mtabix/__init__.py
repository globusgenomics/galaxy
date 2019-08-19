from bgcore.mtabix.iterator import ResultIterator
from bgcore.mtabix.parser import KeyParser
from bgcore.mtabix.reader import BlockReader
from bgcore.mtabix.index import Index


class IndexedMatrix:

    def __init__(self, data, index=None, separator='\t', row=1, column=0, cache_size=100):
        """
            MTabix matrix

        :param data: Data file compressed using bgzip
        :param index: MTabix index file
        :param separator: Field separator to the position in a line of the row identifier (starting at zero)
        :param column: The position in a line of the column identifier (starting at zero)
        :raise IOError: If the data or index files are missing or wrong
        """

        # Default index
        if index is None:
            index = data + ".mtabix"

        # Line parser
        self.parser = KeyParser(separator, row, column)

        # Data and index files
        self.index = Index(index, row, column)
        self.data = BlockReader(data, cache_size=cache_size)

    def close(self):
        self.data.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def filter(self, rows=None, columns=None, header=False):
        """

        Returns a iterator over all the lines with the given rows and/or columns identifiers.
        Each line is a list of values

        :param rows:
        :param columns:
        :return:
        """
        if columns is None:
            columns = self.index.columns

        if rows is None:
            rows = self.index.rows

        # Get all the blocks that we need to read
        blocks = self.index.blocks(rows=rows, columns=columns)

        return ResultIterator(self.data, self.parser, blocks, rows, columns, header)
