import gzip
import os


class Index:
    def __init__(self, index, row, column):

        # Check that the index file exists
        if not os.path.isfile(index):
            raise IOError("Index not found at '{}'".format(index))

        input = gzip.open(index)

        # Read the rows identifiers
        self.rows, self.__row_positions = self.__read_identifiers(input)

        # Read the columns identifiers
        self.columns, self.__column_positions = self.__read_identifiers(input)

        # Read the blocks
        self.__blocks = self.__read_blocks(input)

        # Ordered by rows
        self.sorted_by_rows = row < column

    def blocks(self, rows, columns):

        result = set()

        # It's important to compare the block keys in the sorted order
        firsts = rows if self.sorted_by_rows else columns
        first_positions = self.__row_positions if self.sorted_by_rows else self.__column_positions
        first_block_index = 1 if self.sorted_by_rows else 2

        seconds = columns if self.sorted_by_rows else rows
        second_positions = self.__column_positions if self.sorted_by_rows else self.__row_positions
        second_block_index = 2 if self.sorted_by_rows else 1

        for first in firsts:
            first_position = first_positions[first]

            # Look for the first possible block
            #TODO Improve block search performance using a hash or binary search
            start_i = 0
            for i, block in enumerate(self.__blocks):
                if first_position >= block[first_block_index]:
                    start_i = i - 1 if i != 0 else 0
                    break

            for second in seconds:
                second_position = second_positions[second]

                i = start_i
                while i < len(self.__blocks):

                    # Current block
                    block = self.__blocks[i]

                    # Next block position
                    i += 1

                    if first_position >= block[first_block_index]:

                        # Check if we are at the last block
                        if i != len(self.__blocks):

                            # The next block
                            next_block = self.__blocks[i]

                            if first_position <= next_block[first_block_index]:

                                # It's in a next block
                                if first_position == next_block[first_block_index]:
                                    if second_position >= next_block[second_block_index]:
                                        continue

                                # Otherwise it must be in this block
                                result.add(block[0])
                                break
                            else:
                                continue
                        else:

                            # It's the last block it must be in this block
                            result.add(block[0])
                    else:

                        # The blocks are sorted by the first key
                        # if we are here it means that there is no block
                        # with the given first key
                        break

        return result

    @staticmethod
    def __read_identifiers(input):

        for line in input:
            if line == b'# Keys\n':
                break

        result = []
        positions = {}
        p = 0
        for line in input:

            line = line.strip().decode("utf-8")
            if line == "":
                break

            result.append(line)
            positions[line] = p
            p += 1

        return result, positions

    @staticmethod
    def __read_blocks(input):

        for line in input:
            if line == b'# Blocks\n':
                break

        blocks = []
        for line in input:
            line = line.strip().decode("utf-8")
            if line == "":
                break

            values = line.split('\t')

            offset = int(values[0], 16)
            row = int(values[1])
            column = int(values[2])

            blocks.append([offset, row, column])

        return blocks

