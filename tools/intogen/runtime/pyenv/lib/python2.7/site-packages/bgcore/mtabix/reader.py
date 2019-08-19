import os
import struct
import zlib
import collections

# Block gzip constants
BLOCK_HEADER_LENGTH = 18
BLOCK_LENGTH_OFFSET = 16

# Virtual file pointer constants
SHIFT_AMOUNT = 16
OFFSET_MASK = 0xffff
ADDRESS_MASK = 0xFFFFFFFFFFFF


class BlockReader:
    def __init__(self, filename, cache_size=1000):
        """

        :type index: Index
        :param filename:
        :param index:
        :raise IOError:
        """

        # Check that the data file exists
        if not os.path.isfile(filename):
            raise IOError("The path '{}' is not a file".format(filename))

        self.__filename = filename
        self.__data = open(self.__filename, "rb")
        self.__header = None

        self.__partial_line_ends = {}
        self.__blocks_cache = LRUCache(cache_size)

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def close(self):
        self.__header = None
        if self.__data is not None:
            self.__data.close()

    def read(self, block):

        block_address = (block >> SHIFT_AMOUNT) & ADDRESS_MASK
        block_offset = int(block & OFFSET_MASK)

        block_length, block_content = self.__read_block_bytes(block_address)

        # Remove offset
        block_content = block_content[block_offset:]

        # Cache of partial line ends
        self.__partial_line_ends[block_address] = str(block_content[:block_offset])

        # Check if the last line is partial
        is_partial_last_line = block_content[len(block_content) - 1] != '\n'

        # Split lines
        lines = block_content.splitlines()

        # Complete last line
        if is_partial_last_line:

            next_block_address = block_address + block_length
            if next_block_address not in self.__partial_line_ends:
                # Force to read next block
                next_block_length, next_block_content = self.__read_block_bytes(next_block_address)
                next_block_offset = next_block_content.find('\n')
                self.__partial_line_ends[next_block_address] = str(next_block_content[:next_block_offset])

            complete_line = lines.pop() + self.__partial_line_ends[next_block_address]
            lines.append(complete_line)

        return lines

    def __read_block_bytes(self, block_address):

        # First check cache
        if self.__blocks_cache.has_key(block_address):
            block = self.__blocks_cache.get(block_address)
            return block[0], block[1]

        self.__data.seek(block_address)

        # Read the block header
        header = self.__data.read(BLOCK_HEADER_LENGTH)
        # Extract compressed block length
        block_compressed_length = struct.unpack_from("H", header, offset=BLOCK_LENGTH_OFFSET)[0] + 1

        # Read compressed block
        block_compressed = header + self.__data.read(block_compressed_length - BLOCK_HEADER_LENGTH)

        # Decompress
        lines = zlib.decompress(block_compressed, 15 + 32)

        # Decode string
        lines = lines.decode("utf-8")

        self.__blocks_cache.set(block_address, (block_compressed_length, lines))

        return block_compressed_length, lines

    def header(self):
        if self.__header is None:
            self.__header = self.read(0)[0]

        return self.__header


class LRUCache:
    def __init__(self, capacity):
        self.capacity = capacity
        self.cache = collections.OrderedDict()

    def get(self, key):
        try:
            value = self.cache.pop(key)
            self.cache[key] = value
            return value
        except KeyError:
            return -1

    def set(self, key, value):
        try:
            self.cache.pop(key)
        except KeyError:
            if len(self.cache) >= self.capacity:
                self.cache.popitem(last=False)
        self.cache[key] = value

    def has_key(self, key):
        return key in self.cache